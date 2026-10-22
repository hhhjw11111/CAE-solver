"""Separate source-backed field maps and line comparisons for the fast matrix."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib.pyplot as plt
import plot_f1_main_results as fe


COLORS = {"LQ": "#566B83", "N0": "#B26438", "NT": "#187C72"}
LINE_FIELDS = [
    ("phi", 1, "Potential (V)"), ("E", 1e-5, "Electric field (kV/cm)"),
    ("Q", 1e-12, r"Joule heat ($10^{12}$ W m$^{-3}$)"),
    ("theta", 1, "Temperature rise (K)"),
    ("n_change", 100, "Electron change (%)"),
    ("p_change", 100, "Hole change (%)"),
    ("ux", 1e9, r"$u_x$ (nm)"), ("uy", 1e9, r"$u_y$ (nm)"), ("uz", 1e9, r"$u_z$ (nm)"),
]
MAPS = [
    ("potential_mid", 0, "phi", 1, "Potential, z=0", "V", "viridis"),
    ("electric_field_mid", 0, "E", 1e-5, "Electric field, z=0", "kV/cm", "viridis"),
    ("joule_heat_mid", 0, "Q", 1e-12, "Joule heat, z=0", r"$10^{12}$ W m$^{-3}$", "magma"),
    ("temperature_mid", 0, "theta", 1, "Temperature rise, z=0", "K", "magma"),
    ("uz_mid", 0, "uz", 1e9, "Displacement, z=0", r"$u_z$ (nm)", "RdBu_r"),
    ("ux_mid", 0, "ux", 1e9, "Displacement, z=0", r"$u_x$ (nm)", "RdBu_r"),
    ("uy_mid", 0, "uy", 1e9, "Displacement, z=0", r"$u_y$ (nm)", "RdBu_r"),
    ("electrons_top", 1, "n_change", 100, "Electrons, top surface", "Change (%)", "RdBu_r"),
    ("electrons_bottom", 2, "n_change", 100, "Electrons, bottom surface", "Change (%)", "RdBu_r"),
    ("mobility_mid", 0, "mu", 1, "Electron mobility, z=0", r"$\mu_n/\mu_{n0}$", "viridis"),
    ("holes_top", 1, "p_change", 100, "Hole carriers, top surface", "Change (%)", "RdBu_r"),
    ("stress_mid", 0, "vm", 1e-6, "von Mises stress, z=0", "MPa", "magma"),
]
METRICS = [
    ("joule_W", 1e3, "Total Joule power (mW)"),
    ("theta_max_K", 1, "Peak temperature rise (K)"),
    ("ring_Q_mean_W_m3", 1e-12, r"Annulus mean Joule heat ($10^{12}$ W m$^{-3}$)"),
    ("ring_E_p95_V_m", 1e-5, "Annulus electric-field P95 (kV/cm)"),
    ("ring_abs_n_relative_mean", 100, "Annulus mean absolute electron change (%)"),
    ("right_mean_ux_m", 1e9, r"Right-face mean $u_x$ (nm)"),
]


def values(data, key):
    if key in ("n_change", "p_change"):
        return data[key[0]] - 1
    if key in ("ux", "uy", "uz"):
        return data["u"][:, {"ux": 0, "uy": 1, "uz": 2}[key]]
    return data[key]


def locate_points(xyz, conn, points, hole):
    """Inverse HEX20 map; pick a deterministic one-sided trace at element faces."""
    co = xyz[conn]
    low, high = co.min(axis=1), co.max(axis=1)
    cells = np.full(len(points), -1, dtype=int)
    natural = np.full((len(points), 3), np.nan)
    errors = []
    for i, target in enumerate(points):
        if hole and np.linalg.norm(target[:2]) < 3e-6 - 1e-14:
            continue
        mask = np.all(target >= low-1e-13, axis=1) & np.all(target <= high+1e-13, axis=1)
        if abs(target[2]) < 1e-19:
            mask &= abs(low[:, 2]) < 1e-19
        for cell in np.flatnonzero(mask):
            # Thin curved boundary layers need a near-node starting point.
            q = fe.NATURAL[np.argmin(np.linalg.norm(co[cell]-target, axis=1))].copy()
            for _ in range(18):
                pos = fe.shapes(q[None])[0] @ co[cell]
                delta = pos-target
                if np.linalg.norm(delta) < 1e-16:
                    break
                jac = fe.derivatives(q[None])[0] @ co[cell]
                q -= np.linalg.solve(jac.T, delta)
                if np.max(np.abs(q)) > 5:
                    break
            err = np.linalg.norm(fe.shapes(q[None])[0] @ co[cell]-target)
            if err < 1e-14 and np.max(np.abs(q)) <= 1+1e-7:
                cells[i], natural[i] = cell, np.clip(q, -1, 1)
                errors.append(float(err))
                break
        if cells[i] < 0:
            raise ValueError(f"Point outside represented mesh, not a hole gap: {target}")
    return cells, natural, max(errors, default=0.)


def sample_line(xyz, conn, state, model, located):
    cells, natural, _ = located
    result = {k: np.full(len(cells), np.nan) for k, _, _ in LINE_FIELDS}
    for cell in np.unique(cells[cells >= 0]):
        selected = cells == cell
        data = fe.evaluate(xyz[conn[cell]], state[conn[cell]], natural[selected], model)
        for key in result:
            result[key][selected] = values(data, key)
    return result


def export_gallery(root, delivery, rows, planes):
    lines = delivery/"01_折线对比图"
    maps = delivery/"02_单幅场图"
    source = delivery/"04_数据表与图源"
    for folder in (lines, maps, source):
        folder.mkdir(parents=True, exist_ok=True)
    original_out = fe.OUT
    fe.OUT = maps
    count = 0
    for name, face, key, scale, title, unit, cmap in MAPS:
        arrays = [values(block[face], key)*scale for block in planes.values()]
        bounds = (min(a.min() for a in arrays), max(a.max() for a in arrays))
        if cmap == "RdBu_r":
            limit = max(abs(bounds[0]), abs(bounds[1]))
            bounds = (-limit, limit)
        for case, blocks in planes.items():
            fig, ax = plt.subplots(figsize=(5.2, 4.3), layout="constrained")
            fe.panel(fig, ax, blocks[face], values(blocks[face], key)*scale,
                     case[5:-3].replace("_", "-")+": "+title, unit, cmap, bounds)
            fe.save(fig, case[5:-3]+"_"+name)
            count += 1
        print("MAPS_EXPORTED="+name, flush=True)
    fe.OUT = lines
    for key, scale, label in METRICS:
        fig, ax = plt.subplots(figsize=(5.4, 3.6), layout="constrained")
        for i, (g, color, marker) in enumerate([("H0", "#566B83", "o"), ("H1", "#B26438", "s")]):
            ax.plot(range(3), [r[key]*scale for r in rows[i*3:i*3+3]],
                    color=color, marker=marker, lw=1.5, label=g)
        ax.set(xticks=range(3), xticklabels=["LQ", "N0", "NT"], ylabel=label)
        ax.ticklabel_format(axis="y", useOffset=False)
        ax.grid(axis="y", color=".92")
        ax.legend()
        ax.margins(x=.1, y=.15)
        fe.save(fig, "metric_"+key)
    t = np.linspace(-10, 10, 201)*1e-6
    cuts = {
        "x_y4_z0": (np.column_stack((t, np.full_like(t, 4e-6), np.zeros_like(t))),
                     r"$x$ ($\mu$m)", r"$y=4\,\mu$m, $z=0$"),
        "y_x0_ztop": (np.column_stack((np.zeros_like(t), t, np.full_like(t, .5e-6))),
                       r"$y$ ($\mu$m)", r"$x=0$, $z=0.5\,\mu$m"),
    }
    cache, qa = {}, {"map_count": count, "metric_plot_count": len(METRICS), "cuts": {}}
    for cut, (points, xlabel, title) in cuts.items():
        traces, records = {}, []
        for row in rows:
            case = row["case"]
            geom = case.split("_")[1]
            m, xyz, conn, state = fe.load(case)
            if (geom, cut) not in cache:
                cache[geom, cut] = locate_points(xyz, conn, points, geom == "H1")
            located = cache[geom, cut]
            d = sample_line(xyz, conn, state, m["model"], located)
            traces[case] = d
            qa["cuts"][case+"/"+cut] = {
                "maximum_mapping_error_m": located[2],
                "material_points": int(np.count_nonzero(located[0] >= 0)),
                "hole_gap_points": int(np.count_nonzero(located[0] < 0)),
            }
            for j, point in enumerate(points):
                records.append({"case": row["label"], "s_um": t[j]*1e6,
                                "x_m": point[0], "y_m": point[1], "z_m": point[2],
                                "in_material": int(located[0][j] >= 0),
                                **{k: "" if not np.isfinite(v[j]) else v[j] for k, v in d.items()}})
        with (source/(cut+".csv")).open("w", newline="", encoding="utf-8-sig") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        for key, scale, label in LINE_FIELDS:
            fig, ax = plt.subplots(figsize=(6.0, 3.9), layout="constrained")
            for case, data in traces.items():
                _, g, model, _ = case.split("_")
                ax.plot(t*1e6, data[key]*scale, color=COLORS[model],
                        ls="--" if g == "H0" else "-", lw=1.3,
                        marker="o" if g == "H0" else None, markevery=25, ms=2.8, label=g+"-"+model)
            ax.set(xlabel=xlabel, ylabel=label, title=title, xlim=(-10, 10))
            ax.ticklabel_format(axis="y", useOffset=False)
            ax.grid(color=".92")
            ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(.5, 1.01), fontsize=7)
            ax.margins(y=.22)
            fe.save(fig, cut+"_"+key)
        print("LINE_CUT_EXPORTED="+cut, flush=True)
    qa.update({
        "spatial_plot_count": len(cuts)*len(LINE_FIELDS),
        "method": "Original HEX20 inverse mapping and evaluation; deterministic one-sided gradients; NaN hole gap.",
        "map_color_scales": "Shared across all six cases for each named map; no exaggerated LQ carrier noise.",
        "csv_units": "Coordinates m (s_um in micrometers); phi V; E V/m; Q W/m3; theta K; n_change/p_change dimensionless; ux/uy/uz m.",
        "visual_review": "PENDING",
    })
    (source/"gallery_qa.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    fe.OUT = original_out
    return qa
