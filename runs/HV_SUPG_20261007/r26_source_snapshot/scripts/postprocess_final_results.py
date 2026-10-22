#!/usr/bin/env python3
"""Create the final source-traceable verification figures for TPSolverV2."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update(
    {
        "font.size": 7,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "lines.linewidth": 1.4,
    }
)


BLUE = "#0F4D92"
BLUE_2 = "#5B8FD6"
TEAL = "#42949E"
VIOLET = "#7C6CCF"
RED = "#B64342"
ORANGE = "#E28E2C"
GREEN = "#2E9E44"
GREY = "#767676"
LIGHT = "#D9D9D9"


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "runs" / "20260902T120024Z_G00_G12_L01"
FINAL = ROOT / "final_results"
DATA_DIR = FINAL / "data"
FIG_DIR = FINAL / "figures"
QA_DIR = FINAL / "qa"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_key_values(path: Path) -> dict[str, float | str]:
    result: dict[str, float | str] = {}
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        value = value.strip()
        try:
            result[key.strip()] = float(value)
        except ValueError:
            result[key.strip()] = value
    return result


def parse_solver_result(path: Path):
    metadata: dict[str, float | int | str] = {}
    continuation_header: list[str] | None = None
    continuation: list[list[str]] = []
    node_header: list[str] | None = None
    nodes: list[list[str]] = []
    mode = "metadata"
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("continuation_columns="):
            continuation_header = line.split("=", 1)[1].split(",")
            mode = "continuation"
            continue
        if line.startswith("node_columns="):
            node_header = line.split("=", 1)[1].split(",")
            mode = "nodes"
            continue
        if mode == "metadata":
            if "=" in line:
                key, value = line.split("=", 1)
                value = value.strip()
                try:
                    metadata[key] = float(value)
                except ValueError:
                    metadata[key] = value
        elif mode == "continuation":
            continuation.append([part.strip() for part in line.rstrip(",").split(",")])
        else:
            nodes.append([part.strip() for part in line.split(",")])
    if continuation_header is None or node_header is None:
        raise ValueError(f"Missing table header in {path}")
    continuation_records = []
    for row in continuation:
        record: dict[str, float | int | str] = {}
        for key, value in zip(continuation_header, row):
            if key == "phase":
                record[key] = value
            elif key in {"point", "newton_iterations"}:
                record[key] = int(value)
            else:
                record[key] = float(value)
        continuation_records.append(record)
    node_array = np.asarray(nodes, dtype=float)
    return metadata, continuation_records, node_header, node_array


def write_csv(path: Path, header: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def add_panel_label(ax, label: str, x: float = -0.12, y: float = 1.04, color="black"):
    text_method = ax.text2D if hasattr(ax, "text2D") else ax.text
    text_method(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=color,
    )


def export_figure(fig, base: Path) -> list[Path]:
    base.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix, kwargs in (
        ("svg", {}),
        ("pdf", {}),
        ("png", {"dpi": 300}),
        ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}}),
    ):
        path = base.with_suffix("." + suffix)
        fig.savefig(path, bbox_inches="tight", facecolor="white", **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def build_analytic_figure(g05: dict[str, float | str]):
    elements = np.asarray([1, 2, 4, 8, 16], dtype=float)
    h = 1.0 / elements
    rows = []
    series = {}
    for benchmark in ("v02", "v06"):
        for norm in ("l2", "h1"):
            values = np.asarray(
                [float(g05[f"{benchmark}_ne{int(ne)}_phi_{norm}"]) for ne in elements]
            )
            series[(benchmark, norm)] = values
            for ne, hi, error in zip(elements, h, values):
                rows.append([benchmark.upper(), int(ne), hi, "phi", norm.upper(), error])
    write_csv(
        DATA_DIR / "analytic_convergence.csv",
        ["benchmark", "elements", "normalized_mesh_size", "field", "norm", "relative_error"],
        rows,
    )

    residual_rows = [
        ["V03 electron zero-current", float(g05["v03_normalized_max_current"])],
        ["V04 hole zero-current", float(g05["v04_normalized_max_current"])],
        ["V05 heat balance", float(g05["v05_uniform_source_energy_balance"])],
        ["V06 coupled zero-current", float(g05["v06_normalized_max_current"])],
    ]
    write_csv(DATA_DIR / "analytic_residuals.csv", ["check", "normalized_residual"], residual_rows)

    fig = plt.figure(figsize=(183 / 25.4, 70 / 25.4), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.9])
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    for ax, benchmark, title in zip(axes[:2], ("v02", "v06"), ("V02 electromechanical", "V06 coupled thermo-electro-carrier")):
        ax.loglog(h, series[(benchmark, "l2")], "o-", color=BLUE, ms=4, label=r"$L^2$")
        ax.loglog(h, series[(benchmark, "h1")], "s-", color=ORANGE, ms=3.5, label=r"$H^1$")
        ref3 = series[(benchmark, "l2")][-1] * (h / h[-1]) ** 3
        ref2 = series[(benchmark, "h1")][-1] * (h / h[-1]) ** 2
        ax.loglog(h, ref3, "--", color=BLUE, alpha=0.45, lw=1)
        ax.loglog(h, ref2, "--", color=ORANGE, alpha=0.45, lw=1)
        order_l2 = float(g05[f"{benchmark}_phi_l2_order"])
        order_h1 = float(g05[f"{benchmark}_phi_h1_order"])
        ax.text(0.04, 0.06, f"observed: $p_{{L^2}}={order_l2:.3f}$\n$p_{{H^1}}={order_h1:.3f}$", transform=ax.transAxes)
        ax.set_xlabel("normalized mesh size, $h/L$")
        ax.set_ylabel("relative error")
        ax.set_title(title, fontsize=7.5, pad=4)
        ax.grid(True, which="both", color=LIGHT, lw=0.45)
        ax.legend(loc="upper left", handlelength=1.6)
        ax.invert_xaxis()

    ax = axes[2]
    labels = ["V03\n$J_n$", "V04\n$J_p$", "V05\nheat", "V06\ncoupled"]
    values = np.asarray([row[1] for row in residual_rows])
    colors = [TEAL, VIOLET, GREEN, BLUE_2]
    ax.bar(np.arange(len(values)), values, color=colors, width=0.68, edgecolor="#333333", linewidth=0.5)
    ax.axhline(1e-12, color=RED, ls="--", lw=1, label=r"$10^{-12}$ reference")
    ax.set_yscale("log")
    ax.set_ylim(1e-17, 3e-12)
    ax.set_xticks(np.arange(len(values)), labels, rotation=0, fontsize=6)
    ax.set_ylabel("normalized residual")
    ax.set_title("Analytic residual checks", fontsize=7.5, pad=4)
    ax.grid(True, axis="y", which="both", color=LIGHT, lw=0.45)
    ax.legend(loc="upper right", fontsize=6)

    for label, ax in zip("abc", axes):
        add_panel_label(ax, label)
    fig.suptitle("Analytic benchmarks recover the expected convergence orders", fontsize=9, fontweight="bold")
    return export_figure(fig, FIG_DIR / "figure_1_analytic_validation")


def masked_triangulation(x_nm: np.ndarray, y_nm: np.ndarray, radius_nm: float):
    tri = mtri.Triangulation(x_nm, y_nm)
    xy = np.column_stack([x_nm, y_nm])
    triangles = tri.triangles
    centers = xy[triangles].mean(axis=1)
    center = np.asarray([300.0, 300.0])
    center_radius = np.linalg.norm(centers - center, axis=1)
    max_edge = np.max(
        np.stack(
            [
                np.linalg.norm(xy[triangles[:, 0]] - xy[triangles[:, 1]], axis=1),
                np.linalg.norm(xy[triangles[:, 1]] - xy[triangles[:, 2]], axis=1),
                np.linalg.norm(xy[triangles[:, 2]] - xy[triangles[:, 0]], axis=1),
            ]
        ),
        axis=0,
    )
    tri.set_mask((center_radius < 0.98 * radius_nm) | (max_edge > 95.0))
    return tri


def spatial_panel(fig, ax, tri, values, title, cbar_label, cmap="viridis", symmetric=False):
    if symmetric:
        bound = max(abs(np.nanmin(values)), abs(np.nanmax(values)))
        levels = np.linspace(-bound, bound, 15)
    else:
        levels = np.linspace(np.nanmin(values), np.nanmax(values), 15)
    image = ax.tricontourf(tri, values, levels=levels, cmap=cmap, extend="both")
    ax.set_aspect("equal")
    ax.set_xlim(0, 600)
    ax.set_ylim(0, 600)
    ax.set_xticks([0, 300, 600])
    ax.set_yticks([0, 300, 600])
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_title(title, fontsize=7.5, pad=3)
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(cbar_label, fontsize=6.5)
    cbar.ax.tick_params(labelsize=5.8)


def build_hole_figure(result_path: Path):
    metadata, continuation, header, nodes = parse_solver_result(result_path)
    columns = {name: index for index, name in enumerate(header)}
    top_z = np.max(nodes[:, columns["z"]])
    top = nodes[np.isclose(nodes[:, columns["z"]], top_z, atol=1e-16)]
    x = top[:, columns["x"]] * 1e9
    y = top[:, columns["y"]] * 1e9
    u1 = top[:, columns["u1"]]
    u2 = top[:, columns["u2"]]
    u3 = top[:, columns["u3"]]
    phi = top[:, columns["phi"]]
    eta_n = top[:, columns["eta_n"]]
    eta_p = top[:, columns["eta_p"]]
    theta = top[:, columns["theta"]]
    electron_pct = 100.0 * np.expm1(eta_n)
    hole_pct = 100.0 * np.expm1(eta_p)
    displacement_pm = np.sqrt(u1**2 + u2**2 + u3**2) * 1e12
    radius_nm = float(metadata["hole_radius"]) * 1e9
    tri = masked_triangulation(x, y, radius_nm)

    top_rows = np.column_stack(
        [
            top[:, columns["node_id"]].astype(int),
            x,
            y,
            top[:, columns["z"]] * 1e9,
            u1 * 1e12,
            u2 * 1e12,
            u3 * 1e12,
            displacement_pm,
            phi,
            theta,
            electron_pct,
            hole_pct,
        ]
    )
    write_csv(
        DATA_DIR / "hole_r090_l2_top_nodes.csv",
        ["node_id", "x_nm", "y_nm", "z_nm", "u1_pm", "u2_pm", "u3_pm", "u_magnitude_pm", "phi_V", "temperature_rise_K", "electron_change_percent", "hole_change_percent"],
        top_rows,
    )
    continuation_rows = [
        [
            row["point"],
            row["phase"],
            row["parameter"],
            row["voltage_factor"],
            row["joule_factor"],
            row["newton_iterations"],
            row["residual_rms"],
            row["update_rms"],
            row["charge_balance"],
            row["electron_balance"],
            row["hole_balance"],
            row["heat_balance"],
        ]
        for row in continuation
    ]
    write_csv(
        DATA_DIR / "hole_r090_l2_continuation.csv",
        ["point", "phase", "parameter", "voltage_factor", "joule_factor", "newton_iterations", "residual_rms", "update_rms", "charge_balance", "electron_balance", "hole_balance", "heat_balance"],
        continuation_rows,
    )

    fig = plt.figure(figsize=(183 / 25.4, 140 / 25.4), constrained_layout=True)
    outer = fig.add_gridspec(1, 2, width_ratios=[0.95, 1.05])
    left = outer[0, 0].subgridspec(2, 1, height_ratios=[3.3, 1.0])
    right = outer[0, 1].subgridspec(2, 2)
    ax3d = fig.add_subplot(left[0, 0], projection="3d")
    exaggeration = 1000.0
    xd = x + exaggeration * u1 * 1e9
    yd = y + exaggeration * u2 * 1e9
    zd = top[:, columns["z"]] * 1e9 + exaggeration * u3 * 1e9
    surf = ax3d.plot_trisurf(xd, yd, zd, triangles=tri.get_masked_triangles(), cmap="inferno", linewidth=0.08, antialiased=True, shade=False)
    face_temperature = theta[tri.get_masked_triangles()].mean(axis=1)
    surf.set_array(face_temperature)
    surf.set_clim(np.min(theta), np.max(theta))
    ax3d.set_xlabel("x (nm)", labelpad=2)
    ax3d.set_ylabel("y (nm)", labelpad=2)
    ax3d.set_zlabel("")
    ax3d.set_title(f"R090 top surface, displacement ×{int(exaggeration)}", fontsize=7.5, pad=5)
    ax3d.view_init(elev=32, azim=-55)
    ax3d.set_proj_type("ortho")
    ax3d.set_box_aspect((1, 1, 0.34), zoom=1.05)
    ax3d.tick_params(labelsize=5.5, pad=0)
    cbar = fig.colorbar(surf, ax=ax3d, orientation="horizontal", fraction=0.055, pad=0.03, shrink=0.78, aspect=28)
    cbar.set_label("temperature rise (K)", fontsize=6.5)
    cbar.ax.tick_params(labelsize=5.8)
    add_panel_label(ax3d, "a", x=-0.02, y=0.98)

    axes = [fig.add_subplot(right[0, 0]), fig.add_subplot(right[0, 1]), fig.add_subplot(right[1, 0]), fig.add_subplot(right[1, 1])]
    spatial_panel(fig, axes[0], tri, phi, "Electric potential", "φ (V)", cmap="cividis")
    spatial_panel(fig, axes[1], tri, electron_pct, "Electron redistribution", r"$\Delta n/n_0$ (%)", cmap="coolwarm", symmetric=True)
    spatial_panel(fig, axes[2], tri, hole_pct, "Hole redistribution", r"$\Delta p/p_0$ (%)", cmap="coolwarm", symmetric=True)
    spatial_panel(fig, axes[3], tri, displacement_pm, "Displacement magnitude", "|u| (pm)", cmap="viridis")
    for label, ax in zip("bcde", axes):
        add_panel_label(ax, label, x=-0.34, y=1.10)

    axconv = fig.add_subplot(left[1, 0])
    points = np.asarray([row["point"] for row in continuation[1:]], dtype=int)
    residual = np.maximum(np.asarray([row["residual_rms"] for row in continuation[1:]], dtype=float), 1e-20)
    update = np.maximum(np.asarray([row["update_rms"] for row in continuation[1:]], dtype=float), 1e-20)
    axconv.semilogy(points, residual, "o-", color=BLUE, ms=2.5, label="residual RMS")
    axconv.semilogy(points, update, "s-", color=ORANGE, ms=2.2, label="update RMS")
    axconv.axhline(1e-9, color=GREY, ls="--", lw=0.8)
    axconv.set_xticks(points)
    axconv.set_xlabel("continuation point", fontsize=6)
    axconv.set_ylabel("scaled norm", fontsize=6)
    axconv.tick_params(labelsize=5.5)
    axconv.grid(True, axis="y", which="both", color=LIGHT, lw=0.35)
    axconv.legend(fontsize=5.5, ncol=2, loc="lower left")
    axconv.text(0.98, 0.90, "8 accepted / 0 rejected", transform=axconv.transAxes, fontsize=5.5, ha="right", va="top")
    add_panel_label(axconv, "f", x=-0.08, y=1.02)

    fig.suptitle("Converged three-dimensional nonlinear hole solution", fontsize=9, fontweight="bold")
    outputs = export_figure(fig, FIG_DIR / "figure_2_nonlinear_hole_fields")
    return outputs, metadata


def build_radius_figure(g12: dict[str, float | str]):
    radii = np.asarray([60.0, 90.0, 120.0])
    tags = ["r060", "r090", "r120"]
    metric_names = [
        "total_joule_power",
        "global_joule_mean",
        "annulus_joule_mean",
        "annulus_joule_p95",
        "annulus_electric_mean",
        "annulus_n_relative_mean",
        "annulus_p_relative_mean",
        "annulus_temperature_rise_mean",
        "annulus_displacement_mean",
    ]
    values = {metric: 100.0 * np.asarray([float(g12[f"g12_{tag}_contrast_{metric}"]) for tag in tags]) for metric in metric_names}
    rows = []
    for i, radius in enumerate(radii):
        rows.append([radius] + [values[metric][i] for metric in metric_names])
    write_csv(DATA_DIR / "radius_matrix_contrasts.csv", ["radius_nm"] + [metric + "_percent" for metric in metric_names], rows)

    fig, axes = plt.subplots(2, 2, figsize=(183 / 25.4, 94 / 25.4), constrained_layout=True)
    axes = axes.ravel()
    axes[0].plot(radii, values["total_joule_power"], "o-", color=RED, label="total power")
    axes[0].plot(radii, values["global_joule_mean"], "s-", color=GREY, label="global mean")
    axes[0].axhline(0, color="#333333", lw=0.7)
    axes[0].set_ylabel("contrast (%)")
    axes[0].set_title("Global Joule response", fontsize=7.5)
    axes[0].legend()

    axes[1].plot(radii, values["annulus_joule_mean"], "o-", color=ORANGE, label="annulus mean")
    axes[1].plot(radii, values["annulus_joule_p95"], "s-", color=RED, label="annulus p95")
    axes[1].plot(radii, values["annulus_electric_mean"], "^-", color=BLUE, label="electric field mean")
    axes[1].axhline(0, color="#333333", lw=0.7)
    axes[1].set_ylabel("contrast (%)")
    axes[1].set_title("Local electrical/Joule concentration", fontsize=7.5)
    axes[1].legend(fontsize=6)

    axes[2].plot(radii, values["annulus_n_relative_mean"], "o-", color=BLUE, label="electron")
    axes[2].plot(radii, values["annulus_p_relative_mean"], "s-", color=VIOLET, label="hole")
    axes[2].axhline(0, color="#333333", lw=0.7)
    axes[2].set_ylabel("contrast (%)")
    axes[2].set_title("Carrier redistribution", fontsize=7.5)
    axes[2].legend()

    axes[3].plot(radii, values["annulus_temperature_rise_mean"], "o-", color=TEAL, label="temperature rise")
    axes[3].plot(radii, values["annulus_displacement_mean"], "s-", color=GREEN, label="displacement")
    axes[3].axhline(0, color="#333333", lw=0.7)
    axes[3].set_ylabel("contrast (%)")
    axes[3].set_title("Thermomechanical response", fontsize=7.5)
    axes[3].legend()

    for label, ax in zip("abcd", axes):
        ax.set_xlabel("hole radius (nm)")
        ax.set_xticks(radii)
        ax.grid(True, color=LIGHT, lw=0.45)
        add_panel_label(ax, label)
    fig.suptitle("Radius-dependent contrasts at −5 V (hole versus matched plain control)", fontsize=9, fontweight="bold")
    return export_figure(fig, FIG_DIR / "figure_3_radius_matrix")


def main():
    # This entry point is still bound to the superseded, folded-hole meshes.
    # Keep historical artifacts intact until fresh evidence and nodal data
    # have been explicitly rebound and reviewed together.
    raise RuntimeError(
        "Historical hole geometry was superseded on 2026-09-06. "
        "Do not regenerate final_results from the pinned old evidence. "
        "See docs/geometry_audit_findings_20260906.md."
    )

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    QA_DIR.mkdir(parents=True, exist_ok=True)

    g05_path = EVIDENCE / "G05_run_1.txt"
    g12_path = EVIDENCE / "G12_run_1.txt"
    hole_path = DATA_DIR / "G12_R090_L2_nodal.dat"
    for path in (g05_path, g12_path, hole_path):
        if not path.exists():
            raise FileNotFoundError(path)

    g05 = parse_key_values(g05_path)
    g12 = parse_key_values(g12_path)
    if g05.get("verification_result") != "PASS" or g12.get("verification_result") != "PASS":
        raise RuntimeError("Verification source log is not PASS.")

    analytic_outputs = build_analytic_figure(g05)
    hole_outputs, hole_metadata = build_hole_figure(hole_path)
    radius_outputs = build_radius_figure(g12)

    all_outputs = analytic_outputs + hole_outputs + radius_outputs
    manifest = {
        "backend": "python-matplotlib",
        "source_files": {
            "G05_run_1.txt": sha256(g05_path),
            "G12_run_1.txt": sha256(g12_path),
            "G12_R090_L2_nodal.dat": sha256(hole_path),
        },
        "hole_run": {
            "run_result": hole_metadata["run_result"],
            "node_count": int(float(hole_metadata["node_count"])),
            "element_count": int(float(hole_metadata["element_count"])),
            "accepted_points": int(float(hole_metadata["accepted_points"])),
            "rejected_attempts": int(float(hole_metadata["rejected_attempts"])),
            "maximum_balance_defect": float(hole_metadata["maximum_balance_defect"]),
        },
        "outputs": {path.name: sha256(path) for path in all_outputs},
    }
    manifest_path = QA_DIR / "figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"FIGURE_MANIFEST={manifest_path}")
    for path in all_outputs:
        print(f"FIGURE={path}")


if __name__ == "__main__":
    main()
