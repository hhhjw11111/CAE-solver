"""Assemble the user-authorized reduced-mesh six-case results without a false mesh pass."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil
import numpy as np
import matplotlib.pyplot as plt
import audit_f1_results as audit
import plot_f1_main_results as fe
from summarize_f1_release import supplemental
from plot_fast_f1_gallery import export_gallery


CASES = [f"FAST_{g}_{m}_L0" for g in ["H0", "H1"] for m in ["LQ", "N0", "NT"]]


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--delivery-root", type=Path, required=True)
    args = p.parse_args()
    root, out = args.run_root.resolve(), args.output_dir.resolve()
    delivery = args.delivery_root.resolve()
    source = delivery/"04_数据表与图源"
    source.mkdir(parents=True, exist_ok=True)
    fe.RUNS = audit.ROOT = root
    fe.INPUTS = root/"evidence/inputs"
    fe.OUT = out
    out.mkdir(parents=True, exist_ok=True)
    checks, rows, planes, hashes = {}, [], {}, {}
    for case in CASES:
        m, xyz, conn, state = fe.load(case)
        s = supplemental(case)
        e = read(root/case/"execution.json")
        row = {"case": case, "label": case.replace("FAST_", "").replace("_L0", "").replace("_", "-"),
               **m, **{k: v for k, v in s.items() if not isinstance(v, dict)},
               "elapsed_seconds": e["elapsed_seconds"]}
        rows.append(row)
        checks[case+"/exit_and_acceptance"] = e["exit_code"] == 0 and m["case_acceptance"] == "PASS_NOT_MESH_ACCEPTANCE"
        checks[case+"/positive_fields"] = bool(np.isfinite(state).all() and np.all(state[:, 6] > -300))
        for key, threshold in [
            ("maximum_stage_balance", 1e-8), ("maximum_accepted_residual", 1e-9),
            ("maximum_accepted_update", 1e-9), ("heat_balance_relative", .01),
            ("electron_balance_relative", .01), ("joule_quadrature_4_5_relative", .001),
        ]:
            checks[case+"/"+key] = bool(np.isfinite(m[key]) and 0 <= m[key] < threshold)
        checks[case+"/source_check"] = s["saved_sample_relative_error"] < 2e-6
        for f in ["nodes.csv", "metrics.txt", "field_samples.csv", "input.nml", "stdout.txt"]:
            hashes[case+"/"+f] = sha(root/case/f)
        mid = fe.plane(xyz, conn, state, m["model"], 0.)
        top = fe.plane(xyz, conn, state, m["model"], .5e-6)
        bottom = fe.plane(xyz, conn, state, m["model"], -.5e-6)
        planes[case] = (mid, top, bottom)
        np.savez_compressed(source/(case+"_midplane.npz"), **mid)
        np.savez_compressed(source/(case+"_top.npz"), **top)
        np.savez_compressed(source/(case+"_bottom.npz"), **bottom)
        print("POSTPROCESSED="+case, flush=True)
    checks["same_mesh_within_geometry"] = True
    for g in ["H0", "H1"]:
        coords = [np.loadtxt(root/f"FAST_{g}_{m}_L0/nodes.csv", delimiter=",", skiprows=1, usecols=range(4))
                  for m in ["LQ", "N0", "NT"]]
        checks["same_mesh_within_geometry"] &= all(np.array_equal(coords[0], c) for c in coords[1:])
    fields = ["label", "nodes", "elements", "joule_W", "theta_max_K", "ring_theta_mean_K",
              "thermal_resistance_K_W", "ring_E_mean_V_m", "ring_E_p95_V_m",
              "ring_Q_mean_W_m3", "ring_Q_p95_W_m3", "ring_abs_n_relative_mean",
              "maximum_abs_n_relative", "right_mean_ux_m", "minimum_mu_n_ratio",
              "heat_balance_relative", "electron_balance_relative", "joule_quadrature_4_5_relative",
              "E_p95_V_m", "E_above_1e6_volume_fraction", "elapsed_seconds"]
    with (root/"main_results.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        w = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    (root/"main_results.json").write_text(json.dumps(rows, indent=2, allow_nan=False), encoding="utf-8")
    for name in ["main_results.csv", "main_results.json"]:
        shutil.copy2(root/name, source/name)
    mid, top, bottom = planes["FAST_H1_NT_L0"]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), layout="constrained")
    specs = [
        (mid, mid["phi"], "Potential, midplane", "V", "viridis"),
        (mid, mid["E"]/1e5, "Electric field, midplane", "kV/cm", "viridis"),
        (mid, mid["Q"]/1e12, "Joule heat, midplane", r"$10^{12}$ W m$^{-3}$", "magma"),
        (mid, mid["theta"], "Temperature rise, midplane", "K", "magma"),
        (mid, mid["u"][:, 2]*1e9, "Displacement $u_z$, midplane", "nm", "RdBu_r"),
        (top, 100*(top["n"]-1), "Electrons, top surface", "Change (%)", "RdBu_r"),
    ]
    for ax, spec, letter in zip(axes.flat, specs, "abcdef"):
        fe.panel(fig, ax, *spec)
        ax.text(-.15, 1.13, letter, transform=ax.transAxes, fontweight="bold")
    fe.save(fig, "H1_NT_multifield")
    limit = max(np.max(np.abs(planes[f"FAST_H1_{m}_L0"][i]["n"]-1))*100
                for m in ["LQ", "N0", "NT"] for i in [1, 2])
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), layout="constrained")
    for col, model in enumerate(["LQ", "N0", "NT"]):
        for row, face in enumerate(["top", "bottom"]):
            d = planes[f"FAST_H1_{model}_L0"][row+1]
            fe.panel(fig, axes[row, col], d, 100*(d["n"]-1), model+", "+face,
                     "Change (%)", "RdBu_r", (-limit, limit))
    fe.save(fig, "H1_model_carrier_comparison")
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.5), layout="constrained")
    specs = [("joule_W", 1e3, "Total Joule power", "mW"),
             ("theta_max_K", 1, "Peak temperature rise", "K"),
             ("ring_Q_mean_W_m3", 1e-12, "Annulus mean Joule heat", r"$10^{12}$ W m$^{-3}$"),
             ("ring_E_p95_V_m", 1e-5, "Annulus electric-field P95", "kV/cm"),
             ("ring_abs_n_relative_mean", 100, "Annulus electron deviation", "%"),
             ("right_mean_ux_m", 1e9, "Right-face mean $u_x$", "nm")]
    for ax, (key, scale, title, unit) in zip(axes.flat, specs):
        for i, (label, color, marker) in enumerate([("H0", "#42658A", "o"), ("H1", "#BB713D", "s")]):
            ax.plot(range(3), [r[key]*scale for r in rows[i*3:i*3+3]], marker=marker, color=color, lw=1, ms=4, label=label)
        ax.set(xticks=range(3), xticklabels=["LQ", "N0", "NT"], title=title, ylabel=unit)
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        ax.margins(x=.12, y=.18)
    axes[0, 0].legend(frameon=False)
    fe.save(fig, "six_case_comparison")
    # Compare geometry using the same NT model and shared scales.
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), layout="constrained")
    for col, (key, scale, title, unit, cmap) in enumerate([
        ("E", 1e-5, "Electric field", "kV/cm", "viridis"),
        ("Q", 1e-12, "Joule heat", r"$10^{12}$ W m$^{-3}$", "magma"),
        ("theta", 1, "Temperature rise", "K", "magma"),
    ]):
        bound = (min(planes[f"FAST_{g}_NT_L0"][0][key].min()*scale for g in ["H0", "H1"]),
                 max(planes[f"FAST_{g}_NT_L0"][0][key].max()*scale for g in ["H0", "H1"]))
        for i, g in enumerate(["H0", "H1"]):
            d = planes[f"FAST_{g}_NT_L0"][0]
            fe.panel(fig, axes[i, col], d, d[key]*scale, g+" NT: "+title, unit, cmap, bound)
    fe.save(fig, "NT_geometry_comparison")
    gallery = export_gallery(root, delivery, rows, planes)
    repeat = root/"REPEAT_H1_NT_L0"
    repeat_data = {}
    if (repeat/"nodes.csv").exists():
        original = root/"FAST_H1_NT_L0"
        a = np.loadtxt(original/"nodes.csv", delimiter=",", skiprows=1)
        b = np.loadtxt(repeat/"nodes.csv", delimiter=",", skiprows=1)
        repeat_data = {
            "file_byte_identity": {n: sha(original/n) == sha(repeat/n) for n in ["metrics.txt", "nodes.csv", "field_samples.csv"]},
            "maximum_absolute_field_difference": np.max(np.abs(a-b), axis=0).tolist(),
            "representative_repeat_executed": True,
        }
        checks["representative_repeat_byte_identical"] = all(repeat_data["file_byte_identity"].values())
        re = read(repeat/"execution.json")
        checks["representative_repeat_accepted"] = re.get("exit_code") == 0 and re.get("status") == "PASS_NOT_MESH_ACCEPTANCE"
    else:
        checks["representative_repeat_byte_identical"] = False
    old = root.parent/"F1_20260906_R0144"
    mesh_data = {}
    for name in ["MESH_H1_NT_L0", "MESH_H1_NT_L1"]:
        path = old/name/"metrics.txt"
        if path.exists():
            previous = audit.ROOT
            audit.ROOT = old
            mesh_data[name] = audit.read_metrics(name)
            audit.ROOT = previous
    comparison = {}
    if len(mesh_data) == 2:
        a, b = mesh_data.values()
        comparison = {k: audit.rel(a[k], b[k]) for k in [
            "joule_W", "theta_max_K", "ring_E_mean_V_m", "ring_Q_mean_W_m3",
            "ring_abs_n_relative_mean", "ring_Q_p95_W_m3"]}
    record = {
        "status": "SIX_CASE_REDUCED_MESH_MATRIX_COMPLETE" if all(checks.values()) else "AUDIT_NEEDS_ATTENTION",
        "numerical_checks_passed": all(checks.values()), "checks": checks,
        "mesh_convergence_status": "NOT_FINAL_THREE_LEVEL_MESH_ACCEPTANCE",
        "scope_change": "User authorized reduced meshes on 2026-09-07 for same-day completion. Same physics and numerical tolerances retained.",
        "repeat": repeat_data, "source_hashes": hashes,
        "gallery": gallery,
        "historical_L0_L1_change": comparison,
        "visual_review": "PENDING",
    }
    (root/"fast_audit.json").write_text(json.dumps(record, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({"status": record["status"], "failed": [k for k, v in checks.items() if not v]}, indent=2))


if __name__ == "__main__":
    main()
