"""Audit and visualize the source-backed F3 ZnO hole-radius sweep.

No PDE is solved here. Existing R=0/R=3 F1 evidence is combined with the new
R=2/R=4 runs. All cloud values are evaluated from the original HEX20 fields;
triangles never bridge the hole.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import plot_f1_main_results as fe
from audit_hex20_boundary_geometry import read_mesh


MODELS = ("LQ", "N0", "NT")
RADII_UM = (0, 2, 3, 4)


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_metrics(path: Path) -> dict[str, object]:
    result: dict[str, object] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line or "," in line:
            continue
        key, value = line.split("=", 1)
        key, value = key.strip(), value.strip()
        try:
            result[key] = float(value.replace("D", "E"))
        except ValueError:
            result[key] = value
    return result


def registry(old_root: Path, new_root: Path, old_inputs: Path, new_inputs: Path) -> dict[tuple[int, str], dict[str, Path | str]]:
    entries: dict[tuple[int, str], dict[str, Path | str]] = {}
    for model in MODELS:
        entries[(0, model)] = {
            "case": f"FAST_H0_{model}_L0", "run": old_root / f"FAST_H0_{model}_L0",
            "mesh": old_inputs / "H0_L0.mesh", "source": "F1_REUSED",
        }
        entries[(3, model)] = {
            "case": f"FAST_H1_{model}_L0", "run": old_root / f"FAST_H1_{model}_L0",
            "mesh": old_inputs / "H1_L0.mesh", "source": "F1_REUSED",
        }
        for radius in (2, 4):
            case = f"F3_R{radius:02d}_{model}_L0"
            entries[(radius, model)] = {
                "case": case, "run": new_root / case,
                "mesh": new_inputs / f"F3_R{radius:02d}_L0.mesh", "source": "F3_NEW",
            }
    return entries


def load(entry: dict[str, Path | str]) -> tuple[dict[str, object], np.ndarray, np.ndarray, np.ndarray]:
    run = Path(entry["run"])
    metrics = read_metrics(run / "metrics.txt")
    if metrics.get("case_acceptance") != "PASS_NOT_MESH_ACCEPTANCE":
        raise ValueError(f"Unaccepted case: {entry['case']}")
    xyz, connectivity, _ = read_mesh(Path(entry["mesh"]))
    nodes = np.genfromtxt(run / "nodes.csv", delimiter=",", names=True, usecols=range(13))
    saved_xyz = np.column_stack([nodes[name] for name in ("x_m", "y_m", "z_m")])
    if len(nodes) != len(xyz) or np.max(np.abs(saved_xyz - xyz)) >= 1.0e-19:
        raise ValueError(f"Mesh/result mismatch: {entry['case']}")
    state = np.column_stack([
        nodes[name] for name in ("ux_m", "uy_m", "uz_m", "phi_V", "eta_n", "eta_p", "theta_K")
    ])
    return metrics, xyz, connectivity, state


def sample_check(entry: dict[str, Path | str], metrics: dict[str, object], xyz: np.ndarray, connectivity: np.ndarray, state: np.ndarray) -> float:
    reference = np.loadtxt(Path(entry["run"]) / "field_samples.csv", delimiter=",", skiprows=1, max_rows=1)
    gauss = np.polynomial.legendre.leggauss(5)[0][0]
    sampled = fe.evaluate(xyz[connectivity[0]], state[connectivity[0]], [[gauss, gauss, gauss]], str(metrics["model"]))
    actual = np.r_[
        sampled["Evec"][0], sampled["Q"][0], sampled["theta"][0], sampled["n"][0], sampled["p"][0],
        sampled["u"][0], sampled["mu"][0], sampled["stress"][0],
    ]
    expected = reference[4:]
    scale = np.maximum(np.abs(expected), np.r_[[1e-3] * 3, 1e-3, 1e-8, 1, 1, [1e-18] * 3, 1, [1e-4] * 6])
    scale[:3] = max(np.linalg.norm(expected[:3]), 1e-3)
    return float(np.max(np.abs(actual - expected) / scale))


def save_figure(fig: plt.Figure, directory: Path, stem: str) -> list[str]:
    outputs = []
    for extension, dpi in (("svg", 450), ("pdf", 450), ("png", 450)):
        path = directory / f"{stem}.{extension}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
        outputs.append(str(path))
    tiff = directory / f"{stem}.tiff"
    fig.savefig(tiff, dpi=600, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    outputs.append(str(tiff))
    plt.close(fig)
    return outputs


def line_panel(ax: plt.Axes, rows: list[dict[str, object]], key: str, factor: float, title: str, ylabel: str) -> None:
    styles = {"LQ": ("#42658A", "o"), "N0": ("#BB713D", "s"), "NT": ("#2C9A67", "^")}
    for model in MODELS:
        selected = sorted((row for row in rows if row["model"] == model), key=lambda row: row["radius_um"])
        ax.plot(
            [row["radius_um"] for row in selected], [float(row[key]) * factor for row in selected],
            color=styles[model][0], marker=styles[model][1], linewidth=1.15, markersize=4, label=model,
        )
    ax.set(title=title, xlabel=r"Hole radius $r$ ($\mu$m)", ylabel=ylabel, xticks=RADII_UM)
    ax.tick_params(direction="in")
    ax.margins(x=0.05, y=0.12)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old-run-root", type=Path, required=True)
    parser.add_argument("--new-run-root", type=Path, required=True)
    parser.add_argument("--old-input-root", type=Path, required=True)
    parser.add_argument("--new-input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    old_root, new_root = args.old_run_root.resolve(), args.new_run_root.resolve()
    output = args.output_dir.resolve()
    figures, sources = output / "figures", output / "source_data"
    figures.mkdir(parents=True, exist_ok=False)
    sources.mkdir(parents=True, exist_ok=False)
    entries = registry(old_root, new_root, args.old_input_root.resolve(), args.new_input_root.resolve())

    rows: list[dict[str, object]] = []
    planes: dict[tuple[int, str], dict[str, np.ndarray]] = {}
    checks: dict[str, bool] = {}
    hashes: dict[str, dict[str, str]] = {}
    for radius in RADII_UM:
        for model in MODELS:
            entry = entries[(radius, model)]
            metrics, xyz, connectivity, state = load(entry)
            check = sample_check(entry, metrics, xyz, connectivity, state)
            name = f"R{radius:02d}_{model}"
            checks[name + "/saved_field_reconstruction"] = check < 2.0e-6
            for key, threshold in (
                ("maximum_stage_balance", 1e-8), ("maximum_accepted_residual", 1e-9),
                ("maximum_accepted_update", 1e-9), ("heat_balance_relative", 0.01),
                ("electron_balance_relative", 0.01), ("joule_quadrature_4_5_relative", 0.001),
            ):
                checks[name + "/" + key] = 0.0 <= float(metrics[key]) < threshold
            checks[name + "/low_field_guard"] = float(metrics["E_above_1e6_volume_fraction"]) == 0.0
            checks[name + "/temperature_range_guard"] = float(metrics["theta_max_K"]) < 30.0
            row = {
                **metrics,
                "radius_um": radius, "radius_over_width": radius / 20.0,
                "hole_area_fraction": np.pi * radius**2 / 400.0, "ligament_um": 10 - radius,
                "model": model, "case": entry["case"], "source": entry["source"],
                "nodes": len(xyz), "elements": len(connectivity), "sample_check": check,
            }
            rows.append(row)
            hashes[name] = {
                filename: digest(Path(entry["run"]) / filename)
                for filename in ("input.nml", "metrics.txt", "nodes.csv", "field_samples.csv")
            }
            if model == "NT" and radius > 0:
                midplane = fe.plane(xyz, connectivity, state, model, 0.0)
                top = fe.plane(xyz, connectivity, state, model, 0.5e-6)
                planes[(radius, "mid")] = midplane
                planes[(radius, "top")] = top
                np.savez_compressed(sources / f"R{radius:02d}_NT_midplane.npz", **midplane)
                np.savez_compressed(sources / f"R{radius:02d}_NT_top.npz", **top)

    csv_fields = [
        "radius_um", "radius_over_width", "hole_area_fraction", "ligament_um", "model", "case", "source",
        "nodes", "elements", "joule_W", "terminal_power_W", "theta_max_K", "ring_theta_mean_K",
        "ring_E_mean_V_m", "ring_Q_mean_W_m3", "ring_Q_p95_W_m3",
        "ring_abs_n_relative_mean", "maximum_abs_n_relative", "displacement_max_m", "minimum_mu_n_ratio",
        "ring_displacement_mean_m",
        "heat_balance_relative", "electron_balance_relative", "joule_quadrature_4_5_relative", "E_p95_V_m",
        "E_above_1e6_volume_fraction", "sample_check",
    ]
    with (sources / "radius_sweep_results.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    (sources / "radius_sweep_results.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.7), layout="constrained")
    specifications = (
        ("joule_W", 1e3, "Total Joule power", "mW"),
        ("theta_max_K", 1.0, "Peak temperature rise", "K"),
        ("E_p95_V_m", 1e-5, "Global electric-field P95", "kV cm$^{-1}$"),
        ("maximum_abs_n_relative", 100.0, "Maximum electron deviation", "%"),
        ("displacement_max_m", 1e9, "Maximum displacement", "nm"),
        ("minimum_mu_n_ratio", 1.0, "Minimum electron mobility ratio", r"$\mu_n/\mu_{n0}$"),
    )
    for ax, specification in zip(axes.flat, specifications):
        line_panel(ax, rows, *specification)
    axes[0, 0].legend(frameon=False, ncol=3)
    global_figures = save_figure(fig, figures, "F3_radius_global_metrics")

    hole_rows = [row for row in rows if row["radius_um"] > 0]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.7), layout="constrained")
    local_specs = (
        ("ring_Q_mean_W_m3", 1e-12, "Mean Joule heat", r"$10^{12}$ W m$^{-3}$"),
        ("ring_Q_p95_W_m3", 1e-12, "Joule heat P95", r"$10^{12}$ W m$^{-3}$"),
        ("ring_E_mean_V_m", 1e-5, "Mean electric field", "kV cm$^{-1}$"),
        ("ring_abs_n_relative_mean", 100.0, "Mean electron deviation", "%"),
        ("ring_theta_mean_K", 1.0, "Mean temperature rise", "K"),
        ("ring_displacement_mean_m", 1e9, "Mean displacement", "nm"),
    )
    for ax, specification in zip(axes.flat, local_specs):
        line_panel(ax, hole_rows, *specification)
        ax.set_xticks((2, 3, 4))
    axes[0, 0].legend(frameon=False, ncol=3)
    local_figures = save_figure(fig, figures, "F3_radius_local_metrics")

    q_bounds = (
        min(float(np.min(planes[(radius, "mid")]["Q"])) for radius in (2, 3, 4)) * 1e-12,
        max(float(np.max(planes[(radius, "mid")]["Q"])) for radius in (2, 3, 4)) * 1e-12,
    )
    n_max = max(float(np.max(np.abs(planes[(radius, "mid")]["n"] - 1.0))) for radius in (2, 3, 4)) * 100.0
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), layout="constrained")
    for column, radius in enumerate((2, 3, 4)):
        fe.panel(
            fig, axes[0, column], planes[(radius, "mid")], planes[(radius, "mid")]["Q"] * 1e-12,
            fr"Joule heat, r = {radius} $\mu$m", r"$10^{12}$ W m$^{-3}$", "magma", q_bounds,
        )
        fe.panel(
            fig, axes[1, column], planes[(radius, "mid")], np.abs(planes[(radius, "mid")]["n"] - 1.0) * 100.0,
            fr"Electron disturbance, r = {radius} $\mu$m", "%", "viridis", (0.0, n_max),
        )
    cloud_figures = save_figure(fig, figures, "F3_NT_radius_field_comparison")

    top_limit = max(float(np.max(np.abs(planes[(radius, "top")]["n"] - 1.0))) for radius in (2, 3, 4)) * 100.0
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), layout="constrained")
    for column, radius in enumerate((2, 3, 4)):
        fe.panel(
            fig, axes[column], planes[(radius, "top")], (planes[(radius, "top")]["n"] - 1.0) * 100.0,
            fr"Top-surface electron change, r = {radius} $\mu$m", "%", "RdBu_r", (-top_limit, top_limit),
        )
    surface_figures = save_figure(fig, figures, "F3_NT_top_surface_electron_change")

    audit = {
        "status": "PASS_L0_TREND_ONLY" if all(checks.values()) else "AUDIT_NEEDS_ATTENTION",
        "all_numerical_gates_passed": all(checks.values()), "checks": checks,
        "mesh_acceptance": "NOT_ESTABLISHED; L0 trend matrix only",
        "material_evidence": "ZnO numerical candidate, not same-sample experimental calibration",
        "boundary_interpretation": "Hole wall and non-contact exterior are ideal natural boundaries; exterior dielectric and surface physics omitted.",
        "field_figure_integrity": "Original HEX20 midplane field evaluation; no smoothing across element or hole boundaries; common scales across radii.",
        "source_hashes": hashes,
        "figures": global_figures + local_figures + cloud_figures + surface_figures,
    }
    (output / "radius_sweep_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )
    print(json.dumps({
        "status": audit["status"], "failed_checks": [key for key, value in checks.items() if not value],
        "rows": len(rows), "output": str(output),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
