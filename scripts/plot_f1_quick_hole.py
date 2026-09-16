"""Fast, source-backed H1-L0 hole result figures for an interim presentation."""
from pathlib import Path
import argparse
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_f1_main_results as mainfig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    mainfig.RUNS = args.run_root.resolve()
    mainfig.INPUTS = args.input_root.resolve()
    mainfig.OUT = args.output_dir.resolve()
    mainfig.OUT.mkdir(parents=True, exist_ok=True)
    case = "MESH_H1_NT_L0"
    metrics, xyz, conn, state = mainfig.load(case)
    mid = mainfig.plane(xyz, conn, state, "NT", 0.0)
    top = mainfig.plane(xyz, conn, state, "NT", 0.5e-6)
    specs = [
        (mid, mid["phi"], "Potential, midplane", "V", "viridis"),
        (mid, mid["E"]/1e5, "Electric field, midplane", "kV/cm", "viridis"),
        (mid, mid["Q"]/1e12, "Joule heat, midplane", r"$10^{12}$ W m$^{-3}$", "magma"),
        (mid, mid["theta"], "Temperature rise, midplane", "K", "magma"),
        (mid, mid["u"][:, 2]*1e9, "Axial displacement, midplane", "nm", "RdBu_r"),
        (top, 100*(top["n"]-1), "Electrons, top surface", "Change (%)", "RdBu_r"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), layout="constrained")
    for ax, spec, letter in zip(axes.flat, specs, "abcdef"):
        mainfig.panel(fig, ax, *spec)
        ax.text(-.15, 1.13, letter, transform=ax.transAxes, fontweight="bold")
    for ext in ["png", "svg", "pdf"]:
        fig.savefig(mainfig.OUT/f"H1_NT_L0_multifield.{ext}", dpi=320, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    saved_check = mainfig.check_saved_sample(case, metrics, xyz, conn, state)
    summary = {
        "status": "INTERIM_HOLE_RESULT_NOT_FINAL_MESH_ACCEPTANCE",
        "case": case, "model": metrics["model"], "nodes": len(xyz), "elements": len(conn),
        "metrics": metrics, "saved_sample_relative_error": saved_check,
        "figures": [str(mainfig.OUT/f"H1_NT_L0_multifield.{ext}") for ext in ["png", "svg", "pdf"]],
        "interpretation_boundary": "H1-L0 is a completed stage result; final H1-NT L2 mesh acceptance is still pending.",
        "source_metrics": str(mainfig.RUNS/case/"metrics.txt"),
        "source_nodes": str(mainfig.RUNS/case/"nodes.csv"),
        "source_field_samples": str(mainfig.RUNS/case/"field_samples.csv"),
    }
    (mainfig.OUT/"quick_hole_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({
        "status": summary["status"], "case": case, "nodes": len(xyz), "elements": len(conn),
        "joule_W": metrics["joule_W"], "theta_max_K": metrics["theta_max_K"],
        "ring_Q_mean_W_m3": metrics["ring_Q_mean_W_m3"],
        "saved_sample_relative_error": saved_check,
    }, indent=2))


if __name__ == "__main__":
    main()
