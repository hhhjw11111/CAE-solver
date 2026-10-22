"""Reproduce the existing SUPG analytical verification as a paper figure.

Read-only reuse of all 38 logged cases; no new FEM runs or invented nodal data.
Pe is based on the entire bar, not the largest element. Grids are graded.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from audit_je_paper_reuse_20261018 import validate_acceptance

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "runs/JE_RECONSTRUCTION_20261018/paper_scope"
LOG = PAPER.parent / "stage3_verification/attempt01/SUPG_targeted.txt"
SKILL = Path.home() / ".codex/skills/nature-figure"
sys.path.insert(0, str(SKILL / "scripts"))
from audit_panel_alignment import require_matplotlib_panel_alignment

PE_VALUES = (0.1, 1., 10., 50., 200., 1000.)
EXPECTED_LOG = "a99ff025c1e78f34b40a233d4f54296cd650ec668e17892302e588aff33932aa"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save(p, obj):
    with p.open("x", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write("\n")


def collect():
    validate_acceptance(PAPER / "verification_reuse/acceptance.json")
    if sha(LOG) != EXPECTED_LOG:
        raise ValueError("Frozen SUPG log identity changed")
    rows, pending = [], None
    for line_no, line in enumerate(LOG.read_text(encoding="utf-8").splitlines(), 1):
        if line.startswith("supg_bar_max_element_Pe="):
            pairs = dict(re.findall(r"(\w+)=\s*([^,]+)", line))
            pending = {k: float(v) for k, v in pairs.items()}
            pending["diagnostic_source_line"] = line_no
        if line.startswith("supg_bar_Pe="):
            if pending is None:
                raise ValueError("Missing paired diagnostic")
            fields = dict(re.findall(r"(\w+)=\s*([^,]+)", line))
            row = dict(bar_Pe=float(fields["supg_bar_Pe"]), carrier=int(fields["carrier"]),
                       elements=int(fields["elements"]),
                       density_error=float(fields["density_error"]),
                       current_error=float(fields["current_error"]),
                       max_element_Pe=pending["supg_bar_max_element_Pe"],
                       overshoot=pending["overshoot"],
                       maximum_backward_step=pending["maximum_backward_step"],
                       source_line=line_no, diagnostic_source_line=pending["diagnostic_source_line"],
                       source_file=str(LOG), source_sha256=EXPECTED_LOG, raw_line=line.strip())
            if any(not math.isfinite(row[k]) or row[k] <= 0
                   for k in ("bar_Pe", "density_error", "current_error", "max_element_Pe")):
                raise ValueError("Nonpositive or nonfinite log-axis data")
            rows.append(row)
            pending = None
    if len(rows) != 38 or len({(r["bar_Pe"], r["carrier"], r["elements"]) for r in rows}) != 38:
        raise ValueError("Missing or duplicate benchmark cases")
    finest = []
    for pe in PE_VALUES:
        for carrier in (1, 2):
            group = [r for r in rows if r["bar_Pe"] == pe and r["carrier"] == carrier]
            expected = [16, 32, 64, 128] if pe == 1000 else [16, 32, 64]
            if [r["elements"] for r in group] != expected:
                raise ValueError("Unexpected mesh sequence")
            for a, b in zip(group, group[1:]):
                if b["density_error"] > a["density_error"] * 1.05:
                    raise ValueError("Original refinement condition failed")
            last = group[-1]
            if last["density_error"] > 1e-3 or last["current_error"] > 1e-2:
                raise ValueError("Original finest-grid accuracy failed")
            for r in group:
                if r["elements"] >= 64 and max(r["overshoot"], r["maximum_backward_step"]) > 1e-3:
                    raise ValueError("Original resolved-grid monotonicity failed")
            finest.append(last)
    return rows, finest


def plot(rows, output):
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7,
        "legend.frameon": False, "lines.linewidth": .9,
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })
    fig, axes = plt.subplots(2, 2, figsize=(183 / 25.4, 165 / 25.4))
    fig.subplots_adjust(left=.12, right=.975, bottom=.13, top=.81, wspace=.37, hspace=.48)
    fig.suptitle("S3 | SUPG analytical transport verification", y=.98, fontsize=10, fontweight="bold")
    colors = ["#777777", "#009E73", "#56B4E9", "#CC79A7", "#D55E00", "#0072B2"]
    markers = ["o", "s", "^", "D", "v", "P"]
    for row, carrier in enumerate((1, 2)):
        for col, metric in enumerate(("density_error", "current_error")):
            ax = axes[row, col]
            for pe, color, marker in zip(PE_VALUES, colors, markers):
                g = [r for r in rows if r["bar_Pe"] == pe and r["carrier"] == carrier]
                ax.plot([r["elements"] for r in g], [r[metric] for r in g],
                        color=color, marker=marker, ms=4, mfc="white", label=f"Pe = {pe:g}")
            ax.set_xscale("log", base=2)
            ax.set_yscale("log")
            ax.set_xlim(14, 146)
            ax.set_xticks([16, 32, 64, 128])
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.set_ylim((1e-7, .2) if col == 0 else (3e-5, .4))
            ax.set_ylabel("Normalized density error" if col == 0 else "Scaled current error")
            ax.set_xlabel("Longitudinal elements (graded mesh)")
            ax.set_title(("Electron" if carrier == 1 else "Hole") +
                         (" concentration" if col == 0 else " current"), pad=9)
            ax.tick_params(direction="out", length=3, width=.7)
            ax.axhline(1e-3 if col == 0 else 1e-2, color="#999999", ls="--", lw=.7)
            ax.text(-.20, 1.06, "abcd"[2*row+col], transform=ax.transAxes,
                    fontsize=10, fontweight="bold", va="bottom")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, loc="upper center", bbox_to_anchor=(.54, .935),
               columnspacing=2.5, handlelength=2.5)
    fig.text(.5, .048, "Dashed lines: original finest-grid criteria; coarse grids need not pass.",
             ha="center", fontsize=7.5)
    fig.text(.5, .022, "All 38 logged cases retained. Pe is bar-scale; this is not a 6/8 V device certificate.",
             ha="center", fontsize=7.5)
    fig.canvas.draw()
    stem = output / "S3_SUPG_analytical_verification"
    require_matplotlib_panel_alignment(fig, json_out=stem.with_suffix(".alignment.json"),
        overlay_svg=stem.with_suffix(".alignment.svg"), tolerance_pt=1.5,
        gutter_tolerance_pt=1.5, require_panel_labels=True, strict=True)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".png"), dpi=600)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    rows, finest = collect()
    a.output.mkdir(parents=True, exist_ok=False)
    save(a.output / "figure_contract.json", dict(
        question="Does the existing stabilized carrier discretization converge to its analytical bar solution?",
        claim="Both carrier signs pass the unchanged finest-grid criteria across six bar Pe values.",
        panels=["Electron density", "Electron current", "Hole density", "Hole current"],
        archetype="Paired concentration/current accuracy across transport regimes",
        observations=38, statistical_replicates=False, uncertainty_intervals="Not applicable to deterministic solves",
        data_policy="All logged refinement cases; no fitted order, synthetic FEM curves or omitted failures",
        limits={"density_finest": 1e-3, "current_finest": 1e-2},
        export={"width_mm":183, "height_mm":165, "formats":["pdf","svg","png"], "png_dpi":600},
        style_source="Style-only reuse of existing project supplemental plots; no scientific transformations reused",
        device_mesh_supg_acceptance=False))
    for name, data in [("all_38_cases.csv", rows), ("finest_grid_cases.csv", finest)]:
        with (a.output / name).open("x", encoding="utf-8-sig", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    plot(rows, a.output)
    validate_acceptance(PAPER / "verification_reuse/acceptance.json")
    save(a.output / "provenance.json", dict(
        sources={str(x): sha(x) for x in [LOG, Path(__file__),
            ROOT/"src/verification/supg_check_m.f90", ROOT/"src/fem/transport_supg_m.f90",
            PAPER/"verification_reuse/acceptance.json", SKILL/"scripts/audit_panel_alignment.py"]},
        new_native_runs=0, missing_cases=0, excluded_cases=0, logged_precision="ES12.4 for plotted errors",
        all_finest_pass=True, bar_Pe_not_element_Pe=True,
        largest_element_Pe_at_bar1000_finest=max(r["max_element_Pe"] for r in finest if r["bar_Pe"] == 1000),
        normalization="Density max(abs(c-c_exact))/2; current max(abs(j-j_exact))/max(1,abs(j_exact)); 11 samples per element.",
        limits_unchanged=True, production_acceptance=False))
    print(json.dumps(dict(cases=len(rows), finest_cases=len(finest), original_benchmark_pass=True,
                          production_acceptance=False)))


if __name__ == "__main__":
    main()
