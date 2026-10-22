"""Plot S1/S2 from frozen native logs, without solving or inventing FEM fields.

No statistical uncertainty: each mark is a deterministic accepted solve.
All 13/10 cases are retained. Tables retain exact zeros and original log strings.
Only the new supplemental output directory is writable by this script.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from audit_je_paper_reuse_20261018 import validate_acceptance

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "runs/JE_RECONSTRUCTION_20261018"
PAPER = STUDY / "paper_scope"
ATTEMPT = STUDY / "stage3_verification/attempt01"
SKILL = Path.home() / ".codex/skills/nature-figure"
sys.path.insert(0, str(SKILL / "scripts"))
from audit_panel_alignment import require_matplotlib_panel_alignment

S1_CASES = [f"p{p}_nx{n}" for p in (1, 2) for n in (1, 2, 4, 8)] + [
    "reverse_voltage", "zero_voltage", "half_voltage", "joule_off", "transverse_refinement"]
S2_CASES = [f"s{s}_nx{n}" for s in (-1, 0, 1) for n in (1, 2, 4)] + ["transverse_2x2x2"]
S1_LABELS = [f"P{p}/{n}" for p in (1, 2) for n in (1, 2, 4, 8)] + ["rev", "zero", "half", "off", "3D"]
S2_LABELS = [f"{s}/{n}" for s in ("−", "0", "+") for n in (1, 2, 4)] + ["+/3D"]
COLORS = {"fem": "#0072B2", "analytic": "#111111", "heat": "#D55E00", "power": "#009E73"}
EXPECTED_LOG_SHA = {
    "A1JJE_targeted.txt": "2779ed6760223f75799adc82f182c8934c0de9184bd473ff571aafabf716a945",
    "JEGRAD_targeted.txt": "a4b4cb08d0d7b983b42e98b20432516c95f9771114e1b68c17cbcdc1f6808c60",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def write_csv(path: Path, rows: list[dict]):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("x", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_log(name: str):
    path = ATTEMPT / name
    if sha(path) != EXPECTED_LOG_SHA[name]:
        raise ValueError(f"Frozen source log changed: {path}")
    records = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if "=" not in line:
            continue
        key, raw = (part.strip() for part in line.split("=", 1))
        if key in records:
            raise ValueError(f"Duplicate metric {key}")
        try:
            value = float(raw.replace("D", "E"))
        except ValueError:
            value = raw
        if isinstance(value, float) and not np.isfinite(value):
            raise ValueError(f"Non-finite source metric {key}")
        records[key] = dict(source_file=str(path.relative_to(ROOT)), source_sha256=sha(path),
                            source_line=line_number, key=key, raw=raw, value=value)
    if records["verification_result"]["value"] != "PASS":
        raise ValueError(f"Source verification is not PASS: {name}")
    return records


def extract(records, prefix: str, cases: list[str]):
    result = []
    for i, case in enumerate(cases, 1):
        pre = f"{prefix}_{case}_"
        row = {"index": i, "case": case}
        row.update({key[len(pre):]: rec["value"] for key, rec in records.items() if key.startswith(pre)})
        if row.get("passed") != "T":
            raise ValueError(f"Missing/unaccepted case: {case}")
        result.append(row)
    if int(records[f"{prefix}_solved_cases"]["value"]) != len(cases):
        raise ValueError("Unexpected solve count; do not silently discard cases")
    return result


def field_error(row):
    return max(v for k, v in row.items() if isinstance(v, float) and
               ("_L2_" in k or k in ("max_relative_n_error", "max_relative_p_error")))


def point_and_flux_error(row):
    return max([field_error(row), row["max_flux_error"], row.get("diffusion_work_error", 0.0)] +
               [v for k, v in row.items() if k.startswith("point_")])


def configure():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "legend.frameon": False,
        "lines.linewidth": 0.8, "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def grid(title, labels):
    fig, axes = plt.subplots(2, 2, figsize=(183 / 25.4, 154 / 25.4))
    fig.subplots_adjust(left=0.095, right=0.975, bottom=0.115, top=0.89, wspace=0.35, hspace=0.60)
    fig.suptitle(title, fontsize=10, y=0.98, fontweight="bold")
    for ax, letter in zip(axes.flat, "abcd"):
        ax.text(-0.19, 1.08, letter, transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom")
        ax.set_xticks(np.arange(len(labels)), labels=labels, rotation=90, ha="center")
        ax.set_xlim(-0.65, len(labels) - 0.35)
        ax.tick_params(direction="out", length=3, width=0.7)
        ax.set_axisbelow(True)
    return fig, axes.ravel()


def marks(ax, values, label, color, marker, offset=0.0, size=24):
    x = np.arange(len(values), dtype=float) + offset
    ax.scatter(x, values, s=size, color=color if marker in ("x", "+") else "none",
               edgecolors=color if marker not in ("x", "+") else None,
               marker=marker, linewidths=0.8, label=label, zorder=3)


def log_marks(ax, values, label, color, marker, offset, omissions, scope):
    values = np.asarray(values, dtype=float)
    if (values < 0).any():
        raise ValueError("Error metric cannot be negative")
    good = values > 0
    ax.scatter(np.arange(len(values))[good] + offset, values[good], s=19,
               marker=marker, color=color, label=label, zorder=3)
    omissions[scope] = [int(i + 1) for i in np.flatnonzero(~good)]


def error_axis(ax, ylabel="Recorded error (dimensionless)", lower=1e-18):
    ax.set_yscale("log")
    ax.set_ylim(lower, 1e-5)
    ticks = [1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6]
    ax.set_yticks(ticks, labels=[f"1e{int(np.log10(v))}" for v in ticks])
    ax.set_ylabel(ylabel)


def limit_line(ax, value):
    ax.axhline(value, color="#777777", ls="--", lw=.7)
    ax.text(.98, value * 2.0, f"Limit: 1e{int(np.log10(value))}",
            transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=7.5)


def export(fig, stem: Path):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=stem.with_suffix(".alignment.json"),
        overlay_svg=stem.with_suffix(".alignment.svg"), tolerance_pt=1.5,
        gutter_tolerance_pt=1.5, require_panel_labels=True, strict=True)
    # Fixed physical dimensions, no tight bounding-box cropping or raster editing.
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".png"), dpi=600)
    plt.close(fig)


def plot_s1(rows, out, omissions):
    fig, axs = grid("S1  |  Constant-density Joule benchmark (A1JJE)", S1_LABELS)
    a, b, c, d = axs
    marks(a, [r["numerical_peak_K"] for r in rows], "FEM peak", COLORS["fem"], "o", -.10)
    marks(a, [r["analytic_peak_K"] for r in rows], "Analytical peak", COLORS["analytic"], "x", .10)
    a.set_ylabel("Peak temperature rise (K)")
    a.set_ylim(-0.18, 5.35)
    a.legend(loc="upper left")
    marks(b, [r["joule_W"] * 1e3 for r in rows], "Heat-source integral", COLORS["fem"], "o", -.18, 18)
    marks(b, [(r["heat_left_W"] + r["heat_right_W"]) * 1e3 for r in rows], "Outward heat (left + right)", COLORS["heat"], "s", 0, 18)
    marks(b, [r["heat_factor"] * r["terminal_power_W"] * 1e3 for r in rows], "λ × terminal power", COLORS["power"], "x", .18, 18)
    b.set_ylabel("Thermal power comparison (mW)")
    b.set_ylim(-.09, 2.8)
    b.legend(loc="upper left")
    log_marks(c, [field_error(r) for r in rows], "Maximum field error", COLORS["fem"], "o", -.09, omissions, "S1c_field")
    log_marks(c, [r["max_flux_error"] for r in rows], "Independent flux error", COLORS["heat"], "s", .09, omissions, "S1c_flux")
    error_axis(c)
    limit_line(c, 1e-8)
    c.legend(loc="upper left", fontsize=7.5)
    log_marks(d, [r["quadrature_4_vs_5"] for r in rows], "Order 4 versus order 5", COLORS["power"], "D", 0, omissions, "S1d_quadrature")
    error_axis(d, "Scaled integration difference")
    limit_line(d, 1e-10)
    d.legend(loc="upper left")
    for ax in (c, d):
        ax.set_xlabel("Profile / longitudinal elements; controls")
    fig.text(.5, .018, "All 13 solves retained. Marker offsets aid reading; no spatial FEM curve is reconstructed.",
             ha="center", fontsize=7.5)
    export(fig, out / "S1_A1JJE_summary")


def plot_s2(rows, out, omissions):
    fig, axs = grid("S2  |  Signed diffusion-work benchmark (JEGRAD)", S2_LABELS)
    a, b, c, d = axs
    marks(a, [r["diffusion_work_W"] * 1e9 for r in rows], "FEM diffusion work", COLORS["fem"], "o", -.10)
    marks(a, [r["analytic_diffusion_work_W"] * 1e9 for r in rows], "Analytical diffusion work", COLORS["analytic"], "x", .10)
    a.axhline(0, color="#666666", lw=.65)
    a.set_ylabel("Signed diffusion power (nW)")
    a.set_ylim(-.18, .27)
    a.legend(loc="upper left")
    marks(b, [r["numerical_peak_K"] * 1e6 for r in rows], "FEM peak", COLORS["fem"], "o", -.10)
    marks(b, [r["analytic_peak_K"] * 1e6 for r in rows], "Analytical peak", COLORS["analytic"], "x", .10)
    b.set_ylabel("Peak temperature rise (µK)")
    b.set_ylim(52, 63)
    b.legend(loc="upper left")
    marks(c, [r["joule_W"] * 1e9 for r in rows], "JE volume integral", COLORS["fem"], "o", -.10)
    marks(c, [r["terminal_power_W"] * 1e9 for r in rows], "Terminal power", COLORS["power"], "x", .10)
    c.set_ylabel("Electrical power (nW)")
    c.set_ylim(4.2, 4.95)
    c.legend(loc="upper left")
    log_marks(d, [point_and_flux_error(r) for r in rows], "Max. field / current / flux error", COLORS["fem"], "o", -.09, omissions, "S2d_field_current_flux")
    log_marks(d, [r["quadrature_4_vs_5"] for r in rows], "Order 4 versus order 5", COLORS["heat"], "s", .09, omissions, "S2d_quadrature")
    error_axis(d)
    limit_line(d, 1e-8)
    limit_line(d, 1e-10)
    d.legend(loc="upper left", fontsize=7.5)
    for ax in (c, d):
        ax.set_xlabel("Slope sign / longitudinal elements; 3D check")
    fig.text(.5, .018, "s = −0.05, 0, +0.05; V = 0.02 V. All 10 solves retained; synthetic frozen-transport LQ model.",
             ha="center", fontsize=7.5)
    export(fig, out / "S2_JEGRAD_summary")


def run_qa(out):
    jobs = [("source_preflight", "validate_figure.py", [str(Path(__file__).resolve()), "--json"])]
    for stem in ("S1_A1JJE_summary", "S2_JEGRAD_summary"):
        jobs.extend([
            (stem + ".text_audit", "audit_pdf_text.py", [str(out / (stem + ".pdf")), "--min-pt", "5", "--json"]),
            (stem + ".collision_audit", "audit_figure_collisions.py", [str(out / (stem + ".pdf")),
             "--json-out", str(out / (stem + ".collision.json")), "--overlay-pdf", str(out / (stem + ".collision.pdf"))]),
        ])
    results = []
    for name, script, args in jobs:
        proc = subprocess.run([sys.executable, str(SKILL / "scripts" / script), *args],
                              capture_output=True, text=True, encoding="utf-8", errors="replace")
        (out / (name + ".txt")).write_text(proc.stdout + proc.stderr, encoding="utf-8")
        results.append({"check": name, "exit_code": proc.returncode})
    write_json(out / "automated_qa.json", results)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=PAPER / "supplemental_figures")
    args = parser.parse_args()
    out = args.output.resolve()
    if not out.is_relative_to(PAPER.resolve()):
        raise ValueError("Outputs must remain inside the new paper scope")
    out.mkdir(parents=True, exist_ok=True)
    existing = [p for p in out.iterdir() if p.name != "figure_contract.json"]
    if existing:
        raise FileExistsError("Refuse to overwrite outputs; choose a fresh --output directory")
    contract_path = PAPER / "supplemental_figures/figure_contract.json"
    if not (out / "figure_contract.json").exists():
        (out / "figure_contract.json").write_bytes(contract_path.read_bytes())
    accepted = validate_acceptance(PAPER / "verification_reuse/acceptance.json")
    logs = [read_log(name) for name in EXPECTED_LOG_SHA]
    rows1 = extract(logs[0], "a1j", S1_CASES)
    rows2 = extract(logs[1], "jegrad", S2_CASES)
    for row in rows1:
        row["heat_factor"] = 0.0 if row["case"] == "joule_off" else 1.0
        row["outward_heat_left_plus_right_W"] = row["heat_left_W"] + row["heat_right_W"]
        row["thermal_target_from_terminal_W"] = row["heat_factor"] * row["terminal_power_W"]
        row["max_field_error"] = field_error(row)
        row["absolute_peak_difference_K"] = abs(row["numerical_peak_K"] - row["analytic_peak_K"])
    for row in rows2:
        row["max_field_current_flux_error"] = point_and_flux_error(row)
        row["absolute_peak_difference_K"] = abs(row["numerical_peak_K"] - row["analytic_peak_K"])
        row["heat_left_W"] = "NOT_LOGGED"
        row["heat_right_W"] = "NOT_LOGGED"
    write_csv(out / "S1_all_cases.csv", rows1)
    write_csv(out / "S2_all_cases.csv", rows2)
    raw = [record for log in logs for record in log.values()]
    write_csv(out / "source_metrics_with_lines.csv", raw)
    controls = [r for r in raw if any(s in r["key"] for s in ("negative_", "dual_start", "continuation_difference", "flux_patch", "_NT_"))]
    write_csv(out / "supplemental_controls_and_local_NT.csv", controls)
    omissions = {}
    configure()
    plot_s1(rows1, out, omissions)
    plot_s2(rows2, out, omissions)
    stats = {"S1": {"cases": len(rows1), "maximum_field_error": max(map(field_error, rows1)),
             "maximum_flux_error": max(r["max_flux_error"] for r in rows1),
             "maximum_heat_balance_defect": max(r["heat_balance_defect"] for r in rows1),
             "maximum_quadrature_error": max(r["quadrature_4_vs_5"] for r in rows1)},
             "S2": {"cases": len(rows2), "maximum_field_current_flux_error": max(map(point_and_flux_error, rows2)),
             "maximum_quadrature_error": max(r["quadrature_4_vs_5"] for r in rows2)},
             "logarithmic_plot_exact_zero_omissions_by_case_index": omissions,
             "excluded_solved_cases": 0, "new_native_runs": 0}
    write_json(out / "summary.json", stats)
    qa = run_qa(out)
    sources = [Path(__file__).resolve(), Path(__file__).with_name("audit_je_paper_reuse_20261018.py"),
        PAPER / "verification_reuse/acceptance.json", PAPER / "contract.json", STUDY / "contract.json",
        ROOT / "docs/je_paper_minimal_workflow_20261018.md", contract_path,
        STUDY / "stage3_verification/jegrad_protocol.json",
        ROOT / "docs/a1j_production_joule_flux_audit_spec_v1.md",
        ROOT / "src/verification/production_joule_check_m.f90", ROOT / "src/verification/je_gradient_check_m.f90",
        *[ATTEMPT / name for name in EXPECTED_LOG_SHA],
        *[SKILL / "scripts" / name for name in ("audit_panel_alignment.py", "audit_pdf_text.py", "audit_figure_collisions.py", "validate_figure.py")]]
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(), scope="S1_S2_LOG_REUSE_NO_NEW_NATIVE",
        python=platform.python_version(), matplotlib=matplotlib.__version__, numpy=np.__version__,
        executed_candidate_sha256=accepted["candidate_sha256"],
        source_hashes={str(p.resolve()): sha(p) for p in sources},
        source_attempt=str(ATTEMPT), acceptance_revalidated=True,
        validation_scope="Original frozen candidate execution, not later changed A2P/driver source", qa=qa)
    # Ensure plotting and read-only QA did not change frozen evidence.
    validate_acceptance(PAPER / "verification_reuse/acceptance.json")
    write_json(out / "provenance.json", provenance)
    write_json(out / "generated_file_hashes.json", {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()})
    print(json.dumps(stats, ensure_ascii=False))
    return 1 if any(item["exit_code"] != 0 for item in qa) else 0


if __name__ == "__main__":
    raise SystemExit(main())
