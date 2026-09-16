"""Fail-closed provenance and numeric audit of the completed frozen F1 matrix."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import re
import numpy as np

import audit_f1_results as audit
from summarize_f1_release import MAIN


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--figure-root", type=Path, required=True)
    args = parser.parse_args()
    root, figures = args.run_root.resolve(), args.figure_root.resolve()
    audit.ROOT = root
    checks, details = {}, {}

    def check(name, condition, detail=None):
        checks[name] = bool(condition)
        if detail is not None:
            details[name] = detail

    release = read(root/"evidence/release.json")
    check("release_binary", sha(root/"evidence/TPSolverV2_main.exe") == release["executable_sha256"].lower())
    build = read(root/"evidence/build_manifest.json")
    source_checks = {}
    for entry in build["files"]:
        if entry["path"].startswith("src/"):
            source_checks[entry["path"]] = sha(root/"evidence"/entry["path"]) == entry["sha256"].lower()
    check("source_snapshot_matches_build", bool(source_checks) and all(source_checks.values()), source_checks)
    input_entries = read(root/"evidence/input_hashes.json")
    check("frozen_inputs", all(
        sha(root/"evidence/inputs"/r["path"]) == r["sha256"].lower() for r in input_entries))
    inherited = read(root/"evidence/inherited_evidence.json")
    check("inherited_evidence_unchanged", all(
        sha(root/c["case"]/r["name"]) == r["sha256"].lower()
        for c in inherited["cases"] for r in c["files"]))
    for gate in ["L01", "G07", "G08"]:
        check("targeted_"+gate, "verification_result=PASS" in (root/f"evidence/{gate}_stdout.txt").read_text(encoding="utf-8-sig"))
    for analytic in ["A1", "A2"]:
        check(analytic, "verification_result=PASS" in (root/analytic/f"{analytic}_export_stdout.txt").read_text(encoding="utf-8-sig"))
    small = audit.small("SMALL_H1_LQ_L0", "SMALL_H1_NT_L0_0142")
    check("small_recomputed", small["passed"])
    mesh = audit.mesh([f"MESH_H1_NT_L{i}" for i in range(3)])
    check("mesh_recomputed", mesh["passed"], mesh)
    cases = list(dict.fromkeys([f"MESH_H1_NT_L{i}" for i in range(3)] + MAIN + ["REPEAT_H1_NT_L2"]))
    geometry_coordinates = {}
    for case in cases:
        folder = root/case
        metrics = audit.read_metrics(case)
        supplemental = read(folder/"supplemental_metrics.json") if (folder/"supplemental_metrics.json").exists() else {}
        check(case+"/accepted", metrics.get("case_acceptance") == "PASS_NOT_MESH_ACCEPTANCE")
        for key, limit in [
            ("maximum_stage_balance", 1e-8), ("maximum_accepted_residual", 1e-9),
            ("maximum_accepted_update", 1e-9), ("heat_balance_relative", .01),
            ("electron_balance_relative", .01), ("joule_quadrature_4_5_relative", .001),
        ]:
            check(case+"/"+key, np.isfinite(metrics[key]) and 0 <= metrics[key] < limit)
        check(case+"/jacobian", metrics["minimum_jacobian"] > 0)
        check(case+"/voltage", metrics["voltage_V"] == 3)
        if case != "MESH_H1_NT_L0":
            execution = read(folder/"execution.json")
            check(case+"/execution", execution["exit_code"] == 0 and
                  execution["status"] == "PASS_NOT_MESH_ACCEPTANCE")
            check(case+"/release_identity",
                  execution["executable_sha256"].lower() == release["executable_sha256"].lower()
                  and execution["solver_version"] == release["solver_version"]
                  and execution["storage"] == "OOC"
                  and execution["ooc_max_core_mb"] == release["ooc_max_core_mb"]
                  and execution["ooc_max_swap_mb"] == 0)
            check(case+"/input_identity", sha(folder/"input.nml") == execution["input_sha256"].lower())
            records = read(folder/"output_hashes.json")
            check(case+"/output_hashes", all(sha(folder/r["name"]) == r["sha256"].lower() for r in records))
        nodes = np.loadtxt(folder/"nodes.csv", delimiter=",", skiprows=1, usecols=range(13))
        check(case+"/nodes_finite_positive", np.isfinite(nodes).all()
              and np.all(nodes[:, 10] > -300) and np.all(nodes[:, 11:13] > 0))
        x, y, z = nodes[:, 1:4].T
        left, right = np.abs(x+10e-6) < 1e-15, np.abs(x-10e-6) < 1e-15
        check(case+"/contacts", np.any(left) and np.any(right) and
              np.all(np.abs(nodes[left, 7]) < 1e-12) and
              np.all(np.abs(nodes[right, 7]+3) < 1e-12) and
              np.all(np.abs(nodes[left|right, 8:11]) < 1e-12))
        for anchor, fields in [
            ((-10e-6, -10e-6, -.5e-6), [4, 5, 6]),
            ((-10e-6, 10e-6, -.5e-6), [4, 6]),
            ((-10e-6, -10e-6, .5e-6), [4]),
        ]:
            selected = np.all(np.abs(nodes[:, 1:4]-anchor) < 1e-15, axis=1)
            check(case+"/anchor/"+str(anchor), selected.sum() == 1
                  and np.all(np.abs(nodes[np.ix_(selected, fields)]) < 1e-20))
        if case in MAIN:
            geometry = "H1" if "H1" in case else "H0"
            coords = nodes[:, :4]
            if geometry in geometry_coordinates:
                check(case+"/same_geometry_mesh", np.array_equal(geometry_coordinates[geometry], coords))
            else:
                geometry_coordinates[geometry] = coords.copy()
            check(case+"/table_supplement", bool(supplemental) and
                  supplemental["source_nodes_sha256"] == sha(folder/"nodes.csv") and
                  supplemental["source_metrics_sha256"] == sha(folder/"metrics.txt"))
            # Independent saved volume samples, loaded only after the solve matrix ends.
            sample = np.loadtxt(folder/"field_samples.csv", delimiter=",", skiprows=1)
            check(case+"/samples_finite_positive", np.isfinite(sample).all()
                  and np.all(sample[:, 3] > 0) and np.all(sample[:, 8] > -300)
                  and np.all(sample[:, 9:11] > 0))
            sample_checks = {
                "volume_relative": audit.rel(float(sample[:, 3].sum()), metrics["total_volume_m3"]),
                "joule_relative": audit.rel(float(np.dot(sample[:, 3], sample[:, 7])), metrics["joule_W"]),
                "theta_relative": audit.rel(float(sample[:, 8].max()), metrics["theta_max_K"]),
            }
            check(case+"/saved_volume_integrals", max(sample_checks.values()) < 2e-8, sample_checks)
            del sample
        del nodes
    repeated = {name: sha(root/"MESH_H1_NT_L2"/name) == sha(root/"REPEAT_H1_NT_L2"/name)
                for name in ["metrics.txt", "nodes.csv", "field_samples.csv"]}
    check("repeat_byte_identity", all(repeated.values()), repeated)
    expected_figures = ["A1_joint_analytic", "A2_hole_analytic", "H1_NT_multifield",
                        "H1_model_carrier_comparison", "six_case_comparison"]
    check("figures_present", all((figures/(name+"."+ext)).is_file()
          for name in expected_figures for ext in ["png", "svg", "pdf"]))
    figure_qa = read(figures/"main_figure_qa.json")
    check("figure_FE_source_checks", all(v < 2e-6 for v in figure_qa.values() if isinstance(v, (float, int))))
    check("main_table", len(read(root/"main_results.json")) == 6)
    check("contrast_interpretation", bool(read(root/"contrasts.json")) and bool(read(root/"interpretation.json")["paragraphs"]))
    result = {
        "schema": "F1-release-audit-v1", "passed": all(checks.values()), "checks": checks,
        "details": details, "failed_checks": [k for k, v in checks.items() if not v],
        "date_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "material_passport": {"mode": "run/validate", "verification": "VERIFIED" if all(checks.values()) else "INCOMPLETE",
                             "scope": "Deterministic numerical implementation and conditional candidate model; no experimental calibration."},
        "visual_review": "Separate mandatory manual review of figures and rendered final PDF.",
    }
    (root/"release_audit.json").write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({"passed": result["passed"], "checks": len(checks), "failed_checks": result["failed_checks"]}, indent=2))
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
