"""Report original joint-mesh/SUPG gates by paper observable, without relaxing any.

Missing data stay pending; a group's pass never overrides the original full
gate or a failed local hole comparison. No production acceptance is issued.
"""
from pathlib import Path
import argparse
import csv
import json
import run_je_paper_resumed_20261020 as resume
import run_je_reconstruction_mesh_20261018 as r

GROUPS = {
    "temperature": ["theta_max_K", "ring_theta_mean_K"],
    "power_current": ["joule_W", "electron_right_A"],
    "hole_edge_fields": ["ring_E_mean_V_m", "ring_Q_mean_W_m3", "ring_Q_p95_W_m3", "E_p95_V_m"],
    "electron_ring": ["ring_abs_n_relative_mean"],
    "hole_ring": ["ring_abs_p_relative_mean"],
}


def locate(name):
    for root in (resume.OUTPUT, resume.PRIOR, resume.policy.ORIGINAL):
        if (root/name/"completion.json").exists():
            return root/name
    for root in (resume.OUTPUT, resume.PRIOR, resume.policy.ORIGINAL):
        if (root/name).exists():
            return root/name
    return resume.OUTPUT/name


def snapshot(level):
    reports, rows, states, hashes = {}, [], [], {}

    def preserve(path):
        if path.is_file():
            hashes[str(path)] = r.sha(path)

    def check_point(path):
        item = dict(directory=str(path), state="MISSING_EVIDENCE")
        if not (path/"completion.json").exists():
            if (path/"execution.json").exists():
                execution = r.read(path/"execution.json")
                state, _ = r.classify(path, execution["exit_code"])
                item.update(state="POSTPROCESS_PENDING" if state == r.PASS else state,
                            exit_code=execution["exit_code"])
                preserve(path/"execution.json")
            elif path.exists():
                item["state"] = "INCOMPLETE_NO_EXIT_RECORD"
            states.append(item)
            return
        manifest = r.accepted(path, require_new_completion=True)
        _, metrics = r.classify(path, r.read(path/"execution.json")["exit_code"])
        comparison = r.read(path/"mapped_field_comparison.json")
        for source, digest in comparison["hashes"].items():
            if r.sha(Path(source)) != digest:
                raise ValueError("Changed local evidence: "+source)
            hashes[source] = digest
        power = r.read(path/"power_audit/decomposition.json")
        provenance = r.read(path/"power_audit/provenance.json")
        if provenance["nodes_sha256"] != r.sha(path/"nodes.csv") or provenance["mesh_sha256"] != manifest["mesh_sha256"]:
            raise ValueError("Power source identity mismatch")
        for source, digest in provenance["scripts"].items():
            if r.sha(Path(source)) != digest:
                raise ValueError("Changed power audit: "+source)
            hashes[source] = digest
        item.update(state=r.PASS, local_fields=comparison["comparisons"],
            local_all_fields_passed=comparison["passed"],
            voltage_V=metrics["voltage_V"], supg_scale=metrics["supg_scale"],
            metrics={k: metrics[k] for group in GROUPS.values() for k in group},
            power_gap_all_exterior_percent=100*power["gap_relative"],
            power_gap_electrode_percent=100*power["electrode_gap_relative"],
            power_gap_is_diagnostic_not_new_gate=True)
        for name in ("completion.json", "execution.json", "manifest.json", "input.nml", "nodes.csv",
                     "metrics.txt", "mapped_field_comparison.json", "power_audit/decomposition.json",
                     "power_audit/provenance.json"):
            preserve(path/name)
        states.append(item)

    for voltage in (6, 8):
        mesh_runs = [r.radial18(voltage), locate(r.name((40,20,16), voltage)),
                     locate(r.name((48,22,18), voltage))]
        counts = (40,20,16) if level == 1 else (48,22,18)
        supg_runs = [locate(r.name(counts, voltage, scale)) for scale in (.5,1.,2.)]
        mesh = r.mesh_comparison(mesh_runs, [voltage])
        supg = r.supg_comparison(supg_runs)
        reports[f"{voltage}V"] = dict(mesh=mesh, supg=supg,
            selected_SUPG_counts=list(counts), full_acceptance=False)
        for group, metrics in GROUPS.items():
            for metric in metrics:
                mc = mesh["voltages"][0]["changes"].get(metric)
                sc = supg["comparisons"].get(metric)
                mp = None if mc is None else mc["L1_L2"] is not None and mc["L1_L2"] < mc["limit"]
                sp = None if sc is None else sc["passed"]
                rows.append(dict(voltage_V=voltage, group=group, metric=metric,
                    mesh_last_pair_change=None if mc is None else mc["L1_L2"],
                    original_mesh_limit=None if mc is None else mc["limit"],
                    mesh_metric_status="PENDING_OR_INVALID_EVIDENCE" if mp is None else ("PASS" if mp else "FAIL"),
                    supg_relative_span=None if sc is None else sc["relative_span"],
                    original_supg_limit=None if sc is None else sc["limit"],
                    supg_metric_status="PENDING_OR_INVALID_EVIDENCE" if sp is None else ("PASS" if sp else "FAIL"),
                    joint_and_supg_metric_pass=mp is True and sp is True,
                    production_acceptance=False))
        for path in dict.fromkeys(mesh_runs[1:]+supg_runs):
            check_point(path)
        for path in mesh_runs+supg_runs:
            for name in ("manifest.json", "execution.json", "input.nml", "metrics.txt", "field_audit.json"):
                preserve(path/name)
        base_manifest = r.accepted(mesh_runs[0])
        preserve(Path(base_manifest["mesh"]))
        preserve(Path(base_manifest["executable"]))
    return dict(original_audits=reports, metric_rows=rows, point_states=states,
        source_hashes=hashes, thresholds_changed=False, production_acceptance=False,
        interpretation="Groups organize existing gates only. Local fields, independent power diagnostics and full failed gates remain separate. "
            "A group's PASS is not a certificate for a complete field, other radii, or the full voltage interval.",
        script_sha256=r.sha(Path(__file__)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--joint-level", type=int, choices=[1,2], default=2,
                   help="SUPG mesh to assess; joint-mesh audit always requires all three levels.")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    resume.evidence()
    report = snapshot(a.joint_level)
    a.output.mkdir(parents=True, exist_ok=False)
    r.save(a.output/"observable_audit.json", report)
    with (a.output/"original_metric_gates.csv").open("x", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(report["metric_rows"][0]))
        writer.writeheader()
        writer.writerows(report["metric_rows"])
    lines = ["论文关注量分项审计：原门限不变", "",
        "表内PASS只表示该项原比较通过；不覆盖局部空穴FAIL，不签发生产证书。",
        "PENDING_OR_INVALID_EVIDENCE表示缺失或无效的比较证据，原审计详细错误保留于JSON。",
        "功率差额独立报告，不增加新门限，不使用SUPG补偿凑平。", ""]
    for voltage in (6,8):
        lines.append(f"{voltage} V：")
        for group in GROUPS:
            rr = [v for v in report["metric_rows"] if v["voltage_V"]==voltage and v["group"]==group]
            lines.append(f"  {group}: "+", ".join(
                x["metric"]+" [mesh="+x["mesh_metric_status"]+", SUPG="+x["supg_metric_status"]+"]" for x in rr))
    with (a.output/"REPORT.txt").open("x", encoding="utf-8") as f:
        f.write("\n".join(lines)+"\n")
    print(json.dumps(dict(metric_rows=len(report["metric_rows"]),
        jointly_passed=sum(x["joint_and_supg_metric_pass"] for x in report["metric_rows"]),
        production_acceptance=False)))


if __name__ == "__main__":
    main()
