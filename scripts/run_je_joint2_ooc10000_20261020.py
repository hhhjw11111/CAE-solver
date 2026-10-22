"""Resource-only retry of the documented joint2 PARDISO -9, then SUPG.

Fresh output, original executable/physics/numerical gates/7200 s limit.
Only OOC core budget changes 8000 -> 10000 MB. Stop at the first failed solve.
"""
from pathlib import Path
import json
import subprocess
import sys
import run_je_reconstruction_mesh_20261018 as r
import run_je_paper_resumed_20261020 as prior
import audit_je_paper_observables_20261020 as observables

OUTPUT = r.STUDY/"stage4_paper_mesh_ooc10000_20261020"
FAILED = prior.OUTPUT/"A48_R22_Z18_S1_V060"


def validate():
    evidence = prior.evidence()
    ex = r.read(FAILED/"execution.json")
    if ex["exit_code"] == 0 or ex["timed_out"]:
        raise RuntimeError("Expected resource failure is absent")
    if "PARDISO phase 13 failed with error -9" not in (FAILED/"stdout.txt").read_text():
        raise RuntimeError("Missing original PARDISO -9")
    if ex["input_sha256"] != r.sha(FAILED/"input.nml") or ex["executable_sha256"] != r.EXPECTED_EXE:
        raise RuntimeError("Failed-attempt identity mismatch")
    for v in (6,8):
        prior.resolve_complete(r.name((40,20,16), v))
    return evidence


def complete(name):
    path = OUTPUT/name
    if (path/"completion.json").exists():
        r.accepted(path, require_new_completion=True)
        return path
    return prior.resolve_complete(name)


def report(folder):
    # Explicit study-root binding; original audit functions and limits unchanged.
    original = observables.locate
    observables.locate = lambda name: OUTPUT/name if (OUTPUT/name).exists() else original(name)
    try:
        result = observables.snapshot(2)
    finally:
        observables.locate = original
    result["extra_study_root"] = str(OUTPUT)
    result["resource_retry_runner_sha256"] = r.sha(Path(__file__))
    folder.mkdir(exist_ok=False)
    r.save(folder/"observable_audit.json", result)


def run_point(job, gate):
    prior.policy.preflight()
    source = Path(job["source"])
    dest = Path(job["directory"])
    if dest.exists():
        raise RuntimeError("Existing attempt retained; no overwrite or implicit retry")
    mesh = r.ensure_mesh(prior.OUTPUT, job["counts"])
    old_mapping = prior.OUTPUT/"transfers"/job["name"]
    if old_mapping.exists():
        mapped = old_mapping
        m = r.read(mapped/"transfer_report.json")
        sm = r.accepted(source)
        if (m["source_nodes_sha256"] != r.sha(source/"nodes.csv")
                or m["source_mesh_sha256"] != sm["mesh_sha256"]
                or m["target_mesh_sha256"] != r.sha(mesh) or m["boundary_extension_nodes"] != 0):
            raise RuntimeError("Existing initial-guess transfer identity mismatch")
    else:
        mapped = r.prepare_mapping(OUTPUT, job, mesh)
    if dest.name == FAILED.name:
        template = (r.BASE/"input_V060.nml").read_text(encoding="utf-8-sig")
        changes = dict(case_id="'JE_RECON_"+job["name"]+"'", mesh_file="'"+r.short(mesh)+"'",
                       restart_nodes_file="'"+r.short(mapped/"nodes.csv")+"'",
                       restart_voltage_factor="1.0", supg_scale=str(job["supg_scale"]))
        for key, value in changes.items():
            template = r.replace(template, key, value)
        if template != (FAILED/"input.nml").read_text(encoding="ascii"):
            raise RuntimeError("Retry native input would change; no solve started")
    r.run_native(job, mesh, mapped/"nodes.csv", dest, gate)
    if dest.name == FAILED.name and r.sha(dest/"input.nml") != r.sha(FAILED/"input.nml"):
        raise RuntimeError("Resource-only 6 V retry changed native input")
    subprocess.run([sys.executable,"-X","utf8",str(r.ROOT/"scripts/audit_je_sensitivity_power_20261015.py"),
                    "--runs",str(dest)], check=True)
    fields = r.field_comparison(mapped, source, dest)
    r.save(dest/"mapped_field_comparison.json", fields)
    r.save(dest/"completion.json", dict(power_audit_completed=True,
        field_comparison_passed=fields["passed"], fully_validated=False,
        source_run=str(source), mapping_directory=str(mapped), resource_retry_only=True,
        geometry_and_transfer_only_initial_guess=True))
    print("COMPLETED "+dest.name+"; local fields passed="+str(fields["passed"]), flush=True)


def main():
    prior.policy.preflight()
    evidence = validate()
    OUTPUT.mkdir(exist_ok=False)
    scripts = [Path(__file__).resolve(), Path(r.__file__).resolve(), Path(prior.__file__).resolve(),
               Path(prior.policy.__file__).resolve(), Path(observables.__file__).resolve()]
    hashes = {str(x):r.sha(x) for x in scripts}
    r.save(OUTPUT/"runner_manifest.json", dict(scope="JOINT2_OOC10000_RESOURCE_RETRY",
        original_failure=str(FAILED), failure_hashes={str(FAILED/f):r.sha(FAILED/f)
            for f in ("execution.json","resource_usage.json","input.nml","stdout.txt","manifest.json")},
        evidence=evidence, script_hashes=hashes, changed_runtime_option="OOC 8000 -> 10000 MB",
        threads=r.THREADS, timeout_seconds=r.TIMEOUT, start_disk=r.START_DISK,
        stop_disk=r.STOP_DISK, stop_available_ram=r.MIN_AVAILABLE_RAM,
        thresholds_changed=False, physics_changed=False, production_acceptance=False))
    r.OOC_MB = 10000
    events = []
    try:
        for batch in ("joint","supg"):
            jobs = r.jobs(batch, OUTPUT, 2)
            for job in jobs:
                counts = (40,20,16) if batch == "joint" else (48,22,18)
                job["source"] = str(complete(r.name(counts,job["voltage"])))
            r.save(OUTPUT/(batch+"_plan.json"), dict(jobs=jobs, production_acceptance=False))
            for job in jobs:
                if validate() != evidence:
                    raise RuntimeError("Evidence identity changed")
                for path, h in hashes.items():
                    if r.sha(Path(path)) != h:
                        raise RuntimeError("Batch dependency changed: "+path)
                run_point(job, evidence["gate"])
                events.append(dict(name=job["name"], state=r.PASS))
            report(OUTPUT/("after_"+batch))
        r.save(OUTPUT/"finished.json", dict(events=events, attempts_completed=True, production_acceptance=False))
    except BaseException as exc:
        r.save(OUTPUT/"stopped.json", dict(reason=repr(exc), events=events, production_acceptance=False))
        raise


if __name__ == "__main__":
    main()
