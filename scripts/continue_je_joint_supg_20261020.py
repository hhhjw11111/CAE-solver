"""Finite sequential continuation of the already authorized joint/SUPG study.

Wait for the existing joint1 controller, then run joint2, its report, SUPG2,
and a final observable report. No retries, gate edits, production solves,
schedulers or physical-code changes. Any subprocess failure stops the batch.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ctypes
import json
import subprocess
import sys
import time
import run_je_paper_resumed_20261020 as resume
import run_je_reconstruction_mesh_20261018 as r


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior-controller-pid", type=int, required=True)
    p.add_argument("--prior-finished", type=Path, required=True)
    p.add_argument("--prior-launch", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    out = a.output.resolve()
    if not out.is_relative_to((r.STUDY/"paper_scope").resolve()):
        raise ValueError("Controller records must stay in the paper scope")
    launch = r.read(a.prior_launch)
    if launch["controller_pid"] != a.prior_controller_pid:
        raise ValueError("Prior controller identity mismatch")
    resume.evidence()
    out.mkdir(parents=True, exist_ok=False)
    scripts = [Path(__file__).resolve(), Path(resume.__file__).resolve(),
               Path(r.__file__).resolve(), Path(resume.policy.__file__).resolve(),
               r.ROOT/"scripts/audit_je_paper_observables_20261020.py"]
    frozen = {str(x): r.sha(x) for x in scripts}
    r.save(out/"contract.json", dict(scope="FINITE_JOINT2_THEN_SUPG2",
        prior_controller_pid=a.prior_controller_pid, prior_launch=str(a.prior_launch.resolve()),
        prior_launch_sha256=r.sha(a.prior_launch), prior_finished=str(a.prior_finished.resolve()),
        script_hashes=frozen, threads=8, ooc_core_mb=8000, per_native_timeout_seconds=7200,
        original_resource_guards_unchanged=True, thresholds_changed=False, physics_changed=False,
        retries=0, max_native_processes=1, local_hole_failures_retained=True,
        production_solves_authorized_by_this_controller=False))
    seq = 0

    def event(state, **fields):
        nonlocal seq
        seq += 1
        record = dict(state=state, utc=datetime.now(timezone.utc).isoformat(), **fields)
        r.save(out/f"event_{seq:03d}.json", record)
        print(json.dumps(record, ensure_ascii=False), flush=True)

    def identities():
        for path, h in frozen.items():
            if r.sha(Path(path)) != h:
                raise RuntimeError("Frozen batch dependency changed: "+path)
        resume.evidence()

    def child(label, args):
        identities()
        resume.policy.preflight()
        event("STARTING", step=label, arguments=args)
        with (out/(label+"_stdout.txt")).open("xb") as stdout, (out/(label+"_stderr.txt")).open("xb") as stderr:
            proc = subprocess.Popen([sys.executable, "-X", "utf8"]+args, cwd=r.ROOT, stdout=stdout, stderr=stderr)
            event("CHILD_STARTED", step=label, pid=proc.pid)
            code = proc.wait()
        event("CHILD_EXITED", step=label, exit_code=code)
        if code != 0:
            raise RuntimeError(f"{label} failed with exit {code}; no retry or next step")

    try:
        # A held Windows handle tracks this process object, not a reused PID.
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.OpenProcess.argtypes = [ctypes.c_uint32, ctypes.c_int, ctypes.c_uint32]
        kernel.OpenProcess.restype = ctypes.c_void_p
        kernel.WaitForSingleObject.argtypes = [ctypes.c_void_p, ctypes.c_uint32]
        kernel.WaitForSingleObject.restype = ctypes.c_uint32
        kernel.CloseHandle.argtypes = [ctypes.c_void_p]
        handle = kernel.OpenProcess(0x00100000, False, a.prior_controller_pid)
        event("WAITING_FOR_EXISTING_JOINT1", pid=a.prior_controller_pid, process_handle_opened=bool(handle))
        if handle:
            start = time.monotonic()
            try:
                while True:
                    wait = kernel.WaitForSingleObject(handle, 1000)
                    if wait == 0:
                        break
                    if wait != 258:
                        raise OSError(ctypes.get_last_error(), "Process wait failed")
                    if time.monotonic()-start > 9000:
                        raise RuntimeError("Prior controller wait exceeded 9000s; prior process left untouched")
            finally:
                kernel.CloseHandle(handle)
        if not a.prior_finished.exists() or r.read(a.prior_finished).get("all_requested_attempts_completed") is not True:
            raise RuntimeError("Prior batch has no true completion record; no new solve started")
        for v in (6,8):
            resume.resolve_complete(r.name((40,20,16), v))
        event("PRIOR_JOINT1_COMPLETED", finished_sha256=r.sha(a.prior_finished))
        runner = str(Path(resume.__file__).resolve())
        reporter = str(r.ROOT/"scripts/audit_je_paper_observables_20261020.py")
        child("joint2", [runner, "--batch", "joint", "--joint-level", "2", "--execute"])
        for v in (6,8):
            resume.resolve_complete(r.name((48,22,18),v))
        child("after_joint2_report", [reporter, "--joint-level", "2", "--output", str(out/"after_joint2")])
        child("supg2", [runner, "--batch", "supg", "--joint-level", "2", "--execute"])
        for v in (6,8):
            for s in (.5,1.,2.):
                resume.resolve_complete(r.name((48,22,18),v,s))
        child("final_report", [reporter, "--joint-level", "2", "--output", str(out/"final_observables")])
        event("REQUESTED_BATCHES_COMPLETED_REVIEW_REQUIRED", production_acceptance=False)
    except BaseException as exc:
        event("STOPPED_REVIEW_REQUIRED", reason=repr(exc), retries=0, production_acceptance=False)
        raise


if __name__ == "__main__":
    main()
