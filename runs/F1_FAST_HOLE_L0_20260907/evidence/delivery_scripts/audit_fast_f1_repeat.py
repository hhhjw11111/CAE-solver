"""Close the representative-repeat audit after the main figures are generated."""
from pathlib import Path
import argparse
import json
import numpy as np
from finalize_fast_f1 import read, sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.run_root.resolve()
    record = read(root/"fast_audit.json")
    original, repeat = root/"FAST_H1_NT_L0", root/"REPEAT_H1_NT_L0"
    execution = read(repeat/"execution.json")
    if execution.get("status") != "PASS_NOT_MESH_ACCEPTANCE" or execution.get("exit_code") != 0:
        raise ValueError("Representative repeat has not finished successfully")
    identity = {name: sha(original/name) == sha(repeat/name)
                for name in ("metrics.txt", "nodes.csv", "field_samples.csv")}
    a = np.loadtxt(original/"nodes.csv", delimiter=",", skiprows=1)
    b = np.loadtxt(repeat/"nodes.csv", delimiter=",", skiprows=1)
    record["repeat"] = {
        "representative_repeat_executed": True,
        "file_byte_identity": identity,
        "maximum_absolute_field_difference": np.max(np.abs(a-b), axis=0).tolist(),
        "elapsed_seconds": execution["elapsed_seconds"],
        "original_provenance": read(original/"provenance.json"),
        "repeat_provenance": read(repeat/"provenance.json"),
    }
    p, q = record["repeat"]["original_provenance"], record["repeat"]["repeat_provenance"]
    record["checks"]["repeat_same_configuration"] = all(p[k] == q[k] for k in [
        "executable_sha256", "input_sha256", "threads", "storage", "MKL_DYNAMIC", "KMP_AFFINITY"])
    record["checks"]["representative_repeat_byte_identical"] = all(identity.values())
    record["checks"]["representative_repeat_accepted"] = True
    scales = np.array([1e-10, 1e-10, 1e-10, 3., .1, .1, 10.])
    record["repeat"]["maximum_solver_scaled_state_difference"] = float(
        np.max(np.abs(a[:, 4:11]-b[:, 4:11])/scales))
    record["repeat"]["difference_columns"] = (
        "node,x_m,y_m,z_m,ux_m,uy_m,uz_m,phi_V,eta_n,eta_p,theta_K,n_m3,p_m3").split(",")
    record["case_numerical_checks_passed"] = all(
        v for k, v in record["checks"].items() if "repeat" not in k)
    record["numerical_checks_passed"] = all(record["checks"].values())
    record["status"] = ("SIX_CASE_REDUCED_MESH_MATRIX_COMPLETE" if record["numerical_checks_passed"]
                        else "AUDIT_NEEDS_ATTENTION")
    record["delivery_scope_status"] = (
        "SIX_CASE_L0_COMPLETE_STRICT_REPRODUCTION_NOT_PASSED"
        if record["case_numerical_checks_passed"] and not all(identity.values()) else record["status"])
    record["strict_reproduction_status"] = "PASS" if all(identity.values()) else "FAIL_BYTE_IDENTITY"
    record["reproduction_note"] = (
        "Original byte-identity gate retained and not passed. Small state differences are reported "
        "as diagnostics, not used to replace or relax that gate. All original files are preserved.")
    (root/"fast_audit.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps({"status": record["status"], "repeat": identity}, indent=2))
    if not record["numerical_checks_passed"]:
        raise ValueError("Repeat audit failed; retained all evidence without changing acceptance")


if __name__ == "__main__":
    main()
