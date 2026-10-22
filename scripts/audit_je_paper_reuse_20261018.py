"""Read-only revalidation of existing JE verification evidence for the approved paper scope.

No native executable is started. Existing full-suite failure and all original files stay intact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "runs/JE_RECONSTRUCTION_20261018"
SCOPE = "JE_PAPER_MINIMAL_20261018"
REQUIRED = ["A1", "A2", "A1JJE", "JEGRAD", "G02JE", "F1MJE", "G04JE", "SUPG", "SUPGJ", "SUPGJT"]
LEGACY = ["V01", "V02", "V03", "V04", "V05", "V06", "G02", "G03", "G04", "G05", "G06",
          "G07", "G08", "L01", "F1M", "A1J"]
LATER_VERIFICATION_ONLY = {
    "src/verification/a2_piezo_check_m.f90": "Excluded A2P gauge/Gauss/matrix diagnostics added after this frozen executable was built.",
    "src/verification/verification_driver_m.f90": "Later verification-only A2PMAT registration; not compiled into the reused executable.",
}
EXPECTED_EXES = {
    "baseline": "405e953c8cb0b1eca2a61eccc011945892924b43fca6805e87fd8fec19006d92",
    "candidate": "ad50431beee4c732268fdf384d503cf4491321212aad568a52131457345fe609",
}


def sha(path: Path) -> str:
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_json(path: Path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def relative_key(value: str) -> str:
    return str(value).replace("\\", "/")


def require(condition: bool, message: str):
    if not condition:
        raise ValueError(message)


def audit_reuse(scope_contract: Path) -> dict:
    """Recompute all checks from actual files; return an acceptance payload or raise ValueError."""
    scope_contract = Path(scope_contract).resolve()
    evidence: dict[str, str] = {}

    def record(path: Path, expected: str | None = None) -> str:
        path = Path(path).resolve()
        require(path.is_file(), f"Missing evidence file: {path}")
        digest = sha(path)
        if expected is not None:
            require(digest == expected, f"Evidence SHA-256 mismatch: {path}")
        key = str(path)
        require(key not in evidence or evidence[key] == digest, f"Evidence changed during audit: {path}")
        evidence[key] = digest
        return digest

    record(Path(__file__))
    record(scope_contract)
    scope = read_json(scope_contract)
    require(scope.get("scope") == SCOPE and scope.get("user_approved") is True, "Wrong or unapproved scope")
    require(scope.get("required_ids") == REQUIRED, "Required targeted ID list differs from the approved ten")
    require(scope.get("required_legacy_ids") == LEGACY, "Required legacy ID list differs from the frozen sixteen")
    require(scope.get("excluded_ids") == ["A2P"], "Unexpected excluded verification IDs")
    require(bool(scope.get("authorization")) and bool(scope.get("exclusion_reason")), "Missing scope authorization")
    require(scope.get("required_thresholds_relaxed") is False, "Threshold relaxation is not permitted")
    contract_path = Path(scope["original_contract"]).resolve()
    record(contract_path, scope["original_contract_sha256"])
    record(Path(scope["approved_workflow"]), scope["approved_workflow_sha256"])
    contract = read_json(contract_path)
    record(ROOT / contract["protocol"], contract["protocol_sha256"])
    for key, digest in contract["frozen_gate_scripts"].items():
        record(ROOT / key, digest)
    for name in ("a2p_protocol.json", "jegrad_protocol.json"):
        record(STUDY / "stage3_verification" / name)
    attempt = Path(scope["source_attempt"]).resolve()
    require(attempt == (STUDY / "stage3_verification/attempt01").resolve(), "Unexpected reuse source attempt")
    for name in ("provenance.json", "result.json", "legacy_cli.json", "targeted.json"):
        record(attempt / name)
    provenance = read_json(attempt / "provenance.json")
    original_result = read_json(attempt / "result.json")
    legacy_records = read_json(attempt / "legacy_cli.json")
    targeted_records = read_json(attempt / "targeted.json")
    frozen_sources = {relative_key(p): h for p, h in provenance["sources"].items()}
    require(len(frozen_sources) == len(provenance["sources"]), "Duplicate normalized provenance sources")
    core_records = contract["core_sources_before"]
    require(len(core_records) == 40, "Expected exactly forty frozen core sources")
    core = {}
    for item in core_records:
        key = relative_key(item["path"])
        require(key.startswith("solver_v2/src/") and "/verification/" not in key, "Invalid core source record")
        key = key.removeprefix("solver_v2/")
        require(key not in core, "Duplicate core source record")
        require(frozen_sources.get(key) == item["sha256"], f"Core contract/provenance mismatch: {key}")
        record(ROOT / key, item["sha256"])
        core[key] = item["sha256"]
    require({relative_key(x["path"]).removeprefix("solver_v2/"): x["before"] for x in provenance["core_checks"]}
            == core, "Attempt core-check inventory differs from contract")
    require(all(x["before"] == x["after"] for x in provenance["core_checks"]), "Original core-check mismatch")
    require(provenance.get("core_sources_unchanged") is True, "Attempt did not verify unchanged core")
    current_source_paths = {p.relative_to(ROOT).as_posix(): p for p in (ROOT / "src").rglob("*.f90")}
    require(set(current_source_paths) == set(frozen_sources), "Current source inventory changed; review required")
    source_differences = []
    unchanged_sources = []
    for key, expected in frozen_sources.items():
        actual = record(current_source_paths[key])
        if actual != expected:
            require(key in LATER_VERIFICATION_ONLY, f"Relevant source changed since verification: {key}")
            source_differences.append(dict(path=key, executed_build_sha256=expected, current_sha256=actual,
                                           explanation=LATER_VERIFICATION_ONLY[key]))
        else:
            unchanged_sources.append(key)
    # Required tests must not acquire a dependency on the excluded analytic extension.
    for key, path in current_source_paths.items():
        if key in LATER_VERIFICATION_ONLY:
            continue
        text = path.read_text(encoding="utf-8-sig")
        require(not re.search(r"^\s*use\s+(?:,.*?::\s*)?a2_piezo_check_m\b", text, re.I | re.M),
                f"Required source depends on excluded A2P module: {key}")
    require("src/verification/production_joule_check_m.f90" in unchanged_sources
            and "src/verification/je_gradient_check_m.f90" in unchanged_sources,
            "New heat verification module changed")
    builds = {}
    executables = {}
    for label in ("baseline", "candidate"):
        exe = Path(provenance[label]).resolve()
        require(provenance[label + "_sha256"] == EXPECTED_EXES[label], f"Unexpected {label} frozen executable")
        digest = record(exe, EXPECTED_EXES[label])
        require(original_result[label + "_sha256"] == digest
                and Path(original_result[label]).resolve() == exe, f"{label} result/provenance mismatch")
        manifest_path = exe.parent / "build_manifest.json"
        record(manifest_path)
        manifest = read_json(manifest_path)
        require(manifest.get("schema") == "TPSolverV2-build-manifest-v1"
                and manifest.get("version_verified_against_source_and_executable") is True,
                f"Invalid {label} build manifest")
        entries = {relative_key(r["path"]): r for r in manifest["files"]}
        require(len(entries) == len(manifest["files"]), f"Duplicate {label} manifest entries")
        exe_key = exe.relative_to(ROOT).as_posix()
        require(entries.get(exe_key, {}).get("sha256") == digest
                and entries[exe_key]["bytes"] == exe.stat().st_size, f"{label} executable/manifest mismatch")
        manifest_sources = {p: r["sha256"] for p, r in entries.items() if p.startswith("src/")}
        for key, expected in core.items():
            require(manifest_sources.get(key) == expected, f"{label} core manifest mismatch: {key}")
        if label == "candidate":
            require(manifest_sources == frozen_sources, "Candidate build sources differ from execution provenance")
        record(ROOT / "VERSION", entries["VERSION"]["sha256"])
        builds[label] = dict(manifest=str(manifest_path), manifest_sha256=evidence[str(manifest_path)],
            generated_utc=manifest.get("generated_utc"), source_count=len(manifest_sources),
            verified_core_count=40, executable_verified=True,
            frozen_execution_sources_verified=label == "candidate",
            verification_scope=("all source hashes against execution provenance plus current unchanged core/modules"
                                if label == "candidate" else "forty physical/core source hashes plus executable identity"),
            non_source_document_entries_are_not_current_source_claims=True)
        executables[label] = dict(path=str(exe), sha256=digest)

    def check_log(item: dict, check: str, label: str, executable: Path, passed: bool = True) -> dict:
        log = attempt / f"{check}_{label}.txt"
        require(item.get("check") == check and item.get("label") == label, f"Log record identity: {check}/{label}")
        require(Path(item["log"]).resolve() == log.resolve(), f"Unexpected log path: {check}/{label}")
        command = item.get("command", [])
        require(len(command) >= 3 and Path(command[0]).resolve() == executable
                and command[1:3] == ["verify", check], f"Unexpected command: {check}/{label}")
        digest = record(log, item["sha256"])
        content = log.read_bytes()
        pass_marker = b"verification_result=PASS" in content
        fail_marker = b"verification_result=FAIL" in content
        if passed:
            require(item.get("exit_code") == 0 and item.get("timed_out") is False and item.get("passed") is True,
                    f"Unsuccessful execution record: {check}/{label}")
            require(pass_marker and not fail_marker, f"Log is not a verified PASS: {check}/{label}")
        else:
            require(item.get("exit_code") not in (0, None) and item.get("passed") is False
                    and item.get("timed_out") is False and fail_marker,
                    "Original A2P failure is missing or relabeled")
        return dict(check=check, label=label, exit_code=item["exit_code"], timed_out=item["timed_out"],
                    log=str(log.resolve()), sha256=digest, pass_marker=pass_marker,
                    fail_marker=fail_marker, passed=passed)

    require([r.get("check") for r in legacy_records] == LEGACY, "Legacy inventory missing, duplicate, or reordered")
    legacy = []
    for item in legacy_records:
        check = item["check"]
        baseline = check_log(item["baseline"], check, "baseline", Path(executables["baseline"]["path"]))
        candidate = check_log(item["candidate"], check, "candidate", Path(executables["candidate"]["path"]))
        identical = Path(baseline["log"]).read_bytes() == Path(candidate["log"]).read_bytes()
        require(identical and item.get("byte_identical") is True and item.get("passed") is True,
                f"Legacy bytes differ: {check}")
        legacy.append(dict(check=check, baseline=baseline, candidate=candidate, byte_identical=True, passed=True))
    require(len(targeted_records) == 11 and len({r.get("check") for r in targeted_records}) == 11,
            "Original targeted inventory missing or duplicated")
    target_map = {r["check"]: r for r in targeted_records}
    require(set(target_map) == set(REQUIRED + ["A2P"]), "Original targeted inventory differs")
    targeted = [check_log(target_map[name], name, "targeted", Path(executables["candidate"]["path"]))
                for name in REQUIRED]
    excluded = check_log(target_map["A2P"], "A2P", "targeted", Path(executables["candidate"]["path"]), passed=False)
    require(original_result.get("all_passed") is False
            and original_result.get("targeted_passed") == 10 and original_result.get("targeted_total") == 11
            and original_result.get("legacy_passed") == 16 and original_result.get("legacy_total") == 16,
            "Original full-suite result was changed or does not match source evidence")
    for flag in ("core_sources_unchanged", "frozen_gates_unchanged", "build_manifests_verified", "source_still_same"):
        require(original_result.get(flag) is True, f"Original execution identity check failed: {flag}")
    # A1/A2 exported scientific fields are retained and inventoried; these hashes attest current preservation.
    exports = {}
    for name in ("A1", "A2"):
        files = sorted(p for p in (attempt / f"{name}_fields").rglob("*") if p.is_file())
        require(files, f"Missing original {name} exports")
        exports[name] = {str(p.resolve()): record(p) for p in files}
    result = dict(scope=SCOPE, scope_contract=str(scope_contract), verification_method="REUSE_REVALIDATION_NO_NATIVE_RUN",
        all_required_passed=True, core_sources_unchanged=True, core_source_count=40,
        build_manifests_verified=True, required_ids=REQUIRED, required_legacy_ids=LEGACY,
        legacy_passed=16, legacy_total=16, targeted_passed=10, targeted_total=10,
        source_attempt=str(attempt), baseline=executables["baseline"]["path"],
        baseline_sha256=executables["baseline"]["sha256"], candidate=executables["candidate"]["path"],
        candidate_sha256=executables["candidate"]["sha256"], executables_verified=executables,
        build_identity=builds, unchanged_source_paths=unchanged_sources,
        current_source_differences=source_differences, current_source_equals_frozen_candidate=not source_differences,
        candidate_full_source_identity_scope="Recorded executed snapshot, not all current source files",
        original_full_suite_passed=False, original_failure_retained=True,
        excluded_checks=[dict(**excluded, exclusion_reason=scope["exclusion_reason"], numerical_failure_overridden=False)],
        legacy_checks=legacy, targeted_checks=targeted, retained_A1_A2_exports=exports,
        core_source_hashes=core, evidence_file_hashes=evidence,
        limitations=["No native process was launched by this audit.",
                    "A2P remains failed in its original broader scope; only the user-approved required set changed.",
                    "Later A2P/driver source edits are not certified as a rebuilt executable.",
                    "This verification certificate does not certify mesh/SUPG convergence or production cases.",
                    "Original manifest documentation entries are historical snapshots, not claims that evolving plans are unchanged."])
    return result


def validate_acceptance(path: Path, expected_exe: Path | None = None) -> dict:
    """Fail closed on changed files or claims; no subprocess, file mutation, or native replay."""
    path = Path(path).resolve()
    accepted = read_json(path)
    require(accepted.get("scope") == SCOPE and accepted.get("all_required_passed") is True,
            "Not an accepted JE paper-scope evidence record")
    hashes = accepted.get("evidence_file_hashes")
    require(isinstance(hashes, dict) and bool(hashes), "Missing evidence hash inventory")
    for source, expected in hashes.items():
        require(Path(source).is_file() and sha(Path(source)) == expected, f"Acceptance evidence changed: {source}")
    fresh = audit_reuse(Path(accepted["scope_contract"]))
    for key, value in fresh.items():
        require(accepted.get(key) == value, f"Acceptance claim differs from actual re-audit: {key}")
    if expected_exe is not None:
        exe = Path(expected_exe).resolve()
        candidates = {Path(item["path"]).resolve(): item["sha256"] for item in fresh["executables_verified"].values()}
        require(exe in candidates and sha(exe) == candidates[exe], "Requested executable is not in verified reuse identities")
    return accepted


def write_new_json(path: Path, value: dict):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False)
        stream.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope-contract", type=Path, default=STUDY / "paper_scope/contract.json")
    parser.add_argument("--output", type=Path, default=STUDY / "paper_scope/verification_reuse")
    parser.add_argument("--validate", type=Path, help="Validate an existing certificate without writing anything")
    args = parser.parse_args()
    if args.validate:
        result = validate_acceptance(args.validate)
        print(json.dumps({k: result[k] for k in ("scope", "all_required_passed", "legacy_passed", "targeted_passed")}))
        return 0
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    try:
        result = audit_reuse(args.scope_contract)
    except (ValueError, KeyError, OSError, TypeError) as exc:
        write_new_json(output / "audit_report.json", dict(scope=SCOPE, all_required_passed=False,
            verification_method="REUSE_REVALIDATION_NO_NATIVE_RUN", failure=str(exc)))
        print(str(exc))
        return 1
    result["audit_created_utc"] = datetime.now(timezone.utc).isoformat()
    write_new_json(output / "audit_report.json", result)
    readme = ("# 已有验证证据复用审计\n\n"
        "本目录只读重验已有日志、退出状态、逐字节回归及源码/构建哈希，没有重新运行求解器。\n\n"
        "用户批准范围 JE_PAPER_MINIMAL_20261018：10 项专项全部通过，16 项旧 CLI 基准/候选逐字节一致。"
        "40 个物理/核心源文件与原协议及两份构建清单一致。\n\n"
        "原 11 项专项中的 A2P 仍为 FAIL，原整体 result.json 的 all_passed=false 保留。"
        "A2P 因用户明确缩小论文范围而不再是必做项，并未获得数值豁免或改判。\n\n"
        "执行身份：baseline ARC_JE_FINAL_20261012 与 candidate JE_VERIFY_20261018_FINAL 的可执行文件"
        "均重新核对 SHA-256。候选完整源码清单与当次 provenance 一致；当前 A2P 模块及验证 driver"
        "存在后来修改，逐一记录于 current_source_differences，不能声称当前全部源码等于已运行候选源码。"
        "所有其余源文件（包括 A1/A2、A1JJE、JEGRAD 及物理依赖）重新核对未变。\n\n"
        "两份 build manifest 的核验范围是执行文件、冻结核心，以及候选执行时完整源码清单；"
        "其中历史文档条目不被误称为当前执行计划仍未变化。\n\n"
        "acceptance.json 只接受本范围；网格、局部载流子、SUPG 和9个生产工况仍需后续验收。"
        "下游应调用 validate_acceptance(path, expected_exe=...) 重新检查实际文件，而非只读 PASS 布尔值。\n")
    with (output / "README.md").open("x", encoding="utf-8") as stream:
        stream.write(readme)
    write_new_json(output / "acceptance.json", result)
    validate_acceptance(output / "acceptance.json")
    print(json.dumps(dict(scope=SCOPE, all_required_passed=True, legacy_passed=16, targeted_passed=10,
                         source_differences=len(result["current_source_differences"]),
                         evidence_files=len(result["evidence_file_hashes"])), ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
