"""Read-only audit of printed G12 controls, contrasts and exact-repeat logs.

This independent arithmetic check does not solve the PDE or validate a material.
Exit status is nonzero on any missing, nonfinite or inconsistent reported metric.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

RADII = ("060", "090", "120")
METRICS = (
    "total_joule_power", "global_joule_mean", "annulus_joule_mean",
    "annulus_joule_p95", "annulus_electric_mean", "annulus_n_relative_mean",
    "annulus_p_relative_mean", "annulus_temperature_rise_mean",
    "annulus_displacement_mean",
)


def parse(text: str) -> dict[str, str]:
    values = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        if key in values:
            raise ValueError(f"Duplicate output key: {key}")
        values[key] = value
    return values


def audit_values(values: dict[str, str]) -> dict:
    issues = []
    maximum_error = 0.0
    checked = 0
    annulus_errors = []
    if values.get("verification_result") != "PASS":
        issues.append("Solver did not report PASS")
    if values.get("verification_id") != "G12_radius_matrix":
        issues.append("Unexpected verification identifier")

    def number(key: str) -> float:
        try:
            value = float(values[key].replace("D", "E"))
        except (KeyError, ValueError):
            issues.append(f"Missing or invalid numeric field: {key}")
            return math.nan
        if not math.isfinite(value):
            issues.append(f"Nonfinite value: {key}")
        return value

    for radius in RADII:
        prefix = f"g12_r{radius}_"
        error = number(prefix + "h00_annulus_volume_relative_error")
        annulus_errors.append(error)
        if not (0 <= error < 0.02):
            issues.append(f"Invalid H00 annulus-volume error at r={radius}")
        for metric in METRICS:
            control = number(prefix + "h00_" + metric)
            hole = number(prefix + "h02_" + metric)
            contrast = number(prefix + "contrast_" + metric)
            if not all(map(math.isfinite, (control, hole, contrast))):
                continue
            if control <= 0 or hole < 0:
                issues.append(f"Invalid control or hole metric: r={radius}, {metric}")
                continue
            expected = (hole - control) / control
            discrepancy = abs(contrast - expected) / max(1.0, abs(expected))
            maximum_error = max(maximum_error, discrepancy)
            checked += 1
            if discrepancy > 1e-12:
                issues.append(f"Printed contrast does not match its control: r={radius}, {metric}")
    worst = number("g12_h00_annulus_volume_relative_error")
    if all(map(math.isfinite, annulus_errors)) and (
        not math.isfinite(worst) or abs(worst - max(annulus_errors)) > 1e-14
    ):
        issues.append("H00 summary is not the maximum of all three annulus errors")
    return {
        "status": "FAIL" if issues else "PASS",
        "contrast_checks": checked,
        "maximum_scaled_contrast_discrepancy": maximum_error,
        "issues": issues,
    }


def self_test() -> dict:
    values = {"verification_id": "G12_radius_matrix", "verification_result": "PASS"}
    for index, radius in enumerate(RADII, 1):
        values[f"g12_r{radius}_h00_annulus_volume_relative_error"] = str(index * 0.001)
        for metric in METRICS:
            prefix = f"g12_r{radius}_"
            values[prefix + "h00_" + metric] = str(float(index))
            values[prefix + "h02_" + metric] = str(1.5 * index)
            values[prefix + "contrast_" + metric] = "0.5"
    values["g12_h00_annulus_volume_relative_error"] = "0.003"
    assert audit_values(values)["status"] == "PASS"
    checked = 0
    for radius in RADII:
        for metric in METRICS:
            key = f"g12_r{radius}_h00_{metric}"
            for replacement in ("0", "-1", "NaN", "99"):
                bad = dict(values)
                bad[key] = replacement
                assert audit_values(bad)["status"] == "FAIL", (key, replacement)
                checked += 1
        key = f"g12_r{radius}_h00_annulus_volume_relative_error"
        for replacement in ("0.02", "NaN", "-1"):
            bad = dict(values)
            bad[key] = replacement
            assert audit_values(bad)["status"] == "FAIL", (key, replacement)
            checked += 1
    bad = dict(values)
    bad["g12_h00_annulus_volume_relative_error"] = "0.001"
    assert audit_values(bad)["status"] == "FAIL"
    checked += 1
    # Reproduce the original reporting bug: all printed controls are r=120.
    bad = dict(values)
    for radius in RADII[:-1]:
        for metric in METRICS:
            bad[f"g12_r{radius}_h00_{metric}"] = values[f"g12_r120_h00_{metric}"]
    assert audit_values(bad)["status"] == "FAIL"
    checked += 1
    return {"status": "PASS", "negative_checks": checked}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path, nargs="?")
    parser.add_argument("--repeat", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2))
        return 0
    if args.log is None:
        parser.error("Provide a G12 log or --self-test")
    raw = args.log.read_bytes()
    result = audit_values(parse(raw.decode("utf-8-sig")))
    result["log_sha256"] = hashlib.sha256(raw).hexdigest()
    if args.repeat is not None:
        repeat = args.repeat.read_bytes()
        result["repeat_sha256"] = hashlib.sha256(repeat).hexdigest()
        result["exact_repeat"] = raw == repeat
        if not result["exact_repeat"]:
            result["status"] = "FAIL"
            result["issues"].append("The repeat output is not byte-for-byte identical")
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
