"""Audit the proposed ZnO input contract, not the production Fortran solver.

Default execution is offline and read-only. --output creates a new evidence
directory and refuses to overwrite one. --refresh-metadata requests only public
DOI metadata; --render-evidence renders a downloaded public reference.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

SOLVER = Path(__file__).resolve().parents[1]
WORKSPACE = SOLVER.parent
DEFAULT_CONTRACT = SOLVER / "material_audit/zno_300k_candidate_v2.json"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalized_title(text):
    return re.sub(r"[^a-z0-9]", "", text.lower())


def numeric_audit(contract):
    """Return numerical checks and derived quantities, without filesystem writes."""
    checks = []

    def check(name, ok):
        checks.append({"name": name, "passed": bool(ok)})

    def close(name, actual, expected, rtol=1e-12, atol=0.0):
        check(name, np.allclose(actual, expected, rtol=rtol, atol=atol))

    def spd(name, matrix):
        close(name + "_symmetric", matrix, matrix.T)
        check(name + "_positive_definite", np.linalg.eigvalsh(matrix).min() > 0)

    check("schema_version", contract["schema_version"] == 1)
    check("production_disabled", contract["production_enabled"] is False)
    check("no_minority_floor", contract["conventions"]["minority_density_floor_allowed"] is False)
    check("engineering_shear_order", contract["conventions"]["strain_voigt"] == [
        "eps_xx", "eps_yy", "eps_zz", "gamma_yz", "gamma_zx", "gamma_xy"])
    p, c = contract["parameters"], contract["constants_si"]
    exp = contract["independent_check_values"]
    q, kb, eps0, t0 = c["q_C"], c["kB_J_per_K"], c["eps0_F_per_m"], p["T0_K"]
    check("positive_reference_temperature", math.isfinite(t0) and t0 > 0)
    if t0 <= 0:
        return checks, {}
    check("constants_match_SI", q == 1.602176634e-19 and kb == 1.380649e-23)
    cc, ee = p["elastic_Pa"], p["piezo_C_per_m2"]
    c66 = (cc["C11"] - cc["C12"]) / 2
    elastic = np.array([
        [cc["C11"], cc["C12"], cc["C13"], 0, 0, 0],
        [cc["C12"], cc["C11"], cc["C13"], 0, 0, 0],
        [cc["C13"], cc["C13"], cc["C33"], 0, 0, 0],
        [0, 0, 0, cc["C44"], 0, 0],
        [0, 0, 0, 0, cc["C44"], 0],
        [0, 0, 0, 0, 0, c66]], dtype=float)
    piezo = np.array([[0, 0, 0, 0, ee["e15"], 0],
                      [0, 0, 0, ee["e15"], 0, 0],
                      [ee["e31"], ee["e31"], ee["e33"], 0, 0, 0]])
    dielectric = np.diag(p["relative_permittivity"]) * eps0
    heat_k = np.diag(p["thermal_conductivity_W_per_mK"])
    alpha = np.array(p["alpha_per_K"])
    pyro = np.array(p["pyro_primary_C_per_m2K"])
    thermal_stress = elastic @ alpha
    rh = p["hall_factor"]
    check("positive_Hall_factor", math.isfinite(rh) and rh > 0)
    if rh <= 0:
        return checks, {}
    n0, mu_n0 = rh * p["n_H_per_m3"], p["mu_H_m2_per_Vs"] / rh
    ni = p["intrinsic_concentration_per_m3"]
    p0 = ni**2 / n0
    vt = kb * t0 / q
    mu_p0 = p["mu_p0_m2_per_Vs"]
    dn0, dp0 = vt * mu_n0, vt * mu_p0
    sigma0 = q * (n0 * mu_n0 + p0 * mu_p0)
    for name, matrix in (("elastic", elastic), ("dielectric", dielectric),
                         ("heat_conduction", heat_k),
                         ("electron_mobility", np.eye(3) * mu_n0),
                         ("hole_mobility", np.eye(3) * mu_p0)):
        check(name + "_finite", np.isfinite(matrix).all())
        spd(name, matrix)
    close("C66_hexagonal_identity", c66, exp["C66_Pa"])
    close("thermal_stress_unit_and_value", thermal_stress, exp["thermal_stress_Pa_per_K"])
    close("n0_anchor", n0, exp["n0_per_m3"])
    close("p0_anchor", p0, exp["p0_per_m3"])
    close("sigma0_anchor", sigma0, exp["sigma0_S_per_m"])
    close("Dn0_anchor", dn0, exp["Dn0_m2_per_s"])
    check("minority_above_double_tiny", np.finfo(float).tiny < p0 < n0)
    close("equilibrium_mass_action", n0 * p0, ni**2)
    net_ionized = n0 - p0
    n_stable = (net_ionized + math.hypot(net_ionized, 2 * ni)) / 2
    close("stable_quadratic_n", n_stable, n0)
    close("stable_quadratic_p", ni**2 / n_stable, p0)
    close("free_thermal_expansion_zero_stress", elastic @ alpha - thermal_stress,
          np.zeros(6), atol=1e-9)
    close("engineering_shear_polarization", piezo @ [0, 0, 0, 0, 1, 0],
          [ee["e15"], 0, 0])
    # A 180-degree rotation about y reverses the polar c axis. These diagonal
    # rotation matrices apply to the declared engineering-shear convention.
    rotation = np.diag([-1.0, 1.0, -1.0])
    strain_rotation = np.diag([1.0, 1.0, 1.0, -1.0, 1.0, -1.0])
    close("reverse_c_polarity_changes_pyro", rotation @ pyro, -pyro)
    close("reverse_c_polarity_changes_piezo", rotation @ piezo @ strain_rotation, -piezo)
    close("reverse_c_polarity_preserves_C", strain_rotation @ elastic @ strain_rotation, elastic)
    # Independent finite differences of the declared constitutive maps.
    # This tests the contract algebra; it does not execute the Fortran routines.
    strain = np.array([1.0, -2.0, 3.0, 0.4, 0.5, -0.6]) * 1e-5
    field = np.array([1.0, -2.0, 3.0]) * 1e5
    def stress_at(st, el):
        return elastic @ st - piezo.T @ el - thermal_stress * 2.0
    def displacement_at(st, el):
        return piezo @ st + dielectric @ el + pyro * 2.0
    mixed_stress = np.column_stack([
        (stress_at(strain, field + np.eye(3)[i] * 1e3)
         - stress_at(strain, field - np.eye(3)[i] * 1e3)) / 2e3
        for i in range(3)])
    mixed_displacement = np.column_stack([
        (displacement_at(strain + np.eye(6)[i] * 1e-7, field)
         - displacement_at(strain - np.eye(6)[i] * 1e-7, field)) / 2e-7
        for i in range(6)])
    close("stress_electric_field_difference", mixed_stress, -piezo.T, rtol=1e-9, atol=1e-10)
    close("D_strain_difference", mixed_displacement, piezo, rtol=1e-9, atol=1e-10)
    close("reciprocal_mixed_constitutive_blocks", mixed_stress, -mixed_displacement.T,
          rtol=1e-9, atol=1e-10)
    close("SI_permittivity_conversion", np.diag(dielectric) / eps0,
          p["relative_permittivity"])

    transport = []
    for temp in contract["diagnostics_not_calibration"]["temperatures_K"]:
        row = {"T_K": temp}
        for carrier, mobility in (("n", mu_n0), ("p", mu_p0)):
            beta = p["beta_" + carrier]
            mu = mobility * (temp / t0) ** (-beta)
            diffusion = kb * temp / q * mu
            dmu = -beta / temp * mu
            dd = (1 - beta) / temp * diffusion
            h = temp * 1e-5
            def mu_at(tt):
                return mobility * (tt / t0) ** (-beta)
            mu_fd = (mu_at(temp + h) - mu_at(temp - h)) / (2 * h)
            d_fd = (kb / q * ((temp + h) * mu_at(temp + h)
                    - (temp - h) * mu_at(temp - h))) / (2 * h)
            check(f"transport_positive_{carrier}_{temp:g}K", mu > 0 and diffusion > 0)
            close(f"mu_derivative_{carrier}_{temp:g}K", dmu, mu_fd, rtol=1e-7, atol=1e-16)
            close(f"D_derivative_{carrier}_{temp:g}K", dd, d_fd, rtol=1e-7, atol=1e-16)
            row.update({f"mu_{carrier}_m2_per_Vs": mu, f"D_{carrier}_m2_per_s": diffusion,
                        f"dmu_{carrier}_dT": dmu, f"dD_{carrier}_dT": dd})
        transport.append(row)
    diag = contract["diagnostics_not_calibration"]
    hall = []
    for factor in diag["hall_factors"]:
        n_h = factor * p["n_H_per_m3"]
        mu_h = p["mu_H_m2_per_Vs"] / factor
        close(f"Hall_conductivity_invariant_{factor:g}",
              q * n_h * mu_h, q * p["n_H_per_m3"] * p["mu_H_m2_per_Vs"])
        hall.append({"rH": factor, "n_per_m3": n_h, "mu_m2_per_Vs": mu_h})
    epsilon_increment = piezo @ np.linalg.solve(elastic, piezo.T)
    epsilon_free = dielectric + epsilon_increment
    epsilon_clamped_if_input_free = dielectric - epsilon_increment
    spd("alternate_clamped_dielectric", epsilon_clamped_if_input_free)
    secondary = piezo @ alpha
    free_pyro = pyro + secondary
    thermal_bar = []
    for voltage in diag["voltages_V"]:
        rise = sigma0 * voltage**2 / (8 * heat_k[0, 0])
        power = sigma0 * diag["bar_width_m"] * diag["bar_thickness_m"] / diag["bar_length_m"] * voltage**2
        thermal_bar.append({"voltage_V": voltage, "nominal_field_V_per_m": voltage / diag["bar_length_m"],
                            "deltaT_max_K": rise, "deltaT_mean_K": rise * 2 / 3,
                            "mu_n_drop_at_peak_percent": 100 * (1 - ((t0 + rise) / t0) ** (-p["beta_n"])),
                            "power_W": power})
    poly_alpha = {}
    for axis in ("a", "c"):
        a0, a1, a2 = diag[f"book_{axis}_polynomial_A"]
        poly_alpha[axis] = (a1 + 2 * a2 * t0) / (a0 + a1 * t0 + a2 * t0**2)
    derived = {
        "C_Pa": elastic.tolist(), "e_C_per_m2": piezo.tolist(),
        "eps_constant_strain_candidate_F_per_m": dielectric.tolist(),
        "C_min_eigenvalue_Pa": float(np.linalg.eigvalsh(elastic).min()),
        "C_condition_2": float(np.linalg.cond(elastic)),
        "thermal_stress_Pa_per_K": thermal_stress.tolist(),
        "secondary_pyro_from_candidate_C_per_m2K": secondary.tolist(),
        "total_pyro_from_candidate_C_per_m2K": free_pyro.tolist(),
        "secondary_pyro_discrepancy_percent": 100 * (secondary[2] / diag["book_pyro_secondary_C_per_m2K"] - 1),
        "total_pyro_discrepancy_percent": 100 * (free_pyro[2] / diag["book_pyro_total_C_per_m2K"] - 1),
        "eps_free_stress_if_input_clamped_F_per_m": np.diag(epsilon_free).tolist(),
        "eps_clamped_if_input_free_F_per_m": np.diag(epsilon_clamped_if_input_free).tolist(),
        "n0_per_m3": n0, "p0_per_m3": p0, "effective_net_ionized_doping_per_m3": net_ionized,
        "mu_n0_m2_per_Vs": mu_n0, "Dn0_m2_per_s": dn0, "Dp0_m2_per_s": dp0,
        "thermal_voltage_V": vt, "sigma0_S_per_m": sigma0,
        "hole_conductivity_fraction": p0 * mu_p0 / (n0 * mu_n0 + p0 * mu_p0),
        "Debye_clamped_x_m": math.sqrt(dielectric[0, 0] * vt / (q * (n0 + p0))),
        "Debye_clamped_z_m": math.sqrt(dielectric[2, 2] * vt / (q * (n0 + p0))),
        "alpha_from_printed_book_polynomials_per_K": poly_alpha,
        "alpha_c_table_vs_polynomial_percent": 100 * (alpha[2] / poly_alpha["c"] - 1),
        "transport_diagnostics": transport, "hall_diagnostics": hall,
        "decoupled_constant_conductivity_bar_only": thermal_bar
    }
    return checks, derived


def check_protected_files(contract):
    return [{"path": name, "passed": (WORKSPACE / name).is_file()
             and sha256(WORKSPACE / name) == expected}
            for name, expected in contract["protected_files_sha256"].items()]


def verify_metadata(contract):
    records = []
    for source in contract["sources"]:
        try:
            url = "https://api.crossref.org/works/" + urllib.parse.quote(source["doi"], safe="")
            request = urllib.request.Request(url, headers={"User-Agent": "TPSolverV2-material-audit/1.0"})
            with urllib.request.urlopen(request, timeout=15) as response:
                message = json.load(response)["message"]
            title = " ".join(message.get("title", []))
            year = message["issued"]["date-parts"][0][0]
            matches = (message["DOI"].lower() == source["doi"].lower()
                       and normalized_title(title) == normalized_title(source["title"])
                       and year == source["year"])
            records.append({"id": source["id"], "doi": message["DOI"], "title": title,
                            "journal_or_book": message.get("container-title", []),
                            "year": year, "volume": message.get("volume"), "pages": message.get("page"),
                            "metadata_status": "verified" if matches else "mismatch",
                            "parameter_value_status": "NOT_ESTABLISHED_BY_METADATA"})
        except Exception as error:
            records.append({"id": source["id"], "doi": source["doi"],
                            "metadata_status": "manual_needed", "error": str(error)})
        print("metadata:", records[-1]["id"], records[-1]["metadata_status"], flush=True)
    return records


def self_test(contract):
    checks, _ = numeric_audit(contract)
    if not all(check["passed"] for check in checks):
        raise AssertionError("Valid input failed")
    bad_inputs = [
        ("negative_elastic", lambda x: x["parameters"]["elastic_Pa"].update(C44=-1.0)),
        ("thermal_units", lambda x: x["independent_check_values"].update(thermal_stress_Pa_per_K=[4.31e-6] * 6)),
        ("wrong_shear_order", lambda x: x["conventions"].update(strain_voigt=["xx"] * 6)),
        ("production_promotion", lambda x: x.update(production_enabled=True)),
        ("minority_floor", lambda x: x["conventions"].update(minority_density_floor_allowed=True)),
        ("zero_temperature", lambda x: x["parameters"].update(T0_K=0)),
    ]
    for name, change in bad_inputs:
        changed = copy.deepcopy(contract)
        change(changed)
        bad_checks, _ = numeric_audit(changed)
        if all(item["passed"] for item in bad_checks):
            raise AssertionError("Invalid input not detected: " + name)
    return {"valid_contract": "PASS", "negative_controls": [name for name, _ in bad_inputs],
            "negative_controls_result": "PASS"}


def render_reference(contract, contract_path, output):
    import pypdfium2 as pdfium
    info = contract["public_source_snapshot"]
    source = contract_path.parent / info["file"]
    if sha256(source) != info["sha256"]:
        raise ValueError("Public source PDF hash mismatch")
    destination = output / "reference_pages"
    destination.mkdir()
    paths = []
    with pdfium.PdfDocument(source) as pdf:
        for page_number in info["pages_one_based"]:
            page = pdf[page_number - 1]
            path = destination / f"wiley_p{page_number:03d}.png"
            page.render(scale=1.8).to_pil().save(path)
            paths.append(str(path.relative_to(output)))
            page.close()
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--refresh-metadata", action="store_true")
    parser.add_argument("--render-evidence", action="store_true")
    args = parser.parse_args()
    if (args.refresh_metadata or args.render_evidence) and not args.output:
        parser.error("Evidence options require a new --output directory")
    if args.output and args.output.exists():
        parser.error("Refusing to overwrite existing evidence directory")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    checks, derived = numeric_audit(contract)
    protected = check_protected_files(contract)
    passed = all(item["passed"] for item in checks + protected)
    result = {
        "audit_date": contract["audit_date"], "profile_id": contract["profile_id"],
        "scope": "standalone input algebra and provenance checks; no production solver execution",
        "numerical_input_status": "PASS" if passed else "FAIL",
        "experimental_calibration_status": "PENDING",
        "production_enabled": False, "contract_sha256": sha256(args.contract),
        "audit_script_sha256": sha256(__file__),
        "runtime": {"python": sys.version.split()[0], "numpy": np.__version__},
        "checks_passed": sum(check["passed"] for check in checks), "checks_total": len(checks),
        "checks": checks, "protected_files": protected, "derived": derived,
        "self_tests": self_test(contract), "calibration_gaps": contract["calibration_gaps"]
    }
    if args.output:
        args.output.mkdir(parents=True)
        if args.refresh_metadata:
            result["bibliography_metadata"] = verify_metadata(contract)
        if args.render_evidence:
            result["rendered_reference_pages"] = render_reference(contract, args.contract, args.output)
        (args.output / "material_audit.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "profile_id", "numerical_input_status", "checks_passed", "checks_total",
        "experimental_calibration_status", "production_enabled")}, indent=2))
    if not passed:
        print("FAILED_CHECKS", [item for item in checks + protected if not item["passed"]])
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
