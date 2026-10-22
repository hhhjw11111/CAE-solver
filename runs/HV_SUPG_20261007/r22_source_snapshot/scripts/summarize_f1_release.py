"""Complete F1 table quantities from saved FE fields, without solving a PDE."""
from pathlib import Path
import argparse
import csv
import hashlib
import json

import numpy as np
import audit_f1_results as audit
import plot_f1_main_results as fe

MAIN = [
    "MAIN_H0_LQ_L2", "MAIN_H0_N0_L2", "MAIN_H0_NT_L2",
    "MAIN_H1_LQ_L2", "MAIN_H1_N0_L2", "MESH_H1_NT_L2",
]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def percentile(values, weights, fraction=.95):
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order[np.searchsorted(cumulative, fraction*cumulative[-1])]])


def supplemental(case):
    metrics, xyz, conn, state = fe.load(case)
    check = fe.check_saved_sample(case, metrics, xyz, conn, state)
    gp, gw = np.polynomial.legendre.leggauss(5)
    zp, zw = np.polynomial.legendre.leggauss(3)
    points = np.array([[r, s, t] for r in gp for s in gp for t in zp])
    weights = np.array([a*b*c for a in gw for b in gw for c in zw])
    shape, derivative = fe.shapes(points), fe.derivatives(points)
    ring_E, ring_weights = [], []
    ring_sums = np.zeros(4)
    right_integral = right_area = 0.0
    for nodes in conn:
        coords = xyz[nodes]
        positions = shape @ coords
        radii = np.linalg.norm(positions[:, :2], axis=1)
        selected = (radii >= 3e-6) & (radii <= 3.5e-6)
        if np.any(selected):
            values = fe.evaluate(coords, state[nodes], points, metrics["model"])
            jacobian = np.einsum("pkn,nj->pkj", derivative, coords)
            w = weights * np.linalg.det(jacobian)
            rw = w[selected]
            ring_E.append(values["E"][selected])
            ring_weights.append(rw)
            ring_sums += [
                rw.sum(), np.dot(rw, values["Q"][selected]),
                np.dot(rw, values["E"][selected]),
                np.dot(rw, np.abs(values["n"][selected]-1)),
            ]
        # Integrate the original quadratic contact face, not an unweighted node average.
        for axis in range(3):
            for sign in [-1, 1]:
                face_nodes = np.flatnonzero(fe.NATURAL[:, axis] == sign)
                if not np.all(np.abs(coords[face_nodes, 0]-10e-6) < 1e-15):
                    continue
                axes = [i for i in range(3) if i != axis]
                fp = np.zeros((25, 3))
                fp[:, axis] = sign
                fp[:, axes] = [[r, s] for r in gp for s in gp]
                fw = np.array([a*b for a in gw for b in gw])
                fd = fe.derivatives(fp)
                jac = np.einsum("pkn,nj->pkj", fd, coords)
                area = fw*np.linalg.norm(np.cross(jac[:, axes[0]], jac[:, axes[1]]), axis=1)
                right_area += area.sum()
                right_integral += np.dot(area, (fe.shapes(fp) @ state[nodes])[:, 0])
    if not ring_weights or right_area <= 0:
        raise ValueError("Missing annulus or right-contact samples")
    comparisons = dict(zip(
        ["ring_volume_m3", "ring_Q_mean_W_m3", "ring_E_mean_V_m", "ring_abs_n_relative_mean"],
        [ring_sums[0], *(ring_sums[1:]/ring_sums[0])],
    ))
    errors = {key: audit.rel(value, metrics[key]) for key, value in comparisons.items()}
    # A near-zero LQ carrier perturbation must use the fixed reference density scale.
    errors["ring_abs_n_relative_mean"] = abs(
        comparisons["ring_abs_n_relative_mean"]-metrics["ring_abs_n_relative_mean"])
    if max(errors.values()) >= 2e-7:
        raise ValueError(f"Independent FE ring reconstruction disagrees: {errors}")
    result = {
        "case": case, "saved_sample_relative_error": check,
        "ring_reconstruction_errors": errors,
        "ring_E_p95_V_m": percentile(np.concatenate(ring_E), np.concatenate(ring_weights)),
        "right_mean_ux_m": float(right_integral/right_area),
        "right_contact_area_m2": float(right_area),
        "thermal_resistance_K_W": metrics["theta_max_K"]/metrics["joule_W"],
        "nodes": len(xyz), "elements": len(conn),
        "method": "Original HEX20 fields; 5x5x3 annulus quadrature and 5x5 face quadrature.",
        "source_nodes_sha256": digest(fe.RUNS/case/"nodes.csv"),
        "source_metrics_sha256": digest(fe.RUNS/case/"metrics.txt"),
    }
    (fe.RUNS/case/"supplemental_metrics.json").write_text(
        json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--cases", nargs="+", default=MAIN)
    parser.add_argument("--table", action="store_true")
    args = parser.parse_args()
    fe.RUNS = audit.ROOT = args.run_root.resolve()
    fe.INPUTS = fe.RUNS/"evidence/inputs"
    rows = []
    for case in args.cases:
        result = supplemental(case)
        metrics = audit.read_metrics(case)
        rows.append({**metrics, **{k: v for k, v in result.items() if not isinstance(v, dict)}})
        print(json.dumps(result, allow_nan=False), flush=True)
    if args.table:
        fields = [
            "case", "model", "nodes", "elements", "joule_W", "theta_max_K",
            "ring_theta_mean_K", "thermal_resistance_K_W", "ring_E_mean_V_m",
            "ring_E_p95_V_m", "ring_Q_mean_W_m3", "ring_Q_p95_W_m3",
            "ring_abs_n_relative_mean", "right_mean_ux_m", "minimum_mu_n_ratio",
            "heat_balance_relative", "electron_balance_relative", "joule_quadrature_4_5_relative",
            "E_p95_V_m", "E_above_1e6_volume_fraction",
        ]
        with (fe.RUNS/"main_results.csv").open("w", newline="", encoding="utf-8-sig") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        (fe.RUNS/"main_results.json").write_text(
            json.dumps(rows, indent=2, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main()
