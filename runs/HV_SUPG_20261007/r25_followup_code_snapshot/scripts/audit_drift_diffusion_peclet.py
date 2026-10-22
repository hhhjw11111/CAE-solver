#!/usr/bin/env python3
"""Audit element voltage-drop Peclet indicators for drift--diffusion results.

This is a discretization diagnostic, not an a-posteriori error estimator.  For
each HEX20 element it computes

    Pe_e = (max(phi_a) - min(phi_a)) / (2 k_B T_e / q),

using all element nodes and the nodal mean temperature.  Pe_e > 1 flags an
element whose voltage change exceeds two local thermal voltages and therefore
needs mesh refinement or an exponential-fitting/stabilized transport scheme.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

BOLTZMANN_CONSTANT = 1.380649e-23
ELEMENTARY_CHARGE = 1.602176634e-19


def percentile(sorted_values: list[float], fraction: float) -> float:
    if not sorted_values:
        raise ValueError("Cannot take a percentile of an empty sequence.")
    index = int(fraction * (len(sorted_values) - 1))
    return sorted_values[index]


def read_elements(mesh_path: Path) -> tuple[int, list[list[int]]]:
    with mesh_path.open("r", encoding="utf-8") as stream:
        first = stream.readline().split()
        if len(first) < 2:
            raise ValueError("Mesh header must contain node and element counts.")
        node_count, element_count = int(first[0]), int(first[1])
        for expected_id in range(1, node_count + 1):
            tokens = stream.readline().split()
            if not tokens or int(tokens[0]) != expected_id:
                raise ValueError("Mesh nodes must be contiguous and ordered.")
        elements: list[list[int]] = []
        for expected_id in range(1, element_count + 1):
            tokens = stream.readline().split()
            if len(tokens) != 22 or int(tokens[0]) != expected_id:
                raise ValueError("Expected contiguous HEX20 element records.")
            elements.append([int(value) for value in tokens[2:]])
    return node_count, elements


def read_nodal_fields(nodes_path: Path, node_count: int) -> tuple[list[float], list[float]]:
    phi = [math.nan] * (node_count + 1)
    theta = [math.nan] * (node_count + 1)
    with nodes_path.open("r", encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            node_id = int(row["node"])
            if node_id < 1 or node_id > node_count or math.isfinite(phi[node_id]):
                raise ValueError(f"Invalid or duplicate node id {node_id} in nodes CSV.")
            phi[node_id] = float(row["phi_V"])
            theta[node_id] = float(row["theta_K"])
    if not all(math.isfinite(value) for value in phi[1:] + theta[1:]):
        raise ValueError("Nodes CSV is incomplete or contains a non-finite field.")
    return phi, theta


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    parser.add_argument("nodes", type=Path)
    parser.add_argument("--reference-temperature", type=float, default=300.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.reference_temperature <= 0.0:
        raise ValueError("Reference temperature must be positive.")

    node_count, elements = read_elements(args.mesh)
    phi, theta = read_nodal_fields(args.nodes, node_count)
    peclet: list[float] = []
    drops: list[float] = []
    for connectivity in elements:
        element_phi = [phi[node_id] for node_id in connectivity]
        mean_temperature = args.reference_temperature + sum(
            theta[node_id] for node_id in connectivity
        ) / len(connectivity)
        if mean_temperature <= 0.0:
            raise ValueError("An element has non-positive mean absolute temperature.")
        thermal_voltage = BOLTZMANN_CONSTANT * mean_temperature / ELEMENTARY_CHARGE
        voltage_drop = max(element_phi) - min(element_phi)
        drops.append(voltage_drop)
        peclet.append(voltage_drop / (2.0 * thermal_voltage))

    sorted_pe = sorted(peclet)
    sorted_drop = sorted(drops)
    count_over_one = sum(value > 1.0 for value in peclet)
    count_over_two = sum(value > 2.0 for value in peclet)
    report = {
        "schema": "TPSolverV2-drift-diffusion-peclet-audit-v1",
        "mesh": str(args.mesh.resolve()),
        "nodes": str(args.nodes.resolve()),
        "definition": "Pe_e=(max(phi_node)-min(phi_node))/(2*k_B*T_mean/q)",
        "interpretation": "diagnostic_only_not_an_error_estimator",
        "node_count": node_count,
        "element_count": len(elements),
        "voltage_drop_V": {
            "median": percentile(sorted_drop, 0.5),
            "p95": percentile(sorted_drop, 0.95),
            "maximum": sorted_drop[-1],
        },
        "peclet": {
            "median": percentile(sorted_pe, 0.5),
            "p95": percentile(sorted_pe, 0.95),
            "maximum": sorted_pe[-1],
            "elements_over_1": count_over_one,
            "fraction_over_1": count_over_one / len(elements),
            "elements_over_2": count_over_two,
            "fraction_over_2": count_over_two / len(elements),
        },
        "gate": "FAIL_STANDARD_GALERKIN_HIGH_FIELD" if count_over_one else "PASS_PE_LE_1",
    }
    text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
