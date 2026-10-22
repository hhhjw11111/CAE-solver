"""Reconstruct currents to locate drift/diffusion cancellation in completed NT states.

All quantities are postprocessed from the original HEX20 fields. J.sigma^-1.J
is a transport-dissipation indicator, NOT a proposed replacement heat equation:
thermoelectric heat/energy fluxes would require a separate constitutive model.
No solved field, original drift heat, or acceptance threshold is modified.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from audit_hex20_boundary_geometry import read_mesh, derivatives
from transfer_hv_initial_guess import shapes
from run_hv_sweep import read_metrics, sha

Q = 1.602176634e-19
KB = 1.380649e-23


def indicators(electric, grad_n, grad_p, temperature, n, p):
    mobility_n = .0205*(temperature/300.)**-2
    mobility_p = .005*(temperature/300.)**-2
    sigma_n, sigma_p = Q*n*mobility_n, Q*p*mobility_p
    thermal_voltage = KB*temperature/Q
    drift_n, drift_p = sigma_n[..., None]*electric, sigma_p[..., None]*electric
    diff_n = sigma_n[..., None]*thermal_voltage[..., None]*grad_n
    diff_p = -sigma_p[..., None]*thermal_voltage[..., None]*grad_p
    current_n, current_p = drift_n+diff_n, drift_p+diff_p
    drift_heat = (sigma_n+sigma_p)*np.sum(electric**2, axis=-1)
    drift_heat_z = (sigma_n+sigma_p)*electric[..., 2]**2
    dissipation = np.sum(current_n**2, axis=-1)/sigma_n + np.sum(current_p**2, axis=-1)/sigma_p
    electrical_work = np.sum((current_n+current_p)*electric, axis=-1)
    cancellation_z = abs(current_n[..., 2])/np.maximum(abs(drift_n[..., 2])+abs(diff_n[..., 2]), 1e-300)
    return drift_heat, drift_heat_z, dissipation, electrical_work, cancellation_z


def diagnose(directory):
    manifest = json.loads((directory/'manifest.json').read_text(encoding='utf-8-sig'))
    metrics = read_metrics(directory/'metrics.txt')
    if (metrics.get('model') != 'NT' or not str(metrics.get('source', '')).startswith('ZNO_300K_CANDIDATE_V2;')
            or metrics.get('heat_source') != 'DRIFT_CONDUCTIVITY_NOT_TOTAL_CURRENT_WORK'
            or metrics.get('final_joule_factor') != 1 or metrics.get('final_voltage_factor') != 1
            or not all(math.isfinite(metrics.get(k, math.nan)) and 0 <= metrics[k] < 1e-9
                       for k in ['maximum_accepted_residual', 'maximum_accepted_update'])):
        raise ValueError('A complete frozen-candidate NT state is required')
    mesh = Path(manifest['mesh'])
    if sha(mesh) != manifest['mesh_sha256']:
        raise ValueError('Mesh hash mismatch')
    coords, conn, _ = read_mesh(mesh)
    data = np.genfromtxt(directory/'nodes.csv', delimiter=',', names=True)
    if len(data) != len(coords) or not np.allclose(np.column_stack([data[k] for k in ['x_m','y_m','z_m']]), coords, atol=1e-16, rtol=0):
        raise ValueError('Node/mesh alignment mismatch')
    fields = np.column_stack([data[k] for k in ['phi_V','eta_n','eta_p','theta_K']])
    if not np.isfinite(fields).all():
        raise ValueError('Nonfinite field')
    gp, gw = np.polynomial.legendre.leggauss(5)
    points = np.array([(x,y,z) for z in gp for y in gp for x in gp])
    weights = np.array([wx*wy*wz for wz in gw for wy in gw for wx in gw])
    shape, dn = shapes(points), derivatives(points)
    totals = np.zeros(6)
    for first in range(0, len(conn), 64):
        cells = conn[first:first+64]
        jac = np.einsum('qka,eaj->eqkj', dn, coords[cells])
        det = np.linalg.det(jac)
        if not np.isfinite(det).all() or np.any(det <= 0):
            raise ValueError('Invalid geometry')
        grad = np.einsum('eqjk,qka->eqja', np.linalg.inv(jac), dn)
        values = np.einsum('qa,eaf->eqf', shape, fields[cells])
        gradients = np.einsum('eqja,eaf->eqjf', grad, fields[cells])
        temperature = 300+values[..., 3]
        n, p = 6e22*np.exp(values[..., 1]), ((2.9e-4)**2/6e22)*np.exp(values[..., 2])
        if not np.isfinite(values).all() or np.any(temperature <= 0) or np.any(n <= 0) or np.any(p <= 0):
            raise ValueError('Invalid physical state')
        heat, heat_z, dissipation, work, cancellation_z = indicators(-gradients[..., 0], gradients[..., 1], gradients[..., 2], temperature, n, p)
        w = det*weights
        quantities = [heat, heat_z, dissipation, work, heat_z*(cancellation_z < .1), heat_z*(cancellation_z < .01)]
        if not all(np.isfinite(value).all() for value in quantities):
            raise ValueError('Nonfinite reconstructed current')
        totals += [np.sum(w*value) for value in quantities]
    error = abs(totals[0]-metrics['joule_W'])/max(abs(totals[0]), 1e-300)
    if error >= 1e-10:
        raise ValueError(f'Original drift-heat reconstruction disagrees: {error}')
    return dict(scope='CONSTITUTIVE_DIAGNOSTIC_NOT_NEW_HEAT_MODEL_OR_ACCEPTANCE',
                voltage_V=metrics['voltage_V'], original_case_acceptance=metrics['case_acceptance'],
                mesh_sha256=sha(mesh), nodes_sha256=sha(directory/'nodes.csv'),
                script_sha256=sha(Path(__file__)), original_drift_heat_W=float(totals[0]),
                z_drift_heat_W=float(totals[1]), transport_dissipation_indicator_W=float(totals[2]),
                total_current_electrical_work_W=float(totals[3]),
                z_drift_heat_fraction_with_current_ratio_below_0_1=float(totals[4]/totals[1]) if totals[1] else None,
                z_drift_heat_fraction_with_current_ratio_below_0_01=float(totals[5]/totals[1]) if totals[1] else None,
                original_heat_reconstruction_relative=error,
                current_ratio_definition='abs(Jn_z)/(abs(Jn_drift_z)+abs(Jn_diffusion_z))',
                limitation='Neither J.E nor the dissipation indicator is substituted into the frozen heat equation; thermoelectric terms and energy fluxes have not been derived or calibrated.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = diagnose(args.directory)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
