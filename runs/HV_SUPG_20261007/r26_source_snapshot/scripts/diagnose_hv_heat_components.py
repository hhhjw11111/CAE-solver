"""Explain the frozen drift heat without changing the PDE or its acceptance gates.

This is a postprocessing diagnostic for complete ZNO_F1_NT fields. A difference
between drift heat and terminal electrical power is reported, not redefined as
a discretization error. Assessing the full energy model requires a separate
constitutive derivation; this script does not propose a replacement heat source.
"""
from pathlib import Path
import argparse, csv, json, math
from run_hv_sweep import read_metrics, sha


def diagnose(directory):
    m = read_metrics(directory / 'metrics.txt')
    convergence = [m.get(k) for k in ['maximum_accepted_residual', 'maximum_accepted_update']]
    converged = all(isinstance(v, (int, float)) and math.isfinite(v) and 0 <= v < 1e-9
                    for v in convergence)
    if (m.get('model') != 'NT' or
        m.get('heat_source') != 'DRIFT_CONDUCTIVITY_NOT_TOTAL_CURRENT_WORK' or
        not str(m.get('source', '')).startswith('ZNO_300K_CANDIDATE_V2;') or
        m.get('final_voltage_factor') != 1 or m.get('final_joule_factor') != 1 or not converged):
        raise ValueError('Requires a complete fixed-voltage, full-heat candidate NT state')
    path = directory / 'field_samples.csv'
    volume = heat = recomputed = hot_volume = 0.0
    maximum_field = high_field_heat = high_field_volume = 0.0
    directional = [0.0, 0.0, 0.0]
    t_min, t_max = math.inf, -math.inf
    with path.open(newline='', encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            w, q, theta = (float(row[k]) for k in ['weight_m3', 'Q_W_m3', 'theta_K'])
            e = [float(row[k]) for k in ['Ex_V_m', 'Ey_V_m', 'Ez_V_m']]
            nr, pr = float(row['n_over_n0']), float(row['p_over_p0'])
            temperature = 300.0 + theta
            if (not all(math.isfinite(v) for v in [w, q, temperature, nr, pr, *e]) or
                w <= 0 or temperature <= 0 or min(nr, pr) <= 0):
                raise ValueError('Invalid integration-point field')
            magnitude = math.sqrt(sum(v * v for v in e))
            maximum_field = max(maximum_field, magnitude)
            if magnitude > 1e6:
                high_field_heat += w * q
                high_field_volume += w
            # These constants are the frozen candidate inputs, not fitted values.
            sigma = 1.602176634e-19 * (6e22 * nr * .0205 + (2.9e-4)**2 / 6e22 * pr * .005)
            sigma *= (temperature / 300.0)**(-2)
            components = [w * sigma * v * v for v in e]
            for i in range(3):
                directional[i] += components[i]
            recomputed += sum(components)
            heat += w * q
            volume += w
            hot_volume += w if temperature > 330.0 else 0.0
            t_min, t_max = min(t_min, temperature), max(t_max, temperature)
    if volume <= 0:
        raise ValueError('No quadrature data')
    denom = max(abs(heat), 1e-30)
    discrepancies = dict(
        material_heat_relative=abs(recomputed - heat) / denom,
        volume_integral_relative=abs(heat - m['joule_W']) / denom)
    if max(discrepancies.values()) >= 1e-10:
        raise ValueError(f'Independent drift-heat reconstruction mismatch: {discrepancies}')
    terminal = m['terminal_power_W']
    return dict(
        scope='POSTPROCESSING_DIAGNOSTIC_NOT_MODEL_OR_MESH_ACCEPTANCE',
        voltage_V=m['voltage_V'],
        diagnostic_script_sha256=sha(Path(__file__)),
        field_samples_sha256=sha(path),
        original_case_acceptance=m.get('case_acceptance'),
        heat_source_unchanged=True,
        drift_heat_W=heat, drift_heat_components_W=dict(zip('xyz', directional)),
        drift_heat_fractions=dict(zip('xyz', [v / heat if heat else 0.0 for v in directional])),
        terminal_power_W=terminal, drift_minus_terminal_W=heat - terminal,
        drift_over_terminal=heat / terminal if terminal else None,
        quadrature_temperature_min_K=t_min, quadrature_temperature_max_K=t_max,
        volume_above_330K_fraction=hot_volume / volume,
        maximum_sampled_electric_field_V_m=maximum_field,
        volume_above_1e6_V_m_fraction=high_field_volume / volume,
        drift_heat_above_1e6_V_m_fraction=high_field_heat / heat if heat else 0.0,
        reconstruction=discrepancies,
        limitation='Power definitions differ in the frozen model. This gap alone is not a flux-gate failure or a proof of physical instability.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', type=Path)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    result = diagnose(a.directory)
    output = a.output or a.directory / 'heat_components_diagnostic.json'
    with output.open('x', encoding='utf-8') as f:
        json.dump(result, f, indent=2, allow_nan=False)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
