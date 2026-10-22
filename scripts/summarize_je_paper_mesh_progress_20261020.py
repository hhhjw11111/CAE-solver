"""New immutable progress snapshot; incomplete points never become accepted data."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import run_je_reconstruction_mesh_20261018 as r


def collect(study):
    rows, details, pending = [], [], []
    for run in sorted(study.glob('A*_R*_Z*_S*_V*')):
        if not (run / 'completion.json').exists():
            pending.append(dict(run=str(run), status='INCOMPLETE_NOT_USED'))
            continue
        manifest = r.accepted(run, require_new_completion=True)
        _, m = r.classify(run, r.read(run / 'execution.json')['exit_code'])
        comparison = r.read(run / 'mapped_field_comparison.json')
        for path, digest in comparison['hashes'].items():
            if r.sha(Path(path)) != digest:
                raise ValueError('Changed comparison evidence: ' + path)
        power = r.read(run / 'power_audit/decomposition.json')
        provenance = r.read(run / 'power_audit/provenance.json')
        assert provenance['nodes_sha256'] == r.sha(run / 'nodes.csv')
        assert provenance['mesh_sha256'] == manifest['mesh_sha256']
        for path, digest in provenance['scripts'].items():
            assert r.sha(Path(path)) == digest
        completion = r.read(run / 'completion.json')
        source = Path(completion['source_run'])
        r.accepted(source)
        _, sm = r.classify(source, r.read(source / 'execution.json')['exit_code'])
        mapping = Path(completion.get('mapping_directory', study / 'transfers' / run.name))
        a = np.genfromtxt(mapping / 'nodes.csv', delimiter=',', names=True)
        b = np.genfromtxt(run / 'nodes.csv', delimiter=',', names=True)
        assert len(a) == len(b)
        worst = {}
        for key in ('eta_n', 'eta_p'):
            metric = comparison['comparisons'][key + '_density_perturbation']
            delta = abs(np.expm1(a[key]) - np.expm1(b[key])) / metric['normalization']
            index = int(np.argmax(delta))
            assert abs(float(delta[index]) - metric['maximum_relative']) < 1e-15
            worst[key] = dict(node=int(b['node'][index]),
                xyz_um=[float(b[k][index] * 1e6) for k in ('x_m', 'y_m', 'z_m')],
                source_mapped_perturbation=float(np.expm1(a[key][index])),
                target_solved_perturbation=float(np.expm1(b[key][index])),
                maximum_relative=float(delta[index]), normalization=metric['normalization'],
                passed=metric['passed'], scope='Location only; no nodes removed from original audit')
        row = dict(run=run.name, voltage_V=m['voltage_V'], supg_scale=m['supg_scale'],
            ooc_core_mb=int(manifest['runtime_env']['MKL_PARDISO_OOC_MAX_CORE_SIZE']),
            theta_max_K=m['theta_max_K'], joule_mW=m['joule_W'] * 1e3,
            total_VM_sample_max_MPa=m['stress_von_mises_max_Pa'] * 1e-6,
            total_VM_peak_change_vs_source_percent=100 * abs(m['stress_von_mises_max_Pa']-sm['stress_von_mises_max_Pa']) / max(abs(m['stress_von_mises_max_Pa']),abs(sm['stress_von_mises_max_Pa'])),
            power_gap_all_boundary_percent=100 * power['gap_relative'],
            power_gap_electrode_percent=100 * power['electrode_gap_relative'],
            SUPG_fraction_of_gap_percent=100 * power['SUPG_fraction_of_gap'],
            n_local_change_percent=100 * worst['eta_n']['maximum_relative'],
            p_local_change_percent=100 * worst['eta_p']['maximum_relative'],
            local_field_passed=comparison['passed'], native_seconds=r.read(run/'execution.json')['seconds'],
            status='SOLVE_FLUX_PASS_SPATIAL_SUPG_PENDING')
        rows.append(row)
        files = [run/name for name in ('completion.json','execution.json','manifest.json','metrics.txt','mapped_field_comparison.json','power_audit/decomposition.json')]
        details.append(dict(**row,source=str(source),carrier_worst_locations=worst,
            evidence_hashes={str(p):r.sha(p) for p in files}))
    return rows, details, pending


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study',type=Path,default=r.STUDY/'stage4_paper_mesh')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    rows, details, pending=collect(args.study.resolve())
    args.output.mkdir(parents=True,exist_ok=False)
    r.save(args.output/'snapshot.json',dict(completed=details,incomplete=pending,
        production_acceptance=False,thresholds_changed=False,
        limitation='Completed single-point audits only. Peak stress changes are diagnostics; no stress-field convergence certificate. Incomplete points and old FAIL records are not promoted.',
        script_sha256=r.sha(Path(__file__))))
    if rows:
        with (args.output/'completed_points.csv').open('x',encoding='utf-8-sig',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(completed=len(rows),incomplete=len(pending),production_acceptance=False)))


if __name__=='__main__':main()
