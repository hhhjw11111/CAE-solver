"""Four independent half solves before permitting a symmetric-domain mesh study.

Original FULL results remain immutable. A restricted seed is never a solution.
Each native solve retains the original nonlinear, physical-flux and field gates.
The fixed 1e-8 full/half equivalence gate is not a mesh-convergence certificate.
"""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import traceback
import numpy as np
import run_je_reconstruction_mesh_20261018 as r
import run_je_half_native_20261020 as native
import audit_je_full_field_y_parity_20261019 as parity
from audit_hv_half_y import audit, LIMIT

OUTPUT = r.STUDY/'stage4_paper_half_equivalence_20261020'
PERTURBATION = 0.001
COUNTS = [32, 18, 14]


def fixed_plan():
    return dict(scope='JE_Y_HALF_EQUIVALENCE_BEFORE_MESH_STUDY',
        voltages=[6, 8], full_counts=COUNTS,
        initial_hole_log_perturbations=[0., PERTURBATION],
        equivalence_limit=LIMIT,
        perturbation_shape='delta_eta_p=0.001*max(0,1-(x/1e-5)^2); zero at contacts',
        expected_nonzero_density_change='exp(0.001)-1 approximately 0.10005 percent',
        required='All four native solves and original full/half comparisons must pass',
        integral_convention='Keep actual half integrals; multiply by two only in explicit full reconstruction',
        full_sources={str(v): str(r.radial18(v)) for v in (6,8)},
        resource_budget_mb=native.OOC_MB, timeout_seconds=r.TIMEOUT,
        production_acceptance=False, mesh_validated=False,
        limitations='Symmetric branch only; excludes asymmetric equilibria and dynamical stability')


def seed_check(seed, full, perturb):
    a=np.genfromtxt(full/'nodes.csv', delimiter=',', names=True)
    b=np.genfromtxt(seed/'nodes.csv', delimiter=',', names=True)
    ids=np.load(seed/'full_node_ids.npy')-1
    expected=perturb*np.maximum(0.,1-(b['x_m']/1e-5)**2)
    delta=b['eta_p']-a['eta_p'][ids]
    err=float(np.max(abs(delta-expected)))
    if err>1e-14 or (perturb and float(np.max(abs(delta)))<perturb*.99):
        raise RuntimeError('Perturbed seed was not actually applied')
    return dict(maximum_delta_eta_p=float(np.max(abs(delta))),
        maximum_formula_error=err, initial_guess_only=True,
        seed_sha256=r.sha(seed/'nodes.csv'))


def verify_recovery(run, seed, perturb):
    a=np.genfromtxt(seed/'nodes.csv', delimiter=',', names=True)
    b=np.genfromtxt(run/'nodes.csv', delimiter=',', names=True)
    change=float(np.max(abs(a['eta_p']-b['eta_p'])))
    passed=not perturb or change>PERTURBATION*.9
    return dict(perturbation=perturb, maximum_actual_eta_p_change=change,
        passed=bool(passed),
        note='Nonzero perturbation must actually be removed; full/half audit separately tests every final field')


def run():
    native.preflight()
    gate=r.validate_protocol('paper')
    failed=r.STUDY/'stage4_paper_mesh_ooc12000_20261020/A48_R22_Z18_S1_V060'
    ex=r.read(failed/'execution.json')
    if ex['stop_reason']!='RESOURCE_STOP_AVAILABLE_RAM_BELOW_512_MIB':
        raise RuntimeError('Resource-fallback prerequisite missing')
    if not r.read(failed/'temporary_factor_cleanup.json')['completed']:
        raise RuntimeError('Prior owned-factor cleanup incomplete')
    OUTPUT.mkdir(exist_ok=False)
    deps=[Path(__file__), Path(native.__file__),Path(r.__file__),Path(parity.__file__),
        r.ROOT/'scripts/audit_hv_half_y.py',r.ROOT/'scripts/generate_hv_half_y.py',
        r.ROOT/'scripts/run_hv_sweep.py',r.ROOT/'scripts/hv_runtime.py']
    r.save(OUTPUT/'protocol.json',fixed_plan())
    r.save(OUTPUT/'runner_manifest.json',dict(stage3_gate=gate,
        script_hashes={str(p.resolve()):r.sha(p) for p in deps},
        resource_failure_evidence={str(failed/p):r.sha(failed/p) for p in
            ('execution.json','resource_usage.json','input.nml','temporary_factor_cleanup.json')},
        thresholds_changed=False, core_sources_changed=False))
    prechecks=[parity.audit(r.radial18(v)) for v in (6,8)]
    r.save(OUTPUT/'full_field_parity.json',dict(results=prechecks,
        passed=all(p['full_field_parity_passed'] for p in prechecks)))
    if not all(p['full_field_parity_passed'] for p in prechecks):
        raise RuntimeError('FULL_FIELD_PARITY_FAILED')
    events=[]
    for perturb in (0.,PERTURBATION):
        for v in (6,8):
            full=r.radial18(v); m=r.accepted(full)
            suffix='unperturbed' if not perturb else 'hole_perturbed'
            name=f'A32_R18_Z14_V{v*10:03d}_{suffix}'
            seed=OUTPUT/'seeds'/name
            command=[sys.executable,str(r.ROOT/'scripts/generate_hv_half_y.py'),
                '--mesh',m['mesh'],'--seed-nodes',str(full/'nodes.csv'),
                '--output',str(seed),'--perturb-hole',str(perturb)]
            subprocess.run(command,check=True)
            r.save(seed/'independent_seed_check.json',seed_check(seed,full,perturb))
            dest=OUTPUT/name
            job=dict(name=name, voltage=v,supg_scale=1.,counts=COUNTS)
            print('START_HALF_EQUIVALENCE '+name,flush=True)
            native.run_native(job,seed/'H1_Y_HALF.mesh',seed/'nodes.csv',dest,gate)
            result=audit(full,dest)
            # Save the independent gate even when it fails. Never relabel FAIL.
            evidence=[dest/p for p in ('nodes.csv','metrics.txt','execution.json','input.nml',
                'manifest.json','field_audit.json')]+[full/'nodes.csv',full/'metrics.txt',seed/'manifest.json']
            result['evidence_hashes']={str(p):r.sha(p) for p in evidence}
            r.save(dest/'equivalence.json',result)
            recovery=verify_recovery(dest,seed,perturb)
            r.save(dest/'recovery.json',recovery)
            events.append(dict(run=str(dest),equivalence_passed=result['passed'],
                perturbation_recovery_passed=recovery['passed']))
            print(json.dumps(events[-1]),flush=True)
            if not result['passed'] or not recovery['passed']:
                r.save(OUTPUT/'failed_equivalence.json',dict(events=events,
                    mesh_study_permitted=False, production_acceptance=False))
                raise RuntimeError('HALF_EQUIVALENCE_FAILED: '+name)
    r.save(OUTPUT/'acceptance.json',dict(all_four_equivalence_passed=True,events=events,
        threshold=LIMIT,symmetric_mesh_study_permitted=True,production_acceptance=False,
        mesh_validated=False,limitations=fixed_plan()['limitations']))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    if not args.execute:
        print(json.dumps(fixed_plan(),indent=2))
    else:
        try:run()
        except BaseException as exc:
            if OUTPUT.exists() and not (OUTPUT/'stopped.json').exists():
                r.save(OUTPUT/'stopped.json',dict(error=repr(exc),traceback=traceback.format_exc(),
                    production_acceptance=False))
            raise
