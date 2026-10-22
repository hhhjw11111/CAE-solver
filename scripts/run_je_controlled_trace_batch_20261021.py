"""Controlled SUPG metric comparison; only release after all original gates pass.

Same three joint meshes, same executable, JE, R=0, material, tolerances and
SUPG coefficients as the preserved directional experiment.  Switch the existing
trace metric explicitly; never shrink the coefficient sweep or omit a failure.
"""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import traceback
import run_je_reconstruction_mesh_20261018 as r
import run_je_trace_native_20261021 as native
import run_je_half_equivalence_20261020_v2 as eq
from continue_je_half_mesh_supg_20261020 import physical_power, checked_half
from audit_hv_mesh import audit as mesh_audit
from audit_hv_stabilization import audit as supg_audit
from audit_hv_field_shape import audit as field_audit

ROOT = r.STUDY/'stage5_controlled_numerics_20261021'
OUT = ROOT/'trace_joint_supg'
OLD = r.STUDY/'stage4_paper_half_joint_supg_20261020'
COUNTS = ((32,18,14), (40,20,16), (48,22,18))


def prior(level, voltage):
    if level == 0:
        return eq.OUTPUT/f'A32_R18_Z14_V{voltage*10:03d}_unperturbed'
    return OLD/r.name(COUNTS[level], voltage)


def point(level, voltage, scale=1.):
    return OUT/r.name(COUNTS[level], voltage, scale)


def validate_prerequisites():
    gate = r.validate_protocol('paper')
    path = ROOT/'native_trace_verification_attempt02/acceptance.json'
    acceptance = r.read(path)
    if (acceptance.get('passed') is not True or acceptance.get('exit_code') != 0
            or acceptance.get('original_core_unchanged') is not True
            or acceptance.get('original_thresholds_unchanged') is not True
            or acceptance.get('production_executable_sha256') != r.EXPECTED_EXE):
        raise RuntimeError('Native trace prerequisite failed')
    for inventory in ('verification_source_hashes', 'copied_source_hashes'):
        for source, digest in acceptance[inventory].items():
            if r.sha(Path(source)) != digest:
                raise RuntimeError('Verification source changed: '+source)
    identity = r.read(path.parent/'evidence_identity.json')
    if (r.sha(path) != identity['acceptance_sha256']
            or r.sha(path.parent/'trace_check.exe') != acceptance['verification_executable_sha256']
            or r.sha(path.parent/'trace_check.exe') != identity['verification_executable_sha256']
            or r.sha(path.parent/'verification.txt') != identity['verification_log_sha256']):
        raise RuntimeError('Verification executable/log changed')
    log = (path.parent/'verification.txt').read_text(encoding='utf-8')
    for marker in ('verification_id=TRACE_SUPG_ISO_HIGH_PE',
                   'verification_id=TRACE_SUPG_THIN_HIGH_PE', 'trace_controlled_verification=PASS'):
        if log.count(marker) != 1:
            raise RuntimeError('Incomplete trace high-Pe/low-Pe benchmark')
    acceptance = r.read(ROOT/'checked_block_replay/acceptance.json')
    if (acceptance.get('all_equivalent') is not True or acceptance.get('original_limit') != 1e-8
            or len(acceptance.get('events', [])) != 2
            or {x['voltage'] for x in acceptance['events']} != {6,8}
            or not all(x['passed'] for x in acceptance['events'])):
        raise RuntimeError('Incomplete checked-block equivalence gate')
    from audit_hv_half_y import audit as equivalence
    for v in (6,8):
        folder = ROOT/'checked_block_replay'/f'A32_R18_Z14_V{v*10:03d}'
        checked_half(folder)
        full = Path(r.read(folder/'equivalence.json')['full_run'])
        if not equivalence(full, folder)['passed']:
            raise RuntimeError('Recomputed checked-block equivalence failed')
    for v in (6,8):
        for level in (0,1,2):
            checked_half(prior(level,v))
    return gate


def compare(source, target, label):
    sm, _ = checked_half(source)
    tm, _ = checked_half(target)
    output = OUT/'transfers'/label
    subprocess.run([sys.executable, '-X', 'utf8', str(r.ROOT/'scripts/transfer_hv_initial_guess.py'),
        '--source-mesh', sm['mesh'], '--source-nodes', str(source/'nodes.csv'),
        '--target-mesh', tm['mesh'], '--output', str(output)], check=True)
    mapping = r.read(output/'transfer_report.json')
    if mapping['boundary_extension_nodes'] != 0:
        raise RuntimeError('Extrapolation in local field comparison')
    result = field_audit(output, source, target)
    r.save(target/(label+'_local_comparison.json'), result)
    return result


def execute():
    gate = validate_prerequisites()
    native.preflight()
    native.OOC_MB = 8000
    OUT.mkdir(exist_ok=False)
    (OUT/'audits').mkdir()
    protocol = dict(scope='CONTROLLED_TRACE_SUPG_NUMERICAL_VARIANT',
        counts=COUNTS, voltages=[6,8], scales=[.5,1.,2.],
        existing_input_changed=dict(supg_directional_tau=False),
        resource_mode='EXISTING_CHECKED_BLOCK_SOLVER_FULL_MATRIX_LIMIT_1E-11',
        original_model_and_heat_unchanged=True, thresholds_unchanged=True,
        power_gap_is_diagnostic_no_added_threshold=True,
        production_release_requires_all_original_mesh_supg_and_local_checks=True,
        production_started=False)
    r.save(OUT/'protocol.json', protocol)
    scripts = [Path(__file__), Path(native.__file__), r.ROOT/'scripts/run_je_half_native_20261020_v2.py',
        r.ROOT/'scripts/audit_hv_mesh.py', r.ROOT/'scripts/audit_hv_stabilization.py',
        r.ROOT/'scripts/audit_hv_field_shape.py', r.ROOT/'scripts/transfer_hv_initial_guess.py',
        r.ROOT/'scripts/continue_je_half_mesh_supg_20261020.py']
    prerequisites = [ROOT/'native_trace_verification_attempt02/acceptance.json',
                     ROOT/'native_trace_verification_attempt02/evidence_identity.json',
                     ROOT/'native_trace_verification_attempt02/verification.txt',
                     ROOT/'checked_block_replay/acceptance.json']
    r.save(OUT/'runner_manifest.json', dict(gate=gate, executable_sha256=r.EXPECTED_EXE,
        script_hashes={str(p):r.sha(p) for p in scripts},
        prerequisite_hashes={str(p):r.sha(p) for p in prerequisites}))
    comparisons = []
    # Every geometry is already independently generated and audited.  Its old
    # directional solution is an initial guess only; each trace point is solved.
    for level in (0,1,2):
        for v in (6,8):
            source = prior(level,v)
            sm, _ = checked_half(source)
            destination = point(level,v)
            job = dict(name='TRACE_'+destination.name, voltage=v, supg_scale=1., counts=COUNTS[level])
            native.run_native(job, Path(sm['mesh']), source/'nodes.csv', destination, gate)
            physical_power(destination, destination/'physical_power.json')
            if level:
                result = compare(point(level-1,v), destination, 'mesh_final_pair' if level==2 else 'mesh_first_pair')
                if level == 2:
                    comparisons.append(result)
            r.save(destination/'completion.json', dict(native_solve_flux_fields=True,
                independent_power_audit=True, fully_validated=False))
    mesh_results = [mesh_audit([point(i,v) for i in (0,1,2)], [v]) for v in (6,8)]
    r.save(OUT/'audits/joint_mesh.json', dict(results=mesh_results,
        all_mesh_pass=all(x['all_mesh_pass'] for x in mesh_results)))
    for scale in (.5,2.):
        for v in (6,8):
            source = point(2,v)
            sm, _ = checked_half(source)
            destination = point(2,v,scale)
            native.run_native(dict(name='TRACE_'+destination.name, voltage=v, supg_scale=scale,
                counts=COUNTS[2]), Path(sm['mesh']), source/'nodes.csv', destination, gate)
            comparisons.append(compare(source, destination, 'supg_same_mesh'))
            physical_power(destination, destination/'physical_power.json')
            r.save(destination/'completion.json', dict(native_solve_flux_fields=True,
                independent_power_audit=True, fully_validated=False))
    supg_results = [supg_audit([point(2,v,s) for s in (.5,1.,2.)]) for v in (6,8)]
    passed = (all(x['all_mesh_pass'] for x in mesh_results)
        and all(x['covered_metrics_passed'] for x in supg_results)
        and len(comparisons)==6 and all(x['passed'] for x in comparisons))
    for path, digest in r.read(OUT/'runner_manifest.json')['script_hashes'].items():
        if r.sha(Path(path)) != digest:
            raise RuntimeError('Batch dependency changed during execution: '+path)
    r.validate_protocol('paper')
    r.save(OUT/'final_audit.json', dict(scope=protocol['scope'], mesh_results=mesh_results,
        supg_results=supg_results, final_pair_and_supg_local_pass=all(x['passed'] for x in comparisons),
        all_original_numerical_checks_passed=passed, endpoint_numerical_gate_passed=passed,
        production_started=False,
        limitation='Symmetric r=3um NT endpoints only. Other radii, LQ and intermediate voltages require independent evidence.'))
    print(json.dumps(dict(all_original_numerical_checks_passed=passed, output=str(OUT))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    if not args.execute:
        print(json.dumps(dict(output=str(OUT), counts=COUNTS, voltages=[6,8], scales=[.5,1.,2.],
                              gate_changes=False, physics_changes=False)))
    else:
        try:
            execute()
        except BaseException as exc:
            OUT.mkdir(exist_ok=True)
            r.save(OUT/'stopped.json', dict(error=repr(exc), production_started=False))
            traceback.print_exc()
            raise
