"""Resume the controlled trace batch with voltage-specific comparison mappings.

Retain all first-attempt files. Reuse a completed solve only after identity and
original solve/flux/field gates are recomputed; never overwrite evidence.
"""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import traceback
import run_je_controlled_trace_batch_20261021 as b
import run_je_reconstruction_mesh_20261018 as r
from audit_hv_field_shape import audit as field_audit
from audit_hv_mesh import audit as mesh_audit
from audit_hv_stabilization import audit as supg_audit

ROOT, OUT, COUNTS, point = b.ROOT, b.OUT, b.COUNTS, b.point


def compare(source, target, kind):
    sm, _ = b.checked_half(source)
    tm, _ = b.checked_half(target)
    label = kind+'_'+target.name
    output = OUT/'transfers'/label
    if not output.exists():
        subprocess.run([sys.executable,'-X','utf8',str(r.ROOT/'scripts/transfer_hv_initial_guess.py'),
            '--source-mesh',sm['mesh'],'--source-nodes',str(source/'nodes.csv'),
            '--target-mesh',tm['mesh'],'--output',str(output)],check=True)
    mapping = r.read(output/'transfer_report.json')
    if (mapping['boundary_extension_nodes'] != 0 or mapping['source_nodes_sha256'] != r.sha(source/'nodes.csv')
            or mapping['source_mesh_sha256'] != sm['mesh_sha256']
            or mapping['target_mesh_sha256'] != tm['mesh_sha256']):
        raise RuntimeError('Wrong/mutated mapping identity')
    result = field_audit(output,source,target)
    destination = target/(kind+'_local_comparison_v2.json')
    if destination.exists():
        if r.read(destination) != result:
            raise RuntimeError('Changed preexisting comparison')
    else:
        r.save(destination,result)
    return result


def run_point(level,v,scale,gate):
    destination = point(level,v,scale)
    if destination.exists():
        manifest, metrics = b.checked_half(destination)
        if (metrics.get('supg_directional_tau') != 'F' or metrics['voltage_V'] != v
                or metrics['supg_scale'] != scale or manifest.get('counts') != list(COUNTS[level])
                or manifest['runtime_env']['TPS_PARDISO_SPLIT_HOLES'] != '1'):
            raise RuntimeError('Existing point is not this numerical variant')
    else:
        source = b.prior(level,v) if scale == 1. else point(level,v)
        manifest, _ = b.checked_half(source)
        b.native.run_native(dict(name='TRACE_'+destination.name,voltage=v,supg_scale=scale,counts=COUNTS[level]),
            Path(manifest['mesh']),source/'nodes.csv',destination,gate)
    if not (destination/'physical_power.json').exists():
        b.physical_power(destination,destination/'physical_power.json')
    if not (destination/'completion_v2.json').exists():
        r.save(destination/'completion_v2.json',dict(native_solve_flux_fields=True,
            physical_power_integrals_checked=True,trace_decomposition='Separate trace_power review required',
            fully_validated=False))
    print('CONTROLLED_POINT_COMPLETED '+destination.name,flush=True)
    return destination


def execute():
    gate = b.validate_prerequisites()
    b.native.preflight()
    b.native.OOC_MB = 8000
    if not (OUT/'stopped.json').exists():
        raise RuntimeError('First controller must have stopped before resuming')
    previous = r.read(OUT/'stopped.json')
    if ('FileExistsError' not in previous['error']
            and previous['error'] != 'CONTROLLED_STOP_MAPPING_NAMESPACE_REPAIR'):
        raise RuntimeError('Unexpected original controller failure; inspect before resuming')
    files=[Path(__file__),Path(b.__file__),Path(b.native.__file__),
        r.ROOT/'scripts/audit_hv_field_shape.py',r.ROOT/'scripts/audit_hv_mesh.py',
        r.ROOT/'scripts/audit_hv_stabilization.py',r.ROOT/'scripts/transfer_hv_initial_guess.py',
        OUT/'runner_manifest.json',OUT/'stopped.json']
    r.save(OUT/'runner_manifest_v2.json',dict(original_attempt_retained=True,gate=gate,
        script_hashes={str(p):r.sha(p) for p in files},
        only_orchestration_repaired='Include target identity in comparison mapping directories; validated resume',
        physics_and_thresholds_unchanged=True))
    locals_=[]
    for level in (0,1,2):
        for v in (6,8):
            destination=run_point(level,v,1.,gate)
            if level:
                result=compare(point(level-1,v),destination,'mesh_final_pair' if level==2 else 'mesh_first_pair')
                if level==2:
                    locals_.append(result)
    meshes=[mesh_audit([point(i,v) for i in (0,1,2)],[v]) for v in (6,8)]
    r.save(OUT/'audits/joint_mesh_v2.json',dict(results=meshes,all_mesh_pass=all(x['all_mesh_pass'] for x in meshes)))
    for scale in (.5,2.):
        for v in (6,8):
            destination=run_point(2,v,scale,gate)
            locals_.append(compare(point(2,v),destination,'supg_same_mesh'))
    supgs=[supg_audit([point(2,v,s) for s in (.5,1.,2.)]) for v in (6,8)]
    passed=all(x['all_mesh_pass'] for x in meshes) and all(x['covered_metrics_passed'] for x in supgs) and all(x['passed'] for x in locals_)
    for path,digest in r.read(OUT/'runner_manifest_v2.json')['script_hashes'].items():
        if r.sha(Path(path))!=digest:
            raise RuntimeError('Changed batch dependency: '+path)
    r.validate_protocol('paper')
    r.save(OUT/'final_audit_v2.json',dict(mesh_results=meshes,supg_results=supgs,
        final_pair_and_supg_local_pass=all(x['passed'] for x in locals_),
        all_original_numerical_checks_passed=passed,production_started=False,
        first_attempt_preserved=True,physics_and_thresholds_unchanged=True))
    print(json.dumps(dict(passed=passed,report=str(OUT/'final_audit_v2.json'))),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    if args.execute:
        try:
            execute()
        except BaseException as exc:
            r.save(OUT/'stopped_v2.json',dict(error=repr(exc),production_started=False))
            traceback.print_exc()
            raise
    else:
        print(json.dumps(dict(output=str(OUT),mode='explicit validated resume, no overwrite')))
