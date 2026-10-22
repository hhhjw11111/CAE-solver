"""Finite, sequential joint/SUPG batch conditional on four half-equivalence passes.

No production runs. All physical equations, source definitions and acceptance
auditors stay unchanged. Store raw half integrals and explicit reconstructions.
"""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import time
import traceback
import run_je_reconstruction_mesh_20261018 as r
import run_je_half_native_20261020_v2 as native
import run_je_half_equivalence_20261020_v2 as eq
from audit_hv_half_y import audit as equivalence
from audit_hv_mesh import audit as mesh_audit
from audit_hv_stabilization import audit as supg_audit
from audit_hv_field_shape import audit as field_audit
from audit_je_sensitivity_power_20261015 import run as power_audit

OUTPUT=r.STUDY/'stage4_paper_half_joint_supg_20261020'
COUNTS=((32,18,14),(40,20,16),(48,22,18))
MESH_OOC_MB=8000


def checked_half(folder):
    m=r.read(folder/'manifest.json'); ex=r.read(folder/'execution.json')
    state,metrics=r.classify(folder,ex['exit_code'])
    if state!=r.PASS or r.read(folder/'field_audit.json')['passed'] is not True:
        raise RuntimeError('HALF_SOLVE_FLUX_FIELD_FAILED: '+str(folder))
    if m['domain_symmetry']!='Y_HALF' or metrics['domain_symmetry']!='Y_HALF':
        raise RuntimeError('Half/full study mixing prohibited')
    if (ex['input_sha256']!=r.sha(folder/'input.nml') or
        ex['executable_sha256']!=r.EXPECTED_EXE or m['executable_sha256']!=r.EXPECTED_EXE):
        raise RuntimeError('HALF_EXECUTION_IDENTITY_MISMATCH')
    if r.sha(Path(m['mesh']))!=m['mesh_sha256']:
        raise RuntimeError('HALF_MESH_CHANGED')
    if metrics['heat_source']!='E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE':
        raise RuntimeError('Different heat model')
    return m,metrics


def equivalence_gate():
    certificate=r.read(eq.OUTPUT/'acceptance.json')
    if certificate.get('all_four_equivalence_passed') is not True or len(certificate['events'])!=4:
        raise RuntimeError('Four independent half solves required')
    expected={eq.OUTPUT/f'A32_R18_Z14_V{v*10:03d}_{suffix}'
        for v in (6,8) for suffix in ('unperturbed','hole_perturbed')}
    if {Path(e['run']) for e in certificate['events']}!=expected:
        raise RuntimeError('Duplicate or incomplete equivalence matrix')
    identity=r.read(eq.OUTPUT/'runner_manifest.json')
    for path,digest in identity['script_hashes'].items():
        if r.sha(Path(path))!=digest:raise RuntimeError('Equivalence dependency changed: '+path)
    for event in certificate['events']:
        folder=Path(event['run'])
        checked_half(folder)
        evidence=r.read(folder/'equivalence.json')
        for path,digest in evidence['evidence_hashes'].items():
            if r.sha(Path(path))!=digest:raise RuntimeError('Equivalence evidence changed: '+path)
        repeat=equivalence(Path(evidence['full_run']),folder)
        m=r.read(folder/'manifest.json')
        seed=Path(m['seed_nodes']).parent
        perturb=eq.PERTURBATION if folder.name.endswith('_hole_perturbed') else 0.
        eq.seed_check(seed,Path(evidence['full_run']),perturb)
        recovery=eq.verify_recovery(folder,seed,perturb)
        if not repeat['passed'] or not recovery['passed']:
            raise RuntimeError('Recomputed half equivalence failed')
    return dict(path=str(eq.OUTPUT/'acceptance.json'),
        sha256=r.sha(eq.OUTPUT/'acceptance.json'),all_four_recomputed=True)


def coarse(v):
    return eq.OUTPUT/f'A32_R18_Z14_V{v*10:03d}_unperturbed'


def point(counts,v,scale=1.):
    return OUTPUT/r.name(counts,v,scale)


def full_joint1(v):
    root='stage4_paper_mesh_ooc8000_20261020' if v==6 else 'stage4_paper_mesh_resumed_20261020'
    return r.STUDY/root/r.name(COUNTS[1],v)


def restrict(full_mesh,folder,seed=None):
    command=[sys.executable,str(r.ROOT/'scripts/generate_hv_half_y.py'),
        '--mesh',str(full_mesh),'--output',str(folder)]
    if seed:command+=['--seed-nodes',str(seed)]
    subprocess.run(command,check=True)
    m=r.read(folder/'manifest.json')
    if r.sha(folder/'H1_Y_HALF.mesh')!=m['mesh_sha256']:
        raise RuntimeError('Restricted geometry hash mismatch')
    return folder/'H1_Y_HALF.mesh'


def mapping(source,mesh,label):
    m,_=checked_half(source)
    folder=OUTPUT/'transfers'/label
    subprocess.run([sys.executable,str(r.ROOT/'scripts/transfer_hv_initial_guess.py'),
        '--source-mesh',m['mesh'],'--source-nodes',str(source/'nodes.csv'),
        '--target-mesh',str(mesh),'--output',str(folder)],check=True)
    report=r.read(folder/'transfer_report.json')
    if (report['boundary_extension_nodes']!=0 or
        report['source_mesh_sha256']!=m['mesh_sha256'] or
        report['source_nodes_sha256']!=r.sha(source/'nodes.csv') or
        report['target_mesh_sha256']!=r.sha(mesh)):
        raise RuntimeError('Transfer identity/extrapolation gate failed')
    return folder


def physical_power(folder,destination):
    """Existing all-boundary identity plus separately labelled artificial cut."""
    _,met=checked_half(folder)
    if not (folder/'power_audit').exists():power_audit(folder)
    b=r.read(folder/'power_audit/order5.json')
    cut=b['boundary_sets']['y-']['power_W']
    physical=b['all_exterior_power_W']-cut
    reference=met['physical_exterior_terminal_power_W']
    difference=abs(physical-reference)/max(abs(physical),abs(reference),1e-300)
    if difference>=1e-8:
        raise RuntimeError('Independent physical-exterior power does not match native audit')
    # The cut contribution stays measured, including any nonzero discretization
    # leakage. Nothing is subtracted from the heat or forced to balance.
    result=dict(domain_symmetry='Y_HALF',voltage_V=met['voltage_V'],
        raw_half_JE_W=b['JE_W'],raw_half_all_boundary_power_W=b['all_exterior_power_W'],
        raw_half_symmetry_cut_power_W=cut,raw_half_physical_exterior_power_W=physical,
        physical_exterior_gap_W=b['JE_W']-physical,
        physical_exterior_gap_relative=(b['JE_W']-physical)/b['JE_W'],
        reconstructed_full_JE_W=2*b['JE_W'],
        reconstructed_full_physical_exterior_power_W=2*physical,
        independent_vs_native_relative=difference,comparison_limit=1e-8,
        power_gap_acceptance_threshold=None,
        note='Gap is diagnostic; never heat compensation. Native physical/cut flux gates remain required.',
        evidence_hashes={str(folder/p):r.sha(folder/p) for p in
            ('power_audit/order3.json','power_audit/order5.json','power_audit/decomposition.json',
             'power_audit/provenance.json','metrics.txt','nodes.csv')})
    r.save(destination,result)
    return result


def plan():
    return dict(scope='SYMMETRIC_JOINT_MESH_AND_SUPG_DIAGNOSTIC',
        prerequisite=str(eq.OUTPUT/'acceptance.json'),full_equivalent_counts=COUNTS,
        voltages=[6,8],supg_scales=[.5,1.,2.],
        resource_budget_mb=MESH_OOC_MB,threads=r.THREADS,timeout_seconds=r.TIMEOUT,
        production_started=False,production_auto_release=False,
        order='coarse power audit; joint1 6/8; joint2 6/8; same finest mesh SUPG0.5/2 6/8; original audits',
        no_threshold_changes=True,no_heat_source_changes=True)


def execute():
    gate=r.validate_protocol('paper')
    eqgate=equivalence_gate()
    native.preflight()
    # The same bounded budget already used by the completed FULL joint1 solves.
    # Only the resource setting changes from the smaller equivalence run.
    native.OOC_MB=MESH_OOC_MB
    OUTPUT.mkdir(exist_ok=False)
    (OUTPUT/'audits').mkdir()
    (OUTPUT/'meshes').mkdir()
    r.save(OUTPUT/'protocol.json',plan())
    deps=[Path(__file__),Path(native.__file__),Path(eq.__file__),Path(r.__file__)]
    deps += [r.ROOT/'scripts'/n for n in ('audit_hv_mesh.py','audit_hv_stabilization.py',
        'audit_hv_field_shape.py','audit_je_sensitivity_power_20261015.py',
        'diagnose_je_power_gap_20261015.py','generate_hv_half_y.py','transfer_hv_initial_guess.py')]
    r.save(OUTPUT/'runner_manifest.json',dict(stage3_gate=gate,equivalence_gate=eqgate,
        script_hashes={str(p.resolve()):r.sha(p) for p in deps},
        original_gates_retained=True,production_acceptance=False))
    # Additions go to this fresh study, leaving the four equivalence runs intact.
    for v in (6,8):
        coarse_run=coarse(v)
        auditdir=OUTPUT/'coarse_power'/f'V{v*10:03d}'
        auditdir.mkdir(parents=True)
        # Existing audit writes to its run path; make a source-hashed, explicitly
        # labelled audit view instead of adding files to the frozen solve.
        needed=('input.nml','execution.json','manifest.json','summary.json','field_audit.json','metrics.txt','nodes.csv')
        import shutil
        for name in needed:shutil.copy2(coarse_run/name,auditdir/name)
        r.save(auditdir/'audit_view.json',dict(source_run=str(coarse_run),not_a_new_solve=True,
            hashes={name:r.sha(coarse_run/name) for name in needed}))
        physical_power(auditdir,OUTPUT/'audits'/f'coarse_power_V{v*10:03d}.json')
    locals_={}
    # Joint1 uses the full solution as an initial guess only; its separately
    # generated coarse-half mapping is used exclusively for convergence review.
    for v in (6,8):
        full=full_joint1(v); fm=r.accepted(full,require_new_completion=True)
        name=r.name(COUNTS[1],v); folder=OUTPUT/'meshes'/name
        mesh=restrict(Path(fm['mesh']),folder,full/'nodes.csv')
        source=coarse(v); transfer=mapping(source,mesh,name)
        dest=point(COUNTS[1],v)
        native.run_native(dict(name=name,voltage=v,supg_scale=1.,counts=list(COUNTS[1])),
            mesh,folder/'nodes.csv',dest,gate)
        second_equivalence=equivalence(full,dest)
        r.save(dest/'full_half_equivalence.json',second_equivalence)
        if not second_equivalence['passed']:
            raise RuntimeError('Joint1 full/half equivalence failed')
        comparison=field_audit(transfer,source,dest)
        locals_[name]=comparison
        r.save(dest/'field_comparison.json',comparison)
        physical_power(dest,dest/'physical_power.json')
        r.save(dest/'completion.json',dict(power_audit_completed=True,
            field_comparison_passed=comparison['passed'],fully_validated=False))
    fullmesh=Path(r.read(r.STUDY/'stage4_paper_mesh_ooc12000_20261020/A48_R22_Z18_S1_V060/manifest.json')['mesh'])
    mesh=restrict(fullmesh,OUTPUT/'meshes'/'A48_R22_Z18')
    for v in (6,8):
        name=r.name(COUNTS[2],v);source=point(COUNTS[1],v)
        transfer=mapping(source,mesh,name);dest=point(COUNTS[2],v)
        native.run_native(dict(name=name,voltage=v,supg_scale=1.,counts=list(COUNTS[2])),
            mesh,transfer/'nodes.csv',dest,gate)
        comparison=field_audit(transfer,source,dest);locals_[name]=comparison
        r.save(dest/'field_comparison.json',comparison)
        physical_power(dest,dest/'physical_power.json')
        r.save(dest/'completion.json',dict(power_audit_completed=True,
            field_comparison_passed=comparison['passed'],fully_validated=False))
    mesh_results=[mesh_audit([coarse(v),point(COUNTS[1],v),point(COUNTS[2],v)],[v]) for v in (6,8)]
    r.save(OUTPUT/'joint_mesh_audit.json',dict(results=mesh_results,
        all_mesh_pass=all(x['all_mesh_pass'] for x in mesh_results),production_acceptance=False))
    for scale in (.5,2.):
        for v in (6,8):
            source=point(COUNTS[2],v);name=r.name(COUNTS[2],v,scale);dest=point(COUNTS[2],v,scale)
            transfer=mapping(source,mesh,name)
            native.run_native(dict(name=name,voltage=v,supg_scale=scale,counts=list(COUNTS[2])),
                mesh,transfer/'nodes.csv',dest,gate)
            comparison=field_audit(transfer,source,dest);locals_[name]=comparison
            r.save(dest/'field_comparison.json',comparison)
            physical_power(dest,dest/'physical_power.json')
            r.save(dest/'completion.json',dict(power_audit_completed=True,
                field_comparison_passed=comparison['passed'],fully_validated=False))
    supg_results=[supg_audit([point(COUNTS[2],v,s) for s in (.5,1.,2.)]) for v in (6,8)]
    last_names=[r.name(COUNTS[2],v,s) for s in (1.,.5,2.) for v in (6,8)]
    local_pass=all(locals_[n]['passed'] for n in last_names)
    r.save(OUTPUT/'final_audit.json',dict(mesh_results=mesh_results,supg_results=supg_results,
        final_pair_and_supg_local_pass=local_pass,
        all_original_numerical_checks_passed=local_pass and all(x['all_mesh_pass'] for x in mesh_results)
            and all(x['covered_metrics_passed'] for x in supg_results),
        first_pair_local_failures_retained=True,production_started=False,
        limitation='Only symmetric r=3um endpoints; production voltage range, other radii and LQ still require their own checks'))


def wait_for_equivalence(seconds):
    deadline=time.monotonic()+seconds
    while not (eq.OUTPUT/'acceptance.json').exists():
        if (eq.OUTPUT/'stopped.json').exists():
            raise RuntimeError('Equivalence controller stopped; no mesh child started')
        if time.monotonic()>=deadline:
            raise RuntimeError('Equivalence wait timeout; no mesh child started')
        time.sleep(10)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--wait-for-equivalence',type=float,default=0.)
    args=parser.parse_args()
    if not args.execute:print(json.dumps(plan(),indent=2))
    else:
        try:
            if args.wait_for_equivalence:wait_for_equivalence(args.wait_for_equivalence)
            execute()
        except BaseException as exc:
            if OUTPUT.exists() and not (OUTPUT/'stopped.json').exists():
                r.save(OUTPUT/'stopped.json',dict(error=repr(exc),traceback=traceback.format_exc(),
                    production_acceptance=False))
            raise
