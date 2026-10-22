"""Fixed-physics fine-mesh coefficient and independent directional studies."""
from pathlib import Path
import json, shutil, subprocess, sys, traceback
from run_je_sensitivity_20261015 import ROOT, EXE, sha, run

OLD=ROOT/'runs/JE_SENSITIVITY_20261015'
STUDY=ROOT/'runs/JE_FINE_CROSS_20261016'
SCRIPTS=ROOT/'scripts'

def save(path,value):
    with path.open('x',encoding='utf-8') as f:json.dump(value,f,ensure_ascii=False,indent=2,allow_nan=False)

def accepted(d):
    return (d/'summary.json').exists() and json.loads((d/'summary.json').read_text())['state']=='SOLVE_FLUX_PASS_MESH_PENDING'

def base(v):return OLD/f'A32_S1_V{v*10:03d}'

def jobs():
    result=[]
    for s in [.5,2.]:
        for v in [6,8]:
            result.append(dict(name=f'A32_R14_Z14_S{s:g}_V{v*10:03d}',v=v,scale=s,mesh=str(OLD/'mesh_A32/H1_XY0_Z2.mesh'),source=str(base(v)),mapping=None,counts=[32,14,14],kind='SUPG'))
    for kind,levels in [('radial',[(1,2,16),(2,2,18)]),('thickness',[(0,3,16),(0,4,18)])]:
        for i,(xy,z,n) in enumerate(levels):
            for v in [6,8]:
                counts=[32,n,14] if kind=='radial' else [32,14,n]
                name=f'A32_R{counts[1]}_Z{counts[2]}_S1_V{v*10:03d}'
                previous=base(v) if i==0 else STUDY/f'A32_R{16 if kind=="radial" else 14}_Z{16 if kind=="thickness" else 14}_S1_V{v*10:03d}'
                result.append(dict(name=name,v=v,scale=1.,mesh=str(STUDY/f'mesh_{kind}_{n}/H1_XY{xy}_Z{z}.mesh'),source=str(previous),mapping=str(STUDY/f'transfer_{name}'),counts=counts,kind=kind,xy=xy,z=z))
    # Obtain directional evidence before the potentially expensive upper coefficient.
    return [j for j in result if j['scale']!=2.]+[j for j in result if j['scale']==2.]

def initialize():
    STUDY.mkdir(exist_ok=False)
    manifest=json.loads((EXE.parent/'build_manifest.json').read_text(encoding='utf-8-sig'))
    preserved={str(ROOT/f['path']):f['sha256'] for f in manifest['files'] if f['path'].startswith('src/')}
    preserved[str(EXE)]=sha(EXE)
    # Freeze all previous audit evidence, plus the two referenced baseline fields.
    inventory=json.loads((OLD/'analysis/evidence_hashes.json').read_text(encoding='utf-8'))
    preserved.update({str(OLD/k):v for k,v in inventory.items()})
    preserved[str(OLD/'analysis/evidence_hashes.json')]=sha(OLD/'analysis/evidence_hashes.json')
    for v in [6,8]:
        assert accepted(base(v))
        for name in ['nodes.csv','input.nml','manifest.json','execution.json','metrics.txt','summary.json']:
            preserved[str(base(v)/name)]=sha(base(v)/name)
    assert all(sha(Path(p))==h for p,h in preserved.items()),'Prior evidence/source changed'
    save(STUDY/'preserved_hashes.json',preserved)
    save(STUDY/'contract.json',dict(jobs=jobs(),baseline=[str(base(v)) for v in [6,8]],
        model='JE_CURRENT_WORK_DIAGNOSTIC_ZNO_NT',heat_source_unchanged=True,
        source_and_executable_unchanged=True,threads=8,ooc_core_mb=4000,wall_time_limit_seconds=7200,
        local_carrier_measure='exp(eta)-1; maximum absolute difference / maximum absolute perturbation in either field',
        field_limit=.02,SUPG_metric_limit=.01,mesh_limits_source=str(SCRIPTS/'audit_hv_mesh.py'),
        power_gap='(JE_volume - physical_all_exterior_power)/JE_volume; also report electrode-only',
        no_new_power_gap_acceptance_threshold=True,
        limitations=['Independent directions, not a jointly refined mesh','Fine SUPG uses A32/R14/Z14','JE is a current-work diagnostic, not complete thermodynamic closure']))
    snapshots=STUDY/'script_snapshots';snapshots.mkdir()
    for name in ['run_je_fine_cross_20261016.py','run_je_sensitivity_20261015.py','generate_hv_surface_mesh.py','generate_f1_main_cases.py','transfer_hv_initial_guess.py','audit_je_sensitivity_power_20261015.py','diagnose_je_power_gap_20261015.py','plot_je_8v_clouds_20261013.py','audit_hv_stabilization.py','audit_hv_mesh.py','audit_hv_field_shape.py','summarize_je_sensitivity_20261015.py','run_hv_sweep.py','hv_runtime.py','run_hv_regression.py']:
        shutil.copy2(SCRIPTS/name,snapshots/name)

def main():
    if not STUDY.exists():initialize()
    contract=json.loads((STUDY/'contract.json').read_text(encoding='utf-8'))
    for j in contract['jobs']:
        d=STUDY/j['name'];source=Path(j['source']);mesh=Path(j['mesh'])
        if (d/'summary.json').exists():continue
        if not accepted(source):
            record=STUDY/(j['name']+'_dependency_failure.json')
            if not record.exists():save(record,dict(state='DEPENDENCY_NOT_ACCEPTED',source=str(source)))
            print('DEPENDENCY NOT ACCEPTED '+j['name'],flush=True);continue
        if shutil.disk_usage(STUDY).free<12*1024**3:raise RuntimeError('Less than 12 GiB free; preserve evidence and stop')
        if not mesh.exists():
            subprocess.run([sys.executable,str(SCRIPTS/'generate_hv_surface_mesh.py'),'--xy-level',str(j['xy']),'--z-level',str(j['z']),'--angular-count','32','--output',str(mesh.parent)],check=True)
        mm=json.loads((mesh.parent/'manifest.json').read_text())
        assert list(mm['counts'].values())==j['counts'] and mm['geometry_audit']['minimum_jacobian_3x3x3']>0
        seed=source/'nodes.csv'
        if j['mapping']:
            mapping=Path(j['mapping'])
            if not (mapping/'transfer_report.json').exists():
                sm=json.loads((source/'manifest.json').read_text())
                subprocess.run([sys.executable,str(SCRIPTS/'transfer_hv_initial_guess.py'),'--source-mesh',sm['mesh'],'--source-nodes',str(seed),'--target-mesh',str(mesh),'--output',str(mapping)],check=True)
            report=json.loads((mapping/'transfer_report.json').read_text())
            assert report['source_nodes_sha256']==sha(seed) and report['target_mesh_sha256']==sha(mesh)
            assert report['boundary_extension_nodes']==0,report
            seed=mapping/'nodes.csv'
        run(j['v'],j['scale'],mesh,seed,d,threads=8,timeout_seconds=7200)
        if accepted(d):
            subprocess.run([sys.executable,str(SCRIPTS/'audit_je_sensitivity_power_20261015.py'),'--runs',str(d)],check=True)
        print('STUDY_POINT_COMPLETE '+j['name'],flush=True)
    preserved=json.loads((STUDY/'preserved_hashes.json').read_text(encoding='utf-8'))
    assert all(sha(Path(p))==h for p,h in preserved.items()),'Prior evidence/source changed during study'
    save(STUDY/'execution_complete.json',dict(all_attempts_finished=True,preserved_files=len(preserved),preserved_hashes_match=True))
    print('ALL PLANNED ATTEMPTS FINISHED',flush=True)

if __name__=='__main__':
    try:main()
    except Exception:
        traceback.print_exc();raise
