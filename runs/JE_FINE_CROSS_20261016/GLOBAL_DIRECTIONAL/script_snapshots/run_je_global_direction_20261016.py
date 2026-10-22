"""Independent full-direction JE mesh studies; original model and gates retained."""
from pathlib import Path
import copy,json,shutil,subprocess,sys,time,traceback
from run_je_fine_cross_20261016 import STUDY,ROOT,SCRIPTS,save,sha,accepted
from run_je_sensitivity_resource_20261016 import run

GLOBAL=STUDY/'GLOBAL_DIRECTIONAL'

def initialize():
    assert (GLOBAL/'geometry_validation.json').exists()
    assert json.loads((GLOBAL/'geometry_validation.json').read_text())['baseline_byte_identical']
    merged=copy.deepcopy(json.loads((STUDY/'contract.json').read_text(encoding='utf-8')))
    directional=[]
    for j in merged['jobs']:
        if j['kind'] not in ['radial','thickness']:continue
        n=max(j['counts'][1:]);j['directory']=str(GLOBAL/j['name'])
        j['mesh']=str(GLOBAL/f'mesh_{j["kind"]}_{n}/H1.mesh')
        j['mapping']=str(GLOBAL/f'transfer_{j["name"]}')
        if n==18:j['source']=str(GLOBAL/Path(j['source']).name)
        j.pop('xy',None);j.pop('z',None);directional.append(j)
    merged['full_directional_refinement']=True
    merged['resource_retry']=False
    merged['design_revision']={
        'scope':'Entire radial and thickness logical coordinates refined independently; annulus measurement boundary retained',
        'reason':'Prior 14/16/18 generator refined boundary layers while leaving radial bulk fixed',
        'retired_queue':str(STUDY/'RESOURCE_RETRY/queue_retired_before_solve.json'),
        'initial_four_OOC_failures_preserved':True,'directional_OOC_MB':6000,'threads':8,
        'heat_source_equations_executable_tolerances_unchanged':True,
        'mesh_design_is_not_acceptance':True}
    merged['limitations']=['Independent directions are not a jointly refined mesh','No complete angular convergence certificate',
                           'JE current-work branch is not a complete thermodynamic closure']
    save(GLOBAL/'merged_contract.json',merged)
    save(GLOBAL/'contract.json',{'jobs':directional,'design_revision':merged['design_revision']})
    snapshots=GLOBAL/'script_snapshots';snapshots.mkdir()
    names=['run_je_global_direction_20261016.py','generate_je_global_direction_mesh_20261016.py',
           'run_je_sensitivity_resource_20261016.py','summarize_je_fine_cross_20261016.py','collect_je_global_direction_20261016.ps1',
           'generate_f1_main_cases.py','generate_hv_surface_mesh.py','audit_hex20_boundary_geometry.py','transfer_hv_initial_guess.py',
           'audit_je_sensitivity_power_20261015.py','diagnose_je_power_gap_20261015.py','plot_je_8v_clouds_20261013.py',
           'run_hv_sweep.py','audit_hv_stabilization.py','audit_hv_mesh.py','audit_hv_field_shape.py','summarize_je_sensitivity_20261015.py',
           'run_hv_regression.py','hv_runtime.py','run_je_fine_cross_20261016.py']
    for name in names:shutil.copy2(SCRIPTS/name,snapshots/name)
    save(GLOBAL/'script_hashes.json',{str(SCRIPTS/p.name):sha(p) for p in snapshots.iterdir()})

def prepare_mapping(j):
    source=Path(j['source']);mesh=Path(j['mesh']);mapping=Path(j['mapping'])
    assert accepted(source)
    if not (mapping/'transfer_report.json').exists():
        m=json.loads((source/'manifest.json').read_text())
        subprocess.run([sys.executable,str(SCRIPTS/'transfer_hv_initial_guess.py'),'--source-mesh',m['mesh'],
                        '--source-nodes',str(source/'nodes.csv'),'--target-mesh',str(mesh),'--output',str(mapping)],check=True)
    t=json.loads((mapping/'transfer_report.json').read_text())
    assert t['source_nodes_sha256']==sha(source/'nodes.csv') and t['target_mesh_sha256']==sha(mesh)
    assert t['boundary_extension_nodes']==0

def main():
    if not (GLOBAL/'contract.json').exists():initialize()
    jobs=json.loads((GLOBAL/'contract.json').read_text(encoding='utf-8'))['jobs']
    # Mapping is a small preprocessing task; no second native solver is launched.
    for j in jobs:
        if max(j['counts'][1:])==16:prepare_mapping(j)
    if not (GLOBAL/'preparation_complete.json').exists():save(GLOBAL/'preparation_complete.json',{'first_level_mappings_complete':True})
    print('GLOBAL DIRECTION PREPARATION COMPLETE; waiting for the current SUPG queue',flush=True)
    while not (STUDY/'execution_complete.json').exists() or not (STUDY/'analysis/artifact_hashes.json').exists():time.sleep(10)
    for j in jobs:
        d=Path(j['directory']);source=Path(j['source']);mesh=Path(j['mesh']);mapping=Path(j['mapping'])
        if (d/'summary.json').exists():continue
        if not accepted(source):
            save(GLOBAL/(j['name']+'_dependency_failure.json'),{'state':'DEPENDENCY_NOT_ACCEPTED','source':str(source)})
            continue
        if shutil.disk_usage(GLOBAL).free<16*1024**3:raise RuntimeError('Insufficient disk headroom; preserve evidence and stop')
        prepare_mapping(j)
        run(j['v'],1.,mesh,mapping/'nodes.csv',d,threads=8,timeout_seconds=7200,ooc_core_mb=6000)
        if accepted(d):
            subprocess.run([sys.executable,str(SCRIPTS/'audit_je_sensitivity_power_20261015.py'),'--runs',str(d)],check=True)
        print('GLOBAL DIRECTION POINT COMPLETE '+j['name'],flush=True)
    preserved=json.loads((STUDY/'preserved_hashes.json').read_text(encoding='utf-8'))
    assert all(sha(Path(p))==h for p,h in preserved.items())
    save(GLOBAL/'execution_complete.json',{'planned_directional_attempts_finished':True,'original_evidence_unchanged':True})
    print('ALL GLOBAL DIRECTION ATTEMPTS FINISHED',flush=True)

if __name__=='__main__':
    try:main()
    except Exception:traceback.print_exc();raise
