"""Retry OOC -9 cases with a resource-only 6000 MB budget; retain all attempts."""
from pathlib import Path
import copy,json,shutil,subprocess,sys,time,traceback
from run_je_fine_cross_20261016 import STUDY,ROOT,SCRIPTS,save,sha,accepted
from run_je_sensitivity_resource_20261016 import run

RETRY=STUDY/'RESOURCE_RETRY'

def initialize():
    RETRY.mkdir(exist_ok=False)
    original=json.loads((STUDY/'contract.json').read_text(encoding='utf-8'))
    merged=copy.deepcopy(original);retry_jobs=[]
    for j in merged['jobs']:
        if j['kind'] not in ['radial','thickness']:continue
        j['directory']=str(RETRY/j['name'])
        if max(j['counts'][1:])==18:
            j['source']=str(RETRY/Path(j['source']).name)
            j['mapping']=str(RETRY/Path(j['mapping']).name)
        retry_jobs.append(j)
    merged['resource_retry']=True
    merged['resource_retry_details']={'OOC_MB':6000,'threads':8,'scope':'directional cases only',
        'reason':'PARDISO -9: minimum internal OOC memory exceeded the 4000 MB cap',
        'official_reference':'https://www.intel.com/content/www/us/en/developer/articles/training/how-to-use-ooc-pardiso.html',
        'source_model_heat_mesh_tolerances_unchanged':True,'prior_failed_attempts_preserved':True}
    save(RETRY/'merged_contract.json',merged)
    save(RETRY/'contract.json',dict(jobs=retry_jobs,resource_retry=merged['resource_retry_details']))
    snapshots=RETRY/'script_snapshots';snapshots.mkdir()
    for name in ['run_je_directional_retry_20261016.py','run_je_sensitivity_resource_20261016.py','summarize_je_fine_cross_20261016.py','collect_je_resource_retry_20261016.ps1']:
        shutil.copy2(SCRIPTS/name,snapshots/name)
    save(RETRY/'script_hashes.json',{str(SCRIPTS/p.name):sha(p) for p in snapshots.iterdir()})
    print('RESOURCE RETRY PREPARED; waiting for active fine-grid SUPG queue and its evidence collection',flush=True)

def main():
    if not RETRY.exists():initialize()
    # No concurrent native solves. The first queue remains untouched and is collected first.
    while not (STUDY/'execution_complete.json').exists() or not (STUDY/'analysis/artifact_hashes.json').exists():time.sleep(10)
    contract=json.loads((RETRY/'contract.json').read_text(encoding='utf-8'))
    for j in contract['jobs']:
        d=Path(j['directory']);source=Path(j['source']);mesh=Path(j['mesh']);mapping=Path(j['mapping'])
        if (d/'summary.json').exists():continue
        if not accepted(source):
            save(RETRY/(j['name']+'_dependency_failure.json'),dict(state='DEPENDENCY_NOT_ACCEPTED',source=str(source)))
            continue
        if shutil.disk_usage(RETRY).free<16*1024**3:raise RuntimeError('Insufficient disk headroom; evidence retained')
        if not (mapping/'transfer_report.json').exists():
            m=json.loads((source/'manifest.json').read_text())
            subprocess.run([sys.executable,str(SCRIPTS/'transfer_hv_initial_guess.py'),'--source-mesh',m['mesh'],'--source-nodes',str(source/'nodes.csv'),'--target-mesh',str(mesh),'--output',str(mapping)],check=True)
        transfer=json.loads((mapping/'transfer_report.json').read_text())
        assert transfer['source_nodes_sha256']==sha(source/'nodes.csv') and transfer['target_mesh_sha256']==sha(mesh)
        assert transfer['boundary_extension_nodes']==0
        run(j['v'],1.,mesh,mapping/'nodes.csv',d,threads=8,timeout_seconds=7200,ooc_core_mb=6000)
        if accepted(d):
            subprocess.run([sys.executable,str(SCRIPTS/'audit_je_sensitivity_power_20261015.py'),'--runs',str(d)],check=True)
        print('RESOURCE RETRY POINT COMPLETE '+j['name'],flush=True)
    prior=json.loads((STUDY/'preserved_hashes.json').read_text(encoding='utf-8'))
    assert all(sha(Path(p))==h for p,h in prior.items())
    save(RETRY/'execution_complete.json',dict(planned_retry_attempts_finished=True,original_evidence_unchanged=True))
    print('ALL RESOURCE RETRIES FINISHED',flush=True)

if __name__=='__main__':
    try:main()
    except Exception:traceback.print_exc();raise
