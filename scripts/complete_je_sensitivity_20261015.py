"""Finish only the remaining jobs of this explicitly requested in-session study."""
from pathlib import Path
import json,subprocess,sys,time
from run_je_sensitivity_20261015 import ROOT,STUDY,run,sha
from audit_je_sensitivity_power_20261015 import run as audit_power


def wait_completed(directory,allow_failure=False):
    while not (directory/'summary.json').exists():
        time.sleep(20)
    result=json.loads((directory/'summary.json').read_text())
    if result['state']!='SOLVE_FLUX_PASS_MESH_PENDING' and not allow_failure:raise RuntimeError((directory,result))
    return directory


if __name__=='__main__':
    # Continue the angular study while at most one coarse coefficient job runs.
    directory=wait_completed(STUDY/'A24_S1_V060')
    if not (directory/'power_audit/decomposition.json').exists():audit_power(directory)
    mesh=STUDY/'mesh_A24/H1_XY0_Z2.mesh'
    directory=STUDY/'A24_S1_V080'
    if not (directory/'summary.json').exists():
        run(8,1.,mesh,STUDY/'transfer_A16_A24_V080/nodes.csv',directory,timeout_seconds=7200)
    wait_completed(directory)
    if not (directory/'power_audit/decomposition.json').exists():audit_power(directory)
    # Prepare both fine-mesh initial guesses before waiting for coarse jobs.
    for v in [6,8]:
        source=STUDY/f'A24_S1_V{v*10:03d}';wait_completed(source)
        manifest=json.loads((source/'manifest.json').read_text())
        target_mesh=STUDY/'mesh_A32/H1_XY0_Z2.mesh'
        mapping=STUDY/f'transfer_A24_A32_V{v*10:03d}'
        if not (mapping/'transfer_report.json').exists():
            subprocess.run([sys.executable,str(ROOT/'scripts/transfer_hv_initial_guess.py'),
                            '--source-mesh',manifest['mesh'],'--source-nodes',str(source/'nodes.csv'),
                            '--target-mesh',str(target_mesh),'--output',str(mapping)],check=True)
        report=json.loads((mapping/'transfer_report.json').read_text())
        assert report['source_nodes_sha256']==sha(source/'nodes.csv') and report['target_mesh_sha256']==sha(target_mesh)
    # Finish coarse jobs before starting either largest mesh. Resource-only
    # thread increase uses the same executable, equations and numerical gates.
    for name in ['A16_S2_V060','A16_S2_V080']:
        directory=wait_completed(STUDY/name,allow_failure=True)
        if json.loads((directory/'summary.json').read_text())['state']=='SOLVE_FLUX_PASS_MESH_PENDING':
            print('AUDIT '+name,flush=True)
            if not (directory/'power_audit/decomposition.json').exists():audit_power(directory)
        else:
            print('PRESERVED FAILED COEFFICIENT ATTEMPT '+name,flush=True)
    for v in [6,8]:
        target_mesh=STUDY/'mesh_A32/H1_XY0_Z2.mesh'
        mapping=STUDY/f'transfer_A24_A32_V{v*10:03d}'
        directory=STUDY/f'A32_S1_V{v*10:03d}'
        if not (directory/'summary.json').exists():run(v,1.,target_mesh,mapping/'nodes.csv',directory,threads=8,timeout_seconds=7200)
        wait_completed(directory)
        if not (directory/'power_audit/decomposition.json').exists():audit_power(directory)
    subprocess.run([sys.executable,str(ROOT/'scripts/summarize_je_sensitivity_20261015.py')],check=True)
