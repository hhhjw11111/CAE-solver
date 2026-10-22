"""No-native verification of the lightweight JE input preparer."""
from pathlib import Path
import sys
import json
import hashlib
import subprocess
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
sys.path.insert(0,str(ROOT/'scripts'))
import prepare_je_paper_input_20261018 as prep

def save(path,value):
    with path.open('x',encoding='utf-8') as stream:
        json.dump(value,stream,indent=2,ensure_ascii=False,allow_nan=False)
        stream.write('\n')

def number(text,key):
    return float(re.search(r'(?<![A-Za-z0-9_])'+key+r'\s*=\s*([^,\n/]+)',text)[1])

def expect_reject(call):
    try:call()
    except ValueError as exc:return dict(rejected=True,reason=str(exc))
    raise AssertionError('Expected rejection')

def main():
    sources=[prep.TEMPLATE,ROOT/'scripts/prepare_je_paper_input_20261018.py',
             ROOT/'src/app/f1_case_runner_m.f90',ROOT/'src/nonlinear/continuation_solver_m.f90',
             ROOT/'scripts/generate_je_reconstruction_mesh_20261018.py',ROOT/'scripts/run_hv_point.py']
    frozen={str(p):prep.sha(p) for p in sources}
    save(HERE/'protocol.json',dict(scope='Input-only software verification; no native solve or production acceptance',
        sources=frozen,runner_sha256=prep.sha(Path(__file__)),
        checks=['all nine pure input conditions','6.5/7.5 target and restart factors','three radius manifests',
                'exact coordinate and voltage identity','old 3V drift seed opt-in is initial-guess-only',
                'changed template and wrong geometry rejection','read-only default CLI','fresh prepare CLI',
                'existing-output overwrite rejected']))
    records=[]
    def check(name,call):
        try:
            result=call();records.append(dict(name=name,passed=True,result=result))
            print(name+': PASS',flush=True)
        except Exception as exc:
            records.append(dict(name=name,passed=False,error=repr(exc)))
            save(HERE/'failure.json',dict(records=records))
            raise
    template=prep.TEMPLATE.read_text(encoding='utf-8-sig')
    def pure_case(r,v):
        seed=6. if r==3 else 3.
        text=prep.make_input(template,radius_um=r,voltage=v,mesh_path='C:/mesh/H1.mesh',
                            seed_nodes_path='C:/seed/nodes.csv',seed_voltage=seed)
        assert number(text,'target_voltage')==-v
        assert number(text,'restart_voltage_factor')==seed/v
        assert number(text,'hole_radius')==float(f'{r}e-6')
        assert f'V{round(v*10):03d}' in text
        assert number(text,'residual_tolerance')==number(text,'update_tolerance')==1e-9
        assert number(text,'target_joule_factor')==number(text,'restart_joule_factor')==1.
        assert 'heat_source_current_work=.true.' in text and "material_profile='ZNO_F1_NT'" in text
        return dict(radius_um=r,voltage_V=v,seed_voltage_V=seed,restart_voltage_factor=seed/v,
                    in_memory_input_sha256=hashlib.sha256(text.encode('ascii')).hexdigest())
    for r,v in sorted(prep.CONDITIONS):check(f'pure_r{r}_V{v:g}',lambda r=r,v=v:pure_case(r,v))
    common=dict(radius_um=3,voltage=6.5,mesh_path='C:/mesh/H1.mesh',seed_nodes_path='C:/seed/nodes.csv',seed_voltage=6.)
    check('reject_changed_tolerance',lambda:expect_reject(lambda:prep.make_input(template.replace('residual_tolerance=1e-9','residual_tolerance=1e-8'),**common)))
    check('reject_out_of_scope_radius_voltage',lambda:expect_reject(lambda:prep.make_input(template,**dict(common,radius_um=2))))
    check('reject_seed_above_target',lambda:expect_reject(lambda:prep.make_input(template,**dict(common,seed_voltage=8.))))
    geom=HERE.parent/'geometry_validation_01'
    meshes={r:geom/f'r{r}_A32_R14_Z14/H1.mesh' for r in (2,3,4)}
    coordinates={}
    def geometry(r):
        xyz,m=prep.validate_mesh(meshes[r],meshes[r].parent/'manifest.json',r)
        coordinates[r]=xyz
        return dict(nodes=len(xyz),radius_um=r,mesh_sha256=m['mesh_sha256'])
    for r in (2,3,4):check(f'r{r}_mesh_identity',lambda r=r:geometry(r))
    check('reject_wrong_mesh_radius',lambda:expect_reject(lambda:prep.validate_mesh(meshes[2],meshes[2].parent/'manifest.json',3)))
    source=ROOT/'runs/JE_SENSITIVITY_20261015/A32_S1_V060'
    check('actual_seed_exact_coordinates',lambda:prep.validate_seed_coordinates(source/'nodes.csv',coordinates[3],6.))
    changed=coordinates[3].copy();changed[0,0]+=1e-18
    check('reject_changed_coordinate',lambda:expect_reject(lambda:prep.validate_seed_coordinates(source/'nodes.csv',changed,6.)))
    check('reject_wrong_seed_electrode_voltage',lambda:expect_reject(lambda:prep.validate_seed_coordinates(source/'nodes.csv',coordinates[3],6.5)))
    for r in (2,4):
        old=ROOT/f'runs/F3_RADIUS_L0_20260908/F3_R0{r}_NT_L0/metrics.txt'
        m=prep.metrics(old)
        check(f'old_r{r}_drift_requires_explicit_optin',lambda m=m,r=r:expect_reject(lambda:prep.validate_source_metrics(m,r,3.)))
        def admit(m=m,r=r):
            assert prep.validate_source_metrics(m,r,3.,True)==prep.DRIFT
            return dict(source=str(old),role='INITIAL_GUESS_ONLY',JE_acceptance_evidence=False)
        check(f'old_r{r}_drift_3V_initial_guess_only',admit)
    check('reject_LQ_seed',lambda:expect_reject(lambda:prep.validate_source_metrics(dict(m,model='LQ'),4,3.,True)))
    args=[sys.executable,str(ROOT/'scripts/prepare_je_paper_input_20261018.py'),
          '--radius-um','3','--mesh',str(meshes[3]),'--seed-nodes',str(source/'nodes.csv'),
          '--seed-run',str(source),'--seed-voltage','6']
    def plan_cli(v):
        unused=HERE/f'default_plan_V{round(v*10):03d}_must_not_exist'
        assert not unused.exists()
        proc=subprocess.run(args+['--voltage',str(v),'--output',str(unused)],capture_output=True,timeout=60)
        (HERE/f'plan_V{round(v*10):03d}.stdout.json').write_bytes(proc.stdout)
        (HERE/f'plan_V{round(v*10):03d}.stderr.txt').write_bytes(proc.stderr)
        assert proc.returncode==0,proc.stderr.decode(errors='replace')
        plan=json.loads(proc.stdout)
        assert not unused.exists() and plan['status']=='INPUT_ONLY_NOT_SOLVED'
        assert plan['native_solver_started'] is False and plan['production_acceptance'] is False
        assert number(plan['input_text'],'target_voltage')==-v
        assert plan['restart_voltage_factor']==6./v
        return dict(output_directory_created=False,restart_voltage_factor=plan['restart_voltage_factor'])
    for v in (6.5,7.5):check(f'default_CLI_V{v:g}_read_only',lambda v=v:plan_cli(v))
    prepared=HERE/'prepared_r3_V065_input_only'
    def prepare_cli():
        proc=subprocess.run(args+['--voltage','6.5','--output',str(prepared),'--prepare'],capture_output=True,timeout=60)
        (HERE/'prepare.stdout.json').write_bytes(proc.stdout);(HERE/'prepare.stderr.txt').write_bytes(proc.stderr)
        assert proc.returncode==0,proc.stderr.decode(errors='replace')
        p=prep.read_json(prepared/'preparation_manifest.json')
        assert p['status']=='INPUT_ONLY_NOT_SOLVED' and not p['production_acceptance'] and not p['native_solver_started']
        assert prep.sha(prepared/'input.nml')==p['input_sha256']
        assert set(x.name for x in prepared.iterdir())=={'input.nml','preparation_manifest.json'}
        return dict(directory=str(prepared),input_sha256=p['input_sha256'],status=p['status'])
    check('prepare_CLI_writes_input_only',prepare_cli)
    def reject_overwrite():
        before={str(p):prep.sha(p) for p in prepared.iterdir()}
        proc=subprocess.run(args+['--voltage','6.5','--output',str(prepared),'--prepare'],capture_output=True,timeout=60)
        (HERE/'overwrite_rejection.stderr.txt').write_bytes(proc.stderr)
        assert proc.returncode!=0 and b'FileExistsError' in proc.stderr
        assert all(prep.sha(Path(p))==h for p,h in before.items())
        return dict(rejected=True,existing_files_unchanged=True)
    check('prepare_rejects_existing_output',reject_overwrite)
    unchanged=all(prep.sha(Path(p))==h for p,h in frozen.items())
    save(HERE/'report.json',dict(status='INPUT_PREPARER_CHECKED_NOT_PRODUCTION_ACCEPTANCE',
        all_software_checks_passed=unchanged and all(r['passed'] for r in records),
        checks=len(records),sources_unchanged=unchanged,records=records,
        native_solver_started=False,production_acceptance=False,
        limitation='Three-radius text/geometry identity, actual r3 direct seed and CLI behavior checked. Old drift source eligibility checked separately; no new cross-mesh transfer or re-equilibration has been run.'))
    save(HERE/'artifact_hashes.json',{str(p.relative_to(HERE)):prep.sha(p) for p in HERE.rglob('*') if p.is_file()})

if __name__=='__main__':main()
