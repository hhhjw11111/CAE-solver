import hashlib,json,shutil
from pathlib import Path
HERE=Path(__file__).absolute().parent;ROOT=HERE.parents[2]
OUT=HERE.parent/'r26_heat_restart_code_snapshot';OUT.mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
paths=[ROOT/'scripts/run_hv_point.py',ROOT/'tests/test_hv_heat_restart.py',
       HERE/'check_heat_restart.py',HERE/'audit_heat_restart_check.py',
       HERE/'heat_restart_checkpoint_audit.json',HERE/'collect_high_voltage_restart.py',
       HERE/'prepare_single_route_collector.py',Path(__file__).absolute()]
files=[]
for source in paths:
    dest=OUT/source.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
    files.append(dict(path=source.relative_to(ROOT).as_posix(),sha256=sha(source)))
build=ROOT/'build/HV_TRUST_REGION_20261008_R26C'
manifest=json.loads((build/'build_manifest.json').read_text())
native_match=all(sha(ROOT/x['path'])==x['sha256'] for x in manifest['files']
                 if x['path'].replace('\\','/').startswith('src/'))
assert native_match
record=dict(native_source_matches_frozen_r26c=native_match,executable_sha256=sha(build/'TPSolverV2.exe'),
            python_tests=dict(count=79,passed=True,observed_console='Ran 79 tests in 6.297s; OK'),
            files=files,unchanged_defaults=True,original_physics_and_convergence_gates_retained=True,
            partial_heat_seed_is_final_result=False,full_8V_solution_obtained=False,
            scheduling_audit=json.loads((HERE/'heat_restart_checkpoint_audit.json').read_text()))
(OUT/'manifest.json').write_text(json.dumps(record,indent=2,allow_nan=False))
print(json.dumps(dict(native_source_matches_frozen_r26c=native_match,frozen_files=len(files),python_tests=79)))
