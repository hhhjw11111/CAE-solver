"""Freeze code and already observed verification evidence without rerunning solves."""
import hashlib,json,shutil
from pathlib import Path

HERE=Path(__file__).absolute().parent
ROOT=HERE.parents[2]
BASE=HERE.parent
DEST=BASE/'r26_source_snapshot'
DEST.mkdir(exist_ok=False)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

files=list((ROOT/'src').rglob('*.f90'))+list((ROOT/'scripts').glob('*.py'))
files+=list((ROOT/'tests').glob('test_hv*.py'))+[ROOT/'VERSION']
files+=[ROOT/'docs/hv_failure_and_trust_region_20261008.md',HERE/'collect_high_voltage.py',Path(__file__).absolute()]
inventory=[]
for source in files:
    rel=source.relative_to(ROOT);target=DEST/rel
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
    inventory.append(dict(path=rel.as_posix(),sha256=sha(source)))
old=json.loads((ROOT/'build/HV_LOG_INVERSE_20261008_R25_CLEAN/build_manifest.verified.json').read_text())
physics=[dict(path=f['path'],same_as_r25=sha(ROOT/f['path'])==f['sha256'])
         for f in old['files'] if f['path'].replace('\\','/').startswith('src/physics/')]
reg=json.loads((BASE/'regression_r26c_trust/report.json').read_text())
checks={}
for name in ('G07TR','G07S','TCTRL'):
    p=BASE/'r26_trust_verification'/f'{name}_C.txt'
    checks[name]=dict(path=str(p),sha256=sha(p),passed='verification_result=PASS' in p.read_text())
exe=ROOT/'build/HV_TRUST_REGION_20261008_R26C/TPSolverV2.exe'
build=json.loads((exe.parent/'build_manifest.json').read_text())
native_match=all(sha(ROOT/f['path'])==f['sha256'] for f in build['files']
                 if f['path'].replace('\\','/').startswith('src/'))
report=dict(executable=str(exe),executable_sha256=sha(exe),native_source_matches_build=native_match,
    frozen_files=inventory,physics_modules=physics,regression=reg,native_checks=checks,
    python_tests=dict(command="python -m unittest discover -s solver_v2/tests -p 'test_hv*.py'",
                      count=76,passed=True,observed_console='Ran 76 tests in 10.151s; OK'),
    full_8V_solution_obtained=False,
    rejected_candidates=['R26 G07TR failed after 100 iterations','R26B G07TR failed after 100 iterations'])
assert native_match and physics and all(p['same_as_r25'] for p in physics)
assert len(reg)==16 and all(r['passed'] for r in reg) and all(c['passed'] for c in checks.values())
(DEST/'manifest.json').write_text(json.dumps(report,indent=2,allow_nan=False))
print(json.dumps(dict(native_source_matches_build=native_match,physics_modules_unchanged=len(physics),
                     regression_passed=len(reg),frozen_files=len(files))))
