"""Run frozen JE migration verification; preserve each attempt and fail closed."""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path

from run_hv_regression import runtime_env

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / 'runs/JE_RECONSTRUCTION_20261018'
LEGACY = ['V01','V02','V03','V04','V05','V06','G02','G03','G04','G05','G06',
          'G07','G08','L01','F1M','A1J']

def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()

def write_json(path, data):
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

def check_build(exe, core_only=False):
    manifest = json.loads((exe.parent/'build_manifest.json').read_text(encoding='utf-8-sig'))
    records = [r for r in manifest['files'] if r['path'].startswith('src/')
               and (not core_only or not r['path'].startswith('src/verification/'))]
    executable = [r for r in manifest['files'] if (ROOT/r['path']).resolve()==exe.resolve()]
    return dict(source_count=len(records),
                source_mismatches=[r['path'] for r in records if sha(ROOT/r['path'])!=r['sha256']],
                executable_matches=len(executable)==1 and sha(exe)==executable[0]['sha256'])

def run_check(exe, check, out, label, timeout, exports=False):
    dest = out / f'{check}_{label}.txt'
    command = [str(exe),'verify',check]
    if exports:
        export_dir = out / f'{check}_fields'
        export_dir.mkdir()
        command.append(str(export_dir))
    t0 = time.monotonic()
    with dest.open('xb') as log:
        try:
            completed = subprocess.run(command,cwd=ROOT,env=runtime_env(),stdout=log,
                                       stderr=subprocess.STDOUT,timeout=timeout)
            code = completed.returncode
            timed_out = False
        except subprocess.TimeoutExpired:
            code = None
            timed_out = True
    content = dest.read_bytes()
    export_ok = not exports or any(export_dir.iterdir())
    if check=='A2P' and exports:
        export_ok = ((export_dir/'A2P_summary.csv').is_file()
                     and len(list(export_dir.glob('*_nodes.csv')))==18
                     and len(list(export_dir.glob('*_gauss_rows.csv')))==18
                     and len(list(export_dir.glob('*_gauss_summary.csv')))==18)
    record = dict(check=check,label=label,command=command,exit_code=code,timed_out=timed_out,
                  seconds=time.monotonic()-t0,log=str(dest),sha256=sha(dest),
                  exports_complete=export_ok,
                  passed=code==0 and export_ok and b'verification_result=PASS' in content
                         and b'verification_result=FAIL' not in content)
    print(json.dumps({k:v for k,v in record.items() if k not in ('command','log')}),flush=True)
    return record,content

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--attempt',required=True)
    p.add_argument('--timeout',type=int,default=1800)
    a=p.parse_args()
    if Path(a.attempt).name != a.attempt or a.attempt in ('.','..'):
        p.error('attempt must be a single new directory name')
    candidate = a.candidate.resolve()
    baseline = ROOT / 'build/ARC_JE_FINAL_20261012/TPSolverV2.exe'
    out = STUDY / 'stage3_verification' / a.attempt
    out.mkdir(parents=True,exist_ok=False)
    contract = json.loads((STUDY/'contract.json').read_text(encoding='utf-8'))
    core = []
    for old in contract['core_sources_before']:
        path = ROOT.parent / old['path']
        core.append(dict(path=old['path'],before=old['sha256'],after=sha(path)))
    core_unchanged = all(r['before']==r['after'] for r in core)
    frozen_gates = all(sha(ROOT/p)==h for p,h in contract['frozen_gate_scripts'].items())
    baseline_build = check_build(baseline,core_only=True)
    candidate_build = check_build(candidate)
    builds_match = all(b['source_count']>0 and not b['source_mismatches'] and b['executable_matches']
                       for b in [baseline_build,candidate_build])
    source_manifest = {str(f.relative_to(ROOT)):sha(f) for f in (ROOT/'src').rglob('*.f90')}
    runner_digest = sha(Path(__file__))
    protocol_hashes = {str(f.relative_to(ROOT)):sha(f) for f in (STUDY/'stage3_verification').glob('*protocol.json')}
    for extra in (STUDY/'stage3_verification').glob('*revision.json'):
        protocol_hashes[str(extra.relative_to(ROOT))] = sha(extra)
    write_json(out/'provenance.json',dict(candidate=str(candidate),candidate_sha256=sha(candidate),
               baseline=str(baseline),baseline_sha256=sha(baseline),sources=source_manifest,
               core_checks=core,core_sources_unchanged=core_unchanged,frozen_gates_unchanged=frozen_gates,
               baseline_build=baseline_build,candidate_build=candidate_build,
               runner_sha256=runner_digest,protocol_hashes=protocol_hashes))
    if not core_unchanged or not frozen_gates or not builds_match:
        raise SystemExit('Frozen physical sources or acceptance gates changed; no verification launch')
    legacy=[]
    for check in LEGACY:
        r0,c0=run_check(baseline,check,out,'baseline',a.timeout)
        r1,c1=run_check(candidate,check,out,'candidate',a.timeout)
        legacy.append(dict(check=check,baseline=r0,candidate=r1,byte_identical=c0==c1,
                           passed=r0['passed'] and r1['passed'] and c0==c1))
        write_json(out/'legacy_cli.json',legacy)
    targeted=[]
    for check in contract['stage3']['ids']:
        r,_=run_check(candidate,check,out,'targeted',a.timeout,exports=check in ['A1','A2','A2P'])
        targeted.append(r)
        write_json(out/'targeted.json',targeted)
    source_still_same = all(sha(ROOT/path)==digest for path,digest in source_manifest.items())
    source_still_same = source_still_same and sha(Path(__file__))==runner_digest
    source_still_same = source_still_same and all(sha(ROOT/path)==digest for path,digest in protocol_hashes.items())
    result=dict(all_passed=all(r['passed'] for r in legacy+targeted) and source_still_same,
        core_sources_unchanged=core_unchanged,frozen_gates_unchanged=frozen_gates,
        build_manifests_verified=builds_match,
        source_still_same=source_still_same,legacy_passed=sum(r['passed'] for r in legacy),
        legacy_total=len(legacy),targeted_passed=sum(r['passed'] for r in targeted),
        targeted_total=len(targeted),attempt=str(out),candidate=str(candidate),
        candidate_sha256=sha(candidate),baseline=str(baseline),baseline_sha256=sha(baseline),
        limitation='Numerical verification only; stage4 mesh/SUPG and stage5 production remain separate')
    write_json(out/'result.json',result)
    if result['all_passed']:
        acceptance = out.parent/'acceptance.json'
        if acceptance.exists():
            raise SystemExit('Acceptance already exists; inspect instead of overwriting it')
        write_json(acceptance,result)
    print(json.dumps(result,ensure_ascii=True),flush=True)
    raise SystemExit(0 if result['all_passed'] else 1)

if __name__=='__main__':
    main()
