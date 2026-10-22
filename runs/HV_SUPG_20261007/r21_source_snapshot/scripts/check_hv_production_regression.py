"""Same-input baseline/candidate F1 field regression in separate fresh directories."""
import argparse, hashlib, json, subprocess, time
from pathlib import Path
from run_hv_regression import ROOT,runtime_env

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--baseline',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--case',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    executions=[]
    for label,exe in [('baseline',a.baseline),('candidate',a.candidate)]:
        dest=a.output/label;dest.mkdir();start=time.monotonic()
        with (dest/'stdout.txt').open('wb') as f:
            r=subprocess.run([str(exe.absolute()),'run-f1',str(a.case.absolute()),str(dest.absolute())],
                cwd=ROOT,env=runtime_env(),stdout=f,stderr=subprocess.STDOUT)
        executions.append(dict(label=label,exit_code=r.returncode,seconds=time.monotonic()-start))
        (a.output/'executions.json').write_text(json.dumps(executions,indent=2))
        print(executions[-1],flush=True)
    files={}
    for name in ['metrics.txt','nodes.csv','field_samples.csv']:
        paths=[a.output/x/name for x in ['baseline','candidate']]
        files[name]=dict(exists=all(p.exists() for p in paths))
        if files[name]['exists']:
            hashes=[hashlib.sha256(p.read_bytes()).hexdigest() for p in paths]
            files[name].update(sha256=hashes,byte_identical=hashes[0]==hashes[1])
    passed=all(e['exit_code']==0 for e in executions) and all(v.get('byte_identical') for v in files.values())
    result=dict(passed=passed,files=files,case=str(a.case.absolute()))
    (a.output/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
    return 0 if passed else 1

if __name__=='__main__':raise SystemExit(main())
