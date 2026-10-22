"""Run frozen CLI regressions on both executables and compare stdout bytes."""
import argparse, hashlib, json, os, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def runtime_env():
    env = os.environ.copy()
    env['PATH'] = r'C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin;' + env['PATH']
    env.update(OMP_NUM_THREADS='4', MKL_NUM_THREADS='4', OMP_STACKSIZE='32M', TPS_PARDISO_STORAGE='IC')
    return env

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    report=[]
    for check in ['V01','V02','V03','V04','V05','V06','G02','G03','G04','G05','G06','G07','G08','L01','F1M','A1J']:
        results=[]
        for label,exe in [('baseline',a.baseline),('candidate',a.candidate)]:
            dest=a.output/f'{check}_{label}.txt'
            with dest.open('wb') as f:
                r=subprocess.run([str(exe.resolve()),'verify',check],cwd=ROOT,env=runtime_env(),stdout=f,stderr=subprocess.STDOUT)
            results.append((r.returncode,dest.read_bytes()))
        report.append(dict(check=check,baseline_exit=results[0][0],candidate_exit=results[1][0],
            byte_identical=results[0][1]==results[1][1],
            passed=all(x[0]==0 and b'verification_result=PASS' in x[1] for x in results) and results[0][1]==results[1][1]))
        (a.output/'report.json').write_text(json.dumps(report,indent=2))
        print(json.dumps(report[-1]),flush=True)
    raise SystemExit(0 if all(x['passed'] for x in report) else 1)

if __name__=='__main__':main()
