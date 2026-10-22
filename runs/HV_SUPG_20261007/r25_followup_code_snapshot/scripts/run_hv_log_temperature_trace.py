"""Trace the original equations with temperature parameterization, diagnostic only.

Every saved branch sample requires a separate fixed-voltage solve and physical
flux audit. Neither an inverse solve nor a descending branch proves dynamical
stability, nonexistence of other branches, or complete 0--8 V credibility.
"""
import argparse,json,math,os,subprocess,time
from pathlib import Path
from run_hv_sweep import ROOT,make_input,sha
from run_hv_regression import runtime_env
from hv_runtime import add_storage_arguments,configure_storage,recorded_runtime

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['exe','mesh','seed-nodes','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--seed-voltage',type=float,required=True)
    p.add_argument('--target-voltage',type=float,default=8.)
    p.add_argument('--step-k',type=float,default=4.)
    p.add_argument('--increments',type=int,default=8)
    p.add_argument('--half-y',action='store_true')
    p.add_argument('--threads',type=int,default=2)
    p.add_argument('--max-seconds',type=float,default=7200.)
    add_storage_arguments(p);a=p.parse_args()
    if not all(math.isfinite(x) for x in [a.seed_voltage,a.target_voltage,a.step_k,a.max_seconds]):
        raise ValueError('Finite parameters required')
    if not(0<a.seed_voltage<a.target_voltage<=8 and a.step_k>0 and 1<=a.increments<=100 and a.threads>0 and a.max_seconds>0):
        raise ValueError('Invalid trace parameters')
    for key in ['exe','mesh','seed_nodes','output']:setattr(a,key,Path(os.path.abspath(getattr(a,key))))
    for f in [a.exe,a.mesh,a.seed_nodes]:
        if not f.is_file():raise ValueError('Missing input '+str(f))
    a.output.mkdir(parents=True,exist_ok=False)
    s=make_input(a.mesh,a.target_voltage,a.seed_nodes.parent,a.seed_voltage,1.,False)
    s=s.replace('transport_supg=.true.,','transport_supg=.false.,transport_log_weighted=.true.,')
    s=s.replace((a.seed_nodes.parent/'nodes.csv').as_posix(),a.seed_nodes.as_posix())
    if a.half_y:s=s.replace('&problem_control','&problem_control\n half_y_symmetry=.true.,')
    inp=a.output/'input.nml';inp.write_text(s,encoding='ascii')
    env=runtime_env();env['OMP_NUM_THREADS']=str(a.threads);env['MKL_NUM_THREADS']=str(a.threads)
    env=configure_storage(env,a.output,a.storage,a.ooc_core_mb,a.split_holes)
    manifest=dict(scope='ORIGINAL_MODEL_ALGEBRAIC_BRANCH_DIAGNOSTIC_ONLY',fully_validated=False,
                  executable=str(a.exe),executable_sha256=sha(a.exe),mesh=str(a.mesh),mesh_sha256=sha(a.mesh),
                  seed_nodes=str(a.seed_nodes),seed_sha256=sha(a.seed_nodes),seed_voltage_V=a.seed_voltage,
                  target_voltage_V=a.target_voltage,temperature_step_K=a.step_k,increments=a.increments,
                  domain_symmetry='Y_HALF' if a.half_y else 'FULL',runtime_env=recorded_runtime(env),
                  input_sha256=sha(inp),transport_discretization='HEX20_LOG_WEIGHTED_EXPERIMENTAL')
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    start=time.monotonic();timed_out=False
    with (a.output/'stdout.txt').open('wb') as f:
        proc=subprocess.Popen([str(a.exe),'trace-f1-temperature',str(inp),str(a.output),str(a.step_k),str(a.increments)],
                              cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        (a.output/'process.json').write_text(json.dumps(dict(pid=proc.pid)))
        try:rc=proc.wait(timeout=a.max_seconds)
        except subprocess.TimeoutExpired:proc.kill();rc=proc.wait();timed_out=True
    result=dict(exit_code=rc,seconds=time.monotonic()-start,timed_out=timed_out,input_sha256=sha(inp),
                executable_sha256=sha(a.exe),fully_validated=False,requires_fixed_voltage_flux_replay=True)
    (a.output/'execution.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
    return rc
if __name__=='__main__':raise SystemExit(main())
