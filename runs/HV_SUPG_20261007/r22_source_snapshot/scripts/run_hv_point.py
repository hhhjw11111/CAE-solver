"""Solve and independently audit one voltage, with an optional recorded initial guess.

The output directory must be new. Even a same-voltage transferred guess is
re-equilibrated; no copied or interpolated field is reported as a result.
"""
from pathlib import Path
import argparse,json,math,os,subprocess,time
from run_hv_sweep import ROOT,make_input,sha,classify,audit_nodes
from run_hv_regression import runtime_env
from hv_runtime import add_storage_arguments,configure_storage,recorded_runtime


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',type=Path,required=True);p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--voltage',type=float,required=True)
    p.add_argument('--seed-nodes',type=Path);p.add_argument('--seed-voltage',type=float,default=0.)
    p.add_argument('--supg-scale',type=float,default=1.)
    p.add_argument('--checkpoint-volts',type=float,default=0.)
    p.add_argument('--log-step-limit',type=float,default=1.,help='Step control only; convergence tolerances remain 1e-9')
    p.add_argument('--global-log-limit',action='store_true',help='Experimental global direction scaling')
    p.add_argument('--half-y',action='store_true',help='Experimental y>=0 symmetric branch; requires a half-domain mesh')
    p.add_argument('--log-weighted',action='store_true',help='Experimental density-normalized weak transport, without SUPG')
    p.add_argument('--omp-threads',type=int,default=4)
    p.add_argument('--mkl-threads',type=int,default=4)
    p.add_argument('--max-seconds',type=float,default=7200.)
    add_storage_arguments(p)
    a=p.parse_args()
    if a.log_weighted and a.supg_scale!=1.:
        raise ValueError('SUPG scale is inapplicable to log-weighted transport; leave its placeholder at 1')
    if not all(math.isfinite(x) for x in [a.voltage,a.seed_voltage,a.supg_scale,a.checkpoint_volts,a.max_seconds,a.log_step_limit]):
        raise ValueError('Finite input values are required')
    if (not 0<=a.seed_voltage<=a.voltage<=8 or a.supg_scale<=0 or a.max_seconds<=0 or
        a.mkl_threads<1 or a.omp_threads<1 or a.log_step_limit<=0 or a.ooc_core_mb<=0):
        raise ValueError('Invalid voltage, scale, runtime or thread count')
    if (a.seed_voltage>0 and not a.seed_nodes) or (a.voltage==0 and a.seed_nodes):
        raise ValueError('A nonzero seed voltage requires a file; zero is independently solved without a seed')
    if not 0<=a.checkpoint_volts<=a.voltage:
        raise ValueError('Checkpoint spacing must be zero or at most the target voltage')
    for key in ['exe','mesh','output','seed_nodes']:
        value=getattr(a,key)
        if value is not None:setattr(a,key,Path(os.path.abspath(value)))
    a.output.mkdir(parents=True,exist_ok=False)
    text=make_input(a.mesh,a.voltage,a.seed_nodes.parent if a.seed_nodes else None,a.seed_voltage,a.supg_scale,a.global_log_limit)
    if a.log_weighted:
        text=text.replace('transport_supg=.true.,','transport_supg=.false.,transport_log_weighted=.true.,')
    if a.half_y:text=text.replace('&problem_control','&problem_control\n half_y_symmetry=.true.,')
    text=text.replace('maximum_log_density_update=1.0,',f'maximum_log_density_update={a.log_step_limit:.17g},')
    if a.seed_nodes:
        text=text.replace((a.seed_nodes.parent/'nodes.csv').as_posix(),a.seed_nodes.as_posix())
        (a.output/'seed.json').write_text(json.dumps(dict(path=str(a.seed_nodes),sha256=sha(a.seed_nodes),
            voltage=a.seed_voltage,status='INITIAL_GUESS_REQUIRES_REEQUILIBRATION'),indent=2))
    if a.checkpoint_volts:
        text=text.replace('continuation_secant_predictor=.true.,',
            f'continuation_secant_predictor=.true.,continuation_voltage_checkpoint_spacing={a.checkpoint_volts/a.voltage:.17g},')
    inp=a.output/'input.nml';inp.write_text(text,encoding='ascii')
    env=runtime_env();env['MKL_NUM_THREADS']=str(a.mkl_threads);env['OMP_NUM_THREADS']=str(a.omp_threads)
    env=configure_storage(env,a.output,a.storage,a.ooc_core_mb,a.split_holes)
    manifest=dict(executable=str(a.exe),executable_sha256=sha(a.exe),mesh=str(a.mesh),mesh_sha256=sha(a.mesh),
        voltages=[a.voltage],supg_scale=a.supg_scale,model='FROZEN_DRIFT_HEAT_ZNO_NT',
        runtime_env=recorded_runtime(env),
        log_step_limit=a.log_step_limit,global_log_limit=a.global_log_limit,
        domain_symmetry='Y_HALF' if a.half_y else 'FULL')
    manifest['transport_discretization']='HEX20_LOG_WEIGHTED_EXPERIMENTAL' if a.log_weighted else 'HEX20_CONSISTENT_SUPG_EXPERIMENTAL'
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    start=time.monotonic();timed_out=False
    with (a.output/'stdout.txt').open('wb') as f:
        proc=subprocess.Popen([str(a.exe),'run-f1',str(inp),str(a.output)],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        (a.output/'process.json').write_text(json.dumps(dict(pid=proc.pid)))
        try:rc=proc.wait(timeout=a.max_seconds)
        except subprocess.TimeoutExpired:proc.kill();rc=proc.wait();timed_out=True
    (a.output/'execution.json').write_text(json.dumps(dict(exit_code=rc,seconds=time.monotonic()-start,
        timed_out=timed_out,input_sha256=sha(inp),executable_sha256=sha(a.exe)),indent=2))
    state,_=classify(a.output,rc)
    if state!='FAILED_NO_COMPLETE_STATE':
        fields=audit_nodes(a.output/'nodes.csv',a.mesh,a.voltage,half_y=a.half_y)
        (a.output/'field_audit.json').write_text(json.dumps(fields,indent=2,allow_nan=False))
        if not fields['passed']:state='FAILED_FIELD_AUDIT'
    (a.output/'summary.json').write_text(json.dumps(dict(state=state,mesh_validated=False,
        complete_quantitative_0_8V=False),indent=2))
    print(state,flush=True)
    return 0 if state=='SOLVE_FLUX_PASS_MESH_PENDING' else 1


if __name__=='__main__':raise SystemExit(main())
