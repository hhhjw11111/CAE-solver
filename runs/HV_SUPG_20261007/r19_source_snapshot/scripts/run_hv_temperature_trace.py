"""Trace an algebraic branch using a nodal temperature as the continuation parameter.

Every original PDE is retained; voltage is recovered as an unknown. This runner
does not certify a voltage point, a turning point, or dynamical stability. Each
candidate requires an ordinary fixed-voltage replay and independent flux checks.
"""
from pathlib import Path
import argparse,csv,json,math,os,subprocess,time
from run_hv_sweep import ROOT,make_input,sha,read_metrics
from run_hv_regression import runtime_env
from hv_runtime import add_storage_arguments,configure_storage,recorded_runtime


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',type=Path,required=True);p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--seed-nodes',type=Path,required=True)
    p.add_argument('--seed-voltage',type=float,required=True)
    p.add_argument('--step-k',type=float,default=2.);p.add_argument('--points',type=int,default=8)
    p.add_argument('--supg-scale',type=float,default=1.)
    p.add_argument('--omp-threads',type=int,default=4)
    p.add_argument('--mkl-threads',type=int,default=4)
    p.add_argument('--max-seconds',type=float,default=7200.)
    add_storage_arguments(p)
    a=p.parse_args()
    if (not all(math.isfinite(v) for v in [a.seed_voltage,a.step_k,a.supg_scale,a.max_seconds]) or
        not 0<a.seed_voltage<8 or a.step_k<=0 or not 1<=a.points<=100 or a.supg_scale<=0 or a.max_seconds<=0 or
        min(a.omp_threads,a.mkl_threads)<1 or a.ooc_core_mb<=0):
        raise ValueError('Invalid diagnostic parameters')
    for key in ['exe','mesh','output','seed_nodes']:setattr(a,key,Path(os.path.abspath(getattr(a,key))))
    a.output.mkdir(parents=True,exist_ok=False)
    text=make_input(a.mesh,8.,a.seed_nodes.parent,a.seed_voltage,a.supg_scale,False)
    text=text.replace((a.seed_nodes.parent/'nodes.csv').as_posix(),a.seed_nodes.as_posix())
    inp=a.output/'input.nml';inp.write_text(text,encoding='ascii')
    env=runtime_env();env['OMP_NUM_THREADS']=str(a.omp_threads);env['MKL_NUM_THREADS']=str(a.mkl_threads)
    env=configure_storage(env,a.output,a.storage,a.ooc_core_mb,a.split_holes)
    manifest=dict(executable=str(a.exe),executable_sha256=sha(a.exe),mesh=str(a.mesh),mesh_sha256=sha(a.mesh),
        seed_nodes=str(a.seed_nodes),seed_nodes_sha256=sha(a.seed_nodes),seed_voltage=a.seed_voltage,
        step_K=a.step_k,maximum_increments=a.points,supg_scale=a.supg_scale,
        runtime_env=recorded_runtime(env),
        scope='DIAGNOSTIC_ALGEBRAIC_BRANCH_NOT_DYNAMICAL_STABILITY_OR_VOLTAGE_ACCEPTANCE')
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    timed_out=False;start=time.monotonic()
    with (a.output/'stdout.txt').open('wb') as f:
        proc=subprocess.Popen([str(a.exe),'trace-f1-temperature',str(inp),str(a.output),str(a.step_k),str(a.points)],
            cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        (a.output/'process.json').write_text(json.dumps(dict(pid=proc.pid)))
        try:rc=proc.wait(timeout=a.max_seconds)
        except subprocess.TimeoutExpired:proc.kill();rc=proc.wait();timed_out=True
    (a.output/'execution.json').write_text(json.dumps(dict(exit_code=rc,seconds=time.monotonic()-start,
        timed_out=timed_out,input_sha256=sha(inp),executable_sha256=sha(a.exe)),indent=2))
    rows=[];table=a.output/'inverse_branch.csv'
    if table.exists():
        with table.open(newline='') as f:
            for row in csv.DictReader(f):
                try:values={k:float(v) for k,v in row.items()}
                except (ValueError,TypeError):continue
                if all(math.isfinite(v) for v in values.values()):rows.append(values)
    summary=dict(completed_diagnostic=rc==0,accepted_algebraic_samples=len(rows),
        highest_sampled_voltage_V=max([r['voltage_V'] for r in rows],default=None),
        stop_reason=read_metrics(a.output/'metrics.txt').get('inverse_stop','REQUESTED_STEPS_OR_FAILURE'),
        fixed_voltage_replay_passed=False,mesh_validated=False,dynamical_stability_verified=False,
        complete_quantitative_0_8V=False)
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    print(json.dumps(summary),flush=True)
    return 0 if rc==0 else 1


if __name__=='__main__':raise SystemExit(main())
