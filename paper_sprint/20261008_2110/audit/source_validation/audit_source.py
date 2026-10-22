from pathlib import Path
import hashlib, json, subprocess, datetime, os, re, csv

ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
VAL=ROOT/'paper_sprint/20261008_2110/validation/runtime_retry'
VAL.mkdir(parents=True,exist_ok=True)
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def record(p):
    p=Path(p)
    return {'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}
def run(args,name,cwd=VAL,timeout=180):
    start=datetime.datetime.now(datetime.timezone.utc).isoformat()
    try:
        p=subprocess.run([str(a) for a in args],cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=timeout,env={**os.environ,'PATH':'C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin;'+os.environ.get('PATH',''),'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
        raw=p.stdout;rc=p.returncode
    except subprocess.TimeoutExpired as e:
        raw=e.stdout or b'';rc='TIMEOUT'
    target=VAL/(name+'.txt');target.write_bytes(raw)
    return {'command':[str(a) for a in args],'cwd':str(cwd),'started_utc':start,'ended_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exit_code':rc,'stdout':record(target),'threads':{'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}}

sources=sorted((ROOT/'src').rglob('*.f90'))
source_records=[record(p) for p in sources]
build=ROOT/'build/HV_LOG_INVERSE_20261008_R25_CLEAN'
manifest=json.loads((build/'build_manifest.json').read_text(encoding='utf-8-sig'))
checks=[]
for f in manifest['files']:
    if f['path'].startswith('src/') or f['path'].endswith('TPSolverV2.exe'):
        p=ROOT/f['path']; checks.append({'path':str(p),'recorded_sha256':f['sha256'],'actual_sha256':sha(p),'matches':sha(p)==f['sha256']})
exe=build/'TPSolverV2.exe'
commands=[]
for args,name in [([exe,'version'],'version'),([exe,'verify','V01'],'V01'),([exe,'verify','V02'],'V02'),([exe,'verify','G02'],'G02'),([exe,'verify','F1M'],'F1M'),([exe,'verify','G03'],'G03')]:
    commands.append(run(args,name))
    (VAL/'command_manifest.json').write_text(json.dumps(commands,ensure_ascii=False,indent=2),encoding='utf-8')
arch=ROOT/'runs/HV_SUPG_20261007/log_weighted_production_bar_r22'
commands.append(run([arch/'log_weighted_probe.exe'],'production_bar_replay',timeout=300))
commands[-1]['prior_stdout']=record(arch/'stdout.txt')
commands[-1]['byte_identical_to_prior_stdout']=(VAL/'production_bar_replay.txt').read_bytes()==(arch/'stdout.txt').read_bytes()
(VAL/'command_manifest.json').write_text(json.dumps(commands,ensure_ascii=False,indent=2),encoding='utf-8')
envpaths=[Path('C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/ifx.exe'),Path('C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/mkl_core.lib'),Path('C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/mkl_intel_lp64.lib'),Path('C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/mkl_sequential.lib')]
audit={'audit_time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'repo':str(ROOT),'source_files':source_records,'build_manifest':record(build/'build_manifest.json'),'source_exe_build_matches':checks,'environment_files':[record(p) for p in envpaths if p.exists()],'commands':commands,'no_source_changes':True,'scope':'Read-only source audit and isolated executable replay; existing current build and archived production-kernel harness, no 3D solve.'}
(OUT/'source_audit_runtime_retry.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
for cmd,filename in [(['git','diff','a63ad1db24a5d5813afa4fbbe7eaa6e7235950bc','--','src/physics','src/fem/element_kernel_m.f90','src/app/f1_case_runner_m.f90','src/core/model_types_m.f90'],'baseline_model_diff.patch'),(['git','status','--porcelain=v1'],'git_status_observed.txt')]:
    p=subprocess.run(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (OUT/filename).write_bytes(p.stdout)
rows=[]
for label,p in [('prototype_r21',ROOT/'runs/HV_SUPG_20261007/log_weighted_bar_r21/stdout.txt')]:
    for line in p.read_text(errors='replace').splitlines():
        match=re.match(r'Pe=\s*([\d.]+),carrier=(\d),elements=(\d+),density_and_current_errors=\s*([\d.E+-]+)\s+([\d.E+-]+)',line)
        if match: rows.append(dict(zip(['Pe','carrier','elements','sample_max_density_error','sample_max_current_error'],match.groups()),source=str(p),kind=label))
with (VAL/'archived_bar_errors.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=['Pe','carrier','elements','sample_max_density_error','sample_max_current_error','source','kind']);w.writeheader();w.writerows(rows)
print(json.dumps({'source_count':len(sources),'build_comparison_count':len(checks),'mismatches':[x['path'] for x in checks if not x['matches']],'commands':[(x['command'][-1],x['exit_code']) for x in commands],'archived_bar_rows':len(rows)},indent=2))
