"""Non-overwriting 0--8 V sweep with per-point evidence, never interpolated NT data.

The runner reports algebraic/flux status only. Mesh and model validity remain separate.
Run --prepare-only to inspect all prospective inputs before any solves.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, os, subprocess, time
from run_hv_regression import runtime_env
from hv_runtime import add_storage_arguments,configure_storage,recorded_runtime

ROOT=Path(os.path.abspath(__file__)).parents[1]
FIELDS=['voltage','state','exit_code','theta_max_K','joule_W','electron_right_A',
        'heat_balance_relative','electron_balance_relative','hole_balance_diagnostic',
        'joule_quadrature_4_5_relative','maximum_accepted_residual','maximum_accepted_update']

def sha(path):
    # Field tables can be hundreds of MB; hashing must not duplicate them in RAM
    # while the direct solver is factorizing a large mesh.
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()

def read_metrics(path):
    data={}
    if path.exists():
        for line in path.read_text(encoding='utf-8-sig').splitlines():
            if '=' not in line: continue
            k,v=line.split('=',1)
            try: data[k]=float(v)
            except ValueError: data[k]=v
    return data

def checkpoint_seed(directory,voltage):
    """Locate an actually accepted original-equation state, never an interpolant."""
    matches=[]
    for load_file in sorted(directory.glob('accepted_*_load.txt')):
        m=read_metrics(load_file)
        v=m.get('applied_voltage_V',math.nan)
        nodes=load_file.with_name(load_file.name.replace('_load.txt','_nodes.csv'))
        if (not isinstance(v,(float,int)) or not math.isfinite(v) or
            abs(abs(v)-voltage)>1e-12*max(1.,voltage)):
            continue
        if m.get('joule_factor')!=1:continue
        if (not str(m.get('status','')).startswith('ALGEBRAIC_CHECKPOINT_NOT_FLUX') or not nodes.is_file() or
            not all(isinstance(m.get(k),(float,int)) and math.isfinite(m[k]) and 0<=m[k]<1e-9 for k in ['residual','update'])):
            raise ValueError('Matching checkpoint lacks an accepted full-heat state: '+str(load_file))
        matches.append((nodes,load_file))
    if not matches:raise ValueError(f'No accepted checkpoint at {voltage:g} V; interpolation is not allowed')
    return matches[-1]

def classify(directory,rc):
    m=read_metrics(directory/'metrics.txt')
    required=['final_voltage_factor','final_joule_factor','maximum_accepted_residual',
              'maximum_accepted_update','theta_max_K','joule_W']
    complete=all(isinstance(m.get(k),(float,int)) and math.isfinite(m[k]) for k in required)
    complete=complete and m['final_voltage_factor']==1 and m['final_joule_factor']==1
    complete=complete and 0<=m['maximum_accepted_residual']<1e-9 and 0<=m['maximum_accepted_update']<1e-9
    complete=complete and (directory/'nodes.csv').exists()
    if not complete: return 'FAILED_NO_COMPLETE_STATE',m
    # Recheck the numeric gates independently of the native PASS label. A stale
    # label, a missing metric or a NaN must not promote a state to accepted data.
    limits={'maximum_stage_balance':1e-8,'heat_balance_relative':.01,
            'electron_balance_relative':.01,'hole_balance_diagnostic':.01,
            'joule_quadrature_4_5_relative':.001}
    if m.get('domain_symmetry')=='Y_HALF':
        limits.update({key:.01 for key in ['physical_exterior_heat_balance_relative',
            'physical_exterior_electron_balance_relative','physical_exterior_hole_balance_diagnostic',
            'symmetry_heat_flux_relative','symmetry_electron_flux_relative','symmetry_hole_flux_diagnostic']})
    gates=all(isinstance(m.get(k),(float,int)) and math.isfinite(m[k]) and 0<=m[k]<limit
              for k,limit in limits.items())
    if rc==0 and m.get('case_acceptance')=='PASS_NOT_MESH_ACCEPTANCE' and gates:
        return 'SOLVE_FLUX_PASS_MESH_PENDING',m
    return 'SOLVED_FLUX_GATE_FAILED',m

def audit_nodes(path,mesh,voltage,half_y=False):
    """Independent finite fields, node count, positive density/T and electrode check."""
    with mesh.open(encoding='utf-8-sig') as f:
        count=int(f.readline().split()[0])
        coords=[list(map(float,f.readline().split()[1:4])) for _ in range(count)]
    xmin=min(c[0] for c in coords);xmax=max(c[0] for c in coords)
    errors=[];boundary_error=0.0;minimum_t=math.inf;minimum_n=math.inf;minimum_p=math.inf
    total=0;symmetry_uy=0.0
    if half_y and (abs(min(c[1] for c in coords))>1e-15 or abs(max(c[1] for c in coords)-1e-5)>1e-15):
        errors.append('INVALID_HALF_Y_EXTENT')
    with path.open(encoding='utf-8-sig',newline='') as f:
        for i,row in enumerate(csv.DictReader(f),1):
            total=i
            if i>count: errors.append('EXCESS_NODES');break
            values=[float(row[k]) for k in ['x_m','y_m','z_m','ux_m','uy_m','uz_m','phi_V','eta_n','eta_p','theta_K','n_m3','p_m3']]
            if int(row['node'])!=i or not all(math.isfinite(v) for v in values):
                errors.append('NONFINITE_OR_UNORDERED_NODE');break
            if max(abs(values[j]-coords[i-1][j]) for j in range(3))>1e-16:
                errors.append('MESH_COORDINATE_MISMATCH');break
            minimum_t=min(minimum_t,300+float(row['theta_K']))
            minimum_n=min(minimum_n,float(row['n_m3']));minimum_p=min(minimum_p,float(row['p_m3']))
            if half_y and abs(values[1])<1e-15:symmetry_uy=max(symmetry_uy,abs(float(row['uy_m'])))
            xpos=values[0]
            if abs(xpos-xmin)<1e-12 or abs(xpos-xmax)<1e-12:
                target=0.0 if abs(xpos-xmin)<1e-12 else -voltage
                boundary_error=max(boundary_error,abs(float(row['phi_V'])-target),abs(float(row['theta_K'])))
    if total!=count: errors.append('INCOMPLETE_NODE_COUNT')
    if min(minimum_t,minimum_n,minimum_p)<=0: errors.append('NONPOSITIVE_T_OR_DENSITY')
    if boundary_error>1e-10: errors.append('ELECTRODE_BOUNDARY_ERROR')
    if half_y and symmetry_uy>1e-20:errors.append('SYMMETRY_DISPLACEMENT_ERROR')
    result=dict(passed=not errors,errors=errors,node_count=total,boundary_error=boundary_error,
                minimum_T_K=minimum_t,minimum_n_m3=minimum_n,minimum_p_m3=minimum_p)
    if half_y:result.update(domain_symmetry='Y_HALF',symmetry_uy_m=symmetry_uy)
    return result

def make_input(mesh,voltage,previous,prev_voltage,supg_scale,global_limit):
    restart=''
    if previous is not None:
        restart=(f" restart_nodes_file='{(previous/'nodes.csv').as_posix()}',\n"
                 f" restart_voltage_factor={prev_voltage/voltage:.17g}, restart_joule_factor=1.0,\n")
    denom=max(voltage,1.0)
    step_denom=voltage if voltage>0 else 0.25
    return f"""&case_control
 schema_version=2,case_id='HV_NT_V{int(round(voltage*10)):03d}',model_level='M2_THERMO_TWO_CARRIER',
 mesh_file='{mesh.as_posix()}'
/
&problem_control
 problem_type='CENTER_HOLE_PLATE', material_profile='ZNO_F1_NT',
 target_voltage={-voltage:.17g},target_voltage_factor=1.0,target_joule_factor=1.0,hole_radius=3e-6,
{restart}/
&solver_control
 transport_supg=.true.,supg_scale={supg_scale:.17g},supg_directional_tau=.true.,
 global_log_step_limit={'.true.' if global_limit else '.false.'},
 field_scales=1e-10,1e-10,1e-10,{denom:.17g},0.1,0.1,35.0,
 max_newton=40,newton_stagnation_window=4,residual_tolerance=1e-9,update_tolerance=1e-9,maximum_log_density_update=1.0,
 continuation_full_joule_during_voltage=.true.,
 continuation_secant_predictor=.true.,
 continuation_initial_step={min(1.0,0.25/step_denom):.17g},continuation_maximum_step={min(1.0,0.5/step_denom):.17g},
 continuation_minimum_step={0.005/step_denom:.17g},continuation_growth_factor=1.5,
 continuation_maximum_attempts=200,continuation_maximum_newton_per_attempt=40
/
"""

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--exe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mesh',type=Path,default=ROOT/'cases/F1_20260906/H1_L0.mesh')
    p.add_argument('--supg-scale',type=float,default=1.0)
    p.add_argument('--log-weighted',action='store_true',help='Opt-in original-physics density-normalized weak transport; disables SUPG')
    p.add_argument('--quadrature-order',type=int,choices=[3,5],default=3)
    p.add_argument('--half-y',action='store_true')
    p.add_argument('--start-seed-run',type=Path,help='Initial field for the first requested nonzero voltage; always re-solved.')
    p.add_argument('--seed-checkpoints',type=Path,help='Independently replay accepted checkpoints at the requested voltages; no interpolated seeds.')
    p.add_argument('--omp-threads',type=int,default=4)
    p.add_argument('--mkl-threads',type=int,default=4)
    add_storage_arguments(p)
    p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--seed-sweep',type=Path,help='Completed same-voltage fields used only as re-equilibrated initial guesses.')
    p.add_argument('--voltage',nargs='+',type=float,default=[i/2 for i in range(17)])
    p.add_argument('--component-limit',action='store_true')
    p.add_argument('--global-limit',action='store_true',help='Experimental; can stall on a poorly conditioned carrier block.')
    p.add_argument('--max-seconds-per-point',type=float,default=7200)
    a=p.parse_args()
    a.exe=Path(os.path.abspath(a.exe));a.mesh=Path(os.path.abspath(a.mesh));a.output=Path(os.path.abspath(a.output))
    if not a.voltage or any(not math.isfinite(v) or v<0 or v>8 for v in a.voltage) or a.voltage!=sorted(set(a.voltage)):
        raise ValueError('Voltage list must be unique, ascending and within [0,8].')
    if len({int(round(v*10)) for v in a.voltage})!=len(a.voltage):
        raise ValueError('Requested voltages collide in the sweep folder labels; use run_hv_point for close off-grid points')
    if a.quadrature_order!=3 and not a.log_weighted:
        raise ValueError('Order 5 is restricted to log-weighted transport')
    if a.log_weighted and a.supg_scale!=1:
        raise ValueError('SUPG scale does not apply to log-weighted transport')
    if (not math.isfinite(a.supg_scale) or a.supg_scale<=0 or
        not math.isfinite(a.max_seconds_per_point) or a.max_seconds_per_point<=0 or
        a.omp_threads<1 or a.mkl_threads<1 or a.ooc_core_mb<=0):
        raise ValueError('Invalid runtime or numerical control')
    previous=None;pv=0.0
    checkpoint_sources={}
    if a.seed_checkpoints:
        if a.start_seed_run or a.seed_sweep:
            raise ValueError('Checkpoint replay cannot be combined with other seed modes')
        a.seed_checkpoints=Path(os.path.abspath(a.seed_checkpoints))
        checkpoint_sources={v:checkpoint_seed(a.seed_checkpoints,v) for v in a.voltage if v>0}
    if a.start_seed_run:
        previous=Path(os.path.abspath(a.start_seed_run));seed_metrics=read_metrics(previous/'metrics.txt')
        pv=seed_metrics.get('voltage_V',math.nan)
        if (not math.isfinite(pv) or not 0<=pv<=a.voltage[0] or a.voltage[0]==0 or
            not (previous/'nodes.csv').is_file() or seed_metrics.get('final_voltage_factor')!=1 or
            seed_metrics.get('final_joule_factor')!=1):
            raise ValueError('Starting seed must be a complete field at or below the first nonzero voltage')
    base_env=runtime_env();base_env.update(OMP_NUM_THREADS=str(a.omp_threads),MKL_NUM_THREADS=str(a.mkl_threads))
    for key in list(base_env):
        if key.upper().startswith('MKL_PARDISO_OOC_'):del base_env[key]
    base_env.update(TPS_PARDISO_STORAGE=a.storage,TPS_PARDISO_SPLIT_HOLES='1' if a.split_holes else '0')
    manifest=dict(executable=str(a.exe),executable_sha256=sha(a.exe),mesh=str(a.mesh),mesh_sha256=sha(a.mesh),
        voltages=a.voltage,supg_scale=a.supg_scale,model='FROZEN_DRIFT_HEAT_ZNO_NT',
        acceptance='NO_POINT_IS_MESH_VALIDATED_BY_THIS_SCRIPT',seed_sweep=str(a.seed_sweep) if a.seed_sweep else None,
        domain_symmetry='Y_HALF' if a.half_y else 'FULL',
        transport_discretization='HEX20_LOG_WEIGHTED_EXPERIMENTAL' if a.log_weighted else 'HEX20_CONSISTENT_SUPG_EXPERIMENTAL',
        log_weighted_quadrature_order=a.quadrature_order,
        start_seed_run=str(previous) if previous else None,
        seed_checkpoints=str(a.seed_checkpoints) if a.seed_checkpoints else None,
        seed_checkpoints_manifest_sha256=(sha(a.seed_checkpoints/'manifest.json') if a.seed_checkpoints and
                                         (a.seed_checkpoints/'manifest.json').is_file() else None),
        ooc_core_mb=a.ooc_core_mb,
        runtime_env=recorded_runtime(base_env))
    a.output.mkdir(parents=True,exist_ok=False)
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    columns=FIELDS+(['domain_symmetry','integral_scope','physical_exterior_heat_balance_relative',
                     'symmetry_heat_flux_relative'] if a.half_y else [])
    rows=[dict(voltage=v,state='NOT_RUN',exit_code='') for v in a.voltage]
    if a.half_y:
        for row in rows:row.update(domain_symmetry='Y_HALF',integral_scope='ACTUAL_HALF_DOMAIN_NO_AUTOMATIC_DOUBLING')
    def save():
        with (a.output/'status.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
    save()
    for i,v in enumerate(a.voltage):
        dest=a.output/f'V{int(round(v*10)):03d}';dest.mkdir()
        inp=dest/'input.nml'
        initial,initial_voltage=previous,pv
        checkpoint_nodes=None;checkpoint_load=None
        if v in checkpoint_sources:
            checkpoint_nodes,checkpoint_load=checkpoint_sources[v]
            initial,initial_voltage=checkpoint_nodes.parent,v
        if a.seed_sweep and v>0:
            seed=Path(os.path.abspath(a.seed_sweep))/dest.name
            sm=read_metrics(seed/'metrics.txt')
            if ((seed/'nodes.csv').is_file() and sm.get('voltage_V')==v and
                sm.get('final_voltage_factor')==1 and sm.get('final_joule_factor')==1):
                initial,initial_voltage=seed,v
        text=make_input(a.mesh,v,initial,initial_voltage,a.supg_scale,a.global_limit and not a.component_limit)
        if a.log_weighted:
            text=text.replace('transport_supg=.true.,',
                f'transport_supg=.false.,transport_log_weighted=.true.,log_weighted_quadrature_order={a.quadrature_order},')
        if a.half_y:text=text.replace('&problem_control','&problem_control\n half_y_symmetry=.true.,')
        if checkpoint_nodes:
            text=text.replace((initial/'nodes.csv').as_posix(),checkpoint_nodes.as_posix())
        inp.write_text(text,encoding='ascii')
        env=configure_storage(base_env,dest,a.storage,a.ooc_core_mb,a.split_holes)
        point_manifest=dict(manifest,voltages=[v],runtime_env=recorded_runtime(env))
        (dest/'manifest.json').write_text(json.dumps(point_manifest,indent=2))
        seed_nodes=checkpoint_nodes or (initial/'nodes.csv' if initial else None)
        if seed_nodes is not None and seed_nodes.exists():
            seed_record=dict(path=str(seed_nodes),sha256=sha(seed_nodes),voltage=initial_voltage,
                             status='INITIAL_GUESS_ONLY_REQUIRES_NEW_SOLVE')
            if checkpoint_load:seed_record.update(checkpoint_load=str(checkpoint_load),checkpoint_load_sha256=sha(checkpoint_load))
            (dest/'seed.json').write_text(json.dumps(seed_record,indent=2))
        if a.prepare_only:
            rows[i]['state']='PREPARED_NOT_RUN';previous=dest;pv=v;save();continue
        rows[i]['state']='RUNNING';save();print(f'START {dest.name}',flush=True)
        t=time.monotonic();rc=None;timed_out=False
        with (dest/'stdout.txt').open('wb') as f:
            process=subprocess.Popen([str(a.exe),'run-f1',str(inp),str(dest)],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
            (dest/'process.json').write_text(json.dumps(dict(pid=process.pid)))
            try: rc=process.wait(timeout=a.max_seconds_per_point)
            except subprocess.TimeoutExpired:
                process.kill();rc=process.wait();timed_out=True
        (dest/'execution.json').write_text(json.dumps(dict(exit_code=rc,seconds=time.monotonic()-t,
            timed_out=timed_out,input_sha256=sha(inp),executable_sha256=sha(a.exe)),indent=2))
        state,m=classify(dest,rc)
        if state!='FAILED_NO_COMPLETE_STATE':
            fields=audit_nodes(dest/'nodes.csv',a.mesh,v,half_y=a.half_y)
            (dest/'field_audit.json').write_text(json.dumps(fields,indent=2,allow_nan=False))
            if not fields['passed']: state='FAILED_FIELD_AUDIT'
        rows[i].update({k:m[k] for k in columns if k in m})
        rows[i].update(state=state,exit_code=rc);save();print(f'END {dest.name}: {state}',flush=True)
        if state.startswith('FAILED_'): break
        # Failed flux points may seed exploratory solves, but NEVER count as accepted results.
        previous=dest;pv=v
    all_flux=all(r['state']=='SOLVE_FLUX_PASS_MESH_PENDING' for r in rows)
    (a.output/'summary.json').write_text(json.dumps(dict(all_requested_solve_flux_pass=all_flux,
        mesh_validated=False,complete_quantitative_0_8V=False),indent=2))
    return 0 if all_flux or a.prepare_only else 1

if __name__=='__main__':raise SystemExit(main())
