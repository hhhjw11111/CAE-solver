"""Strict three-level audit of requested voltage points; missing evidence fails closed."""
import argparse,csv,json,math
from pathlib import Path
from run_hv_sweep import read_metrics,classify,sha

LIMITS={'joule_W':.01,'theta_max_K':.01,'ring_theta_mean_K':.01,'electron_right_A':.01,
        'ring_E_mean_V_m':.02,'ring_Q_mean_W_m3':.02,'ring_abs_n_relative_mean':.02,
        'ring_abs_p_relative_mean':.02,
        'ring_Q_p95_W_m3':.05,'E_p95_V_m':.05}

def relative(a,b):
    if not math.isfinite(a) or not math.isfinite(b):return None
    return abs(a-b)/max(abs(a),abs(b),1e-300)

def audit(levels,voltages):
    # Three copies of one mesh, or different physics/executables, are not a
    # mesh-convergence study. Check the original manifests before comparisons.
    provenance_errors=[];manifests=[]
    if len(levels)!=3:provenance_errors.append('THREE_LEVELS_REQUIRED')
    if not voltages:provenance_errors.append('VOLTAGES_REQUIRED')
    for root in levels:
        path=root/'manifest.json'
        if not path.exists():
            provenance_errors.append('MISSING_MANIFEST');continue
        try:
            m=json.loads(path.read_text(encoding='utf-8-sig'));manifests.append(m)
            for key in ['mesh','executable']:
                source=Path(m.get(key) or '')
                if not source.is_file() or sha(source)!=m.get(key+'_sha256'):
                    provenance_errors.append(key.upper()+'_HASH_MISMATCH')
        except (OSError,ValueError,TypeError):
            provenance_errors.append('INVALID_MANIFEST')
    if len(manifests)==3:
        if len({m.get('domain_symmetry','FULL') for m in manifests})!=1:
            provenance_errors.append('INCONSISTENT_DOMAIN_SYMMETRY')
        if len({m.get('mesh_sha256') for m in manifests})!=3:
            provenance_errors.append('THREE_DISTINCT_MESHES_REQUIRED')
        for key in ['executable_sha256','supg_scale','model']:
            if any(key not in m for m in manifests) or len({m.get(key) for m in manifests})!=1:
                provenance_errors.append('INCONSISTENT_'+key.upper())
    records=[]
    for voltage in voltages:
        row={'voltage_V':voltage,'complete_solve_flux_all_levels':False,'mesh_pass':False,
             'status':'MISSING_EVIDENCE','changes':{}}
        data=[]
        for root in levels:
            # Accept either a sweep root or an actual single-point run. Never
            # rename/copy fields merely to make a completed study layout.
            d=root if (root/'metrics.txt').is_file() else root/f'V{round(voltage*10):03d}'
            if not(d/'execution.json').exists():break
            try:
                ex=json.loads((d/'execution.json').read_text())
                state,m=classify(d,ex['exit_code'])
                if state!='SOLVE_FLUX_PASS_MESH_PENDING':row['status']='SOLVE_OR_FLUX_FAIL';break
                if not(d/'field_audit.json').exists() or json.loads((d/'field_audit.json').read_text()).get('passed') is not True:
                    row['status']='FIELD_AUDIT_MISSING_OR_FAILED';break
                if m.get('voltage_V')!=voltage:
                    row['status']='VOLTAGE_MISMATCH';break
                manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8-sig'))
                if m.get('supg_scale')!=manifest.get('supg_scale'):
                    row['status']='SUPG_SCALE_MISMATCH';break
                if m.get('domain_symmetry','FULL')!=manifest.get('domain_symmetry','FULL'):
                    row['status']='DOMAIN_SYMMETRY_MISMATCH';break
                if not all(k in m and isinstance(m[k],float) and math.isfinite(m[k]) for k in LIMITS):
                    row['status']='MISSING_OR_NONFINITE_METRIC';break
                data.append(m)
            except (OSError,ValueError,TypeError,KeyError):
                row['status']='INVALID_POINT_EVIDENCE';break
        if len(data)==3:
            for key in ['model','heat_source','source','supg_directional_tau']:
                values=[m.get(key) for m in data]
                if None in values or len(set(values))!=1:
                    provenance_errors.append('INCONSISTENT_'+key.upper())
        if len(data)==3 and not provenance_errors:
            row['complete_solve_flux_all_levels']=True
            row['changes']={k:dict(L0_L1=relative(data[0][k],data[1][k]),
                L1_L2=relative(data[1][k],data[2][k]),limit=v) for k,v in LIMITS.items()}
            row['mesh_pass']=all(d['L1_L2'] is not None and d['L1_L2']<d['limit'] for d in row['changes'].values())
            row['status']='MESH_PASS_STABILIZATION_SENSITIVITY_PENDING' if row['mesh_pass'] else 'MESH_FAIL'
        elif len(data)==3:
            row['status']='INVALID_STUDY_PROVENANCE'
        records.append(row)
    return dict(voltages=records,all_mesh_pass=bool(records) and not provenance_errors and all(r['mesh_pass'] for r in records),
                provenance_errors=sorted(set(provenance_errors)),
                fully_validated=False,limitation='Requires stabilization-sensitivity and model-scope review; no experimental calibration implied.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--levels',nargs=3,type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--voltage',nargs='+',type=float,default=[i/2 for i in range(17)])
    a=p.parse_args();r=audit(a.levels,a.voltage)
    with a.output.open('x',encoding='utf-8') as f:json.dump(r,f,indent=2,allow_nan=False)
    print(json.dumps(r,indent=2));return 0 if r['all_mesh_pass'] else 1

if __name__=='__main__':raise SystemExit(main())
