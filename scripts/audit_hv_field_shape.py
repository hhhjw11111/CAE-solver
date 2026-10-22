"""Compare independently solved fields after an explicitly recorded mesh mapping.

Mapping data is interpolation for comparison, never a solved field. Carrier
errors are scaled by the density perturbation, not the large reference density.
The additional predeclared pointwise tolerance is 2%; no interval certificate.
"""
import argparse,json
from pathlib import Path
import numpy as np
from run_hv_sweep import sha,classify

TOLERANCE=.02

def difference(a,b,scale=None):
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Nonfinite field')
    denominator=max(float(np.max(abs(a))),float(np.max(abs(b))),1e-300) if scale is None else scale
    delta=abs(a-b)/denominator
    return dict(maximum_relative=float(np.max(delta)),rms_relative=float(np.sqrt(np.mean(delta**2))),
                normalization=denominator,limit=TOLERANCE,passed=bool(np.max(delta)<TOLERANCE))

def audit(mapping,source,target):
    records=[]; manifests=[]
    for directory in [source,target]:
        execution=json.loads((directory/'execution.json').read_text())
        manifest_check=json.loads((directory/'manifest.json').read_text())
        state,metrics=classify(directory,execution['exit_code'])
        if state!='SOLVE_FLUX_PASS_MESH_PENDING' or json.loads((directory/'field_audit.json').read_text()).get('passed') is not True:
            raise ValueError('Both source and target require completed solve/flux/field checks')
        if execution.get('input_sha256')!=sha(directory/'input.nml'):
            raise ValueError('Changed solver input')
        if execution.get('executable_sha256')!=manifest_check['executable_sha256']:
            raise ValueError('Execution identity mismatch')
        if sha(Path(manifest_check['executable']))!=manifest_check['executable_sha256']:
            raise ValueError('Changed solver executable')
        if sha(Path(manifest_check['mesh']))!=manifest_check['mesh_sha256']:
            raise ValueError('Changed solver mesh')
        if metrics.get('log_weighted_quadrature_order',3)!=manifest_check.get('log_weighted_quadrature_order',3):
            raise ValueError('Quadrature provenance mismatch')
        records.append(metrics)
        manifests.append(manifest_check)
    for key in ['voltage_V','model','heat_source','source','transport_discretization','domain_symmetry']:
        if records[0].get(key)!=records[1].get(key):
            raise ValueError('Different problem definition: '+key)
    if records[0].get('log_weighted_quadrature_order',3)!=records[1].get('log_weighted_quadrature_order',3):
        raise ValueError('Different quadrature orders in a mesh-field comparison')
    report=json.loads((mapping/'transfer_report.json').read_text())
    manifest=json.loads((target/'manifest.json').read_text())
    if report['source_nodes_sha256']!=sha(source/'nodes.csv') or report['target_mesh_sha256']!=manifest['mesh_sha256']:
        raise ValueError('Mapping provenance mismatch')
    if report['source_mesh_sha256']!=manifests[0]['mesh_sha256']:
        raise ValueError('Mapping source mesh mismatch')
    if sha(Path(manifest['mesh']))!=manifest['mesh_sha256']:
        raise ValueError('Changed target mesh')
    a=np.genfromtxt(mapping/'nodes.csv',delimiter=',',names=True)
    b=np.genfromtxt(target/'nodes.csv',delimiter=',',names=True)
    if len(a)!=len(b) or any(not np.array_equal(a[k],b[k]) for k in ['node','x_m','y_m','z_m']):
        raise ValueError('Different node coordinates or order')
    metrics={k:difference(a[k],b[k]) for k in ['phi_V','theta_K']}
    displacement_scale=max(max(float(np.max(abs(v[k]))) for k in ['ux_m','uy_m','uz_m']) for v in [a,b])
    for k in ['ux_m','uy_m','uz_m']:
        metrics[k]=difference(a[k],b[k],max(displacement_scale,1e-300))
    for k in ['eta_n','eta_p']:
        metrics[k+'_density_perturbation']=difference(np.expm1(a[k]),np.expm1(b[k]))
    files=[mapping/'nodes.csv',mapping/'transfer_report.json']
    for directory in [source,target]:
        files += [directory/name for name in ['nodes.csv','manifest.json','execution.json','input.nml','metrics.txt','field_audit.json']]
    return dict(scope='POINTWISE_MAPPED_FIELD_SHAPE_ONLY',comparisons=metrics,
                passed=all(v['passed'] for v in metrics.values()),fully_validated=False,
                hashes={str(f.absolute()):sha(f) for f in files})

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['mapping','source','target','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();result=audit(a.mapping,a.source,a.target)
    with a.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(result,indent=2));return 0 if result['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
