"""Read-only NT/JE transport diagnostics on complete saved HEX20 fields.

No native solver is started. All Gauss samples are retained; the statistics and
endpoint envelopes are descriptive, with no new Pe/gradient acceptance limit.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import re
import shutil
import numpy as np
import plot_je_8v_clouds_20261013 as fe

ROOT=Path(__file__).resolve().parents[1]
JE='E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE'
SOURCES=['src/core/constants_m.f90','src/physics/zno_candidate_m.f90',
         'src/physics/transport_m.f90','src/fem/transport_supg_m.f90','src/mesh/hex20_shape_m.f90']
UNITS={'xyz_m':'m','natural_coordinates':'dimensionless','weight_m3':'m^3','detJ_m3':'m^3',
       'E_V_m':'V/m','E_magnitude_V_m':'V/m','temperature_K':'K','grad_temperature_K_m':'K/m',
       'h_direction_m':'m','g_direction_m2':'1/m^2','Pe':'dimensionless',
       'eta_n':'dimensionless','eta_p':'dimensionless','n_over_n0':'dimensionless','p_over_p0':'dimensionless',
       'n_relative_perturbation':'n/n0-1','p_relative_perturbation':'p/p0-1',
       'n_abs_relative_perturbation':'abs(n/n0-1)','p_abs_relative_perturbation':'abs(p/p0-1)',
       'grad_eta_n_per_m':'1/m','grad_eta_p_per_m':'1/m',
       'grad_n_relative_per_m':'1/m','grad_p_relative_per_m':'1/m',
       'grad_eta_n_magnitude_per_m':'1/m','grad_eta_p_magnitude_per_m':'1/m',
       'grad_n_relative_magnitude_per_m':'1/m','grad_p_relative_magnitude_per_m':'1/m'}


def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def save(path,value):
    with Path(path).open('x',encoding='utf-8') as stream:
        json.dump(value,stream,indent=2,ensure_ascii=False,allow_nan=False);stream.write('\n')


def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def input_value(text,key):
    matches=re.findall(r'(?<![A-Za-z0-9_])'+key+r"\s*=\s*('[^']*'|[^,\n/]+)",text)
    if len(matches)!=1:raise ValueError('Expected exactly one input '+key)
    return matches[0].strip().strip("'")


def rule(order):
    gp,gw=np.polynomial.legendre.leggauss(order)
    points=np.array([[r,s,t] for t in gp for s in gp for r in gp])
    weights=np.array([a*b*c for c in gw for b in gw for a in gw])
    return points,weights


def direction(E,T,G):
    """Exact scalar-Einstein reduction of native regularized gdir, not |cell|."""
    trace=np.trace(G,axis1=1,axis2=2)
    e2=np.sum(E*E,axis=1)
    vt=fe.KB*T/fe.CHARGE
    epsilon_E=1e-12*vt*vt*trace
    directional=(np.einsum('pi,pij,pj->p',E,G,E)+epsilon_E*trace/3)/(e2+epsilon_E)
    h=2/np.sqrt(directional)
    magnitude=np.sqrt(e2)
    pe=magnitude*h/(2*vt)
    return directional,h,pe,magnitude


def sample(coords,state,points,weights):
    n,dn=fe.shapes(points),fe.derivatives(points)
    jac=np.einsum('pka,aj->pkj',dn,coords)
    det=np.linalg.det(jac)
    if np.any(det<=0):raise ValueError('Nonpositive sampled Jacobian')
    ji=np.linalg.inv(jac)
    gradient=np.einsum('pik,pka->pia',ji,dn)
    G=np.einsum('pik,pjk->pij',ji,ji)
    value=n@state
    grad=np.einsum('pia,af->pif',gradient,state)
    E=-grad[:,:,3];T=fe.T0+value[:,6]
    if np.any(T<=0):raise ValueError('Nonpositive Gauss-point temperature')
    gdir,h,pe,emag=direction(E,T,G)
    with np.errstate(over='raise',invalid='raise'):
        nr,pr=np.exp(value[:,4]),np.exp(value[:,5])
        un,up=np.expm1(value[:,4]),np.expm1(value[:,5])
        gun=nr[:,None]*grad[:,:,4];gup=pr[:,None]*grad[:,:,5]
    if np.any(nr<=0) or np.any(pr<=0):raise ValueError('Unrepresentable positive carrier ratio')
    result=dict(xyz_m=n@coords,natural_coordinates=points.copy(),weight_m3=det*weights,detJ_m3=det,
        E_V_m=E,E_magnitude_V_m=emag,temperature_K=T,grad_temperature_K_m=grad[:,:,6],
        h_direction_m=h,g_direction_m2=gdir,Pe=pe,eta_n=value[:,4],eta_p=value[:,5],
        n_over_n0=nr,p_over_p0=pr,n_relative_perturbation=un,p_relative_perturbation=up,
        n_abs_relative_perturbation=np.abs(un),p_abs_relative_perturbation=np.abs(up),
        grad_eta_n_per_m=grad[:,:,4],grad_eta_p_per_m=grad[:,:,5],
        grad_n_relative_per_m=gun,grad_p_relative_per_m=gup,
        grad_eta_n_magnitude_per_m=np.linalg.norm(grad[:,:,4],axis=1),
        grad_eta_p_magnitude_per_m=np.linalg.norm(grad[:,:,5],axis=1),
        grad_n_relative_magnitude_per_m=np.linalg.norm(gun,axis=1),
        grad_p_relative_magnitude_per_m=np.linalg.norm(gup,axis=1))
    if not all(np.isfinite(x).all() for x in result.values()):
        raise ValueError('Nonfinite diagnostic sample; no filtering or clipping permitted')
    return result


def closed_form_check():
    """Independent affine-box expressions; tolerances check arithmetic only."""
    lengths=np.array([2e-6,4e-6,1e-6]);center=np.array([1e-6,-2e-6,.3e-6])
    coords=center+fe.NATURAL*lengths/2
    gn=np.array([1e4,-2e4,3e4]);gp=np.array([-3e4,1e4,2e4])
    cases=[]
    for order in (3,5):
        pts,w=rule(order);xyz=center+pts*lengths/2
        for electric in (np.zeros(3),np.array([1e5,0.,0.]),np.array([1e5,-2e5,3e4])):
            T=350.;state=np.zeros((20,7));state[:,3]=-coords@electric
            state[:,4]=.13+coords@gn;state[:,5]=-.17+coords@gp;state[:,6]=T-300.
            d=sample(coords,state,pts,w)
            gdiag=4/lengths**2;tr=float(gdiag.sum());vt=fe.KB*T/fe.CHARGE
            epsilon=1e-12*vt*vt*tr
            exact_g=(float(np.sum(electric**2*gdiag))+epsilon*tr/3)/(float(electric@electric)+epsilon)
            exact_h=2/math.sqrt(exact_g);exact_pe=float(np.linalg.norm(electric))*exact_h/(2*vt)
            en=.13+xyz@gn;ep=-.17+xyz@gp
            np.testing.assert_allclose(d['xyz_m'],xyz,rtol=2e-12,atol=1e-20)
            np.testing.assert_allclose(d['E_V_m'],np.broadcast_to(electric,d['E_V_m'].shape),rtol=2e-12,atol=1e-7)
            np.testing.assert_allclose(d['temperature_K'],T,rtol=2e-12,atol=1e-10)
            np.testing.assert_allclose(d['h_direction_m'],exact_h,rtol=2e-12,atol=1e-20)
            np.testing.assert_allclose(d['Pe'],exact_pe,rtol=2e-12,atol=1e-12)
            np.testing.assert_allclose(d['n_relative_perturbation'],np.expm1(en),rtol=2e-12,atol=1e-13)
            np.testing.assert_allclose(d['p_relative_perturbation'],np.expm1(ep),rtol=2e-12,atol=1e-13)
            np.testing.assert_allclose(d['grad_n_relative_per_m'],np.exp(en)[:,None]*gn,rtol=2e-12,atol=1e-7)
            np.testing.assert_allclose(d['grad_p_relative_per_m'],np.exp(ep)[:,None]*gp,rtol=2e-12,atol=1e-7)
            np.testing.assert_allclose(d['weight_m3'].sum(),np.prod(lengths),rtol=2e-12,atol=0)
            # Separate check with the native velocity/diffusion form for both carriers.
            errors=[]
            for mu0,beta in [(fe.MU_N0,fe.BETA_N),(fe.MU_P0,fe.BETA_P)]:
                mu=mu0*(T/300.)**(-beta);diff=vt*mu;v=mu*electric
                eps=1e-12*diff*diff*tr
                native_g=(float(np.sum(v*v*gdiag))+eps*tr/3)/(float(v@v)+eps)
                native_pe=float(np.linalg.norm(v))/(diff*math.sqrt(native_g))
                errors.append(abs(native_pe-exact_pe)/max(abs(exact_pe),1.))
            if max(errors)>2e-12:raise AssertionError('Einstein/native directional identity failed')
            cases.append(dict(order=order,E_V_m=electric.tolist(),temperature_K=T,h_exact_m=exact_h,
                Pe_exact=exact_pe,Pe_max_absolute_error=float(np.max(abs(d['Pe']-exact_pe))),
                native_n_p_Pe_relative_errors=errors,zero_field=(bool(np.all(electric==0)))))
    return dict(status='CLOSED_FORM_IMPLEMENTATION_CHECKS_PASSED',cases=cases,
        scope='Affine rectangular HEX20, uniform E/T and affine eta with exact exp(eta) gradients; numerical arithmetic checks only, no production threshold.')


def evidence(directory):
    directory=Path(directory).resolve(strict=True)
    # Completion checks precede any node-field read, so incomplete runs are excluded explicitly.
    manifest=read(directory/'manifest.json');execution=read(directory/'execution.json')
    summary=read(directory/'summary.json');metrics=fe.read_metrics(directory/'metrics.txt')
    if execution.get('exit_code')!=0 or execution.get('timed_out') or summary.get('state')!='SOLVE_FLUX_PASS_MESH_PENDING':
        raise ValueError('Require a completed saved solve; never sample an active or failed native run')
    if not read(directory/'field_audit.json').get('passed'):
        raise ValueError('Saved field audit missing or failed')
    mesh=Path(manifest['mesh']).resolve(strict=True);exe=Path(manifest['executable']).resolve(strict=True)
    if sha(mesh)!=manifest['mesh_sha256'] or sha(exe)!=manifest['executable_sha256']:
        raise ValueError('Mesh or executable identity mismatch')
    if execution['executable_sha256']!=manifest['executable_sha256'] or execution['input_sha256']!=sha(directory/'input.nml'):
        raise ValueError('Input/execution identity mismatch')
    text=(directory/'input.nml').read_text(encoding='utf-8-sig')
    if input_value(text,'material_profile')!='ZNO_F1_NT' or input_value(text,'heat_source_current_work').lower()!='.true.':
        raise ValueError('Only actual NT JE fields are supported')
    if input_value(text,'transport_supg').lower()!='.true.' or input_value(text,'supg_directional_tau').lower()!='.true.':
        raise ValueError('Requires the existing directional SUPG branch')
    stated=Path(input_value(text,'mesh_file'))
    if not stated.is_absolute():stated=directory/stated
    if not os.path.samefile(stated,mesh):raise ValueError('Input points to a different mesh')
    if (metrics.get('model')!='NT' or metrics.get('heat_source')!=JE or metrics.get('supg_directional_tau')!='T'
        or metrics.get('transport_discretization')!='HEX20_CONSISTENT_SUPG_EXPERIMENTAL'
        or metrics.get('final_voltage_factor')!=1. or metrics.get('final_joule_factor')!=1.):
        raise ValueError('Wrong saved model or incomplete load')
    for key in ('target_voltage_factor','target_joule_factor'):
        if float(input_value(text,key))!=1.:raise ValueError('Nonunit target load')
    if metrics['voltage_V']!=abs(float(input_value(text,'target_voltage'))):raise ValueError('Voltage identity mismatch')
    if metrics['supg_scale']!=float(input_value(text,'supg_scale')) or metrics['supg_scale']!=manifest['supg_scale']:
        raise ValueError('SUPG scale identity mismatch')
    if not math.isclose(metrics['hole_radius_m'],float(input_value(text,'hole_radius')),rel_tol=1e-13):
        raise ValueError('Radius identity mismatch')
    build_path=exe.parent/'build_manifest.json';build=read(build_path)
    entries={r['path'].replace('\\','/'):r['sha256'].lower() for r in build['files']}
    if entries.get(exe.relative_to(ROOT).as_posix())!=sha(exe):raise ValueError('Build/executable identity mismatch')
    for source in SOURCES:
        if entries.get(source)!=sha(ROOT/source):raise ValueError('Native definition changed since this build: '+source)
    files=[directory/name for name in ['manifest.json','execution.json','summary.json','field_audit.json','input.nml','metrics.txt','nodes.csv']]
    files += [mesh,exe,build_path,Path(__file__),Path(fe.__file__)]+[ROOT/source for source in SOURCES]
    return manifest,metrics,mesh,{str(p):sha(p) for p in files}


def statistics(data,key):
    values=data[key];weights=data['weight_m3'];total=float(weights.sum())
    order=np.argsort(values,kind='stable');cumulative=np.cumsum(weights[order])/total
    def location(index):
        return dict(value=float(values[index]),sample_index_0based=int(index),
            element_1based=int(data['element_1based'][index]),gauss_point_1based=int(data['gauss_point_1based'][index]),
            xyz_m=data['xyz_m'][index].tolist(),natural_coordinates=data['natural_coordinates'][index].tolist())
    return dict(units=UNITS[key],minimum=location(int(np.argmin(values))),maximum=location(int(np.argmax(values))),
        volume_weighted_mean=float(weights@values/total),volume_weighted_rms=float(np.sqrt(weights@(values*values)/total)),
        volume_weighted_percentiles={str(q):float(np.interp(q/100,cumulative,values[order])) for q in (50,95,99)})


def audit(directory, output, order):
    manifest,metrics,mesh,hashes=evidence(directory)
    output.mkdir(parents=True,exist_ok=False)
    save(output/'provenance.json',dict(run=str(directory),order=order,hashes=hashes,
        original_scope='Saved solve/flux accepted but mesh and stabilization certification remains separate.'))
    xyz,conn=fe.read_mesh(mesh)
    nodes=np.loadtxt(Path(directory)/'nodes.csv',delimiter=',',skiprows=1)
    if nodes.shape!=(len(xyz),13) or not np.array_equal(nodes[:,0],np.arange(1,len(xyz)+1)):
        raise ValueError('Wrong saved node table dimensions/order')
    if not np.array_equal(nodes[:,1:4],xyz):raise ValueError('Node/mesh coordinates differ')
    if not np.isfinite(nodes).all():raise ValueError('Nonfinite saved state')
    state=nodes[:,4:11];points,weights=rule(order)
    pieces={key:[] for key in UNITS};ne=len(conn);nq=len(points)
    for element,numbers in enumerate(conn):
        item=sample(xyz[numbers],state[numbers],points,weights)
        for key in pieces:pieces[key].append(item[key])
    data={key:np.concatenate(values) for key,values in pieces.items()}
    data['element_1based']=np.repeat(np.arange(1,ne+1,dtype=np.int32),nq)
    data['gauss_point_1based']=np.tile(np.arange(1,nq+1,dtype=np.int32),ne)
    if len(data['Pe'])!=ne*order**3:raise AssertionError('Incomplete sampling')
    np.savez_compressed(output/'all_gauss_samples.npz',**data)
    scalar_keys=[key for key in UNITS if data[key].ndim==1 and key not in ('weight_m3','detJ_m3','g_direction_m2')]
    stats={key:statistics(data,key) for key in scalar_keys}
    report=dict(status='DIAGNOSTIC_ONLY_NO_ACCEPTANCE_THRESHOLD',run=str(directory),voltage_V=metrics['voltage_V'],
        radius_m=metrics['hole_radius_m'],supg_scale=metrics['supg_scale'],order=order,
        mesh_sha256=manifest['mesh_sha256'],executable_sha256=manifest['executable_sha256'],
        elements=ne,sample_count=len(data['Pe']),filtered_sample_count=0,clipped_sample_count=0,
        sampled_volume_m3=float(data['weight_m3'].sum()),minimum_jacobian_m3=float(data['detJ_m3'].min()),
        statistics=stats,E_component_ranges_V_m=[[float(data['E_V_m'][:,j].min()),float(data['E_V_m'][:,j].max())] for j in range(3)],
        zero_field_sample_count=int(np.count_nonzero(data['E_magnitude_V_m']==0.)),
        source_hashes_unchanged=all(sha(Path(p))==h for p,h in hashes.items()),
        all_samples_sha256=sha(output/'all_gauss_samples.npz'),
        all_samples_schema={key:dict(shape=list(value.shape),dtype=str(value.dtype),units=UNITS.get(key,'1-based index')) for key,value in data.items()},
        limitations=['Pe depends on mesh metric and local temperature; compare voltage effects on the same mesh and model.',
            'Gauss extrema are sampled diagnostics, not certified continuous extrema or monotonicity.',
            'The common scalar-Einstein Pe is identical for electrons/holes; their perturbations and gradients differ.',
            'This is not an error estimator, stability certificate, energy-model closure or production acceptance.'])
    if not report['source_hashes_unchanged']:raise ValueError('Source changed during sampling')
    save(output/'report.json',report)
    print(json.dumps(dict(run=Path(directory).name,order=order,samples=report['sample_count'],
        Pe_max=stats['Pe']['maximum']['value'],T_range_K=[stats['temperature_K'][k]['value'] for k in ('minimum','maximum')],
        n_perturbation_max=stats['n_abs_relative_perturbation']['maximum']['value'],
        p_perturbation_max=stats['p_abs_relative_perturbation']['maximum']['value'])),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',nargs='+',type=Path,required=True)
    p.add_argument('--order',type=int,choices=(3,5),default=3)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();output=a.output.resolve();runs=[d.resolve(strict=True) for d in a.runs]
    if len(set(d.name for d in runs))!=len(runs):p.error('Run directory names must be unique')
    if any(output==d or output.is_relative_to(d) for d in runs):p.error('New output must be outside source runs')
    output.mkdir(parents=True,exist_ok=False)
    save(output/'protocol.json',dict(scope='Read-only diagnostic; no new production acceptance threshold',
        order=a.order,runs=[str(d) for d in runs],script_sha256=sha(Path(__file__)),
        definition={'J':'dx/dxi stored as natural rows / physical columns; G=J^{-1} J^{-T}',
                    'epsilon_E':'1e-12*(kB*T/q)^2*trace(G)',
                    'gdir':'(E.G.E+epsilon_E*trace(G)/3)/(|E|^2+epsilon_E)',
                    'h':'2/sqrt(gdir), meters; matches native directional tau regularization',
                    'Pe':'q*|E|*h/(2*kB*T), common to n/p for scalar Einstein ZnO transport',
                    'zero_field':'E=0 gives Pe=0 and h=2*sqrt(3/trace(G)); retain zero-field points',
                    'carrier_perturbation':'expm1(eta); signed values and absolute amplitude retained',
                    'carrier_gradient':'grad(n/n0-1)=exp(eta_n)*grad(eta_n), similarly p; units 1/m',
                    'statistics':'all positive quadrature weights, volume weighted, no exclusions or clipping'},
        completeness='One row per element per tensor Gauss point; full vectors and scalar samples retained in NPZ',
        native_solver_started=False,production_acceptance=False))
    shutil.copy2(__file__,output/Path(__file__).name)
    save(output/'closed_form_check.json',closed_form_check())
    reports=[audit(d,output/d.name,a.order) for d in runs]
    keys=['Pe','n_abs_relative_perturbation','p_abs_relative_perturbation',
          'grad_n_relative_magnitude_per_m','grad_p_relative_magnitude_per_m','E_magnitude_V_m','temperature_K']
    envelope={key:max((dict(run=r['run'],voltage_V=r['voltage_V'],**r['statistics'][key]['maximum']) for r in reports),key=lambda x:x['value']) for key in keys}
    save(output/'summary.json',dict(status='DIAGNOSTIC_REFERENCE_ONLY_NOT_ACCEPTANCE',
        reports=[str(output/Path(r['run']).name/'report.json') for r in reports],
        descriptive_maxima=envelope,all_sources_unchanged=all(r['source_hashes_unchanged'] for r in reports),
        same_mesh=len({r['mesh_sha256'] for r in reports})==1,same_order=True,native_solver_started=False,
        production_acceptance=False,limitation='Endpoint sample maxima are descriptive references, not PASS limits for intermediate voltages.'))
    save(output/'artifact_hashes.json',{str(f.relative_to(output)):sha(f) for f in output.rglob('*') if f.is_file()})


if __name__=='__main__':main()
