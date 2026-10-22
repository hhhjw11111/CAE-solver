"""Read-only, bounded-memory HEX20 diagnostics of completed 3 V evidence.
Run: python -B solver_v2/paper_sprint/20261008_2110/scripts/postprocess_3v.py
No solver invocation. Outputs belong solely to this sprint.
"""
from pathlib import Path
import sys, os, json, csv, hashlib, itertools, re, time
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parents[1]; SOLVER=OUT.parents[1]
sys.path.insert(0,str(OUT/'scripts'))
from hex20_geometry import read_mesh, derivatives, NATURAL
Q=1.602176634e-19; KB=1.380649e-23; N0=6e22; P0=(2.9e-4)**2/N0
SIGMA=Q*(N0*.0205+P0*.005); K=60.
for d in ['data','data/profiles','figures','analysis','audit/postprocess']: (OUT/d).mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'], 'font.size':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def dump(p,x): p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def csvwrite(p,rows):
    if not rows:return
    with p.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def metrics(p):
    result={}
    for line in p.read_text(encoding='utf-8-sig').splitlines():
        if '=' in line:
            k,v=line.split('=',1)
            try: result[k]=float(v)
            except ValueError: result[k]=v
    return result
def shapes(points):
    points=np.asarray(points);f=1+points[:,None,:]*NATURAL[None,:,:]
    n=np.zeros((len(points),20))
    n[:,:8]=np.prod(f[:,:8],axis=2)*(np.sum(points[:,None,:]*NATURAL[None,:8,:],axis=2)-2)/8
    for a in range(8,20):
        z=int(np.flatnonzero(NATURAL[a]==0)[0]);other=[k for k in range(3) if k!=z]
        n[:,a]=(1-points[:,z]**2)*f[:,a,other[0]]*f[:,a,other[1]]/4
    return n

elastic=np.zeros((6,6));elastic[0,0]=elastic[1,1]=209.7e9;elastic[2,2]=210.9e9
elastic[0,1]=elastic[1,0]=121.1e9
elastic[0,2]=elastic[2,0]=elastic[1,2]=elastic[2,1]=105.1e9
elastic[3,3]=elastic[4,4]=42.47e9;elastic[5,5]=44.3e9
piezo=np.zeros((3,6));piezo[0,4]=piezo[1,3]=-.37;piezo[2,:3]=[-.62,-.62,.96]
thermal=elastic@np.array([4.31e-6,4.31e-6,2.49e-6,0,0,0])

def fields(co,ue,points,model):
    # co: elements x 20 x 3; ue: elements x 20 x 7.
    n=shapes(points);dn=derivatives(points)
    jac=np.einsum('qka,eaj->eqkj',dn,co);det=np.linalg.det(jac)
    if np.any(det<=0) or not np.isfinite(det).all():raise ValueError('Nonpositive sampled Jacobian')
    grad=np.einsum('eqjk,qka->eqja',np.linalg.inv(jac),dn)
    val=np.einsum('qa,eaf->eqf',n,ue);g=np.einsum('eqja,eaf->eqjf',grad,ue)
    xyz=np.einsum('qa,eaj->eqj',n,co);E=-g[...,3];T=300+val[...,6]
    nr=1+val[...,4] if model=='LQ' else np.exp(val[...,4]);pr=1+val[...,5] if model=='LQ' else np.exp(val[...,5])
    if np.any(T<=0) or min(nr.min(),pr.min())<=0:raise ValueError('Invalid physical field')
    mu=(T/300)**-2 if model=='NT' else np.ones_like(T)
    sn=Q*N0*.0205*mu*(1 if model=='LQ' else nr);sp=Q*P0*.005*mu*(1 if model=='LQ' else pr)
    vt=KB*(T if model=='NT' else np.full_like(T,300.))/Q
    diff=vt[...,None]*(sn[...,None]*g[...,4]-sp[...,None]*g[...,5])
    drift=(sn+sp)[...,None]*E;pd=(sn+sp)*np.sum(E*E,axis=-1);pdiff=np.sum(diff*E,axis=-1)
    strain=np.stack([g[...,0,0],g[...,1,1],g[...,2,2],g[...,2,1]+g[...,1,2],g[...,0,2]+g[...,2,0],g[...,1,0]+g[...,0,1]],axis=-1)
    stress=strain@elastic.T-E@piezo-val[...,6,None]*thermal
    vm=np.sqrt(.5*((stress[...,0]-stress[...,1])**2+(stress[...,1]-stress[...,2])**2+(stress[...,2]-stress[...,0])**2)+3*np.sum(stress[...,3:]**2,axis=-1))
    return dict(val=val,g=g,xyz=xyz,E=E,T=T,nr=nr,pr=pr,mu=mu,pd=pd,pdiff=pdiff,J=drift+diff,stress=stress,vm=vm,det=det)

def inverse_samples(xyz,conn,state,targets,model):
    co=xyz[conn];lo=co.min(1)-1e-18;hi=co.max(1)+1e-18;rows=[]
    for target in targets:
        ids=np.flatnonzero(np.all((target>=lo)&(target<=hi),axis=1));xs=co[ids];xi=np.zeros((len(ids),3))
        if not len(ids):raise ValueError('No element at material target '+str(target))
        for _ in range(15):
            rr=np.einsum('pa,paj->pj',shapes(xi),xs)-target
            jj=np.einsum('pka,paj->pjk',derivatives(xi),xs)
            xi-=np.clip(np.linalg.solve(jj,rr[...,None])[...,0],-1.,1.)
        error=np.linalg.norm(np.einsum('pa,paj->pj',shapes(xi),xs)-target,axis=1)
        valid=(error<1e-17)&(np.max(abs(xi),axis=1)<=1+1e-8)
        if not valid.any():raise ValueError('No non-extrapolated inverse map at '+str(target))
        j=int(np.flatnonzero(valid)[0]);eid=int(ids[j]);f=fields(co[eid:eid+1],state[conn[eid:eid+1]],xi[j:j+1],model)
        v=f['val'][0,0];s=f['stress'][0,0];e=f['E'][0,0];r=np.hypot(target[0],target[1]);c=target[0]/r;ss=target[1]/r
        rows.append(dict(x_m=float(target[0]),y_m=float(target[1]),z_m=float(target[2]),element=eid+1,inverse_residual_m=float(error[j]),
            phi_V=float(v[3]),eta_n=float(v[4]),eta_p=float(v[5]),n_relative=float(f['nr'][0,0]-1),p_relative=float(f['pr'][0,0]-1),theta_K=float(v[6]),
            Ex_V_m=float(e[0]),Ey_V_m=float(e[1]),Ez_V_m=float(e[2]),ux_m=float(v[0]),uy_m=float(v[1]),uz_m=float(v[2]),
            srr_Pa=float(c*c*s[0]+ss*ss*s[1]+2*c*ss*s[5]),stt_Pa=float(ss*ss*s[0]+c*c*s[1]-2*c*ss*s[5]),vm_Pa=float(f['vm'][0,0])))
    return rows

F1=SOLVER/'runs/F1_FAST_HOLE_L0_20260907';HV=SOLVER/'runs/HV_SUPG_20261007'
CASES=[]
for geom,model in itertools.product(['H0','H1'],['LQ','N0','NT']):
    name=f'FAST_{geom}_{model}_L0'; CASES.append(dict(name=name,root=F1/name,mesh=F1/f'evidence/inputs/{geom}_L0.mesh',model=model,geom=geom,family='F1_L0',half=False,exe=F1/'evidence/TPSolverV2.exe',source=F1/'evidence/src'))
for level in ['L0','L1','L2']:
    name=f'r23_log_half_{level}_V030'; rr=HV/name; mm=json.loads((rr/'manifest.json').read_text())
    CASES.append(dict(name=name,root=rr,mesh=Path(mm['mesh']),model='NT',geom='H1',family='R23_'+level,half=True,exe=Path(mm['executable']),source=HV/'r23_source_snapshot/src'))
rr=HV/'r25_followthrough/V030';mm=json.loads((rr/'manifest.json').read_text())
CASES.append(dict(name='r25_A64_V030',root=rr,mesh=Path(mm['mesh']),model='NT',geom='H1',family='R25_A64',half=True,exe=Path(mm['executable']),source=HV/'r25_source_snapshot/src'))

def analyze(case):
    root=case['root'];name=case['name'];m=metrics(root/'metrics.txt');print('START',name,flush=True)
    xyz,conn,_=read_mesh(case['mesh']);d=np.genfromtxt(root/'nodes.csv',delimiter=',',names=True)
    if len(d)!=len(xyz) or np.max(abs(xyz-np.column_stack([d[k] for k in ['x_m','y_m','z_m']])))>1e-19:raise ValueError('mesh/field coordinate mismatch')
    if not np.array_equal(d['node'],np.arange(1,len(d)+1)):raise ValueError('nodes ID mismatch')
    state=np.column_stack([d[k] for k in ['ux_m','uy_m','uz_m','phi_V','eta_n','eta_p','theta_K']]);del d
    if m['voltage_V']!=3 or m['case_acceptance']!='PASS_NOT_MESH_ACCEPTANCE':raise ValueError('Not accepted 3 V algebraic/flux case')
    hashes={str(p.resolve()):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in [root/'input.nml',case['mesh'],root/'nodes.csv',root/'metrics.txt',root/'field_samples.csv',case['exe']]}
    manifest=json.loads((root/'manifest.json').read_text()) if (root/'manifest.json').exists() else None
    if manifest:
        assert hashes[str(case['exe'].resolve())]['sha256']==manifest['executable_sha256']
        assert hashes[str(case['mesh'].resolve())]['sha256']==manifest['mesh_sha256']
    else:
        build=json.loads((F1/'evidence/build_manifest.json').read_text())
        assert hashes[str(case['exe'].resolve())]['sha256']==next(x['sha256'] for x in build['files'] if x['path'].endswith('TPSolverV2.exe'))
    source_records=[]
    for p in sorted(case['source'].rglob('*.f90')):
        source_records.append(dict(path=str(p.resolve()),sha256=sha(p)))
    if not source_records:raise ValueError('Source snapshot missing')
    gp,gw=np.polynomial.legendre.leggauss(5);points=np.array([(x,y,z) for z in gp for y in gp for x in gp]);weights=np.array([x*y*z for z in gw for y in gw for x in gw])
    sample=np.loadtxt(root/'field_samples.csv',delimiter=',',skiprows=1,max_rows=1)
    ff=fields(xyz[conn[:1]],state[conn[:1]],points[:1],case['model']);vv=ff['val'][0,0]
    actual=np.r_[ff['E'][0,0],ff['pd'][0,0],vv[6],ff['nr'][0,0],ff['pr'][0,0],vv[:3],ff['mu'][0,0],ff['stress'][0,0]]
    scale=np.maximum(abs(sample[4:]),np.r_[[1e-3]*3,1e-3,1e-8,1,1,[1e-18]*3,1,[1e-4]*6]);scale[:3]=max(np.linalg.norm(sample[4:7]),1e-3)
    reconstruction=float(np.max(abs(actual-sample[4:])/scale))
    if reconstruction>2e-6:raise ValueError('Saved material sample mismatch '+str(reconstruction))
    phiL=float(np.mean(state[np.isclose(xyz[:,0],-1e-5,atol=1e-16,rtol=0),3]));phiR=float(np.mean(state[np.isclose(xyz[:,0],1e-5,atol=1e-16,rtol=0),3]))
    vol=np.zeros(len(conn));sumv=np.zeros(7);refmax=0.;Elist=[];Wlist=[]
    for start in range(0,len(conn),24):
        cs=conn[start:start+24];f=fields(xyz[cs],state[cs],points,case['model']);w=f['det']*weights
        reference=SIGMA/(2*K)*(f['val'][...,3]-phiL)*(phiR-f['val'][...,3]);resid=f['val'][...,6]-reference
        quantities=[f['pd'],f['pdiff'],abs(f['pdiff']),f['pd']+f['pdiff'],resid**2,resid,np.ones_like(resid)]
        sumv += [float(np.sum(w*v)) for v in quantities];vol[start:start+len(cs)]=w.sum(1)
        refmax=max(refmax,float(abs(resid).max()));Elist.append(np.linalg.norm(f['E'],axis=-1).ravel());Wlist.append(w.ravel())
    reconheat=abs(sumv[0]-m['joule_W'])/abs(m['joule_W'])
    if reconheat>1e-9:raise ValueError('Heat integral differs '+str(reconheat))
    Es=np.concatenate(Elist);Ws=np.concatenate(Wlist);order=np.argsort(Es);E95=float(Es[order][np.searchsorted(np.cumsum(Ws[order]),.95*Ws.sum())]);del Elist,Wlist,Es,Ws
    center=shapes([[0,0,0]])[0]@xyz[conn];drop=np.ptp(state[conn,3],axis=1);te=300+state[conn,6].mean(1)
    pe_local=drop/(2*KB*te/Q);pe=pe_local if case['model']=='NT' else drop/(2*KB*300/Q)
    region=np.full(len(conn),'bulk',dtype=object)
    region[np.linalg.norm(center[:,:2],axis=1)<=3.5e-6]='hole_neighbourhood' if case['geom']=='H1' else 'same_radius_control'
    region[abs(center[:,2])>=.45e-6]='surface_layers';region[np.max(abs(xyz[conn,0]),axis=1)>=1e-5-1e-16]='contact'
    per=[]
    for r in ['all',*sorted(set(region))]:
        ix=np.ones(len(conn),bool) if r=='all' else region==r
        per.append(dict(case=name,region=r,count=int(ix.sum()),volume_m3=float(vol[ix].sum()),Pe_median=float(np.median(pe[ix])),Pe_p95=float(np.quantile(pe[ix],.95)),Pe_max=float(pe[ix].max()),count_fraction_Pe_gt_1=float(np.mean(pe[ix]>1)),volume_fraction_Pe_gt_1=float(vol[ix&(pe>1)].sum()/vol[ix].sum()),count_fraction_Pe_gt_2=float(np.mean(pe[ix]>2)),volume_fraction_Pe_gt_2=float(vol[ix&(pe>2)].sum()/vol[ix].sum())))
    erows=[dict(element=i+1,x_m=float(center[i,0]),y_m=float(center[i,1]),z_m=float(center[i,2]),volume_m3=float(vol[i]),region=region[i],Pe_model_ratio=float(pe[i]),Pe_localT_background=float(pe_local[i])) for i in range(len(conn))]
    csvwrite(OUT/f'data/{name}_element_peclet.csv',erows)
    z=np.unique(np.r_[np.linspace(-.5e-6,.5e-6,401),xyz[:,2]])
    prof=inverse_samples(xyz,conn,state,np.column_stack([np.zeros(len(z)),np.full(len(z),3.25e-6),z]),case['model'])
    mid=min(prof,key=lambda p:abs(p['z_m']));midu=np.array([mid[k] for k in ['ux_m','uy_m','uz_m']])
    for row in prof:
        row['u_relative_midplane_m']=float(np.linalg.norm(np.array([row[k] for k in ['ux_m','uy_m','uz_m']])-midu))
        row['theta_reference_K']=SIGMA/(2*K)*(row['phi_V']-phiL)*(phiR-row['phi_V']);row['rT_K']=row['theta_K']-row['theta_reference_K']
    csvwrite(OUT/f'data/profiles/{name}_thickness.csv',prof)
    angles=np.linspace(.25,179.75,121) if case['half'] else np.linspace(.25,359.75,241)
    ang=np.radians(angles);ring=inverse_samples(xyz,conn,state,np.c_[3.25e-6*np.cos(ang),3.25e-6*np.sin(ang),np.zeros(len(ang))],case['model'])
    for row,angle in zip(ring,angles):row['angle_degree']=float(angle)
    csvwrite(OUT/f'data/profiles/{name}_ring_stress.csv',ring)
    result=dict(case=name,family=case['family'],model=case['model'],geometry=case['geom'],voltage_V=3.,integral_scope='ACTUAL_HALF_DOMAIN' if case['half'] else 'FULL_DOMAIN',nodes=len(xyz),elements=len(conn),
        volume_m3=float(sumv[6]),Pd_W=float(sumv[0]),Pdiff_W=float(sumv[1]),Pdiff_abs_W=float(sumv[2]),P_JdotE_W=float(sumv[3]),Pterm_all_boundary_W=m['terminal_power_W'],Pterm_physical_exterior_W=m.get('physical_exterior_terminal_power_W',m['terminal_power_W']),heat_out_all_boundary_W=m['heat_out_W'],heat_out_physical_exterior_W=m.get('physical_exterior_heat_out_W',m['heat_out_W']),
        power_identity_gap_relative=float(abs(sumv[3]-m['terminal_power_W'])/abs(m['terminal_power_W'])),Pdiff_abs_over_Pd=float(sumv[2]/sumv[0]),reference_rms_K=float(np.sqrt(sumv[4]/sumv[6])),reference_mean_K=float(sumv[5]/sumv[6]),reference_maxabs_K=refmax,theta_max_saved_K=m['theta_max_K'],E95_volume_weighted_V_m=E95,
        thickness_change_m=prof[-1]['uz_m']-prof[0]['uz_m'],relative_displacement_line_max_m=max(x['u_relative_midplane_m'] for x in prof),ring_midplane_max_abs_stt_Pa=max(abs(x['stt_Pa']) for x in ring),ring_midplane_max_vm_Pa=max(x['vm_Pa'] for x in ring),
        saved_sample_scaled_error=reconstruction,drift_heat_reconstruction_relative=reconheat,heat_flux_gap=m.get('physical_exterior_heat_balance_relative',m['heat_balance_relative']),Pe_max=float(pe.max()),Pe_count_gt1=float(np.mean(pe>1)),Pe_volume_gt1=float(vol[pe>1].sum()/vol.sum()),
        status='EXPLORATORY_L0' if case['family']=='F1_L0' else 'MESH_PENDING_LOCAL_FAIL_3_PERCENT' if case['family'].startswith('R23') else 'SOLVE_FLUX_PASS_MESH_PENDING')
    dump(OUT/f'audit/postprocess/{name}.json',dict(metrics=result,original_metrics=m,paths={k:str(case[k].resolve()) for k in ['root','mesh','exe','source']},hashes=hashes,source_hashes=source_records,pe_definition='Node range / (2 D/mu); NT Tmean, N0/LQ T0. LQ is background diagnostic only.',quadrature_order=5))
    print('DONE',name,'Pd',sumv[0],'diff',sumv[1],flush=True)
    return result,per

def savefig(fig,name):
    for ext in ['png','svg','pdf']:fig.savefig(OUT/f'figures/{name}.{ext}',dpi=250,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def readcsv(p):return np.genfromtxt(p,delimiter=',',names=True,encoding='utf-8-sig')
def figures(results,regions):
    old=[r for r in results if r['family']=='F1_L0'];colors={'LQ':'#888888','N0':'#3974a6','NT':'#c76732'}
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.3),layout='constrained')
    for r in old:
        p=readcsv(OUT/f"data/profiles/{r['case']}_thickness.csv");label=r['geometry']+' '+r['model'];sty='--' if r['geometry']=='H0' else '-'
        for ax,key,fac,ylabel in zip(axs.flat,['n_relative','theta_K','Ez_V_m','rT_K'],[100,1,1e-5,1],[r'Electron perturbation (%)',r'Temperature rise (K)',r'$E_z$ ($10^5$ V/m)',r'$\theta-\theta_{ref}$ (K)']):
            ax.plot(p['z_m']*1e9,p[key]*fac,sty,color=colors[r['model']],lw=1,label=label);ax.set(xlabel='z (nm)',ylabel=ylabel)
    axs[0,0].legend(ncol=2,fontsize=7);fig.suptitle('3 V F1 L0 exploratory; mesh acceptance incomplete\nMaterial line: x=0, y=3.25 µm (a+d; d=0.25 µm)',fontsize=9);savefig(fig,'fig06_thickness_profiles')
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.2),layout='constrained')
    names=[r['geometry']+' '+r['model'] for r in old];xx=np.arange(len(old))
    for j,(key,label) in enumerate([('Pd_W','Drift heat'),('Pdiff_W','Signed diffusion work'),('Pdiff_abs_W','Absolute diffusion work'),('Pterm_all_boundary_W','Terminal power')]):axs[0].bar(xx+(j-1.5)*.19,[r[key]*1e3 for r in old],.18,label=label,color=['#3974a6','#c76732','#c6a567','#555555'][j])
    axs[0].set(xticks=xx,xticklabels=names,ylabel='Power (mW)',ylim=(-.2,2.65));axs[0].tick_params(axis='x',rotation=45);axs[0].legend(fontsize=6)
    axs[1].bar(xx,[r['reference_rms_K'] for r in old],color=[colors[r['model']] for r in old]);axs[1].set(xticks=xx,xticklabels=names,ylabel='Reference residual RMS (K)');axs[1].tick_params(axis='x',rotation=45)
    fig.suptitle('3 V F1 L0 exploratory; mesh pending / reduced drift-heat model',fontsize=9);savefig(fig,'fig08_power_and_reference')
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.2),layout='constrained')
    for r in old:
        p=readcsv(OUT/f"data/profiles/{r['case']}_ring_stress.csv");pr=readcsv(OUT/f"data/profiles/{r['case']}_thickness.csv");style='--' if r['geometry']=='H0' else '-';lab=r['geometry']+' '+r['model']
        axs[0].plot(p['angle_degree'],p['stt_Pa']/1e6,style,color=colors[r['model']],label=lab,lw=1)
        axs[1].plot(pr['z_m']*1e9,pr['u_relative_midplane_m']*1e9,style,color=colors[r['model']],lw=1)
    axs[0].set(xlabel='Angle (degree)',ylabel=r'Unsmoothed $\sigma_{\theta\theta}$ (MPa)');axs[1].set(xlabel='z (nm)',ylabel='Displacement relative to line midpoint (nm)');axs[0].legend(ncol=2,fontsize=6)
    fig.suptitle('3 V F1 L0 exploratory; stress at r=3.25 µm, z=0\nElement-side traces; no nodal stress smoothing; mesh pending',fontsize=9);savefig(fig,'fig07_mechanics')
    selected=['FAST_H1_NT_L0','r23_log_half_L2_V030','r25_A64_V030'];regionsel=['contact','surface_layers','hole_neighbourhood','bulk'];fig,axs=plt.subplots(1,2,figsize=(7.2,3.3),layout='constrained')
    for j,c in enumerate(selected):
        rs={r['region']:r for r in regions if r['case']==c};x=np.arange(4)+(j-1)*.25
        for ax,key in zip(axs,['count_fraction_Pe_gt_1','volume_fraction_Pe_gt_1']):ax.bar(x,[100*rs.get(k,{}).get(key,np.nan) for k in regionsel],.24,label=c.replace('FAST_H1_','').replace('log_half_',''),color=['#3974a6','#999999','#c76732'][j])
    for ax in axs:ax.set(xticks=np.arange(4),xticklabels=['contact','surface','near hole','bulk'],ylim=(0,105));ax.tick_params(axis='x',rotation=30)
    axs[0].set_ylabel('Element count with Pe > 1 (%)');axs[1].set_ylabel('Volume in elements with Pe > 1 (%)');axs[0].legend(fontsize=6)
    fig.suptitle('3 V Pe diagnostic (node voltage range); not an error bound\nF1 L0 exploratory; R23/R25 mesh pending; R23 local carrier FAIL retained',fontsize=9);savefig(fig,'fig03_peclet_regions')
    fshape=json.loads((HV/'field_shape_r23_L1_L2_V030.json').read_text())['comparisons'];fig,ax=plt.subplots(figsize=(6.4,3.1),layout='constrained')
    labels=list(fshape);vals=[100*fshape[k]['maximum_relative'] for k in labels];ax.bar(np.arange(len(vals)),vals,color=['#c76732' if v>2 else '#3974a6' for v in vals]);ax.axhline(2,color='#a22',ls='--',label='Frozen limit 2%');ax.set(xticks=np.arange(len(vals)),xticklabels=['phi','theta','ux','uy','uz','n perturb.','p perturb.'],ylabel='Mapped maximum relative difference (%)');ax.legend(fontsize=7);ax.set_title('3 V R23 L1 → L2: local carrier FAIL; mesh certificate incomplete',fontsize=9);savefig(fig,'fig03b_spatial_failure')

def main():
    start=time.time();results=[];regions=[];errors=[]
    for case in CASES:
        try:r,p=analyze(case);results.append(r);regions.extend(p)
        except Exception as e:
            errors.append(dict(case=case['name'],error=repr(e)));print('BLOCKED',case['name'],repr(e),flush=True)
        csvwrite(OUT/'data/metrics_summary.csv',results);csvwrite(OUT/'data/peclet_regions.csv',regions);dump(OUT/'audit/postprocess/errors.json',errors)
    figures(results,regions)
    registry=[]
    for r in results:
        registry.append(dict(case=r['case'],model=r['model'],voltage_V=3,geometry=r['geometry'],integral_scope=r['integral_scope'],L1_data='HASHED_AND_COORDINATES_MATCH',L2_algebraic='EXISTING_PASS_RECHECKED_METRICS',L3_flux='EXISTING_PASS_RECHECKED_METRICS',L4_quadrature='5_POINT_POSTPROCESS_RECONSTRUCTION_PASS',L5_spatial=r['status'],L6_model='REDUCED_DRIFT_HEAT_ONLY',L7_experiment='NOT_CALIBRATED',source_path=str(next(c['root'].resolve() for c in CASES if c['name']==r['case']))))
    csvwrite(OUT/'data/case_registry.csv',registry);dump(OUT/'data/metrics_summary.json',dict(scope='Postprocessing only, no new 3D solves',elapsed_seconds=time.time()-start,results=results,blocked=errors));print('COMPLETE',time.time()-start,flush=True)
if __name__=='__main__':main()
