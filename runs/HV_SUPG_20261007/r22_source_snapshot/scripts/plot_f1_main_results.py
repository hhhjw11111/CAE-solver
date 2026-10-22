"""Exact FE plane evaluation of saved F1 solutions, then source-backed figures.

No PDE solves, scattered-data smoothing, or z-mirror assumptions are used.
Each plotted triangle belongs to one original HEX20 face/plane patch.
"""
from pathlib import Path
import argparse,json,itertools
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from audit_hex20_boundary_geometry import NATURAL,derivatives,read_mesh
from audit_f1_results import read_metrics,ROOT as RUNS

SOLVER=Path(__file__).resolve().parents[1]
OUT=SOLVER.parent/'output/figures/F1_20260906'
INPUTS=SOLVER/'cases/F1_20260906'
PASSPORT=json.loads((SOLVER/'material_audit/zno_300k_candidate_v2.json').read_text(encoding='utf-8-sig'))
P=PASSPORT['parameters'];C=np.zeros((6,6));e=np.zeros((3,6))
C[0,0]=C[1,1]=P['elastic_Pa']['C11'];C[2,2]=P['elastic_Pa']['C33']
C[0,1]=C[1,0]=P['elastic_Pa']['C12']
C[0,2]=C[2,0]=C[1,2]=C[2,1]=P['elastic_Pa']['C13']
C[3,3]=C[4,4]=P['elastic_Pa']['C44'];C[5,5]=(C[0,0]-C[0,1])/2
e[0,4]=e[1,3]=P['piezo_C_per_m2']['e15'];e[2,:]=[P['piezo_C_per_m2']['e31']]*2+[P['piezo_C_per_m2']['e33']]+[0]*3
stress_T=C@np.array(P['alpha_per_K']);N0=P['n_H_per_m3'];P0=P['intrinsic_concentration_per_m3']**2/N0
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':7,'axes.titlesize':8,'legend.fontsize':7,'svg.fonttype':'none','pdf.fonttype':42,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'legend.frameon':False})

def shapes(points):
    pts=np.asarray(points);fac=1+pts[:,None,:]*NATURAL[None,:,:]
    out=np.empty((len(pts),20))
    out[:,:8]=np.prod(fac[:,:8,:],axis=2)*(np.sum(pts[:,None,:]*NATURAL[None,:8,:],axis=2)-2)/8
    for node in range(8,20):
        zero=int(np.flatnonzero(NATURAL[node]==0)[0]);others=[a for a in range(3) if a!=zero]
        out[:,node]=(1-pts[:,zero]**2)*fac[:,node,others[0]]*fac[:,node,others[1]]/4
    return out

def evaluate(coords,ue,points,model):
    shape=shapes(points);dn=derivatives(points)
    jac=np.einsum('pkn,nj->pkj',dn,coords)
    if np.any(np.linalg.det(jac)<=0):raise ValueError('Non-positive sampled geometry')
    grad=np.linalg.solve(jac,dn)
    vals=shape@ue;dvals=np.einsum('pkn,ns->pks',grad,ue)
    E=-dvals[:,:,3];theta=vals[:,6];T=300+theta
    n=1+vals[:,4] if model=='LQ' else np.exp(vals[:,4])
    p=1+vals[:,5] if model=='LQ' else np.exp(vals[:,5])
    mu_ratio=(T/300)**(-P['beta_n']) if model=='NT' else np.ones(len(T))
    mup_ratio=(T/300)**(-P['beta_p']) if model=='NT' else np.ones(len(T))
    Q=1.602176634e-19*((N0*(1 if model=='LQ' else n)*P['mu_H_m2_per_Vs']*mu_ratio)+
        (P0*(1 if model=='LQ' else p)*P['mu_p0_m2_per_Vs']*mup_ratio))*np.sum(E**2,axis=1)
    strain=np.column_stack([dvals[:,0,0],dvals[:,1,1],dvals[:,2,2],
        dvals[:,2,1]+dvals[:,1,2],dvals[:,0,2]+dvals[:,2,0],dvals[:,1,0]+dvals[:,0,1]])
    stress=strain@C.T-E@e-theta[:,None]*stress_T
    vm=np.sqrt(.5*((stress[:,0]-stress[:,1])**2+(stress[:,1]-stress[:,2])**2+
                  (stress[:,2]-stress[:,0])**2)+3*np.sum(stress[:,3:]**2,axis=1))
    return {'xyz':shape@coords,'phi':vals[:,3],'E':np.linalg.norm(E,axis=1),'Evec':E,
        'Q':Q,'theta':theta,'n':n,'p':p,'u':vals[:,:3],'mu':mu_ratio,'stress':stress,'vm':vm}

def load(case):
    m=read_metrics(case)
    if m.get('case_acceptance')!='PASS_NOT_MESH_ACCEPTANCE':raise ValueError('Case not accepted: '+case)
    geom='H1' if 'H1' in case else 'H0';level=next(t for t in case.split('_') if t in ['L0','L1','L2'])
    xyz,conn,_=read_mesh(INPUTS/f'{geom}_{level}.mesh')
    d=np.genfromtxt(RUNS/case/'nodes.csv',delimiter=',',names=True,usecols=range(13))
    assert len(d)==len(xyz)
    assert np.max(np.abs(xyz-np.column_stack([d[c] for c in ['x_m','y_m','z_m']])))<1e-19
    state=np.column_stack([d[c] for c in ['ux_m','uy_m','uz_m','phi_V','eta_n','eta_p','theta_K']])
    return m,xyz,conn,state

def check_saved_sample(case,m,xyz,conn,state):
    ref=np.loadtxt(RUNS/case/'field_samples.csv',delimiter=',',skiprows=1,max_rows=1)
    gp=np.polynomial.legendre.leggauss(5)[0][0]
    sampled=evaluate(xyz[conn[0]],state[conn[0]],[[gp,gp,gp]],m['model'])
    actual=np.r_[sampled['Evec'][0],sampled['Q'][0],sampled['theta'][0],sampled['n'][0],sampled['p'][0],
                  sampled['u'][0],sampled['mu'][0],sampled['stress'][0]]
    expected=ref[4:]
    # A vector scale avoids dividing a near-zero electric component by cancellation noise.
    scale=np.maximum(np.abs(expected),np.r_[[1e-3]*3,1e-3,1e-8,1,1,[1e-18]*3,1,[1e-4]*6])
    scale[:3]=max(np.linalg.norm(expected[:3]),1e-3)
    err=float(np.max(np.abs(actual-expected)/scale))
    if not err<2e-6:raise ValueError(f'FE plotting material/coordinate mismatch: {err:g}')
    return err

def plane(xyz,conn,state,model,z):
    positions=[];field_blocks={};triangles=[];offset=0
    t=np.linspace(-1,1,5)
    for nodes in conn:
        co=xyz[nodes];lo=co[:,2].min();hi=co[:,2].max()
        # Select exactly one side of the midplane, not an average of gradients.
        selected=(abs(lo-z)<1e-19) if z<0.5e-6 else (abs(hi-z)<1e-19)
        if not selected:continue
        zz=1 if z==0.5e-6 else -1
        points=np.array([[r,s,zz] for s in t for r in t])
        out=evaluate(co,state[nodes],points,model);positions.append(out.pop('xyz'))
        for k,v in out.items():field_blocks.setdefault(k,[]).append(v)
        for j in range(4):
            for i in range(4):
                a=offset+5*j+i;triangles += [[a,a+1,a+6],[a,a+6,a+5]]
        offset+=25
    if not positions:raise ValueError('No original mesh patch intersects requested plane')
    pos=np.concatenate(positions)
    out={k:np.concatenate(v) for k,v in field_blocks.items()}
    out['xyz']=pos;out['triangles']=np.array(triangles)
    return out

def save(fig,name):
    for ext in ['svg','pdf','png']:fig.savefig(OUT/(name+'.'+ext),dpi=260,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def panel(fig,ax,data,values,title,label,cmap='viridis',bounds=None):
    xy=data['xyz'][:,:2]*1e6;tri=mtri.Triangulation(xy[:,0],xy[:,1],data['triangles'])
    vmin,vmax=(float(np.min(values)),float(np.max(values))) if bounds is None else bounds
    if vmax-vmin<1e-14:maxabs=max(abs(vmax),1);vmin-=maxabs*1e-8;vmax+=maxabs*1e-8
    # Gouraud uses only within-original-element linear triangles, with no cross-hole links.
    artist=ax.tripcolor(tri,values,shading='gouraud',cmap=cmap,vmin=vmin,vmax=vmax,rasterized=True)
    ax.set(aspect='equal',xlim=(-10,10),ylim=(-10,10),xticks=[-10,0,10],yticks=[-10,0,10],
           xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)',title=title)
    cb=fig.colorbar(artist,ax=ax,fraction=.05,pad=.025);cb.set_label(label)

def main():
    global RUNS,OUT,INPUTS
    parser=argparse.ArgumentParser();parser.add_argument('--sample-check',nargs=1)
    parser.add_argument('--run-root',type=Path,default=RUNS)
    parser.add_argument('--output-dir',type=Path,default=OUT)
    parser.add_argument('--input-root',type=Path,default=INPUTS)
    args=parser.parse_args()
    RUNS=args.run_root.resolve();OUT=args.output_dir.resolve();INPUTS=args.input_root.resolve()
    import audit_f1_results
    audit_f1_results.ROOT=RUNS
    OUT.mkdir(parents=True,exist_ok=True)
    if args.sample_check:
        case=args.sample_check[0];m,xyz,conn,state=load(case)
        print(json.dumps({'case':case,'saved_sample_relative_error':check_saved_sample(case,m,xyz,conn,state)}));return
    mesh_check=json.loads((RUNS/'mesh_comparison.json').read_text())
    if not mesh_check['passed']:raise ValueError('Mesh comparison has not passed')
    case='MESH_H1_NT_L2';m,xyz,conn,state=load(case)
    qa={case:check_saved_sample(case,m,xyz,conn,state)}
    mid=plane(xyz,conn,state,'NT',0);top=plane(xyz,conn,state,'NT',0.5e-6)
    bottom=plane(xyz,conn,state,'NT',-0.5e-6)
    fig,axes=plt.subplots(2,3,figsize=(7.2,4.8),layout='constrained')
    specs=[(mid,mid['phi'],'Potential, midplane','V','viridis'),
       (mid,mid['E']/1e5,'Electric field, midplane','kV/cm','viridis'),
       (mid,mid['Q']/1e12,'Joule heat, midplane',r'$10^{12}$ W m$^{-3}$','magma'),
       (mid,mid['theta'],'Temperature rise, midplane','K','magma'),
       (mid,mid['u'][:,2]*1e9,'Axial displacement, midplane','nm','RdBu_r'),
       (top,100*(top['n']-1),'Electrons, top surface','Change (%)','RdBu_r')]
    for ax,spec,lbl in zip(axes.flat,specs,'abcdef'):
        panel(fig,ax,*spec);ax.text(-.15,1.13,lbl,transform=ax.transAxes,fontweight='bold')
    save(fig,'H1_NT_multifield')
    np.savez_compressed(OUT/'H1_NT_midplane_source.npz',**mid)
    np.savez_compressed(OUT/'H1_NT_top_source.npz',**top)
    np.savez_compressed(OUT/'H1_NT_bottom_source.npz',**bottom)
    comparisons=[]
    for model,case in [('LQ','MAIN_H1_LQ_L2'),('N0','MAIN_H1_N0_L2'),('NT','MESH_H1_NT_L2')]:
        mm,xx,cc,ss=load(case);qa[case]=check_saved_sample(case,mm,xx,cc,ss)
        comparisons.append((model,plane(xx,cc,ss,model,0.5e-6),plane(xx,cc,ss,model,-0.5e-6)))
    limit=max(np.max(np.abs(d['n']-1))*100 for _,top,bottom in comparisons for d in [top,bottom])
    fig,axes=plt.subplots(2,3,figsize=(7.2,4.8),layout='constrained')
    for col,(model,top,bottom) in enumerate(comparisons):
        for row,(label,data) in enumerate([('top',top),('bottom',bottom)]):
            panel(fig,axes[row,col],data,100*(data['n']-1),f'{model}, {label} surface','Change (%)','RdBu_r',(-limit,limit))
    for ax,lbl in zip(axes.flat,'abcdef'):ax.text(-.15,1.13,lbl,transform=ax.transAxes,fontweight='bold')
    save(fig,'H1_model_carrier_comparison')
    qa['conclusion']='Main nonlinear fields and carrier model contrasts, conditioned on the frozen synthetic/effective material assumptions.'
    qa['planes']='Exact FE evaluations at z=0 and z=+/-0.5 um. Gradients use one-sided element traces; no z symmetry is imposed.'
    qa['integrity']='Field triangles remain within original HEX20 patches; no cross-hole interpolation; common carrier color limits across models and both faces.'
    (OUT/'main_figure_qa.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
    print(OUT)

if __name__=='__main__':main()
