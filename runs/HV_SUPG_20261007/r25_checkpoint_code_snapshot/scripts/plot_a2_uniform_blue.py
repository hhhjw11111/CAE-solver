"""Uniform, interleaved analytic/FEM markers from saved A2 HEX20 fields."""
from pathlib import Path
import csv,json,shutil
import numpy as np
import matplotlib.pyplot as plt
import plot_a2_skill_delivery as b

OUT=b.ROOT/'output/A2_解析解FEM_蓝色统一间距_20260909'
COLOR='#2676D2'
OUT.mkdir(parents=True,exist_ok=True)
b.OUT=OUT
nodes=np.genfromtxt(b.RUN/'A2_L2_nodes.csv',delimiter=',',names=True)
coords,conn,err=b.reconstruct_l2_connectivity(nodes)
phi=nodes['phi_V']; a=b.A_UM*1e-6; half=b.L_UM*1e-6
records=[]
def evaluate(element,r,s):
    n,d=b.hex20_shape_numpy(r,s,0.)
    co=coords[conn[element]]; v=phi[conn[element]]
    field=-np.linalg.solve(d@co,d)@v
    return float(n@v*1e3),float(field@field/b.E0**2)

def center(xi):
    x=(2*xi-1)*half; et=0 if x>0 else 16
    for er in range(4):
        element=et*4+er; co=coords[conn[element]]
        low=b.hex20_shape_numpy(-1.,-1.,0.)[0]@co
        high=b.hex20_shape_numpy(1.,-1.,0.)[0]@co
        r=2*(x-low[0])/(high[0]-low[0])-1
        if -1-1e-9<=r<=1+1e-9:
            return evaluate(element,float(np.clip(r,-1,1)),-1.)
    raise ValueError(xi)

def wall(xi):
    # Invert the actual curved HEX20 face map to the requested polar angle.
    angle=(xi%1)*2*np.pi
    et=int(np.floor((xi%1)*32))%32
    co=coords[conn[et*4]]
    lower,upper=-1.,1.
    for _ in range(55):
        s=(lower+upper)/2
        xyz=b.hex20_shape_numpy(-1.,s,0.)[0]@co
        theta=float(np.arctan2(xyz[1],xyz[0])%(2*np.pi))
        if theta<angle: lower=s
        else: upper=s
    s=-1. if abs(angle-et*2*np.pi/32)<1e-12 else (lower+upper)/2
    xyz=b.hex20_shape_numpy(-1.,s,0.)[0]@co
    angle_error=np.angle(np.exp(1j*(np.arctan2(xyz[1],xyz[0])-angle)))
    assert abs(angle_error)<1e-12
    return evaluate(et*4,-1.,s)

fig,axs=plt.subplots(2,2,figsize=(183/25.4,143/25.4))
fig.subplots_adjust(left=.095,right=.98,bottom=.105,top=.95,wspace=.32,hspace=.40)
titles=['(a) 中心线电势','(b) 孔壁电势','(c) 中心线焦耳热源','(d) 孔壁焦耳热源']
labels=[r'中心线电势 $\phi$（mV）',r'孔壁电势 $\phi$（mV）',r'中心线焦耳热源 $Q_J/Q_\infty$',r'孔壁焦耳热源 $Q_J/Q_\infty$']
for i,ax in enumerate(axs.ravel()):
    is_center=i%2==0; quantity=0 if i<2 else 1
    intervals=[(0.,.35,8),(.65,1.,8)] if is_center else [(0.,1.,21)]
    for segment,(lo,hi,count) in enumerate(intervals):
        xf=np.linspace(lo,hi,count); xm=(xf[:-1]+xf[1:])/2; dense=np.linspace(lo,hi,401)
        def exact(t):
            if is_center:
                x=(2*t-1)*half
                return -b.E0*(x+a*a/x)*1e3 if quantity==0 else (1-a*a/x**2)**2
            return -2*b.E0*a*np.cos(2*np.pi*t)*1e3 if quantity==0 else 4*np.sin(2*np.pi*t)**2
        yf=np.array([(center(t) if is_center else wall(t))[quantity] for t in xf])
        ax.plot(dense,exact(dense),color=COLOR,lw=1)
        ax.plot(xm,exact(xm),'o',color=COLOR,ms=2.8,label='解析解' if segment==0 else None)
        ax.plot(xf,yf,linestyle='none',marker='^',mfc='white',mec=COLOR,mew=.85,ms=4.2,label='FEM（L2）' if segment==0 else None)
        for role,xx,yy in [('analytic_curve',dense,exact(dense)),('analytic_marker',xm,exact(xm)),('FEM_marker',xf,yf)]:
            records.extend([i+1,role,float(x),float(y)] for x,y in zip(xx,yy))
        assert np.min(np.abs(xf[:,None]-xm[None,:]))>.02
    ax.set(ylabel=labels[i],xlim=(0,1),xticks=np.linspace(0,1,6),
           xlabel=r'中心线位置 $(x+L/2)/L$' if is_center else r'孔壁角度 $\psi/(2\pi)$')
    ax.tick_params(direction='in',top=True,right=True,labelsize=7,length=3,width=.65)
    for sp in ax.spines.values(): sp.set_visible(True);sp.set_color('black')
    if is_center:
        ax.axvspan(.35,.65,color='#e8e8e8',zorder=-10)
        ax.text(.5,.1,'孔洞',transform=ax.transAxes,ha='center',fontsize=7,color='.48')
    ax.legend(loc=['upper right','upper right','lower left','upper center'][i],fontsize=6.6,frameon=False)
    ax.set_ylim([(-120,120),(-70,70),(-.04,.92),(-.15,4.25)][i])
fig.text(.5,.014,r'$L=20\,\mu\mathrm{m}$，$a=3\,\mu\mathrm{m}$，$E_0=10^4\,\mathrm{V/m}$；$Q_\infty=\sigma E_0^2$',ha='center',fontsize=7)
b.save_all(fig,'A2_blue_uniform_four_panel')
with (OUT/'plotted_data.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f);w.writerow(['panel','role','normalized_coordinate','value']);w.writerows(records)
audit={'FEM':'Uniform samples evaluated from saved L2 potential using HEX20 shape functions and gradients; no PDE rerun, fitting, or value shifts.',
       'analytic':'Independent closed forms; filled markers halfway between FEM marker coordinates.',
       'centerline_FEM_step':.05,'wall_FEM_step':.05,'analytic_marker_offset':.025,'coordinate_reconstruction_error_m':err,
       'sources':{p.name:b.sha256(p) for p in [b.RUN/'A2_L2_nodes.csv',b.RUN/'A2_L2_wall.csv']},
       'all_uniform_sample_max_abs_errors':{}}
for panel in range(1,5):
    subset=[r for r in records if r[0]==panel and r[1]=='FEM_marker']
    xx=np.array([r[2] for r in subset]); yy=np.array([r[3] for r in subset]);x=(2*xx-1)*half
    if panel==1: expected=-b.E0*(x+a*a/x)*1e3
    elif panel==2: expected=-2*b.E0*a*np.cos(2*np.pi*xx)*1e3
    elif panel==3: expected=(1-a*a/x**2)**2
    else: expected=4*np.sin(2*np.pi*xx)**2
    audit['all_uniform_sample_max_abs_errors'][str(panel)]=float(np.max(np.abs(yy-expected)))
(OUT/'sampling_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
assert '<text' in (OUT/'A2_blue_uniform_four_panel.svg').read_text(encoding='utf-8')
print(json.dumps(audit['all_uniform_sample_max_abs_errors']))
