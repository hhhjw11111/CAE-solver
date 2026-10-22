"""True F3 geometry and representative saved HEX20 boundary mesh, reference layout."""
from pathlib import Path
import json, csv, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.collections import PolyCollection

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/figures/F2_voltage_geometry_vertical_bold_20260909'
OUT.mkdir(parents=True,exist_ok=True)
MESH=ROOT/'solver_v2/cases/F1_20260906/H1_L2.mesh'
plt.rcParams.update({'font.family':'Times New Roman','font.weight':'bold',
 'mathtext.fontset':'custom','mathtext.rm':'Times New Roman:weight=bold',
 'mathtext.it':'Times New Roman:style=italic:weight=bold',
 'mathtext.bf':'Times New Roman:weight=bold',
 'font.size':12,'svg.fonttype':'none','pdf.fonttype':42})
# Contract: define radius-only study; true mesh on left, exact-plan geometry on right.
fig=plt.figure(figsize=(6.0,8.4),facecolor='white')
left=fig.add_axes([.055,.54,.89,.445]); right=fig.add_axes([.02,.025,.96,.51])
red='#ed1717'
def text(ax,x,y,s,**kw):
    kw['fontsize']=kw.pop('fontsize',9)*1.4
    kw['fontweight']='bold'
    ax.text(x,y,s,ha=kw.pop('ha','center'),va=kw.pop('va','center'),**kw)
def arrow(ax,a,b,color='black',style='-|>',lw=.7,scale=8):
    ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle=style,color=color,lw=lw,mutation_scale=scale,shrinkA=0,shrinkB=0))
def project(x):
    x=np.asarray(x)
    return np.stack((.985*x[...,0]-.32*x[...,1],.15*x[...,0]+.68*x[...,1]+.72*x[...,2]),axis=-1)
lines=MESH.read_text().splitlines(); nn,ne=map(int,lines[0].split())
xyz=np.array([[float(v) for v in s.split()[1:4]] for s in lines[1:nn+1]])*1e6
conn=np.array([[int(v)-1 for v in s.split()[2:22]] for s in lines[nn+1:nn+ne+1]])
faces=[[0,8,1,9,2,10,3,11],[4,12,5,13,6,14,7,15],
       [0,8,1,17,5,12,4,16],[1,9,2,18,6,13,5,17],
       [2,10,3,19,7,14,6,18],[3,11,0,16,4,15,7,19]]
exterior={}
for el in conn:
    for face in faces:
        ids=el[face]; key=tuple(sorted(ids[::2]))
        if key in exterior: exterior[key]=None
        else: exterior[key]=ids
polys=[]
normal=np.cross([.985,-.32,0],[.15,.68,.72])
def surface_point(pts,u,v):
    # Quadratic eight-node face interpolation; drawing refinement only.
    weights=np.array([
        -(1-u)*(1-v)*(1+u+v)/4, (1-u*u)*(1-v)/2,
        -(1+u)*(1-v)*(1-u+v)/4, (1+u)*(1-v*v)/2,
        -(1+u)*(1+v)*(1-u-v)/4, (1-u*u)*(1+v)/2,
        -(1-u)*(1+v)*(1+u-v)/4, (1-u)*(1-v*v)/2])
    return weights@pts
grid=np.linspace(-1,1,4)
for ids in exterior.values():
    if ids is None: continue
    pts=xyz[ids]; centre=pts.mean(axis=0)
    shade='#d0d0d0' if np.allclose(pts[:,2],.5) else '#adadad'
    for u0,u1 in zip(grid[:-1],grid[1:]):
        for v0,v1 in zip(grid[:-1],grid[1:]):
            refined=np.array([surface_point(pts,u,v) for u,v in
                [(u0,v0),(u1,v0),(u1,v1),(u0,v1)]])
            polys.append((refined.mean(axis=0)@normal,project(refined),shade))
polys.sort(key=lambda p:p[0])
left.add_collection(PolyCollection([p[1] for p in polys],facecolors=[p[2] for p in polys],edgecolors='#535353',linewidths=.18))
left.set(xlim=(-15.5,17.5),ylim=(-12,10.5),aspect='equal');left.axis('off')
q=project([10,-10,.5]); qb=project([10,-10,-.5]); xx=q[0]-.7
left.plot([q[0],xx+.5],[q[1],q[1]],'k-',lw=.5)
left.plot([qb[0],xx+.5],[qb[1],qb[1]],'k-',lw=.5)
arrow(left,(xx,q[1]+1.3),(xx,q[1]),scale=6)
arrow(left,(xx,qb[1]-1.3),(xx,qb[1]),scale=6)
text(left,xx,qb[1]-2.3,r'$t=1\ \mathrm{\mu m}$',fontsize=9)
text(left,0,-11.3,'(a)',fontsize=11)

right.set(xlim=(-20.5,22),ylim=(-16.5,12),aspect='equal');right.axis('off')
right.add_patch(Rectangle((-10,-10),20,20,facecolor='#d4d4d4',edgecolor='black',lw=.7))
for x in (-10,10):right.plot([x,x],[-10,10],color=red,lw=1)
right.add_patch(Circle((0,0),3,facecolor='white',edgecolor='black',lw=.85))
right.plot(0,0,'ko',ms=2.3);arrow(right,(0,0),(2.82,1.026));text(right,1.1,1.1,r'$r$',fontsize=10)
text(right,3.8,7.1,'Natural\nboundary',fontsize=9,linespacing=1.15)
arrow(right,(3.3,5.1),(1.7,2.47),scale=7)
arrow(right,(-8.2,5.8),(-5.8,5.8),scale=7);arrow(right,(-8.2,5.8),(-8.2,7.9),scale=7)
text(right,-4.8,5.8,r'$x_1$');text(right,-8.2,8.9,r'$x_2$')
text(right,-15.3,.6,r'$\varphi=0\ \mathrm{V}$'+'\n'+r'$T=300\ \mathrm{K}$'+'\n'+r'$\Delta n=\Delta p=0$',color=red,fontsize=8.5,linespacing=1.5)
text(right,15,.6,r'$\varphi=-V_{\mathrm{app}}$'+'\n'+r'$T=300\ \mathrm{K}$'+'\n'+r'$\Delta n=\Delta p=0$',color=red,fontsize=8.5,linespacing=1.5)
text(right,0,-5.5,r'$r=3\ \mathrm{\mu m}$',fontsize=10)
text(right,0,-8.0,r'$V_{\mathrm{app}}=0$–$10\ \mathrm{V}$',color=red,fontsize=9)
for x in (-10,10):right.plot([x,x],[-10.5,-12.4],color='black',lw=.5)
arrow(right,(-10,-11.8),(10,-11.8),style='<->',scale=7)
text(right,0,-13.0,r'$L=20\ \mathrm{\mu m}$',fontsize=10)
for y in (-10,10):right.plot([19.5,20.6],[y,y],color='black',lw=.5)
arrow(right,(20,-10),(20,10),style='<->',scale=7)
text(right,21.4,0,r'$L=20\ \mathrm{\mu m}$',rotation=90,fontsize=9)
text(right,0,-15,'(b)',fontsize=11)
stem=OUT/'F2_voltage_geometry_vertical'
for ext in ['png','svg','pdf','tiff']:
    kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
    fig.savefig(stem.with_suffix('.'+ext),dpi=600 if ext=='tiff' else 450,facecolor='white',**kw)
plt.close(fig)
contract={'claim':'Applied-voltage study geometry with fixed hole radius and thermal contacts',
 'archetype':'schematic-led composite','size_inches':[6.0,8.4],'backend':'python',
 'layout':'a above b; no thermal-load annotations in a; compact panel labels',
 'typography':'All text 1.4 times prior size; bold Times New Roman including math; expanded side annotation spacing',
 'mesh_display':'Each saved L2 boundary face subdivided 3x3 by quadratic interpolation for a dense schematic; not a newly solved FEM mesh and not evidence of mesh convergence',
 'geometry_um':[20,20,1],'radii_um':[3],'mesh':str(MESH.relative_to(ROOT)),
 'applied_voltage_V':[0,.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10],
 'left_panel':'Refined drawing of HEX20 boundary geometry at r=3 um; physical thickness retained',
 'right_panel':'True plan aspect ratio; fixed r=3 um; variable applied voltage',
 'contacts':{'left_phi_V':0,'right_phi_V':'-V_app','both_T_K':300,'both_density_perturbation':0},
 'natural':'Zero traction, normal electric displacement, carrier currents and heat flux on unelectroded surfaces and hole wall',
 'anchors_um':{'A':[-10,-10,-.5],'B':[-10,10,-.5],'C':[-10,-10,.5]},
 'constraints':'u1(A,B,C)=0; u2(A)=0; u3(A,B)=0',
 'source':['solver_v2/docs/f2_voltage_sweep_notes_20260909.md','solver_v2/src/app/f1_case_runner_m.f90'],
 'scope':'Geometry only; applied-voltage range includes attempted cases, not an assertion of convergence. 7-10 V have no complete NT solution in current records; no additional solve.'}
stem.with_suffix('.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')
with stem.with_suffix('.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f);w.writerow(['L_um','t_um','r_um','left_phi_V','right_phi_V','contact_T_K'])
    w.writerows([[20,1,3,0,-v,300] for v in contract['applied_voltage_V']])
shutil.copy2(__file__,OUT/Path(__file__).name)
assert '<text' in stem.with_suffix('.svg').read_text(encoding='utf-8')
print(stem.with_suffix('.png'))
