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
OUT=ROOT/'output/figures/F3_radius_geometry_20260909'
OUT.mkdir(parents=True,exist_ok=True)
MESH=ROOT/'solver_v2/cases/F1_20260906/H1_L0.mesh'
plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix',
 'font.size':9,'svg.fonttype':'none','pdf.fonttype':42})
# Contract: define radius-only study; true mesh on left, exact-plan geometry on right.
fig=plt.figure(figsize=(10.6,5.1),facecolor='white')
left=fig.add_axes([.015,.19,.455,.73]); right=fig.add_axes([.47,.19,.52,.73])
red='#ed1717'
def text(ax,x,y,s,**kw):
    ax.text(x,y,s,ha=kw.pop('ha','center'),va=kw.pop('va','center'),**kw)
def arrow(ax,a,b,color='black',style='-|>',lw=.7,scale=8):
    ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle=style,color=color,lw=lw,mutation_scale=scale,shrinkA=0,shrinkB=0))
def project(x):
    x=np.asarray(x)
    return np.stack((.94*x[...,0]+.342*x[...,1],-.22*x[...,0]+.605*x[...,1]+.765*x[...,2]),axis=-1)
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
normal=np.cross([.94,.342,0],[-.22,.605,.765])
for ids in exterior.values():
    if ids is None: continue
    pts=xyz[ids]; centre=pts.mean(axis=0)
    shade='#d0d0d0' if np.allclose(pts[:,2],.5) else '#adadad'
    polys.append((centre@normal,project(pts),shade))
polys.sort(key=lambda p:p[0])
left.add_collection(PolyCollection([p[1] for p in polys],facecolors=[p[2] for p in polys],edgecolors='#535353',linewidths=.22))
left.set(xlim=(-15,16),ylim=(-12,12),aspect='equal');left.axis('off')
# Electrodes are opposing x faces. Both are at 300 K; top/bottom are adiabatic.
for x in (-10,10):
    q=project([[x,-10,.5],[x,10,.5]])
    left.plot(q[:,0],q[:,1],color=red,lw=.9)
text(left,-7.6,10.2,r'$T=300\ \mathrm{K}$',color=red,fontsize=10)
arrow(left,(-7.6,9.2),project([-10,2,.5]),color=red)
text(left,12,-5.8,r'$T=300\ \mathrm{K}$',color=red,fontsize=10)
arrow(left,(12,-4.9),project([10,1,.5]),color=red)
q=project([10,-10,.5]); qb=project([10,-10,-.5]); xx=q[0]+2.0
left.plot([q[0],xx+.5],[q[1],q[1]],'k-',lw=.5)
left.plot([qb[0],xx+.5],[qb[1],qb[1]],'k-',lw=.5)
arrow(left,(xx,q[1]+1.3),(xx,q[1]),scale=6)
arrow(left,(xx,qb[1]-1.3),(xx,qb[1]),scale=6)
text(left,xx+1.0,qb[1]-.7,r'$t=1\ \mathrm{\mu m}$',ha='left',fontsize=9)
origin=np.array([-11,-8.8]);
for vec,label,off in [([2.5,0,0],r'$x_1$',(.6,-.1)),([0,2.5,0],r'$x_2$',(-.3,.5)),([0,0,2.5],r'$x_3\parallel +c$',(-.7,.7))]:
    tip=origin+project(vec);arrow(left,origin,tip,scale=6);text(left,*(tip+off),label,fontsize=8)
text(left,0,-12.3,r'(a)  3D mesh, $r=3\ \mathrm{\mu m}$',fontsize=10)

right.set(xlim=(-17.8,18.8),ylim=(-15,13),aspect='equal');right.axis('off')
right.add_patch(Rectangle((-10,-10),20,20,facecolor='#d4d4d4',edgecolor='black',lw=.7))
for x in (-10,10):right.plot([x,x],[-10,10],color=red,lw=1)
right.add_patch(Circle((0,0),3,facecolor='white',edgecolor='black',lw=.85))
right.plot(0,0,'ko',ms=2.3);arrow(right,(0,0),(2.82,1.026));text(right,1.1,1.1,r'$r$',fontsize=10)
text(right,3.8,7.1,'Natural\nboundary',fontsize=9,linespacing=1.15)
arrow(right,(3.3,5.1),(1.7,2.47),scale=7)
arrow(right,(-8.2,6.6),(-5.8,6.6),scale=7);arrow(right,(-8.2,6.6),(-8.2,8.8),scale=7)
text(right,-5.25,6.6,r'$x_1$');text(right,-8.2,9.4,r'$x_2$')
text(right,-13.5,.6,r'$\varphi=0\ \mathrm{V}$'+'\n'+r'$T=300\ \mathrm{K}$'+'\n'+r'$\Delta n=\Delta p=0$',color=red,fontsize=8.5,linespacing=1.7)
text(right,13.35,.6,r'$\varphi=-3\ \mathrm{V}$'+'\n'+r'$T=300\ \mathrm{K}$'+'\n'+r'$\Delta n=\Delta p=0$',color=red,fontsize=8.5,linespacing=1.7)
text(right,0,-5.9,r'$r=0\,/\,2\,/\,3\,/\,4\ \mathrm{\mu m}$',fontsize=10)
text(right,0,-8.1,r'$r=0$: no-hole control',fontsize=8)
for x in (-10,10):right.plot([x,x],[-10.5,-12.4],color='black',lw=.5)
arrow(right,(-10,-11.8),(10,-11.8),style='<->',scale=7)
text(right,0,-13.0,r'$L=20\ \mathrm{\mu m}$',fontsize=10)
for y in (-10,10):right.plot([16.7,17.8],[y,y],color='black',lw=.5)
arrow(right,(17.2,-10),(17.2,10),style='<->',scale=7)
text(right,18.25,0,r'$L=20\ \mathrm{\mu m}$',rotation=90,fontsize=9)
text(right,0,-15,'(b)  Geometry and boundary conditions',fontsize=10)
fig.text(.5,.106,'Unelectroded surfaces and hole wall: traction-free, electrically insulating and adiabatic.',ha='center',fontsize=9)
fig.text(.5,.061,'Mechanical constraints: $u_1(A,B,C)=0$, $u_2(A)=0$, $u_3(A,B)=0$.',ha='center',fontsize=9)
fig.text(.5,.025,r'$A=(-10,-10,-0.5)$, $B=(-10,10,-0.5)$, $C=(-10,-10,0.5)$ $\mathrm{\mu m}$.',ha='center',fontsize=8)
stem=OUT/'F3_radius_geometry'
for ext in ['png','svg','pdf','tiff']:
    kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
    fig.savefig(stem.with_suffix('.'+ext),dpi=600 if ext=='tiff' else 450,facecolor='white',**kw)
plt.close(fig)
contract={'claim':'Radius-only study geometry with fixed electrical and thermal contacts',
 'archetype':'schematic-led composite','size_inches':[10.6,5.1],'backend':'python',
 'geometry_um':[20,20,1],'radii_um':[0,2,3,4],'mesh':str(MESH.relative_to(ROOT)),
 'left_panel':'Actual exterior HEX20 face mesh at r=3 um; physical thickness retained',
 'right_panel':'True plan aspect ratio; representative r=3 um; listed independent radius cases',
 'contacts':{'left_phi_V':0,'right_phi_V':-3,'both_T_K':300,'both_density_perturbation':0},
 'natural':'Zero traction, normal electric displacement, carrier currents and heat flux on unelectroded surfaces and hole wall',
 'anchors_um':{'A':[-10,-10,-.5],'B':[-10,10,-.5],'C':[-10,-10,.5]},
 'constraints':'u1(A,B,C)=0; u2(A)=0; u3(A,B)=0',
 'source':['solver_v2/docs/f3_radius_sweep_protocol_20260908.md','solver_v2/src/app/f1_case_runner_m.f90'],
 'scope':'Geometry only; no additional numerical solve or field-data claims'}
stem.with_suffix('.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')
with stem.with_suffix('.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f);w.writerow(['L_um','t_um','r_um','left_phi_V','right_phi_V','contact_T_K'])
    w.writerows([[20,1,r,0,-3,300] for r in [0,2,3,4]])
shutil.copy2(__file__,OUT/Path(__file__).name)
assert '<text' in stem.with_suffix('.svg').read_text(encoding='utf-8')
print(stem.with_suffix('.png'))
