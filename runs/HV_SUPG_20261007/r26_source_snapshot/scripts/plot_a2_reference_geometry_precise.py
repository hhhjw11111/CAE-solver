"""A2 geometry drawn to scale in the supplied compact reference style."""
from pathlib import Path
import json, csv, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/figures/A2_geometry_reference_20260909'
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family':'Times New Roman', 'mathtext.fontset':'stix',
                     'font.size':9, 'svg.fonttype':'none', 'pdf.fonttype':42})
fig, ax = plt.subplots(figsize=(4.9,4.65))
fig.subplots_adjust(0,0,1,1)
ax.set(xlim=(-17,17), ylim=(-17,15.3), aspect='equal')
ax.axis('off')
blue='#ff0000'
ax.add_patch(Rectangle((-10,-10),20,20, facecolor='#d4d4d4',edgecolor=blue,lw=.85))
ax.add_patch(Circle((0,0),3,facecolor='white',edgecolor='black',lw=.8))
ax.plot(0,0,'ko',ms=2)
def arrow(a,b,color='black',style='-|>',lw=.65,scale=7):
    ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle=style,color=color,lw=lw,mutation_scale=scale,shrinkA=0,shrinkB=0))
def txt(x,y,s,**kw):
    ax.text(x,y,s,ha=kw.pop('ha','center'),va=kw.pop('va','center'),**kw)
# Direction triad uses the actual solver plane x1-x2; x3 is out of plane.
arrow((-8.1,6.1),(-5.7,6.1)); arrow((-8.1,6.1),(-8.1,8.5))
txt(-5.15,6.05,r'$x_1$');txt(-8.1,9.15,r'$x_2$')
# Radius and insulating wall annotation.
arrow((0,0),(2.819,1.026));txt(1.2,1.05,r'$r$')
txt(3.85,6.9,'Natural\nboundary',fontsize=9,linespacing=1.1)
txt(3.85,4.8,r'$\mathbf{E}\!\cdot\!\mathbf{n}=0$',fontsize=9)
arrow((3.5,3.95),(1.72,2.46),scale=6)
txt(0,-6.0,r'$r=3\ \mathrm{\mu m}$',fontsize=10)
txt(0,-8.1,r'$t=1\ \mathrm{\mu m}$  (along $x_3$)',fontsize=8.5)
# One blue condition applies to the entire outer perimeter.
txt(0,13.75,r'$\varphi=-E_0x_1\!\left(1+\dfrac{r^2}{x_1^2+x_2^2}\right)$',color=blue,fontsize=10)
txt(0,11.45,r'$\Delta n=\Delta p=0$  on all four outer edges',color=blue,fontsize=8.5)
# True dimension lines with extensions and inward arrowheads.
for x in (-10,10): ax.plot([x,x],[-10.6,-13.0],color='black',lw=.55)
arrow((-10,-12.3),(10,-12.3),style='<->',scale=7)
txt(0,-13.7,r'$L=20\ \mathrm{\mu m}$',fontsize=10)
for y in (-10,10): ax.plot([12.7,14.1],[y,y],color='black',lw=.55)
arrow((13.35,-10),(13.35,10),style='<->',scale=7)
txt(14.8,0,r'$L=20\ \mathrm{\mu m}$',rotation=90,fontsize=10)
txt(-13.2,0,r'$\mathbf{u}=0$'+'\n'+r'$\Delta T=0$'+'\nthroughout',fontsize=8.5,linespacing=1.5)
txt(0,-16,r'$E_0=10^4\ \mathrm{V\,m^{-1}}$',fontsize=9)
stem=OUT/'A2_geometry_reference'
for ext in ('svg','pdf','png','tiff'):
    opts={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
    fig.savefig(stem.with_suffix('.'+ext),dpi=600 if ext=='tiff' else 450,facecolor='white',**opts)
plt.close(fig)
contract={'model':'A2_F1_insulated_hole','dimensions_um':[20,20,1],'radius_um':3,
 'plane':['x1','x2'],'E0_V_per_m':10000,'outer_boundary':'phi=-E0*x1*(1+r^2/(x1^2+x2^2)); delta_n=delta_p=0',
 'hole_wall':'E dot n = 0; natural carrier boundary','throughout':'u=0; delta_T=0',
 'source':'solver_v2/src/verification/a2_hole_check_m.f90',
 'run':'solver_v2/runs/F1_20260906_R0144/A2/A2_stdout.txt',
 'geometry':'Exact physical in-plane aspect ratio; hole diameter / plate width = 0.3',
 'scope':'Geometry schematic only; synthetic isothermal A2 benchmark; no new solver run'}
stem.with_suffix('.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')
with stem.with_suffix('.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.writer(f);w.writerow(['parameter','value','unit'])
    w.writerows([['L',20,'um'],['r',3,'um'],['t',1,'um'],['E0',10000,'V/m']])
shutil.copy2(__file__,OUT/Path(__file__).name)
assert '<text' in stem.with_suffix('.svg').read_text(encoding='utf-8')
print(stem.with_suffix('.png'))
