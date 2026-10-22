"""Higher-contrast, discrete rainbow rendering of verified H1-L0 Joule data."""
from pathlib import Path
import json
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'output/H1_L0_焦耳热_JET_20260909'
OUT = ROOT/'output/H1_L0_焦耳热_增强对比_20260909'
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':12,
                     'svg.fonttype':'none','pdf.fonttype':42})
# Quantitative grid: expose hole-edge hotspots and low-field lateral regions.
# Same linear 0--16 scale for every model. Discrete one-unit bands; no clipping,
# data smoothing, geometric changes or cross-hole triangulation.
levels = np.arange(17, dtype=float)
cmap = LinearSegmentedColormap.from_list('contrast_rainbow',[
    (0.,'#000080'),(.10,'#0000ff'),(.23,'#00cfff'),(.34,'#00ffd4'),
    (.47,'#4cff43'),(.62,'#ffff00'),(.78,'#ff8500'),(.92,'#ff0000'),(1.,'#b80000')],N=16)
norm = BoundaryNorm(levels,cmap.N)
data={}
audit={'scale_W_m3':[0,16e12],'band_width_W_m3':1e12,'normalization':'linear; discrete bands',
       'data_modified':False,'clipped':False,'models':{}}
for model in ['LQ','N0','NT']:
    path=SOURCE/f'H1_{model}_L0_midplane.npz'
    with np.load(path) as f: data[model]={k:f[k] for k in f.files}
    q=data[model]['Q']*1e-12
    assert np.all(np.isfinite(q)) and q.min()>=0 and q.max()<=16
    audit['models'][model]={'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                             'Q_max_W_m3':float(q.max()*1e12)}

def draw(ax,model):
    d=data[model]; xy=d['xyz'][:,:2]*1e6
    tri=mtri.Triangulation(xy[:,0],xy[:,1],d['triangles'])
    im=ax.tricontourf(tri,d['Q']*1e-12,levels=levels,cmap=cmap,norm=norm)
    t=np.linspace(0,2*np.pi,400)
    ax.plot(3*np.cos(t),3*np.sin(t),'k',lw=.8)
    for s in ax.spines.values(): s.set_linewidth(.85)
    ax.set(aspect='equal',xlim=(-10,10),ylim=(-10,10),xticks=[-10,-5,0,5,10],
           yticks=[-10,-5,0,5,10],xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)',title=f'H1-{model} | L0')
    return im

def save(fig,name):
    for ext in ['png','pdf','svg']:
        fig.savefig(OUT/f'{name}.{ext}',dpi=320,bbox_inches='tight',facecolor='white')
    plt.close(fig)

label=r'$Q_J$ ($10^{12}$ W m$^{-3}$)'
for model in data:
    fig,ax=plt.subplots(figsize=(6.2,5.4),layout='constrained')
    im=draw(ax,model)
    fig.colorbar(im,ax=ax,shrink=.9,pad=.035,ticks=[0,2,4,6,8,10,12,14,16],label=label)
    save(fig,f'H1_{model}_L0_joule_high_contrast')
fig,axs=plt.subplots(1,3,figsize=(12,4.5),layout='constrained')
for ax,model in zip(axs,data): im=draw(ax,model)
fig.colorbar(im,ax=axs,shrink=.80,pad=.025,ticks=[0,4,8,12,16],label=label)
fig.suptitle(r'Joule heat density | $z=0$, $|\Delta V|=3$ V',fontsize=14)
save(fig,'H1_L0_three_models_high_contrast')
(OUT/'render_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print('Exported three individual maps and one comparison in PNG/PDF/SVG.')
