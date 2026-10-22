"""F3 radius comparison: physical plane maps and field-height surfaces.

Figure contract: compare spatial fields at the three frozen radii. Linear shared
normalization per quantity, expanded green palette, original element patches.
The third plotting coordinate is the displayed scalar, not displacement.
"""
from pathlib import Path
import hashlib,json,csv,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'output/F3_孔洞半径扫描_20260909_v3/source_data'
OUT=ROOT/'output/F3_焦耳热_黄色扩展红色收窄_20260909'
OUT.mkdir(parents=True,exist_ok=True)
for folder in ['pairs','comparisons','source_data']: (OUT/folder).mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Microsoft YaHei','Arial','DejaVu Sans'],
 'font.size':9,'axes.titlesize':11,'axes.labelsize':9,'svg.fonttype':'none','pdf.fonttype':42,
 'axes.unicode_minus':False,'figure.facecolor':'white','savefig.facecolor':'white'})
stops=[(0,'#0000b5'),(.09,'#003dff'),(.16,'#00cee9'),(.24,'#00db35'),
       (.32,'#42eb12'),(.40,'#d4f300'),(.46,'#fff000'),(.63,'#ffdf00'),
       (.84,'#ff9800'),(.94,'#ff5500'),(.97,'#ff0800'),(1,'#bb0000')]
cmap=LinearSegmentedColormap.from_list('expanded_yellow_reduced_red',stops,N=256)
# Signed carrier map: zero is green; negative/positive endpoints blue/red.
signed=LinearSegmentedColormap.from_list('signed_blue_green_red',[
 (0,'#0000bd'),(.18,'#006bff'),(.32,'#00cfcd'),(.43,'#00db48'),(.5,'#30e526'),
 (.60,'#7ae517'),(.70,'#e5ed00'),(.82,'#ff9c00'),(1,'#ed0000')],N=256)
specs=[
 {'key':'Q','face':'midplane','title':'焦耳热密度','label':r'$Q_J$ ($10^{12}$ W m$^{-3}$)', 'scale':1e-12,'limits':(0,17),'ticks':[0,4,8,12,16]},
 {'key':'dn','face':'midplane','title':'电子相对扰动（中面）','label':r'$\Delta n/n_0$ ($10^{-3}$ %)', 'scale':1e5,'limits':(-1.8,1.8),'ticks':[-1.8,-.9,0,.9,1.8]},
 {'key':'phi','face':'midplane','title':'电势分布','label':r'$\phi$ (V)', 'scale':1,'limits':(-3,0),'ticks':[-3,-2.25,-1.5,-.75,0]},
 {'key':'E','face':'midplane','title':'电场强度','label':r'$|\mathbf{E}|$ (kV cm$^{-1}$)', 'scale':1e-5,'limits':(0,3),'ticks':[0,.75,1.5,2.25,3]},
 {'key':'dn','face':'top','title':'电子相对扰动（顶面补充）','label':r'$\Delta n/n_0$ (%)', 'scale':100,'limits':(-30,30),'ticks':[-30,-15,0,15,30]},
]
specs=specs[:1]
blocks={}; metadata={'model':'NT','voltage_V':3,'radii_um':[2,3,4],
 'geometry_um':[20,20,1],'xy_limits_um':[-10,10],'status':'PASS_L0_TREND_ONLY',
 'normalization':'Linear, fixed limits shared by all radii for each field; no clipping, per-map rescaling or gamma.',
 'palette_stops':stops,'height_semantics':'Scalar field values; not geometric thickness or displacement.',
 'carrier_definition':'Signed relative electron density change, n/n0-1; not log-density eta_n; not hole carriers.',
 'integrity':'Original HEX20 patch triangulation; no interpolation across the hole; no averaging between element patches.',
 'sources':{},'fields':[]}
for radius in [2,3,4]:
    for face in ['midplane']:
        p=SRC/f'R0{radius}_NT_{face}.npz'
        with np.load(p) as f: blocks[radius,face]={k:f[k] for k in f.files}
        d=blocks[radius,face]
        assert np.allclose(d['xyz'][:,2],0 if face=='midplane' else .5e-6,atol=1e-18)
        metadata['sources'][p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        np.savez_compressed(OUT/'source_data'/p.name,**{k:d[k] for k in ['xyz','triangles','Q','n','phi','E']})

def values(d,s): return (d['n']-1 if s['key']=='dn' else d[s['key']])*s['scale']
def palette(s): return signed if s['key']=='dn' else cmap
def triangulation(d): return mtri.Triangulation(d['xyz'][:,0]*1e6,d['xyz'][:,1]*1e6,d['triangles'])
def flat(ax,d,v,radius,s):
    ax.tripcolor(triangulation(d),v,cmap=palette(s),norm=Normalize(*s['limits']),shading='gouraud',rasterized=True)
    t=np.linspace(0,2*np.pi,400)
    ax.plot(radius*np.cos(t),radius*np.sin(t),color='black',lw=.7)
    ax.set(aspect='equal',xlim=(-10,10),ylim=(-10,10),xticks=[-10,-5,0,5,10],yticks=[-10,-5,0,5,10],
           xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)')
    for edge in ax.spines.values(): edge.set_linewidth(.75)

def surface(ax,d,v,s,compact=False):
    ax.plot_trisurf(triangulation(d),v,cmap=palette(s),norm=Normalize(*s['limits']),
                   linewidth=0,antialiased=False,shade=False,rasterized=True)
    ax.set(xlim=(-10,10),ylim=(-10,10),zlim=s['limits'],xticks=[-10,0,10],yticks=[-10,0,10],
           zticks=s['ticks'][::2],xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)')
    if not compact: ax.set_zlabel(s['label'],labelpad=9)
    ax.set_box_aspect((1,1,.65));ax.view_init(elev=29,azim=-60)
    ax.tick_params(labelsize=7 if compact else 8,pad=1)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
        axis.pane.fill=False;axis._axinfo['grid']['color']=(.6,.6,.6,.3)

def save(fig,path):
    for ext in ['png','svg','pdf']:
        fig.savefig(path.with_suffix('.'+ext),dpi=320,bbox_inches='tight',facecolor='white')
    plt.close(fig)

for s in specs:
    stats=[]
    for r in [2,3,4]:
        d=blocks[r,s['face']];v=values(d,s)
        assert np.isfinite(v).all() and v.min()>=s['limits'][0] and v.max()<=s['limits'][1],(s,r,v.min(),v.max())
        stats.append({'radius_um':r,'display_min':float(v.min()),'display_max':float(v.max())})
        fig=plt.figure(figsize=(11.5,5.2))
        left=fig.add_axes([.07,.17,.34,.70]);flat(left,d,v,r,s)
        left.set_title('(a) 平面云图',pad=10)
        right=fig.add_axes([.44,.09,.43,.81],projection='3d');surface(right,d,v,s)
        right.set_title('(b) 场值曲面',pad=10)
        cax=fig.add_axes([.94,.23,.015,.52])
        fig.colorbar(ScalarMappable(norm=Normalize(*s['limits']),cmap=palette(s)),cax=cax,ticks=s['ticks'],label=s['label'])
        fig.suptitle(f"{s['title']} | NT, r = {r} μm",y=.99,fontsize=14)
        ztext='0' if s['face']=='midplane' else '0.5'
        fig.text(.5,.025,f'z = {ztext} μm；3 V；L0 趋势结果；曲面高度表示场值',ha='center',fontsize=9)
        save(fig,OUT/'pairs'/f"R0{r}_{s['key']}_{s['face']}_map_height")
    fig=plt.figure(figsize=(12.6,8.3))
    gs=fig.add_gridspec(2,3,left=.055,right=.89,bottom=.10,top=.91,wspace=.28,hspace=.22)
    for j,r in enumerate([2,3,4]):
        d=blocks[r,s['face']];v=values(d,s)
        ax=fig.add_subplot(gs[0,j]);flat(ax,d,v,r,s);ax.set_title(f'r = {r} μm')
        ax=fig.add_subplot(gs[1,j],projection='3d');surface(ax,d,v,s,compact=True)
    cax=fig.add_axes([.94,.24,.012,.49])
    fig.colorbar(ScalarMappable(norm=Normalize(*s['limits']),cmap=palette(s)),cax=cax,ticks=s['ticks'],label=s['label'])
    fig.suptitle(s['title']+' | NT · 三半径统一坐标与色标',fontsize=15,y=.98)
    fig.text(.5,.035,('z = 0 μm' if s['face']=='midplane' else 'z = 0.5 μm')+'；3 V；L0 趋势结果；下排高度表示场值',ha='center',fontsize=10)
    save(fig,OUT/'comparisons'/f"compare_{s['key']}_{s['face']}")
    metadata['fields'].append({**s,'statistics':stats})
    print('EXPORTED '+s['key']+' '+s['face'],flush=True)

(OUT/'figure_contract.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copyfile(__file__,OUT/Path(__file__).name)
links=[]
for s in specs:
    name=f"comparisons/compare_{s['key']}_{s['face']}"
    pairs=' · '.join(f'<a href="pairs/R0{r}_Q_midplane_map_height.png">r={r} μm 平面＋曲面</a>' for r in [2,3,4])
    links.append(f'<h2>{s["title"]}</h2><p>{pairs}</p><p><a href="{name}.pdf">PDF</a> · <a href="{name}.svg">SVG</a></p><img src="{name}.png">')
(OUT/'图像索引.html').write_text('<!doctype html><meta charset="utf-8"><title>F3 三半径场图</title><style>body{max-width:1200px;margin:30px auto;font-family:Arial,sans-serif}img{width:100%}h2{margin-top:45px}</style><h1>F3 三半径场图</h1><p>同物理量共享线性色标与坐标范围。主图为中面，顶面电子图为补充。三维高度为场值，非变形。</p>'+''.join(links),encoding='utf-8')
assert len(list((OUT/'pairs').glob('*.png')))==3
assert len(list((OUT/'comparisons').glob('*.png')))==1
for f in OUT.rglob('*.svg'): assert '<text' in f.read_text(encoding='utf-8')
print('COMPLETE: 3 Joule paired figures + 1 radius comparison, PNG/SVG/PDF; editable text verified.')
