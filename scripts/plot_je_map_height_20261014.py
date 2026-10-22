"""Verified 6/7/8 V current-work fields, using the requested map/height layout."""
from pathlib import Path
import argparse
import csv
import json
import shutil
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
import plot_je_8v_clouds_20261013 as fe

QA = Path('C:/Users/Administrator/.codex/skills/nature-figure/scripts')
sys.path.insert(0, str(QA))
from audit_panel_alignment import require_matplotlib_panel_alignment

ROOT = fe.SOLVER.parent
RUNS = fe.SOLVER/'runs/JE_6_8_HALF_V_20261013'
CACHE = RUNS/'clouds_6_7_8V_template_20261013'
OUT = ROOT/'output/JE_绘图风格修复_20261014'
STYLE_PATH = fe.SOLVER/'scripts/styles/joule_map_height_larger_surfaces.json'
VOLTAGES = [6,7,8]
FACES = {'midplane':0.,'top':.5e-6,'bottom':-.5e-6}
# All surfaces show physical scalar values, not geometric displacement.
SPECS = [
 {'id':'01_joule_mid','key':'Q','factor':1e-12,'face':'midplane','title':'焦耳热密度（J·E）',
  'label':r'$Q_{J\cdot E}$ ($10^{12}$ W m$^{-3}$)','limits':[0,110],'ticks':[0,25,50,75,100],'cmap':'heat'},
 {'id':'02_potential_mid','key':'phi','factor':1.,'face':'midplane','title':'电势分布',
  'label':r'$\phi$ (V)','limits':[-8,0],'ticks':[-8,-6,-4,-2,0],'cmap':'heat'},
 {'id':'03_electron_mid','key':'dn_percent','factor':1.,'face':'midplane','title':'电子相对浓度扰动（中面）',
  'label':r'$\Delta n/n_0$ (%)','limits':[-.2,.2],'ticks':[-.2,-.1,0,.1,.2],'cmap':'signed'},
 {'id':'04_hole_mid','key':'dp_percent','factor':1.,'face':'midplane','title':'空穴相对浓度扰动（中面）',
  'label':r'$\Delta p/p_0$ (%)','limits':[-35,35],'ticks':[-35,-17.5,0,17.5,35],'cmap':'signed'},
 {'id':'05_total_stress_mid','key':'vm','factor':1e-6,'face':'midplane','title':'耦合总等效应力（含热效应）',
  'label':r'$\sigma_{\mathrm{VM}}$ (MPa)','limits':[0,6],'ticks':[0,1.5,3,4.5,6],'cmap':'heat'},
 {'id':'06_thermal_term_mid','key':'thermal_xx','factor':1e-6,'face':'midplane','title':'热膨胀本构应力项（xx 分量）',
  'label':r'$-\beta_{xx}\theta$ (MPa)','limits':[-45,0],'ticks':[-45,-30,-15,0],'cmap':'heat'},
 {'id':'07_electron_bottom','key':'dn_percent','factor':1.,'face':'bottom','title':'电子相对浓度扰动（下表面）',
  'label':r'$\Delta n/n_0$ (%)','limits':[-300,300],'ticks':[-300,-150,0,150,300],'cmap':'signed'},
 {'id':'08_hole_top','key':'dp_percent','factor':1.,'face':'top','title':'空穴相对浓度扰动（上表面）',
  'label':r'$\Delta p/p_0$ (%)','limits':[-800,800],'ticks':[-800,-400,0,400,800],'cmap':'signed'},
]

def write_json(path, obj):
    fe.write_json(path,obj)

def verify_and_load(only=None):
    for sub in ['figures','source_data','qa']:
        (OUT/sub).mkdir(parents=True,exist_ok=True)
    contract = {
      'question':'How do actual 6/7/8 V current-work-branch spatial fields change?',
      'archetype':'paired 2D field maps and enlarged scalar-height surfaces',
      'evidence_roles':{'top':'hotspot location and spatial footprint','bottom':'height, signs and peak amplitudes on the identical fields'},
      'reference_adaptation':'style and layout only; voltage replaces radius',
      'selection':'All requested physical quantities included. Three representative integer voltages selected for the supplied three-column template; available 6.5/7.5 V states remain in the verification table.',
      'geometry_um':[20,20,1],'radius_um':3,'voltages_V':VOLTAGES,
      'branch':'E_DOT_TOTAL_CURRENT_DIAGNOSTIC_NOT_THERMODYNAMIC_CLOSURE',
      'scope':'Full-load solved fields and physical-flux gates passed; mesh/stabilization acceptance pending. Not the historical drift-only branch.',
      'stress':'sigma=C:epsilon-e^T E-beta*theta; vm uses total stress. -beta_xx*theta is a constitutive term, not an independently equilibrated thermal-only solution.',
      'carriers':{'n0_m3':fe.N0,'p0_m3':fe.P0,'definition':'signed 100*(exp(eta)-1); relative hole changes do not imply appreciable absolute hole density'},
      'scales':'fixed linear shared color and height limits within each figure; no per-voltage normalization; signed carrier maps have zero at green',
      'integrity':'HEX20 element-local evaluation; original 16-subdivision face triangulation; all cached slice samples retained; no smoothing across elements or hole.',
      'height':'scalar field value, not deformation','statistics':'deterministic single-mesh fields, no statistical interval',
      'export':'Frozen F3 larger-surfaces layout with tight bounding box; PNG 320 dpi plus PDF and editable SVG; no journal submission format claimed',
      'style':{'path':str(STYLE_PATH.relative_to(ROOT)), 'sha256':fe.sha256(STYLE_PATH),
               'reference':'F3_焦耳热_下排再次放大_20260909/comparisons/compare_Q_midplane.png'},
      'panels':[s for s in SPECS if not only or s['id'] in only],
    }
    write_json(OUT/'figure_contract.json',contract)
    old = json.loads((CACHE/'provenance.json').read_text(encoding='utf-8'))
    mesh = fe.SOLVER/'cases/HV_SURFACE_20261008/XY0_Z2/H1_XY0_Z2.mesh'
    assert fe.sha256(mesh)==old['mesh']['sha256']
    assert fe.sha256(Path(fe.__file__))==old['evaluator']['sha256']
    xyz,conn = fe.read_mesh(mesh)
    assert len(xyz)==13920 and len(conn)==3136
    verified=[]; data={}; source_hashes={}
    for v in [6,6.5,7,7.5,8]:
        run=RUNS/f'V{int(v*10):03d}'
        m=fe.read_metrics(run/'metrics.txt')
        inp=(RUNS/f'input_V{int(v*10):03d}.nml').read_text(encoding='utf-8-sig')
        assert m['voltage_V']==v and m['final_voltage_factor']==m['final_joule_factor']==1
        assert m['case_acceptance']=='PASS_NOT_MESH_ACCEPTANCE'
        assert m['heat_source']==contract['branch'] and 'heat_source_current_work=.true.' in inp
        assert max(m['maximum_accepted_residual'],m['maximum_accepted_update'])<1e-9
        assert m['heat_balance_relative']<.01 and m['electron_balance_relative']<.01
        assert m['joule_quadrature_4_5_relative']<.001
        verified.append({'voltage_V':v,'temperature_rise_max_K':m['theta_max_K'],
          'joule_power_mW':m['joule_W']*1e3,'total_VM_sample_max_MPa':m['stress_von_mises_max_Pa']*1e-6,
          'heat_balance_error_percent':100*m['heat_balance_relative'],
          'final_voltage_factor':m['final_voltage_factor'],'final_joule_factor':m['final_joule_factor'],
          'residual':m['maximum_accepted_residual'],'update':m['maximum_accepted_update'],
          'status':m['case_acceptance'],'plotted':v in VOLTAGES})
        if v not in VOLTAGES:
            continue
        for name in ['nodes.csv','metrics.txt','field_samples.csv']:
            digest=fe.sha256(run/name)
            assert digest==old[str(v)][name]['sha256'],(v,name)
            source_hashes[str((run/name).relative_to(ROOT))]=digest
        native=json.loads((CACHE/'qa'/f'V{v}_native_audit.json').read_text(encoding='utf-8'))
        assert native['status']=='PASS_POSTPROCESSING_RECONSTRUCTION'
        assert max(native['maximum_error_normalized_by_column_peak'])<5e-7
        shutil.copy2(CACHE/'qa'/f'V{v}_native_audit.json',OUT/'qa'/f'V{v}_prior_native_audit.json')
        nodes=np.loadtxt(run/'nodes.csv',delimiter=',',skiprows=1)
        assert np.isfinite(nodes).all() and np.max(np.abs(nodes[:,1:4]-xyz))<1e-19
        assert np.max(abs(nodes[np.isclose(xyz[:,0],-1e-5,atol=1e-18,rtol=0),7]))<1e-10
        assert np.max(abs(nodes[np.isclose(xyz[:,0],1e-5,atol=1e-18,rtol=0),7]+v))<1e-10
        state=nodes[:,4:11]
        audit=[]
        for face,z in FACES.items():
            path=CACHE/'source_data'/f'V{v}_{face}.npz'
            d=dict(np.load(path))
            assert np.isfinite(d['xyz']).all() and np.max(abs(d['xyz'][:,2]-z))<1e-18
            assert np.all(d['nr']>0) and np.all(d['pr']>0)
            elements=np.unique(d['element_id'])
            chosen=elements[np.linspace(0,len(elements)-1,11,dtype=int)]
            side=1. if z>0 else -1.
            t=np.linspace(-1,1,17)
            points=np.array([[r,s,side] for s in t for r in t])
            max_err=0.
            for element in chosen:
                e=int(element)-1
                reconstructed=fe.evaluate(xyz[conn[e]],state[conn[e]],points)
                mask=d['element_id']==element
                for key in ['xyz','phi','Q','theta','nr','pr','dn_percent','dp_percent','vm','stress_thermal']:
                    a=d[key][mask]; b=reconstructed[key]
                    err=float(np.max(abs(a-b))/max(float(np.max(abs(b))),1e-20))
                    assert err<1e-9,(v,face,element,key,err)
                    max_err=max(max_err,err)
            assert np.all(d['element_id'][d['triangles']]==d['element_id'][d['triangles']][:,:1])
            audit.append({'face':face,'node_samples':len(d['xyz']),'triangles':len(d['triangles']),
              'sampled_elements':chosen.tolist(),'maximum_cache_reconstruction_error':max_err})
            shutil.copy2(path,OUT/'source_data'/path.name)
            source_hashes[str(path.relative_to(ROOT))]=fe.sha256(path)
            data[v,face]=d
        write_json(OUT/'qa'/f'V{v}_cache_validation.json',audit)
        print(f'VERIFIED {v} V: actual full load; native provenance and cached slice evaluation match',flush=True)
    write_json(OUT/'verified_results.json',verified)
    write_json(OUT/'provenance.json',{'sources':source_hashes,'mesh_sha256':fe.sha256(mesh),
       'evaluator_sha256':fe.sha256(Path(fe.__file__)),'mesh_nodes':len(xyz),'HEX20_elements':len(conn)})
    with (OUT/'verified_results.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(verified[0]));w.writeheader();w.writerows(verified)
    return data

def value(d,s):
    return (d['stress_thermal'][:,0] if s['key']=='thermal_xx' else d[s['key']])*s['factor']

def render(data, only=None):
    style=json.loads(STYLE_PATH.read_text(encoding='utf-8-sig'))
    heat=LinearSegmentedColormap.from_list('reference_heat',style['palette_stops'],N=style['colormap_levels'])
    signed=LinearSegmentedColormap.from_list('signed_blue_green_red',[
      (0,'#0000bd'),(.18,'#006bff'),(.32,'#00cfcd'),(.43,'#00db48'),(.5,'#30e526'),
      (.60,'#7ae517'),(.70,'#e5ed00'),(.82,'#ff9c00'),(1,'#ed0000')],N=256)
    mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':style['fonts'],
      'font.size':9,'axes.titlesize':11,'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9,
      'axes.linewidth':.8,'axes.labelpad':4,'xtick.major.pad':3.5,'ytick.major.pad':3.5,
      'xtick.major.size':3.5,'ytick.major.size':3.5,
      'svg.fonttype':'none','pdf.fonttype':42,'axes.unicode_minus':False,
      'xtick.direction':'out','ytick.direction':'out','figure.facecolor':'white','savefig.facecolor':'white'})
    extrema=[]
    for s in SPECS:
        if only and s['id'] not in only:
            continue
        cmap=heat if s['cmap']=='heat' else signed
        norm=Normalize(*s['limits'])
        fig=plt.figure(figsize=style['comparison_figsize_inches'])
        layout=style['comparison_layout'].copy()
        rows=layout.pop('rows');columns=layout.pop('columns')
        gs=fig.add_gridspec(rows,columns,**layout)
        upper=[];lower=[]
        for j,v in enumerate(VOLTAGES):
            d=data[v,s['face']]; values=value(d,s)
            assert np.isfinite(values).all()
            assert values.min()>=s['limits'][0]-1e-9 and values.max()<=s['limits'][1]+1e-9
            xy=d['xyz'][:,:2]*1e6
            tri=mtri.Triangulation(xy[:,0],xy[:,1],d['triangles'])
            ax=fig.add_subplot(gs[0,j]);upper.append(ax)
            ax.tripcolor(tri,values,shading='gouraud',cmap=cmap,norm=norm,rasterized=True)
            angle=np.linspace(0,2*np.pi,400)
            ax.plot(3*np.cos(angle),3*np.sin(angle),color='black',lw=.7)
            ax.set(aspect='equal',xlim=style['xy_limits_um'],ylim=style['xy_limits_um'],
                   xticks=style['xy_ticks_um'],yticks=style['xy_ticks_um'],
                   xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)',title=f'{v} V')
            for edge in ax.spines.values():edge.set_linewidth(.75)
            bx=fig.add_subplot(gs[1,j],projection='3d');lower.append(bx)
            bx.set_position(style['lower_row_enlargement']['axes_positions'][j])
            bx.patch.set_alpha(0)
            bx.plot_trisurf(tri,values,cmap=cmap,norm=norm,linewidth=0,
                           antialiased=False,shade=False,rasterized=True)
            bx.set(xlim=(-10,10),ylim=(-10,10),zlim=s['limits'],xticks=[-10,0,10],yticks=[-10,0,10],
                   zticks=s['ticks'][::2],xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)')
            view=style['surface_view'];enlargement=style['lower_row_enlargement']
            bx.set_box_aspect(view['box_aspect'],zoom=enlargement['projection_zoom'])
            bx.view_init(elev=view['elev'],azim=view['azim'])
            bx.tick_params(labelsize=enlargement['tick_fontsize'],pad=1)
            for axis in [bx.xaxis,bx.yaxis,bx.zaxis]:
                axis.pane.fill=False;axis._axinfo['grid']['color']=(.6,.6,.6,.3)
            for kind,index in [('min',np.argmin(values)),('max',np.argmax(values))]:
                extrema.append({'figure':s['id'],'voltage_V':v,'face':s['face'],'field':s['key'],
                    'extreme':kind,'display_value':float(values[index]),'x_um':float(xy[index,0]),
                    'y_um':float(xy[index,1]),'z_um':FACES[s['face']]*1e6})
        cax=fig.add_axes(style['colorbar_axes'])
        cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,ticks=s['ticks'])
        cb.set_label(s['label'],fontsize=9,labelpad=4);cb.ax.tick_params(labelsize=9)
        fig.suptitle(s['title'],fontsize=15,y=.98)
        # Upper and lower rows deliberately have different sizes, as requested.
        # Equal physical widths, heights and gutters are checked within each row.
        top_ids=[f'map_{v}V' for v in VOLTAGES]
        bottom_ids=[f'surface_{v}V' for v in VOLTAGES]
        ids=top_ids+bottom_ids
        require_matplotlib_panel_alignment(fig,axes=upper+lower,panel_ids=ids,
            row_groups=[top_ids,bottom_ids],column_groups=[[a,b] for a,b in zip(top_ids,bottom_ids)],
            exemptions=[{'panels':bottom_ids,'checks':['column'],
              'reason':'User requests enlarged lower-row 3D surfaces; column widths intentionally differ from the upper 2D maps. Both rows retain equal widths and gutters internally.'}],
            json_out=OUT/'qa'/f'{s["id"]}_alignment.json',strict=True)
        stem=OUT/'figures'/s['id']
        for ext in style['export']['formats']:
            fig.savefig(stem.with_suffix('.'+ext),dpi=style['export']['dpi'],
                        bbox_inches='tight',facecolor='white')
        plt.close(fig)
        assert '<text' in stem.with_suffix('.svg').read_text(encoding='utf-8')
        print('EXPORTED '+s['id'],flush=True)
    extrema_name='plane_extrema.csv' if not only else 'plane_extrema_'+only[0]+'.csv'
    with (OUT/extrema_name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(extrema[0]));w.writeheader();w.writerows(extrema)

def write_index(only=None):
    text='''# 新版高电压场图：6、7、8 V

已核对五个半伏点的实际输入和运行指标；全部达到完整电压与热源因子 1，代数求解及独立通量门限通过，状态为 PASS_NOT_MESH_ACCEPTANCE。

本图采用 6、7、8 V 三个整数电压以匹配三列模板，统一孔半径 3 μm、20 × 20 × 1 μm 几何和 XY0_Z2 网格。不是孔径扫描，也不是旧外推曲线。

热源采用 J·E 总电流功诊断分支，不是此前漂移电导热源分支，不等同于完整热电能量闭合。当前是单网格候选结果，尚未取得该分支的网格/稳定化认证。

每张图上排为二维云图，下排为同一场量的放大三维曲面，曲面高度不是变形。同一物理量在三个电压下使用相同色标和高度范围；坐标统一为 −10 至 10 μm。不调整单张图的数据或颜色范围。

主图均取中面 z=0；另补充下表面 z=−0.5 μm 的电子及上表面 z=+0.5 μm 的空穴图。载流子保留正负号，使用 Δn/n0、Δp/p0 的百分比，不是对数变量 eta，也不取绝对值。

## 应力解释

“耦合总等效应力（含热效应）”是完整应力 σ=C:ε−eᵀE−βθ 的 von Mises 值；不能改称独立热应力。
“热膨胀本构应力项（xx 分量）”是 −βxxθ，负号表示该本构项的压缩方向；它不含位移松弛产生的弹性抵消，不是单独求解的热应力场。

## 载流子解释

n0=6×10²² m⁻³，p0≈1.40×10⁻³⁰ m⁻³。空穴相对扰动较大不表示绝对空穴浓度或其导电贡献显著。表面采样峰值与全域高斯采样最大值的采样位置不同，数值可能不同。

## 数据与质量

原始节点场、指标、积分采样文件、网格及重构程序的哈希与已验收后处理记录逐一核对。复用既有全积分点重构审核，并重新在每个电压的每个截面抽查 11 个单元的 HEX20 计算，检查缓存数据一致性。
保留全部截面三角面片，无跨孔插值，无跨单元平滑；导数场的接缝与中面微弱载流子波纹原样保留，不能仅凭这些局部纹理认定物理振荡。
颜色为用户指定蓝—绿—黄—红模板；颜色映射线性且三电压共用，带符号载流子场零值为绿色。

导出 PNG、PDF、SVG。verified_results.csv 为五个电压的核验摘要，plane_extrema.csv 为图中采样极值；source_data 为未改写的截面数据副本，qa 为图像检查记录。
'''
    (OUT/'说明.md').write_text(text,encoding='utf-8')
    links=[]
    for s in SPECS:
        if only and s['id'] not in only:
            continue
        name='figures/'+s['id']
        links.append(f'<section><h2>{s["title"]}</h2><p><a href="{name}.png">高清 PNG</a> · <a href="{name}.pdf">PDF</a> · <a href="{name}.svg">SVG</a></p><img src="{name}.png"></section>')
    (OUT/'图像索引.html').write_text('<!doctype html><meta charset="utf-8"><title>6–8 V 场云图与放大曲面</title><style>body{max-width:1400px;margin:30px auto;font-family:Arial,"Microsoft YaHei",sans-serif;color:#222}img{width:100%}section{margin:40px 0}p{line-height:1.6}</style><h1>6、7、8 V：场云图与放大曲面</h1><p>J·E 总电流功诊断分支，r=3 μm。同物理量共用坐标、色标和高度范围。完整载荷求解及通量门限通过；单网格结果，未作网格认证。下排高度表示场值。</p><p>总应力与热膨胀本构项分别标示；电子与空穴均保留扰动正负号。</p>'+''.join(links),encoding='utf-8')
    shutil.copy2(__file__,OUT/Path(__file__).name)
    shutil.copy2(fe.__file__,OUT/Path(fe.__file__).name)
    shutil.copy2(STYLE_PATH,OUT/STYLE_PATH.name)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--render-only',action='store_true')
    parser.add_argument('--only',nargs='+',choices=[s['id'] for s in SPECS])
    parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args()
    OUT=args.output.resolve()
    if args.render_only:
        data={(v,face):dict(np.load(OUT/'source_data'/f'V{v}_{face}.npz')) for v in VOLTAGES for face in FACES}
    else:
        if OUT.exists():
            raise FileExistsError('Choose a new --output directory to preserve previous exports: '+str(OUT))
        data=verify_and_load(args.only)
    render(data,args.only)
    write_index(args.only)
