"""Four voltage-response panels: exact LQ scaling dashed line and complete NT solid line."""
from pathlib import Path
import csv
import hashlib
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / 'output/F2_孔洞电压扩展_0_10V_20260909/F2_voltage_sweep_0_10V_source.csv'
OUT = ROOT / 'output/F2_电压0_8V_蓝虚线红实线_20260909'
OUT.mkdir(parents=True, exist_ok=True)
with SRC.open(encoding='utf-8-sig', newline='') as f:
    rows = [r for r in csv.DictReader(f) if float(r['voltage_V']) <= 8]
anchor = next(r for r in rows if float(r['voltage_V']) == 3)
complete = [r for r in rows if r['nt_solve_complete'].lower() == 'true']
vx = np.arange(0, 8.01, .5)
vn = np.array([float(r['voltage_V']) for r in complete])
assert vn[-1] == 6 and all(float(r['last_accepted_joule_factor']) == 1 for r in complete)
specs = [
    ('theta_max_K', 1, 2, r'最大温升 $\theta_{\max}$ (K)', 'global_5x5x5_quadrature_sample_maximum'),
    ('joule_W', 1e3, 2, r'总焦耳功率 $P_J$ (mW)', 'solid_domain_volume_integral_of_drift_joule_source'),
    ('ring_E_mean_V_m', 1e-5, 1, r'孔周平均电场 $\langle |E|\rangle$ (kV/cm)', 'volume_weighted_mean_in_hole_annulus'),
    ('ring_Q_mean_W_m3', 1e-12, 2, r'孔周平均焦耳热密度 $\langle Q_J\rangle$ ($10^{12}$ W m$^{-3}$)', 'volume_weighted_mean_in_hole_annulus'),
]
plt.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['Microsoft YaHei','Arial','DejaVu Sans'],
    'font.size':11,'axes.labelsize':12,'legend.fontsize':10,'axes.linewidth':1,
    'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
    'axes.unicode_minus':False,'svg.fonttype':'none','pdf.fonttype':42,
    'figure.facecolor':'white','savefig.facecolor':'white'})
fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.0))
fig.subplots_adjust(left=.095, right=.975, bottom=.085, top=.97, wspace=.29, hspace=.34)
records = []
for ax, (key, scale, exponent, ylabel, sampling) in zip(axes.flat, specs):
    linear = float(anchor['lq_' + key]) * (vx / 3)**exponent
    nonlinear = np.array([float(r['nt_' + key]) for r in complete])
    assert np.isfinite(linear).all() and np.isfinite(nonlinear).all()
    # All actual complete NT states are included, irrespective of heat-flux audit.
    vdense = np.linspace(0, 8, 321)
    ldense = float(anchor['lq_' + key]) * (vdense / 3)**exponent
    ax.plot(vdense, ldense*scale, linestyle='--', lw=1.8,
            color='#2676D2', label='线性化 LQ', zorder=3)
    ax.plot(vn, nonlinear*scale, linestyle='-', lw=1.8, color='#D94B40',
            label='非线性 NT', zorder=2)
    ax.set(xlim=(0,8), xlabel=r'电压幅值 $V$ (V)', ylabel=ylabel)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    upper = max(max(linear*scale), max(nonlinear*scale))
    ax.set_ylim(-.035*upper, 1.10*upper)
    ax.tick_params(length=4, width=1)
    ax.legend(loc='upper left', frameon=False, handlelength=2.6)
    for v, val in zip(vx, linear):
        records.append({'model':'LQ','voltage_V':v,'metric':key,'value_SI':val,
                        'display_value':val*scale,'sampling':sampling,'status':'exact_scaling_of_3V_FEM'})
    for row, val in zip(complete, nonlinear):
        records.append({'model':'NT','voltage_V':row['voltage_V'],'metric':key,'value_SI':val,
                        'display_value':val*scale,'sampling':sampling,'status':row['nt_status']})
stem = OUT / 'F2_0_8V_蓝虚线_NT实线_四指标'
for ext in ['png','pdf','svg','tiff']:
    options = {'pil_kwargs':{'compression':'tiff_lzw'}} if ext == 'tiff' else {}
    fig.savefig(stem.with_suffix('.'+ext), dpi=600 if ext=='tiff' else 450,
                bbox_inches='tight', **options)
plt.close(fig)
with (OUT/'plot_data.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
with (OUT/'source_status_0_8V.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
contract={'source':str(SRC.relative_to(ROOT)),'sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),
    'geometry_um':[20,20,1],'hole_radius_um':3,'voltage_axis_V':[0,8],
    'LQ':'Blue dashed line, exact scaling of 3 V reference; sampled data exported every 0.5 V; not independent solves.',
    'NT':'Red solid piecewise line through all complete states 0 to 6 V; includes heat-flux audit failures.',
    'missing_NT_V':[7,8],'missing_reason':'No complete full-Joule solution; no extrapolation.',
    'annulus_um':{'radial_range':[3,3.5],'z_range':[-.5,.5]},
    'ring_quadrature':'5x5x3 element quadrature, physical-radius selection and volume weights',
    'temperature_maximum':'maximum over 5x5x5 sample points in every element; location may vary with voltage',
    'total_power':'Integral of solver drift-conduction Joule source; not assumed equal to terminal VI.',
    'styling':'Reference-style four-panel layout; white background, inward ticks, no grid, blue dashed line/red solid line.'}
(OUT/'figure_contract.json').write_text(json.dumps(contract,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'采集位置与图注.md').write_text('''# 0–8 V 四指标图

20 × 20 × 1 μm 中心贯穿圆孔 ZnO 薄板，孔半径 3 μm；左端 0 V、右端 −V，两端温升为零。

蓝色虚线：线性化 LQ，以 3 V FEM 结果精确缩放，覆盖 0–8 V；电场随 V 缩放，热源与温升随 V² 缩放。
红色实线：NT 真实完整解，0–5 V 每隔 0.5 V，另有 6 V；4、4.5、5、6 V 虽未通过热流审计，仍全部纳入实线。没有 5.5 V 的独立 NT 点。
7、8 V 无完整全焦耳热解，所以红线到 6 V；未补点、外推或把部分加载状态当作最终解。图上颜色/线型仅区分模型，不代表审计通过。

## 采集区域

| 子图 | 采集方式 |
|---|---|
| 左上：最大温升 | 全部实体单元 5×5×5 独立积分采样点中的最大温升；不是固定坐标点，最大位置可随电压改变。 |
| 右上：总焦耳功率 | 对整个实体体积积分程序定义的漂移电导焦耳热密度，孔洞排除；不是孔壁单点或端口 VI。 |
| 左下：孔周平均电场 | 孔外 0.5 μm 宽、贯穿整个厚度的环带内，对电场模长作体积加权平均。 |
| 右下：孔周平均焦耳热密度 | 同一环带内焦耳热密度的体积加权平均。 |

环带范围：3 ≤ sqrt(x²+y²) ≤ 3.5 μm，−0.5 ≤ z ≤ 0.5 μm。环带统计采用 5×5×3 单元积分点及其体积权重，不是中面圆周平均。

当前为 L0 数据。详细求解/审计状态保留于 source_status_0_8V.csv。
''',encoding='utf-8')
shutil.copyfile(__file__,OUT/Path(__file__).name)
assert '<text' in stem.with_suffix('.svg').read_text(encoding='utf-8')
print('DONE: four panels, PNG/PDF/SVG/TIFF, blue dashed LQ and solid NT, data CSV and sampling definitions.')
