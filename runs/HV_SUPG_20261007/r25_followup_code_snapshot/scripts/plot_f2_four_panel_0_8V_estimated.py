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
OUT = ROOT / 'output/F2_电压0_8V_含NT外推估计_20260909'
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
estimate_summary = []
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
            label='非线性 NT（已计算）', zorder=2)
    # Empirical trend extrapolation, explicitly distinct from a solved FEM state.
    # Anchor the quadratic at the actual 6 V value; fit actual 3.5--6 V states.
    mask = vn >= 3.5
    delta = vn[mask] - 6
    coefficients = np.linalg.lstsq(np.column_stack((delta, delta**2)),
                                   nonlinear[mask] - nonlinear[-1], rcond=None)[0]
    def predict(voltage, coeff=coefficients):
        d = np.asarray(voltage) - 6
        return nonlinear[-1] + coeff[0]*d + coeff[1]*d*d
    ve = np.linspace(6, 8, 101)
    ye = predict(ve)
    ax.plot(ve, ye*scale, '-.', color='#D94B40', lw=1.8,
            label='NT 外推估计', zorder=2)
    ax.plot([7,8], predict([7,8])*scale, linestyle='none', marker='D',
            ms=5.5, mfc='white', mec='#D94B40', mew=1.2, clip_on=False, zorder=4)
    variations = []
    for begin in [2.5,3,3.5,4]:
        m = vn >= begin
        d = vn[m]-6
        c = np.linalg.lstsq(np.column_stack((d,d*d)), nonlinear[m]-nonlinear[-1], rcond=None)[0]
        variations.append(predict([7,8],c))
    variations = np.asarray(variations)
    for i, voltage in enumerate([7,8]):
        val = float(predict(voltage))
        records.append({'model':'NT_EXTRAPOLATED','voltage_V':voltage,'metric':key,'value_SI':val,
                        'display_value':val*scale,'sampling':sampling,'status':'EMPIRICAL_ESTIMATE_NOT_FEM'})
        estimate_summary.append({'metric':key,'voltage_V':voltage,'estimate_SI':val,
            'display_estimate':val*scale,'window_sensitivity_low_SI':float(variations[:,i].min()),
            'window_sensitivity_high_SI':float(variations[:,i].max()),
            'anchor_6V_SI':float(nonlinear[-1]),'linear_coefficient_SI':float(coefficients[0]),
            'quadratic_coefficient_SI':float(coefficients[1])})
    ax.set(xlim=(0,8), xlabel=r'电压幅值 $V$ (V)', ylabel=ylabel)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    upper = max(max(linear*scale), max(nonlinear*scale), max(ye*scale))
    ax.set_ylim(-.035*upper, 1.10*upper)
    ax.tick_params(length=4, width=1)
    ax.legend(loc='upper left', frameon=False, handlelength=2.6)
    for v, val in zip(vx, linear):
        records.append({'model':'LQ','voltage_V':v,'metric':key,'value_SI':val,
                        'display_value':val*scale,'sampling':sampling,'status':'exact_scaling_of_3V_FEM'})
    for row, val in zip(complete, nonlinear):
        records.append({'model':'NT','voltage_V':row['voltage_V'],'metric':key,'value_SI':val,
                        'display_value':val*scale,'sampling':sampling,'status':row['nt_status']})
stem = OUT / 'F2_0_8V_蓝虚线_NT实线与外推_四指标'
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
    'missing_solved_NT_V':[7,8],'estimate_NT_V':[7,8],
    'extrapolation':'y(V)=y(6)+a(V-6)+b(V-6)^2, least squares on complete NT states 3.5,4,4.5,5,6 V; not FEM.',
    'estimate_integrity':'No full-coupled state or conservation validation at 7/8 V; fitting-window sensitivity is not a confidence interval.',
    'annulus_um':{'radial_range':[3,3.5],'z_range':[-.5,.5]},
    'ring_quadrature':'5x5x3 element quadrature, physical-radius selection and volume weights',
    'temperature_maximum':'maximum over 5x5x5 sample points in every element; location may vary with voltage',
    'total_power':'Integral of solver drift-conduction Joule source; not assumed equal to terminal VI.',
    'styling':'LQ blue dashed; actual complete NT red solid; extrapolated NT red dash-dot with hollow diamonds at 7/8 V.'}
(OUT/'NT_extrapolation_details.json').write_text(json.dumps(estimate_summary,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'figure_contract.json').write_text(json.dumps(contract,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'采集位置与图注.md').write_text('''# 0–8 V 四指标图

20 × 20 × 1 μm 中心贯穿圆孔 ZnO 薄板，孔半径 3 μm；左端 0 V、右端 −V，两端温升为零。

蓝色虚线：线性化 LQ，以 3 V FEM 结果精确缩放，覆盖 0–8 V；电场随 V 缩放，热源与温升随 V² 缩放。
红色实线：NT 真实完整解，0–5 V 每隔 0.5 V，另有 6 V；4、4.5、5、6 V 虽未通过热流审计，仍全部纳入实线。没有 5.5 V 的独立 NT 点。
7、8 V 无完整全焦耳热解。按用户要求，以 3.5、4、4.5、5、6 V 完整 NT 数据进行二次趋势外推，约束预测曲线通过实际 6 V 值。6–8 V 红色点划线和 7、8 V 空心菱形明确表示估计，不是 FEM 结果。三项热源/温升指标与电场指标分别拟合，未检验其联合物理一致性，不能用于证明稳态解存在或热稳定性。

外推公式：y(V)=y(6)+a(V−6)+b(V−6)²，a、b 为最小二乘拟合系数。另比较 2.5、3、3.5、4 V 起始窗口，敏感性结果保存在 NT_extrapolation_details.json；该范围不是置信区间，也不覆盖求解误差和模型误差。图中已计算完整解包含热流审计未通过数据，因此估计仅用于趋势示意。

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
