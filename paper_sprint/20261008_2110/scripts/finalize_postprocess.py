"""Provenance-backed report, geometry checks, analytical evidence and schematic."""
import postprocess_3v as p
import json, csv, numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Circle
from hex20_geometry import sample
import itertools

def main():
    results=json.loads((p.OUT/'data/metrics_summary.json').read_text(encoding='utf-8'))['results'];geometry={};registry=[]
    for c in p.CASES:
        key=str(c['mesh'])
        if key not in geometry:
            x,conn,_=p.read_mesh(c['mesh']);gp=np.polynomial.legendre.leggauss(5)[0];points=[]
            for axis,side in itertools.product(range(3),[-1.,1.]):
                other=[k for k in range(3) if k!=axis]
                for a,b in itertools.product(gp,repeat=2):
                    pt=np.zeros(3);pt[axis]=side;pt[other]=[a,b];points.append(pt)
            geometry[key]=sample(x,conn,points);geometry[key]['mesh_sha256']=p.sha(c['mesh'])
            if geometry[key]['status']!='PASS_SAMPLED':raise ValueError('Boundary geometry sample failed '+key)
        r=next(r for r in results if r['case']==c['name'])
        registry.append(dict(case=r['case'],model=r['model'],voltage_V=3,geometry=r['geometry'],integral_scope=r['integral_scope'],
            provenance_status='PASS',geometry_status='PASS',solve_status='PASS',flux_status='PASS',spatial_status='FAIL' if r['family'].startswith('R23') else 'PENDING',model_scope_status='PENDING',claim_status='PENDING',
            geometry_detail='5x5 samples on all 6 HEX20 faces plus 5x5x5 volume; not global positivity proof',spatial_detail=r['status'],claim_detail='Exploratory numerical response only; no device validation or complete local-field certificate',source_path=str(c['root'].resolve())))
    p.dump(p.OUT/'audit/postprocess/boundary_geometry_samples.json',geometry);p.csvwrite(p.OUT/'data/case_registry.csv',registry)
    # Saved analytical tests remain distinct from 3 V production cases.
    an=p.SOLVER/'runs/F1_20260906';a1p=an/'A1/A1_export_stdout.txt';a2p=an/'A2/A2_export_stdout.txt';a1=p.metrics(a1p);a2=p.metrics(a2p)
    rows1=[dict(elements=n,phi_L2=a1[f'a1_ne{n}_phi_l2'],phi_H1=a1[f'a1_ne{n}_phi_h1'],eta_n_L2=a1[f'a1_ne{n}_eta_n_l2'],eta_p_L2=a1[f'a1_ne{n}_eta_p_l2']) for n in [1,2,4,8,16]]
    rows2=[dict(level=i,phi_L2=a2[f'a2_L{i}_phi_relative_L2'],E_L2=a2[f'a2_L{i}_E_relative_L2'],wall_normal_RMS=a2[f'a2_L{i}_wall_normal_E_RMS_over_E0']) for i in range(3)]
    p.csvwrite(p.OUT/'data/A1_error_source.csv',rows1);p.csvwrite(p.OUT/'data/A2_error_source.csv',rows2)
    p.dump(p.OUT/'audit/postprocess/analytic_sources.json',dict(paths={str(f.resolve()):p.sha(f) for f in [a1p,a2p,an/'A1/A1_stdout.txt',an/'A2/A2_stdout.txt']},evidence='Saved historical analytical logs; no new solve. Synthetic A1/A2, not full ZnO device validation.'))
    fig,ax=plt.subplots(1,2,figsize=(7.2,3.2),layout='constrained')
    for key,lab in [('phi_L2',r'$\phi$: $L^2$'),('phi_H1',r'$\phi$: $H^1$ seminorm')]:ax[0].loglog([r['elements'] for r in rows1],[r[key] for r in rows1],'o-',label=lab)
    ax[0].set(xlabel='Elements along x',ylabel='Relative error',title='A1 saved coupled analytical test');ax[0].legend(fontsize=7)
    for key,lab in [('phi_L2',r'$\phi$: $L^2$'),('E_L2',r'$E$: $L^2$'),('wall_normal_RMS',r'$E_n$: wall RMS / $E_0$')]:ax[1].semilogy([r['level'] for r in rows2],[r[key] for r in rows2],'o-',label=lab)
    ax[1].set(xlabel='Mesh level',ylabel='Normalized error',xticks=[0,1,2],title='A2 saved insulated-hole test');ax[1].legend(fontsize=7)
    fig.suptitle('Historical synthetic analytical evidence; actual errors\nDoes not certify 3 V production fields (L0 exploratory / mesh pending)',fontsize=9);p.savefig(fig,'fig02_analytic_errors')
    fig=plt.figure(figsize=(7.2,4.9),layout='constrained');gs=fig.add_gridspec(2,2,width_ratios=[1,1.1]);ax=fig.add_subplot(gs[:,0]);info=fig.add_subplot(gs[0,1]);side=fig.add_subplot(gs[1,1])
    ax.add_patch(Rectangle((-10,-10),20,20,fc='#e9eff3',ec='#344b60',lw=1));ax.add_patch(Circle((0,0),3,fc='white',ec='#344b60',lw=1));ax.plot([0],[3.25],'o',color='#c76732',ms=4);ax.annotate('Thickness line\n(0, 3.25 µm)',xy=(0,3.25),xytext=(1,7),fontsize=8,arrowprops={'arrowstyle':'->'});ax.annotate('a = 3 µm',xy=(3,0),xytext=(4,-3),fontsize=8,arrowprops={'arrowstyle':'->'});ax.text(-10.8,0,'φ = 0\nθ = 0\nηn = ηp = 0',ha='right',va='center',fontsize=8);ax.text(10.8,0,'φ = −3 V\nθ = 0\nηn = ηp = 0',ha='left',va='center',fontsize=8)
    ax.set(xlim=(-17,18),ylim=(-12,13),aspect='equal',xlabel='x (µm)',ylabel='y (µm)',title='H1 plan view; H0 omits the hole');ax.set_xticks([-10,0,10]);ax.set_yticks([-10,0,10])
    info.axis('off');info.text(0,1,'Schematic, not a computed field\n20 × 20 × 1 µm plate; c-axis +z\n\nOther exterior surfaces and hole wall:\ntraction = 0; D·ν = 0\nJn·ν = Jp·ν = 0; heat flux·ν = 0\n\nFull-domain rigid-body anchors (µm):\nA=(−10,−10,−0.5): ux=uy=uz=0\nB=(−10,+10,−0.5): ux=uz=0\nC=(−10,−10,+0.5): ux=0',va='top',fontsize=8,linespacing=1.35)
    side.add_patch(Rectangle((-10,-.5),20,1,fc='#e9eff3',ec='#344b60'));side.plot([0,0],[-.5,.5],color='#c76732',lw=2);side.annotate('c-axis +z',xy=(7,.9),xytext=(7,.1),arrowprops={'arrowstyle':'->'},fontsize=8);side.set(xlim=(-11,11),ylim=(-1,1.25),xlabel='x (µm)',ylabel='z (µm)',title='Section y=3.25 µm (outside the hole)');side.set_yticks([-.5,0,.5])
    fig.suptitle('3 V model and sampling schematic; production mesh validation remains pending',fontsize=9);p.savefig(fig,'fig01_geometry_boundary')
    # Tables and explanatory prose are generated from the delivered numeric data.
    lookup={r['case']:r for r in results};n0=lookup['FAST_H1_N0_L0'];nt=lookup['FAST_H1_NT_L0'];h0=lookup['FAST_H0_NT_L0']
    contrasts=dict(H1_NT_vs_N0_theta_percent=100*(nt['theta_max_saved_K']/n0['theta_max_saved_K']-1),H1_NT_vs_N0_Pd_percent=100*(nt['Pd_W']/n0['Pd_W']-1),NT_H1_vs_H0_theta_percent=100*(nt['theta_max_saved_K']/h0['theta_max_saved_K']-1),NT_H1_vs_H0_Pd_percent=100*(nt['Pd_W']/h0['Pd_W']-1),constant_reference_peak_K=p.SIGMA*9/(8*p.K))
    p.dump(p.OUT/'data/contrasts_exploratory.json',contrasts)
    lines=['# 3 V 零新增求解后处理报告','','本轮重构10个已完成状态：6个F1全域L0（H0/H1×LQ/N0/NT）、3个R23半域NT网格层级及1个R25 A64半域NT。没有调用求解器，没有增加电压、几何或物理。F1的同模型几何对照与同几何模型对照均为历史粗网格探索，未完成空间验收。R23/R25采用试验LOG_WEIGHTED格式，不能当作与F1相同离散格式的复算。','','## 证据与重构','','input、mesh、nodes、field_samples、metrics及exe逐文件SHA-256保存在audit/postprocess/各case.json；源码快照逐文件哈希与对应归档manifest复核，详见provenance_gate_recheck.json。原始结果所在runs是工作区链接，审计保存resolve后的实际路径；原文件仅只读。10例原exit code为0，残差/增量/阶段平衡/独立通量/4–5阶热源积分按原阈值重新数值核查。本次归档身份检查不等价于从源码重新编译。', '',f"全体第一5阶采样点电场/热源/温度/浓度/位移/迁移率/应力重构最大尺度误差{max(r['saved_sample_scaled_error'] for r in results):.3g}；独立5×5×5积分漂移热相对差最大{max(r['drift_heat_reconstruction_relative'] for r in results):.3g}。形函数先插值η再指数；LQ使用1+η和参考密度电流/热源，N0冻结μ与D，NT使用当地T。",'', '几何通过5阶体点及每单元六面5×5点正Jacobian抽样；这是一项必要的局部映射检查，不是全域正性数学证明。废弃旧graded-hole不在输入清单。厚度和中面重构对每个物理取样点求逆等参映射，不作域外外推。','','## 功率与参考偏离（F1 L0，探索性，空间待验收）','','| Case | Pd / mW | Pdiff / mW | ∫abs(Jdiff·E) / mW | Pterm / mW | Qout / mW | 采样峰温升 / K | rT RMS / K |','|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in results[:6]:lines.append(f"|{r['case']}|{r['Pd_W']*1e3:.8f}|{r['Pdiff_W']*1e3:.8f}|{r['Pdiff_abs_W']*1e3:.8f}|{r['Pterm_all_boundary_W']*1e3:.8f}|{r['heat_out_all_boundary_W']*1e3:.8f}|{r['theta_max_saved_K']:.8f}|{r['reference_rms_K']:.8f}|")
    lines += ['',f"H1 NT相对N0的采样峰温升变化{contrasts['H1_NT_vs_N0_theta_percent']:.6f}%，漂移功率变化{contrasts['H1_NT_vs_N0_Pd_percent']:.6f}%；同为NT的H1相对H0采样峰温升变化{contrasts['NT_H1_vs_H0_theta_percent']:.6f}%。这些是新核对原始场后的历史粗网格趋势，不进入正式摘要的定量结论。H1 NT的绝对扩散做功/Pd={nt['Pdiff_abs_over_Pd']*100:.6f}%，高于上述微小反馈量级；不能用约化漂移热源直接证明真实器件的精细热反馈机制。",'',f"常系数参考θref=σ0(φ−φL)(φR−φ)/(2k)，其3V理论峰值为{contrasts['constant_reference_peak_K']:.9f} K。rT以5阶体积权重计算；该偏离包含模型与离散差异，不是纯半导体贡献。原metrics中的theta_max是保存的积分采样最大值（不等价连续场严格最大值）。",'', '半域功率保留真实半域值，表中字段分别列全部边界和物理外边界，不自动乘2。P_JdotE与Pterm的差单独保存；未将它强行设为零。漂移/扩散重构没有修改热方程。','','## Pe与空间证据','','Pe使用20节点电势范围，非二次场内部精确极差。N0使用固定T0对应D/μ；NT使用节点平均T；LQ的Pe仅为电势背景指标。区域的最终定义为：接触区优先选择接触x=±10μm的第一层单元；其余按中心|z|≥0.45μm为表层，再按r≤3.5μm为孔周，其余为体内粗区。上下表层合并；每区计数与体积比例分别输出。初始|x中心|≥9μm出现空区，已在region_definition_correction.json记录诊断修正；全域Pe和物理指标没有变化。','','| Case | Pe max | 元素计数Pe>1 | 体积权重Pe>1 | 空间状态 |','|---|---:|---:|---:|---|']
    for r in results:lines.append(f"|{r['case']}|{r['Pe_max']:.5f}|{100*r['Pe_count_gt1']:.3f}%|{100*r['Pe_volume_gt1']:.3f}%|{r['status']}|")
    lines += ['', 'R23三层全局/环带指标原通过不能覆盖局部场失败。L1→L2映射电子与空穴扰动最大差分别3.0065437%与3.0097661%，超过原2%门槛，保留FAIL。Pe>1也是风险诊断，不是任意HEX20的误差界；LOG_WEIGHTED结果不能只凭Pe获得PASS。','','## 厚度、力学与图源','','厚度线为(x,y)=(0,3.25μm)，距r=3μm孔壁0.25μm；H0使用相同坐标。每条线含401个均匀z点并加入真实网格层面，保存ηn/ηp、n/p相对扰动、φ、Ez、θ及参考偏离。场梯度为所选单元内迹，接缝保留。孔周应力在r=3.25μm、z=0取样并变换成σrr和σθθ，不作节点平滑。厚度变化uz(top)−uz(bottom)保留真实均匀应变；相对位移仅减厚度线中点平移。半域与全域刚体约束不同，原始位移绝对峰值不直接作版本差异。', '',f"H1 NT厚度变化{nt['thickness_change_m']*1e12:.6f} pm，线上相对位移最大{nt['relative_displacement_line_max_m']*1e12:.6f} pm；孔外中面σθθ绝对采样最大{nt['ring_midplane_max_abs_stt_Pa']/1e6:.6f} MPa。它们均未做独立空间验收，不据此声称ZnO屈服或失效。",'', '图件清单（均有PNG/SVG/PDF）：', '', '- fig01_geometry_boundary：模型/边界/孔外厚度线；示意与结果分开。','- fig02_analytic_errors：从保存A1/A2日志读取真实误差，无新求解；合成解析例不替代生产验证。','- fig03_peclet_regions：同3V不同已存网格的风险分布。','- fig03b_spatial_failure：R23 3V局部载流子3.01% FAIL。','- fig04_midplane_fields：F1 H0/H1 NT中面E/J/θ共同色标，三角片不跨原单元和孔洞。','- fig06_thickness_profiles：F1六例真实材料线。','- fig07_mechanics：未平滑应力与相对位移。','- fig08_power_and_reference：漂移/扩散/端口功率与电热参考偏离。', '', 'CSV图源位于data/metrics_summary.csv、peclet_regions.csv、*_element_peclet.csv、*_midplane.csv及独立triangles.csv、profiles/*.csv、A1_error_source.csv与A2_error_source.csv。fig03b原始来源为runs/HV_SUPG_20261007/field_shape_r23_L1_L2_V030.json，已在图件清单记录哈希。','','## 阻塞与允许表述','','1. F1六例没有完整三层空间证据；R23局部载流子门槛明确FAIL。正式局部场与小比例反馈结论BLOCKED_SPATIAL；图用于定位与方法展示。','2. H0旧环带积分欠采样约28.48%，本轮删除H0/H1环带增强比，不拿色图替代区域积分验证。','3. NT/N0在新R25网格缺少同条件N0状态，不能从旧L0模型差值推断新网格模型差已收敛。','4. 材料参数未实验标定，热源未包含完整漂移–扩散能量闭合；真实器件机制BLOCKED_MODEL_SCOPE。','5. 本轮没有求解独立常导电生产例，解析参考只是解释/诊断工具；LQ不得冒充严格常导电验证。','6. 资源与现有活动任务限制下未新增三维求解，未中断既有高压；未扩新物理/几何/电压。','','## 实际执行与改动','','运行postprocess_3v.py（退出0；153.77秒）、supplemental_3v.py（退出0）、refine_pe_regions.py（退出0）、本脚本finalize_postprocess.py。全部为新目录脚本；没有修改src、历史输入/结果、exe、分支或标签。首次草稿脚本N0/LQ热电压数组广播报错已修复并完整重跑；该后处理实现错误不改写原算例状态。可重复命令如下：','','```text','python -B solver_v2/paper_sprint/20261008_2110/scripts/postprocess_3v.py','python -B solver_v2/paper_sprint/20261008_2110/scripts/supplemental_3v.py','python -B solver_v2/paper_sprint/20261008_2110/scripts/refine_pe_regions.py','python -B solver_v2/paper_sprint/20261008_2110/scripts/finalize_postprocess.py','```','']
    (p.OUT/'analysis/postprocess_report.md').write_text('\n'.join(lines),encoding='utf-8')
    p.dump(p.OUT/'audit/postprocess/figure_sources.json',{'spatial_failure':dict(path=str((p.HV/'field_shape_r23_L1_L2_V030.json').resolve()),sha256=p.sha(p.HV/'field_shape_r23_L1_L2_V030.json')),'geometry_contract':dict(path=str((p.SOLVER/'src/app/f1_case_runner_m.f90').resolve()),sha256=p.sha(p.SOLVER/'src/app/f1_case_runner_m.f90'))})
    manifest=[]
    for folder in ['data','figures','scripts','analysis','audit/postprocess']:
        for file in sorted((p.OUT/folder).rglob('*')):
            if file.is_file() and file.name!='hash_manifest.json' and '__pycache__' not in str(file):manifest.append(dict(path=file.relative_to(p.OUT).as_posix(),bytes=file.stat().st_size,sha256=p.sha(file)))
    p.dump(p.OUT/'audit/postprocess/hash_manifest.json',manifest)
    print('FINALIZED',len(manifest),'files',flush=True)
if __name__=='__main__':main()
