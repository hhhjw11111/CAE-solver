"""Create the final F1 PDF from actual saved evidence; incomplete runs require --partial."""
from pathlib import Path
import argparse,hashlib,json,sys
from reportlab.platypus import BaseDocTemplate,Frame,PageTemplate,PageBreak,Spacer,Image,Table,TableStyle
from reportlab.lib import colors
from pypdf import PdfReader
from build_case_freeze_pdf import p,STYLES,WIDTH,HEIGHT,MARGIN,AVAILABLE,NAVY,TEAL,GREY,LIGHT,LINE

ROOT=Path(__file__).resolve().parents[1];RUNS=ROOT/'solver_v2/runs/F1_20260906'
FIG=ROOT/'output/figures/F1_20260906'
CASES=[('H0-LQ','MAIN_H0_LQ_L2'),('H0-N0','MAIN_H0_N0_L2'),('H0-NT','MAIN_H0_NT_L2'),
       ('H1-LQ','MAIN_H1_LQ_L2'),('H1-N0','MAIN_H1_N0_L2'),('H1-NT','MESH_H1_NT_L2')]

def kv(path):
    out={}
    if not path.exists():return out
    for line in path.read_text(encoding='utf-8-sig').splitlines():
        if '=' not in line or line.startswith(('stage=','balance_detail_point=')):continue
        k,v=line.split('=',1)
        try:out[k]=float(v)
        except ValueError:out[k]=v.strip()
    return out

def js(path):return json.loads(path.read_text(encoding='utf-8-sig')) if path.exists() else {}
def fmt(v,scale=1):return '未完成' if v is None else f'{v*scale:.5g}'
def yes(v):return '通过' if v else '未通过 / 未完成'
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def table(rows,widths=None):
    widths=widths or [1/len(rows[0])]*len(rows[0])
    tab=Table([[p(str(c),'headcell' if i==0 else 'cell') for c in row] for i,row in enumerate(rows)],
       colWidths=[w*AVAILABLE for w in widths],repeatRows=1,hAlign='LEFT')
    tab.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),NAVY),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,LIGHT]),
        ('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),
        ('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('LINEBELOW',(0,-1),(-1,-1),.5,LINE)]))
    return tab

def fig(name):
    path=FIG/(name+'.png')
    if not path.exists():return p('该图尚未生成；不以预期图替代计算结果。','small')
    from PIL import Image as PILImage
    w,h=PILImage.open(path).size
    return Image(str(path),width=AVAILABLE,height=AVAILABLE*h/w)

def header(canvas,doc):
    canvas.saveState();canvas.setFillColor(TEAL);canvas.rect(MARGIN,HEIGHT-36,26,3,fill=1,stroke=0)
    canvas.setFont('Hei',8);canvas.setFillColor(GREY)
    canvas.drawString(MARGIN+35,HEIGHT-35,'TPSolverV2 / F1 计算结果与证据')
    canvas.setStrokeColor(LINE);canvas.line(MARGIN,39,WIDTH-MARGIN,39)
    canvas.setFont('Song',8);canvas.drawString(MARGIN,25,'候选有效材料数值研究；不等同于实验标定或真实器件验证')
    canvas.setFont('Hei',9);canvas.drawRightString(WIDTH-MARGIN,25,f'{doc.page:02d}');canvas.restoreState()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--partial',action='store_true');args=parser.parse_args()
    a1=kv(RUNS/'A1/A1_export_stdout.txt');a2=kv(RUNS/'A2/A2_export_stdout.txt')
    small=js(RUNS/'small_comparison.json');mesh=js(RUNS/'mesh_comparison.json')
    main_data={label:kv(RUNS/case/'metrics.txt') for label,case in CASES}
    repeat_ok=False;repeat_records={}
    if all((RUNS/'REPEAT_H1_NT_L2'/name).exists() for name in ['metrics.txt','nodes.csv','field_samples.csv']):
        for name in ['metrics.txt','nodes.csv','field_samples.csv']:
            repeat_records[name]=[digest(RUNS/'MESH_H1_NT_L2'/name),digest(RUNS/'REPEAT_H1_NT_L2'/name)]
        repeat_ok=all(a==b for a,b in repeat_records.values())
    complete=(a1.get('verification_result')=='PASS' and a2.get('verification_result')=='PASS' and
        small.get('passed',False) and mesh.get('passed',False) and repeat_ok and
        all(m.get('case_acceptance')=='PASS_NOT_MESH_ACCEPTANCE' for m in main_data.values()))
    if not complete and not args.partial:raise SystemExit('Frozen execution is incomplete; final report refused (use --partial only for a labelled progress report).')
    story=[p('孔洞多场耦合研究','title'),p('F1 计算结果与证据报告','subtitle'),
        p('状态：'+('冻结计算计划已完成' if complete else '阶段结果；完整计划尚未完成'),'h3'),
        p('本报告仅汇总真实保存的数值结果。两类解析解验证数值实现；候选 ZnO 主算例用于比较孔洞、非线性及温度相关输运的共同影响。材料尚未完成同一试样的实验标定。'),
        table([['阶段','本轮范围','状态'],['A1','一维热-压电-双载流子，5 层网格',yes(a1.get('verification_result')=='PASS')],
            ['A2','圆孔绝缘边界解析，3 层网格',yes(a2.get('verification_result')=='PASS')],
            ['小载荷','H1-LQ / NT，0.1 V',yes(small.get('passed',False))],
            ['主网格','H1-NT 三层边界层网格',yes(mesh.get('passed',False))],
            ['主矩阵','H0/H1 × LQ/N0/NT，3 V',yes(all(m.get('case_acceptance')=='PASS_NOT_MESH_ACCEPTANCE' for m in main_data.values()))],
            ['代表性复现','最终 H1-NT 一次复算、文件逐字一致',yes(repeat_ok)]],[.18,.56,.26]),Spacer(1,12),
        p('已完成的旧材料 G12/G13 和核心门禁继续作为内部历史证据；它们不替代本轮 ZnO 计算，也不追加为论文的第三类解析解。'),
        p('执行中只处理直接阻碍当前算例的问题：解析数据导出、零热源时的量纲归一化，以及细网格所需的内存管理。未追加半径、材料、热导率或迁移率指数扫描。')]
    story += [PageBreak(),p('01 / A1 联合解析对比','title'),
       p('1 × 1 × 1 μm 合成一维棱柱，焦耳热与复合关闭；左端位移、电势、两类扰动和温升为零；右端零牵引、零两类电流、电势为零、温升 1 K。横向位移固定只用于这一维验证。'),
       fig('A1_joint_analytic'),Spacer(1,6),
       table([['x 向单元数','电势相对 L² 误差','电势 H¹ 半范数误差']]+[[n,fmt(a1.get(f'a1_ne{n}_phi_l2')),fmt(a1.get(f'a1_ne{n}_phi_h1'))] for n in [1,2,4,8,16]],[.2,.4,.4]),
       p('最终收敛阶：L²='+fmt(a1.get('a1_phi_l2_order'))+'，H¹='+fmt(a1.get('a1_phi_h1_order'))+'。电子、空穴误差与电势一致到浮点精度；最大归一化零电流残差为 '+fmt(a1.get('a1_maximum_normalized_current'))+'。','small')]
    story += [PageBreak(),p('02 / A2 圆孔解析对比','title'),
       p('20 × 20 × 1 μm 方板，贯通孔半径 3 μm；四条外侧面施加 φ=−E₀x(1+a²/r²)，E₀=10⁴ V/m；孔壁及上下表面为零电学/载流子法向通量。温度和位移固定，压电关闭。'),
       fig('A2_hole_analytic'),Spacer(1,6),
       table([['网格','φ-L²','E-L²','孔壁 En-RMS/E₀']]+[[f'L{i}',fmt(a2.get(f'a2_L{i}_phi_relative_L2')),fmt(a2.get(f'a2_L{i}_E_relative_L2')),fmt(a2.get(f'a2_L{i}_wall_normal_E_RMS_over_E0'))] for i in range(3)]),
       p('最细网格孔壁电场、焦耳热峰值倍率分别为 '+fmt(a2.get('a2_L2_wall_Et_peak_over_E0'))+' 与 '+fmt(a2.get('a2_L2_wall_Q_peak_over_Qinf'))+'，解析值为 2 和 4。这些倍率只适用于本约化解析算例，不强加给真实非线性主模型。','small')]
    story += [PageBreak(),p('03 / 三维模型与材料','title'),
       p('H0：20 × 20 × 1 μm 无孔方板。H1：同尺寸，中心贯通圆孔 a=3 μm。环带为 3≤r≤3.5 μm；H0 使用同位置虚拟环带。晶体 c 轴沿 +z，电压沿 x；不施加 z 镜像对称。'),
       table([['位置','冻结边界条件'],['x=−10 μm','φ=0，T=300 K，n=n₀，p=p₀'],['x=+10 μm','φ=−3 V，T=300 K，n=n₀，p=p₀'],
         ['其余外表面及孔壁','D·ν=Jn·ν=Jp·ν=0；热绝缘；零牵引'],
         ['最小刚体约束','A(−10,−10,−0.5)：ux=uy=uz=0；B(−10,+10,−0.5)：ux=uz=0；C(−10,−10,+0.5)：ux=0。坐标单位 μm。']],[.24,.76]),
       p('候选材料核心参数','h3'),
       table([['参数','数值 / 假设'],['弹性与压电','C11=209.7、C12=121.1、C13=105.1、C33=210.9、C44=42.47 GPa；e31=−0.62、e33=0.96、e15=−0.37 C/m²'],
         ['介电与热学','εr=(8.54,8.54,10.2)；α=(4.31,4.31,2.49)×10⁻⁶/K；一次热释电 pz=−6.9×10⁻⁶ C/(m² K)；k=60 W/(m K)'],
         ['载流子','n₀=6×10²² m⁻³；p₀=1.40167×10⁻³⁰ m⁻³；不人为抬高空穴浓度'],
         ['输运','μn0=0.0205、μp0=0.005 m²/(V s)；βn=βp=2；NT 中 D=kBTμ/q']],[.24,.76]),
       p('LQ：线性相对浓度、μ/D 固定；先等温机电载流子，再二阶热源与热场，最后仅更新热位移。N0：完整非线性，但 μ 和 D 同时固定。NT：完整非线性且 μ(T)、D(T) 均更新。'),
       p('热源冻结为漂移电导型 QJ；完整漂移-扩散电流含扩散分量，因此本模型不要求 IV 与该约化热源积分处处相等。')]
    story += [PageBreak(),p('04 / 小载荷与网格验收','title'),
       table([['小载荷项目','实际差值 / 指标']]+[[label,fmt(small.get(k),100)+'%'] for label,k in [
         ('电势相对差','phi_relative_L2_nodes'),('位移相对差','u_relative_L2_nodes'),('最大载流子偏移','maximum_density_shift'),
         ('Q/V² 相对差','Q_over_V2_relative_difference'),('θ/V² 相对差','theta_over_V2_relative_difference')]] ,[.58,.42]),
       p('本几何的等温 LQ 载流子变化接近零：以 n₀、p₀ 作为物理浓度比较尺度，不用接近零的线性浓度扰动作分母。'),
       p('H1-NT 主网格：最细两层比较','h3'),
       table([['指标','L1→L2 变化','阈值']]+[[k,fmt(v.get('L1_L2'),100)+'%',fmt(v.get('limit'),100)+'%'] for k,v in mesh.get('metric_changes',{}).items()],[.5,.25,.25]),
       p('三层孔洞网格单元数为 2240、4608、8064，首层厚度 14、7、3.5 nm。输入细化尺度不是已解析屏蔽层的证明，以上验收结果才决定是否使用最终网格。'),
       p('零热源阶段的守恒指标采用预先可计算的无孔恒电导参考功率归一化，并保存实际残差（W）与原始分母；开启热源后仍使用原守恒判据。这不修改 PDE 或放宽 10⁻⁸ 阈值。','small')]
    rows=[['结构 / 模型','功率 mW','峰值温升 K','环带 Q 均值 10¹² W/m³','环带电子偏移 %','最小 μ/μ₀']]
    for label,_ in CASES:
        m=main_data[label];rows.append([label,fmt(m.get('joule_W'),1e3),fmt(m.get('theta_max_K')),
          fmt(m.get('ring_Q_mean_W_m3'),1e-12),fmt(m.get('ring_abs_n_relative_mean'),100),fmt(m.get('minimum_mu_n_ratio'))])
    story += [PageBreak(),p('05 / 六个主算例结果','title'),table(rows,[.17,.14,.15,.2,.19,.15]),Spacer(1,12),
       p('比较时保持几何、边界和材料一致，只改变冻结方案中指定的模型。H1-NT 直接复用网格验收的最终层，不重复计算同一点。'),
       p('如何解释模型差值','h3'),
       p('N0−LQ 是浓度非线性与双向耦合相对于一次弱耦合的总体影响；不能全部归为某一个漂移乘积项。NT−N0 是迁移率和扩散系数温度相关性的共同影响，不能称为“纯迁移率效应”。'),
       p('任何物理差值都应与网格变化量一起报告。差值小于约两倍网格变化量时，本轮只能判为“当前精度下不可分辨”，不通过提高电压来制造明显结果。'),
       p('独立守恒与适用性记录','h3')]
    rows=[['算例','热平衡相对差','电子电流相对差','E-P95 kV/cm']]
    for label,_ in CASES:
        m=main_data[label];rows.append([label,fmt(m.get('heat_balance_relative')),fmt(m.get('electron_balance_relative')),fmt(m.get('E_p95_V_m'),1e-5)])
    story += [table(rows),p('热平衡阈值为 1%，独立电子电流为 1%，4/5 阶热源积分差为 0.1%。10 kV/cm 的 E-P95 仅为工作性筛选，不是实验适用区证明。','small')]
    story += [PageBreak(),p('06 / 完整模型多场分布','title'),fig('H1_NT_multifield'),
       p('电势、电场、焦耳热、温升与位移采用 z=0 的原有限元插值；电子浓度图采用 z=+0.5 μm。图像没有跨孔插值，也没有人为平滑或镜像补齐。'),
       p('解释顺序应为：先比较同模型有孔/无孔的局部电场与热源，再看温度扩散和边界约束带来的空间分布，最后检查热-机电反馈形成的载流子重分布。局部热源峰值不必与温度峰值重合。')]
    story += [PageBreak(),p('07 / 载流子重分布与模型差异','title'),fig('H1_model_carrier_comparison'),
       p('上下表面分别从实际有限元场计算，所有模型使用相同色标。极性晶体不预设上下镜像关系；等温线性浓度近乎不变的情形也保留，不把零结果删去。'),
       p('空穴方程确实参与求解，但该候选 n 型材料的 p₀ 极小。可检查空穴方程、正性及相对变化，却不能把其绝对贡献夸大为主要传热或导电机制。')]
    hashes={}
    for name in ['TPSolverV2_analytic.exe','TPSolverV2_main_0142.exe']:
        path=RUNS/'evidence'/name
        if path.exists():hashes[name]=digest(path)
    story += [PageBreak(),p('08 / 证据与结论边界','title'),
       p('可信层次','h3'),
       p('解析验证通过，说明指定约化问题中的数值实现达到预设精度；网格与独立通量通过，说明这些主算例在当前离散与闭合模型下可重复计算。二者均不自动证明真实材料参数或全部能量物理已经实验验证。'),
       p('仍需明确的材料与物理假设','h3'),
       p('介电测量的机械约束、同一试样的 n(T)/迁移率与 Hall 因子、一次热释电参数、热参数混用、表面补偿和外部电场等问题仍然开放。屏蔽层里漂移与扩散显著抵消时，局部发热解释尤其需要谨慎。'),
       p('保存的复现信息','h3'),p('输入网格与参数、源码快照、可执行文件、逐步日志、节点数据、独立积分场样本及图源均有对应。初次失败记录不被覆盖。代表性复算一致性：'+yes(repeat_ok)+'。'),
       p('当前主计算使用磁盘辅助的直接线性求解，仅改变因子存储方式，不启用近似迭代线性求解。实现依据为 Intel oneMKL PARDISO 官方 Fortran 文档的 Storage of Matrices 小节。','small'),
       p('https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2026-0/onemkl-pardiso-parallel-direct-sparse-solver-iface.html','small')]
    for name,h in hashes.items():story += [p(name+'：','small'),p(h,'small')]
    out=ROOT/'output/pdf'/('F1_完整计算结果报告.pdf' if complete else 'F1_阶段计算结果报告.pdf')
    out.parent.mkdir(parents=True,exist_ok=True)
    doc=BaseDocTemplate(str(out),pagesize=(WIDTH,HEIGHT),title='F1 多场耦合计算结果',author='TPSolverV2 项目')
    frame=Frame(MARGIN,51,AVAILABLE,HEIGHT-109,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
    doc.addPageTemplates(PageTemplate(id='f1',frames=[frame],onPage=header));doc.build(story)
    reader=PdfReader(out);qa=ROOT/'tmp/pdfs/F1_results';qa.mkdir(parents=True,exist_ok=True)
    info={'complete':complete,'pdf':str(out),'pages':len(reader.pages),'pdf_sha256':digest(out),
      'text_characters_per_page':[len(pg.extract_text() or '') for pg in reader.pages],
      'visual_review':'PENDING','binary_sha256':hashes,'representative_repeat':repeat_records}
    (qa/'text_check.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(info,ensure_ascii=True,indent=2))

if __name__=='__main__':main()
