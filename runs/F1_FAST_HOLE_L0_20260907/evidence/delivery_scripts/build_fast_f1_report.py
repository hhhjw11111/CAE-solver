"""Chinese presentation report of the completed, explicitly reduced F1 matrix."""
from pathlib import Path
import argparse
import hashlib
import json
from reportlab.platypus import BaseDocTemplate, Frame, PageTemplate, PageBreak, Spacer, Image, Table, TableStyle
from reportlab.lib import colors
from pypdf import PdfReader
from build_case_freeze_pdf import p, WIDTH, HEIGHT, MARGIN, AVAILABLE, NAVY, TEAL, GREY, LIGHT, LINE


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def table(rows, ratios=None):
    ratios = ratios or [1/len(rows[0])]*len(rows[0])
    result = Table([[p(str(c), "headcell" if i == 0 else "cell") for c in row] for i, row in enumerate(rows)],
                   colWidths=[r*AVAILABLE for r in ratios], repeatRows=1, hAlign="LEFT")
    result.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6), ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LINEBELOW", (0, -1), (-1, -1), .5, LINE),
    ]))
    return result


def picture(path, maxheight=400):
    from PIL import Image as PILImage
    with PILImage.open(path) as im:
        w, h = im.size
    scale = min(AVAILABLE/w, maxheight/h)
    return Image(str(path), width=w*scale, height=h*scale)


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Hei", 8)
    canvas.setFillColor(GREY)
    canvas.drawString(MARGIN, HEIGHT-32, "TPSolverV2 / 2026-09-07 / 六算例快速汇报")
    canvas.setStrokeColor(LINE)
    canvas.line(MARGIN, 40, WIDTH-MARGIN, 40)
    canvas.setFont("Song", 8)
    canvas.drawString(MARGIN, 25, "L0 缩减网格矩阵；候选材料数值研究，非实验标定")
    canvas.setFont("Hei", 9)
    canvas.drawRightString(WIDTH-MARGIN, 25, str(doc.page))
    canvas.restoreState()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--figure-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root, figs, output = args.run_root.resolve(), args.figure_root.resolve(), args.output
    rows, audit = read(root/"main_results.json"), read(root/"fast_audit.json")
    if len(rows) != 6 or not audit["case_numerical_checks_passed"] or not audit["repeat"]["representative_repeat_executed"]:
        raise ValueError("All six cases must pass and the representative repeat must be executed")
    lookup = {r["label"]: r for r in rows}
    h0, h1 = lookup["H0-NT"], lookup["H1-NT"]
    def fmt(x, scale=1):
        return f"{x*scale:.5g}"
    changes = audit["historical_L0_L1_change"]
    source_figs = root.parents[2]/"output/figures/F1_20260906"
    # Resolve the project root independently of junction aliases.
    source_figs = Path(__file__).resolve().parents[1]/"output/figures/F1_20260906"
    story = [
        p("孔洞多场耦合研究", "title"), p("六算例结果与汇报要点", "subtitle"),
        p("交付状态：H0/H1 × LQ/N0/NT 六个 3 V 算例已完成；代表性 H1-NT 已复算。",
          "h3"),
        p("六算例收敛与守恒检查通过。复算存在浮点末位差异，原逐字节一致门槛未通过；"
          "本报告保留该未通过项，不将近似一致冒充严格复现通过。", "small"),
        p("根据 2026 年 9 月 7 日的加速要求，六算例统一采用 L0 级网格。每个几何内三个模型使用完全相同的网格。"
          "不改变材料、尺寸、孔径、边界和收敛容差，不新增参数扫描。原 L1 验证通过并保留，原 L2 异常退出记录保留。"
          "本交付不是原 L2 三层网格方案的最终验收。"),
        table([["算例", "总功率 mW", "最高温升 K", "环带温升 K", "用时 min"]]+[
            [r["label"], fmt(r["joule_W"], 1e3), fmt(r["theta_max_K"]),
             fmt(r["ring_theta_mean_K"]), fmt(r["elapsed_seconds"], 1/60)] for r in rows]),
        p("可以直接汇报的观察", "h3"),
        p(f"完整 NT 模型中，有孔/无孔总焦耳功率比为 {h1['joule_W']/h0['joule_W']:.5g}，"
          f"最高温升差为 {h1['theta_max_K']-h0['theta_max_K']:.5g} K；"
          f"同一 3–3.5 μm 环带的平均热源比为 {h1['ring_Q_mean_W_m3']/h0['ring_Q_mean_W_m3']:.5g}。"
          "局部热源集中和整体功率变化必须分开描述。"),
        p("NT−N0 代表温度相关迁移率与扩散系数的共同影响；N0−LQ 包含浓度非线性与双向反馈。"
          "这些差异不能归因于单一机制。环带局部指标另有粗网格采样偏差，见第 7 页；"
          "细小差异仍需结合网格变化量判断。"),
        PageBreak(), p("01 / 孔洞完整模型多场图", "title"),
        picture(figs/"H1_NT_multifield.png"),
        p("图中电势、电场、焦耳热、温升与 z 向位移取 z=0 中面；电子偏移取上表面。"
          "全部由原 HEX20 场直接评价，没有跨孔插值、人为平滑或上下镜像补齐。"),
        p(f"H1-NT：总功率 {fmt(h1['joule_W'],1e3)} mW，最高温升 {fmt(h1['theta_max_K'])} K，"
          f"环带平均电子绝对相对偏移 {fmt(h1['ring_abs_n_relative_mean'],100)}%。"
          "表面颜色峰值与体积积分点峰值的采样位置不同，不能混为同一统计量。"),
        PageBreak(), p("02 / 有孔与无孔同模型对照", "title"),
        picture(figs/"NT_geometry_comparison.png"),
        p("上下行分别为无孔和有孔，均采用 NT 模型；每列共用同一色标。"
          "孔壁附近的电场/热源变化与温度场的扩散响应是不同的空间现象。"),
        p("L0 孔洞网格 2,240 个 HEX20 单元、10,144 节点；无孔网格 1,960 个单元、9,345 节点。"
          "外尺寸和边界条件相同，孔洞有无使网格拓扑与单元数不同，不能称为两几何节点一一匹配。"),
        PageBreak(), p("03 / 六算例整体与局部指标", "title"),
        picture(figs/"six_case_comparison.png", 330),
        table([["模型", "热阻 K/W", "环带 E-P95 kV/cm", "环带 Q-P95 10¹² W/m³", "右端 ux nm"]]+[
            [r["label"], fmt(r["thermal_resistance_K_W"]), fmt(r["ring_E_p95_V_m"],1e-5),
             fmt(r["ring_Q_p95_W_m3"],1e-12), fmt(r["right_mean_ux_m"],1e9)] for r in rows],
             [.15,.19,.22,.25,.19]),
        p("热阻=最高温升/总焦耳功率；端面位移是面积积分均值；P95 是体积加权分位数。"
          "曲线连线只辅助对照，不表示额外参数扫描或统计置信区间。", "small"),
        PageBreak(), p("04 / 电子重分布与反馈", "title"),
        picture(figs/"H1_model_carrier_comparison.png"),
        p("上、下表面分别计算，所有模型采用统一色标。LQ 等温载流子变化接近舍入误差，"
          "不能把图中极小起伏解释为一阶载流子集中；N0/NT 才包含完整浓度非线性与双向热反馈。"),
        p("候选 n 型 ZnO 的平衡空穴浓度极小。空穴方程参与求解，但不把它的相对变化夸大为导电或传热的主因。"),
        PageBreak(), p("05 / 典型空间剖面折线", "title"),
        picture(figs.parent/"01_折线对比图/x_y4_z0_theta.png", 250),
        picture(figs.parent/"01_折线对比图/y_x0_ztop_E.png", 250),
        p("上图沿 y=4 μm 的中面直线取温升，下图沿 x=0 的上表面直线取电场模。"
          "虚线为无孔、实线为有孔；同一模型保持同色。有孔曲线在孔内断开，无孔曲线连续。"
          "HEX20 梯度在单元交界处取确定的一侧值，不做跨单元梯度平滑。", "small"),
        PageBreak(), p("06 / 独立验收与网格限制", "title"),
        table([["算例", "热流误差 %", "电子电流误差 %", "4/5 阶热源差 %"]]+[
            [r["label"], fmt(r["heat_balance_relative"],100), fmt(r["electron_balance_relative"],100),
             fmt(r["joule_quadrature_4_5_relative"],100)] for r in rows]),
        p("全部保留原阈值：缩放残差与更新 <10⁻⁹，离散守恒 <10⁻⁸，独立热流与电子电流 <1%，"
          "4/5 阶焦耳积分差 <0.1%。缩减网格不意味着放宽收敛或守恒。"),
        p("已有 H1-NT：原 L0→L1 变化", "h3"),
        table([["指标", "相对变化 %"]]+[[k, fmt(v,100)] for k,v in changes.items()],[.72,.28]),
        p("这只有两层证据，不能替代原三层最终网格验收，也不能将 H1-NT 的网格变化量当作所有模型差值的严格误差界。"
          "小于约两倍相应网格变化量的效应，仅作当前精度下不可分辨的趋势。"),
        p(f"额外发现：解析环带体积约 1.02102×10⁻¹⁷ m³，H0 积分点筛选得到 "
          f"{h0['ring_volume_m3']:.5g} m³，约低估 "
          f"{(1-h0['ring_volume_m3']/(3.141592653589793*(3.5e-6**2-3e-6**2)*1e-6))*100:.2f}%。"
          "H0 环带不贴合单元边界，粗网格积分点覆盖不足；因此环带平均值和 P95 仅作当前离散趋势，"
          "不能把约 1.45 倍的局部热源比当成已收敛定量结论。后续可先改善局部后处理积分，再决定是否重算。",
          "small"),
        p("本轮采用缩减网格、线程化直接求解和内存内存储，并对后续启动增加可用提交内存检查。"
          "孔洞 N0 采用更大初始延续步，保留失败自动回退与原容差。不同网格和并发条件的耗时不可用于宣称纯线程加速比。", "small"),
        PageBreak(), p("07 / 已完成的两类解析验证", "title"),
        picture(source_figs/"A1_joint_analytic.png", 245),
        p("A1：一维联合解，五层网格；最终 L² 约 4.94325×10⁻⁶，H¹ 半范数约 1.61919×10⁻⁴，"
          "观测阶约 2.99008/1.99146。零电流误差约 9.48×10⁻¹⁴。", "small"),
        picture(source_figs/"A2_hole_analytic.png", 245),
        p("A2：绝缘圆孔解析解，三层网格；最终电场 L² 约 0.541%，孔壁切向电场/焦耳热峰值倍率 "
          "1.999857/3.999428。该 2 倍/4 倍只属于约化解析问题。", "small"),
        PageBreak(), p("08 / 复现、边界与后续修正", "title"),
        p("一次 H1-NT 代表性复算已完成，原程序、输入、4 线程和运行配置保持一致，"
          "但三个核心结果文件的 SHA-256 不同，原字节一致门槛未通过。"
          "输入、网格、源码、可执行文件、运行日志、节点场、积分场样本、图源与结果表均保留。"),
        table([["复现文件", "字节一致"]]+[
            [k, "通过" if v else "未通过"] for k,v in audit["repeat"]["file_byte_identity"].items()],[.65,.35]),
        p(f"复算诊断：最大电势差 {audit['repeat']['maximum_absolute_field_difference'][7]:.4g} V，"
          f"最大温升差 {audit['repeat']['maximum_absolute_field_difference'][10]:.4g} K，"
          f"按原求解场尺度计算的最大节点状态差 "
          f"{audit['repeat']['maximum_solver_scaled_state_difference']:.4g}。"
          "这些差异仅作为数值对照记录，不替代原逐字节一致验收。"),
        p("研究条件", "h3"),
        p("20×20×1 μm 板，贯通孔半径 3 μm，晶体 +c 轴沿 +z。左右电极 0/−3 V、300 K、"
          "n=n₀、p=p₀；其他表面零电学/载流子法向通量、绝热、零牵引；只有六个最小刚体约束。"),
        p("LQ：参考输运与单向二阶焦耳热；N0：指数浓度和完整漂移扩散，但 μ/D 均冻结；"
          "NT：μ(T)=μ₀(T/T₀)⁻²，D=kBTμ/q。焦耳热采用漂移电导型闭合，不是全部漂移扩散能量项。"),
        p("今天可以交付什么", "h3"),
        p("已交付六算例的确定性数值结果和场分布，可用于检查量级、热点位置与模型差异，并据此讨论是否修改算例。"
          "尚未通过原 L2 三层网格最终验收和本次严格字节复现；不声称具体 ZnO 样品的实验预测。"),
        p("下一轮只修正由这些结果暴露的直接问题，保留当前结果，不通过改变材料或提高电压制造更明显的现象。"),
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    doc = BaseDocTemplate(str(output), pagesize=(WIDTH, HEIGHT), title="F1 六算例结果与汇报", author="TPSolverV2")
    doc.addPageTemplates(PageTemplate(id="fast", frames=[Frame(
        MARGIN, 52, AVAILABLE, HEIGHT-108, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)],
        onPage=footer))
    doc.build(story)
    pdf = PdfReader(output)
    info = {"pdf": str(output), "pages": len(pdf.pages),
            "characters": [len(p.extract_text() or "") for p in pdf.pages],
            "sha256": hashlib.sha256(output.read_bytes()).hexdigest(), "visual_review": "PENDING"}
    (root/"report_qa.json").write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(info, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
