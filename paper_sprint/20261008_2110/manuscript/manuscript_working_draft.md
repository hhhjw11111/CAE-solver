# 含圆孔ZnO压电半导体板的非线性热—电—力耦合响应

**WORKING_DRAFT｜2026年10月8日｜中文论文工作初稿 v0.2**

目标期刊：《固体力学学报》。作者、单位、基金与通信作者信息由作者提供；此处不代填。本文为可继续修改的研究初稿，不是已达到投稿条件的定稿。正文中的“待核验”说明所需证据，绝不代表零值或已通过。第5节保留的L0结果均为**历史粗网格探索性数据，未完成空间验收**。摘要与结论不采用这些数值作为经验证的定量发现。

## 摘要

电压驱动的含孔压电半导体板同时涉及载流子重分布、自发热及热压电变形，其局部响应的可信预测需要区分非线性方程求解收敛与空间离散精度。针对两端等温的ZnO板，本文建立三维二十节点六面体有限元模型，联合描述小变形力学、电位移、双载流子漂移—扩散及含漂移电导热源的稳态热传导。以无孔板和半径3 μm贯穿圆孔板为对象，在3 V电压下比较顺序线性参考模型、常输运系数非线性模型以及温度相关输运模型。采用独立端面电流与热流核对、漂移和扩散功率分解、常系数电热参考关系和材料内部厚度剖面约束结果解释。现有证据支持开展指定约化模型下的探索性比较，但局部载流子场及模型差值仍需补充分辨率证据。本文据此给出模型适用边界与可复核的分析方法，当前不提出圆孔增强幅度、独立迁移率贡献或真实器件性能的定量结论。

**关键词：** 压电半导体；圆孔板；热电力耦合；漂移—扩散；有限元

## Nonlinear thermo-electro-mechanical response of a ZnO piezoelectric semiconductor plate with a circular hole

### Abstract

Voltage-driven perforated piezoelectric semiconductor plates involve coupled carrier redistribution, self-heating and thermo-piezoelectric deformation. Reliable prediction of their local response requires nonlinear convergence and spatial resolution to be assessed separately. This study formulates a three-dimensional twenty-node hexahedral finite element model coupling small-strain mechanics, electric displacement, two-carrier drift–diffusion and steady heat conduction with a drift-conductivity heat source. A plain ZnO plate and a plate containing a through-hole of radius 3 μm are considered under a voltage of 3 V with isothermal end contacts. A sequential linear reference, a nonlinear model with frozen transport coefficients and a model with temperature-dependent mobility and diffusivity are distinguished. Independent current and heat-flux audits, drift and diffusion power decomposition, a constant-coefficient electrothermal reference and through-thickness profiles are used to constrain interpretation. Available evidence permits exploratory comparisons within the specified reduced model; spatial resolution of local carrier fields and inter-model differences remains to be established. The present working draft therefore provides the governing model and a reproducible interpretation framework without claiming a validated enhancement magnitude, an isolated mobility contribution or experimentally calibrated device performance.

**Keywords:** piezoelectric semiconductor; perforated plate; thermo-electro-mechanical coupling; drift–diffusion; finite element method

## 1 引言

压电半导体中的极化电场能够驱动自由载流子迁移，载流子形成的空间电荷又会改变电场及机械响应。这一耦合使几何缺陷的作用不能仅用弹性应力集中描述。对于施加电压的含孔薄板，圆孔同时改变电流通道与导热路径；电流产生的温升进一步通过热膨胀、热释电和输运系数变化反馈到载流子与变形场。因此，需要在一致的电、热和机械边界下比较几何效应，并判断观测到的变化是否超出普通电热过程能够解释的范围。

已有研究为压电半导体板的有限元分析提供了直接基础。Zhao等[1]采用一阶Mindlin展开构造压电半导体板模型，通过Galerkin弱式建立有限元方程，并与COMSOL自定义偏微分方程模块进行比较。该工作说明了载流子与板弯曲响应的耦合及其数值处理途径，但其板理论和验证条件与本文的三维含孔、自发热问题不同，不能作为当前程序已经通过外部软件验证的依据。

热载荷下的多场耦合则具有不同的驱动方式。Kong等[2]从三维热压电半导体方程出发建立包含热释电的夹层梁模型，研究给定热载荷下温度、变形、极化和载流子之间的耦合，并给出特定电学边界下的伸长、弯曲和屈曲分析。本文关注外加电压产生热源后的稳态含孔板响应，采用的是接触电势、等温端面及绝热孔壁条件。二者可用于比较耦合结构和模型假设，但不能直接互换解析解、参数标定或边界条件。

这一问题还存在数值层面的区分。较小的Newton残差表示给定离散方程得到充分求解，并不自动保证薄载流子层、沿电场方向的输运和孔缘局部场已被网格解析。漂移与扩散可能局部抵消，使总电流较小而分项较大；仅依据总电流或端面热平衡也不能判断局部载流子曲线正确。若温度相关模型与常系数模型之间的差异较小，其解释还需要模型差值本身的空间敏感性证据。

本文采用三维HEX20离散，在相同材料、3 V电压和温度边界下比较无孔与圆孔板。研究分为三项相互关联的工作：明确顺序线性参考、常输运非线性及温变输运模型的差别；以独立功率、通量和局部场检查评估既有结果的适用范围；借助常系数电热关系判断圆孔与温度反馈的解释边界。研究范围限于给定约化热源及材料输入下的稳态数值响应，不涉及新增实体实验、高场输运、挠曲电或新的缺陷形状。

## 2 控制方程、材料与边界条件

### 2.1 场变量及热压电本构

设材料域为Ω，位移、电势与相对温升分别为u、φ和θ，总温度为T=T₀+θ。采用SI单位、小应变和准静态假设，定义

$$
\boldsymbol\varepsilon=\operatorname{sym}\nabla\mathbf u,\qquad
\mathbf E=-\nabla\phi. \tag{1}
$$

为与程序中的工程剪应变一致，Voigt应变向量采用[εxx, εyy, εzz, γyz, γzx, γxy]ᵀ，其中γij=2εij。晶体正c轴与全局正z轴重合。本构关系写为

$$
\boldsymbol\sigma=\mathbf C\boldsymbol\varepsilon-\mathbf e^T\mathbf E-\mathbf b\theta,
\qquad
\mathbf D=\mathbf e\boldsymbol\varepsilon+\boldsymbol\kappa\mathbf E+\boldsymbol\pi\theta,
\qquad \mathbf b=\mathbf C\boldsymbol\alpha. \tag{2}
$$

其中C、e、κ、α分别表示弹性、压电、介电及热膨胀参数；π为热释电系数，以避免与空穴浓度p混淆。b的单位为Pa/K。静力平衡和Gauss方程为

$$
\nabla\!\cdot\boldsymbol\sigma=\mathbf0,\qquad
\nabla\!\cdot\mathbf D=\rho,\qquad
\rho=q(p-n+N_D^+-N_A^-). \tag{3}
$$

本模型描述参考态附近的增量热压电响应。常量自发极化及其表面补偿没有独立标定；孔内空气电场也未作为独立区域求解。因此，式(2)—(3)连同下述自然边界构成本文采用的模型契约，不被解释为已经重现真实器件的全部界面物理。

### 2.2 漂移—扩散与约化热源

电子和空穴常规电流密度分别定义为

$$
\mathbf J_n=q n\boldsymbol\mu_n\mathbf E+q\mathbf D_n\nabla n,\qquad
\mathbf J_p=q p\boldsymbol\mu_p\mathbf E-q\mathbf D_p\nabla p. \tag{4}
$$

本轮模型中产生复合项为零，连续性方程为∇·Jn=0和∇·Jp=0。N0与NT模型采用对数浓度变量

$$
n=n_0\exp(\eta_n),\qquad p=p_0\exp(\eta_p), \tag{5}
$$

从而∇n=n∇ηn、∇p=p∇ηp。在有限变量值下这一表示保证浓度为正，但不提供空间精度保证。NT采用

$$
\boldsymbol\mu_s(T)=\boldsymbol\mu_{s0}(T/T_0)^{-\beta_s},\qquad
\mathbf D_s(T)=\frac{k_BT}{q}\boldsymbol\mu_s(T),\qquad s=n,p. \tag{6}
$$

N0同时将μs和Ds冻结在T₀；它并非仅冻结迁移率而保留扩散系数中的T。该模型未加入温度梯度驱动的独立热扩散项、高场迁移率、陷阱、非完全电离或完整载流子能量方程。

热流和热平衡方程为

$$
\mathbf q_h=-\mathbf k\nabla T,\qquad
\nabla\!\cdot\mathbf q_h=Q_d,
\qquad Q_d=q\left(n\mathbf E\!\cdot\boldsymbol\mu_n\mathbf E+
p\mathbf E\!\cdot\boldsymbol\mu_p\mathbf E\right). \tag{7}
$$

式(7)采用漂移电导型热源；在最终生产状态中热源延续系数等于1。它不等同于(Jn+Jp)·E，也不包含完整漂移—扩散能量传递。后续功率审计用于定量说明这一约化的后果，并不静默将热源更换为另一物理定义。

### 2.3 三种模型的比较含义

LQ是一阶相对浓度及顺序反馈参考。其浓度写为n=n₀(1+ηn)、p=p₀(1+ηp)，电流采用参考浓度乘以漂移—扩散括号项，去除ηE和η∇η的高阶乘积；保留的二阶漂移热也以n₀、p₀计算。该模型依次执行三个阶段：固定θ=0并关闭热源，求解电—力—载流子场；固定前六个字段，恢复热源求解θ；固定φ、ηn、ηp和θ，仅更新位移。第三阶段之后不返回求解电场和载流子，因此最终组合状态不能默认满足NT模型式样的七场同时平衡。

N0和NT均保留式(5)及完整密度非线性，通过同时耦合求解形成稳态。二者之差反映本模型内迁移率和扩散系数共同随温度变化后的反馈；它不构成“纯迁移率贡献”。同样，N0与LQ的差异同时包括密度非线性和求解反馈结构变化，不能作为载流子非线性的单项因果分解。本文比较时保留这两项限制。

### 2.4 几何与边界

板的外尺寸为20 μm×20 μm×1 μm，即x、y∈[−10,10] μm，z∈[−0.5,0.5] μm。H0为无孔板，H1为以z轴为中心、半径a=3 μm的贯穿圆孔板。规定V*=3 V，左端x=−10 μm取φ=0，右端x=10 μm取φ=−V*。左右端面均取ηn=ηp=0和θ=0，以固定参考浓度并保持T₀。其余外表面和孔壁采用

$$
\boldsymbol\sigma\boldsymbol\nu=\mathbf0,\quad
\mathbf D\!\cdot\boldsymbol\nu=0,\quad
\mathbf J_n\!\cdot\boldsymbol\nu=\mathbf J_p\!\cdot\boldsymbol\nu=0,\quad
\mathbf q_h\!\cdot\boldsymbol\nu=0, \tag{8}
$$

其中ν为外法向。零法向电位移只规定材料侧的边界条件，不表示孔内电场已被求解；边界中也未包含孔壁陷阱、接触热阻或表面复合。

全域模型用三个角点消除刚体模态。A=(−10,−10,−0.5) μm处固定u₁、u₂、u₃；B=(−10,10,−0.5) μm处固定u₁、u₃；C=(−10,−10,0.5) μm处固定u₁。这里没有新增外部机械载荷。不同算例的位移比较保持相同约束和参考系，并同时使用相对位移、厚度变化与应力，避免将刚体参考差异解释为圆孔效应。已有半域结果若作为补充证据，必须单列对称面条件和实际半域积分，不能与全域数值直接混用。

![图1 几何、晶向、边界与材料内部厚度取样线；示意而非计算场](C:/Users/Administrator/Desktop/codex-workspaces/Console3/solver_v2/paper_sprint/20261008_2110/figures/fig01_geometry_boundary.png)

**图1** 几何、晶向、端面电热边界、A/B/C机械约束与材料内部厚度取样线。Fig. 1 Geometry, crystal orientation, boundary conditions, rigid-body constraints and through-thickness sampling lines. 图为几何示意，不含计算结果；(0,3.25 μm)取样线位于孔外材料内。

### 2.5 材料输入与适用性

表1给出当前程序实际输入，不将其等同于同一样品的完整实验标定。参数来自既有材料审计和源码契约，弹性、压电、介电、热学及输运参数存在来源与约束不同的问题；例如将绝热测量弹性常数用于当前稳态本构、将Hall参数用于漂移描述，以及μ(T)幂律指数的选取，都应视为有效输入假设。空穴参考浓度极小，保留空穴自由度并不表示当前ZnO算例已呈现可观测的双极输运贡献。

**表1 当前ZnO候选材料输入（非同样品实验定标）**  
Table 1 Candidate ZnO inputs used in the model; not a same-specimen calibration.

| 参数 | 输入值 | 单位或限制 |
|---|---:|---|
| C₁₁=C₂₂；C₁₂；C₁₃=C₂₃；C₃₃ | 209.7；121.1；105.1；210.9 | GPa |
| C₄₄=C₅₅；C₆₆ | 42.47；44.3 | GPa；工程剪应变约定 |
| e₃₁=e₃₂；e₃₃；e₁₅=e₂₄ | −0.62；0.96；−0.37 | C/m² |
| κ₁₁/ε₀=κ₂₂/ε₀；κ₃₃/ε₀ | 8.54；10.2 | 定应变形式的候选解释 |
| α₁=α₂；α₃ | 4.31×10⁻⁶；2.49×10⁻⁶ | K⁻¹ |
| π₃ | −6.9×10⁻⁶ | C/(m²·K)；π₁=π₂=0 |
| k | 60 | W/(m·K)，各向同性有效假设 |
| T₀；n₀ | 300；6.0×10²² | K；m⁻³ |
| nᵢ；p₀=nᵢ²/n₀ | 2.9×10⁻⁴；1.4016667×10⁻³⁰ | m⁻³；p₀为由输入计算值 |
| μn₀；μp₀ | 0.0205；0.005 | m²/(V·s)，各向同性 |
| βn；βp | 2；2 | 局部有效模型假设 |
| N_D⁺；N_A⁻ | n₀−p₀；0 | m⁻³，参考电中性 |

## 3 有限元离散与非线性求解

### 3.1 弱式与空间离散

对满足本质边界条件的试函数v、wφ、wn、wp和wT，标准Galerkin弱式可写为

$$
\int_\Omega\boldsymbol\varepsilon(\mathbf v)^T\boldsymbol\sigma\,d\Omega
-\int_{\Gamma_t}\mathbf v\!\cdot\bar{\mathbf t}\,d\Gamma=0, \tag{9}
$$
$$
\int_\Omega\nabla w_\phi\!\cdot\mathbf D\,d\Omega
+\int_\Omega w_\phi\rho\,d\Omega
-\int_{\Gamma_D}w_\phi\bar D_\nu\,d\Gamma=0, \tag{10}
$$
$$
\int_\Omega\nabla w_s\!\cdot\mathbf J_s\,d\Omega
-\int_{\Gamma_{J_s}}w_s\bar J_{s\nu}\,d\Gamma=0,\quad s=n,p, \tag{11}
$$
$$
\int_\Omega\nabla w_T\!\cdot\mathbf k\nabla T\,d\Omega
-\int_\Omega w_T Q_d\,d\Omega
+\int_{\Gamma_q}w_T\bar q_\nu\,d\Gamma=0. \tag{12}
$$

式(10)的ρ号由∇·D=ρ及E=−∇φ共同确定；式(12)使用外向热流约定。求解器可对整行乘以−1而不改变零残差问题，但不能只改变源项或边界项符号。

每个HEX20单元具有20个节点，每节点变量顺序为(u₁,u₂,u₃,φ,ηn,ηp,θ)，单元未知量共140个。各场采用同一等参形函数插值，几何映射和场梯度使用一致的Jacobian。组装后的残差向量按五类方程组织，即力学三行、电势一行、电子一行、空穴一行和温度一行。默认采用三阶张量积Gauss积分；独立通量与体热源使用各自记录的积分规则，其提高阶次检查与网格加密分开报告。

当前程序还包含默认关闭的log-weighted输运离散。该形式在载流子方程中使用w=exp(−ηs)N及对应一致线性化；物理电流与热源保持式(4)和式(7)，但测试函数空间及离散代数平衡含义发生变化。若复用显式启用该选项的结果，其离散形式在表2中单列，物理电流仍以独立边界通量验收，不能把加权自由行平衡直接称为物理电流守恒。本文不新增稳定化物理或将不同离散形式未经说明合并比较。

### 3.2 一致线性化与迭代控制

设所有自由度构成向量x，非线性残差为R(x)。Newton步满足K(xᵏ)Δx=−R(xᵏ)，其中K=∂R/∂x包含场间交叉块。对式(5)—(6)，关键导数包括

$$
\frac{\partial n}{\partial\eta_n}=n,\quad
\frac{\partial p}{\partial\eta_p}=p,\quad
\frac{\partial\boldsymbol\mu_s}{\partial T}=-\frac{\beta_s}{T}\boldsymbol\mu_s,
\quad
\frac{\partial\mathbf D_s}{\partial T}=\frac{1-\beta_s}{T}\mathbf D_s. \tag{13}
$$

例如电子电流关于温度的导数为qn[(∂μn/∂T)E+(∂Dn/∂T)∇ηn]。热源的浓度、电场和温度导数也纳入同一Jacobian，从而使温度反馈与残差所用本构一致。N0中所有上述输运温度导数为零；LQ使用其一阶模型自身的导数，而非沿用非线性浓度导数。

程序采用字段尺度化、稀疏直接求解、对数密度增量控制和基于真实残差的回溯线搜索。电压与热源延续仅改变到达指定最终状态的求解路径；中间电压或部分热源系数的解不被登记为完整工况。历史标准F1状态采用最大阶段平衡量小于10⁻⁸、残差和增量小于10⁻⁹等原门槛，并同时检查独立通量。本文不放宽既有门槛来追求表面上的完成状态。

### 3.3 比较量与空间敏感性

主要整体量为采样峰温升、端面电流、漂移热总功率、外边界热流、相对位移及厚度变化。局部量包括孔外材料的电场、浓度比、热源和应力。积分、均值与分位数必须使用真实体积权重；应力点最大值和节点峰值单独保留，不能与加权95%分位数混称。若相同环带在两种几何中包含的材料不同，需明确公共有效材料区域及积分方法。

对同一观测量M，模型差异定义为δM=(MNT−MN0)/MN0，但只在分母远离零且尺度有明确意义时使用。相邻网格差衡量敏感性，不构成严格误差上界。对于非均匀、不同方向同时变化的网格，不无依据拟合收敛阶或套用Richardson外推。局部载流子扰动、NT−N0差值及几何差值均需要各自的分辨率检查，整体功率稳定不替代这些要求。

## 4 验证方法与电热参考关系

### 4.1 求解、通量与空间证据分层

本文分别检查来源、几何、代数求解、物理通量、空间分辨率与模型适用范围。只有输入、网格、源码或可执行文件及结果场可以对应，数值才进入比较；全域与半域、Galerkin与加权离散、完整热源与延续中间态分别登记。独立热流检查以|∮qh·νdΓ−∫Qd dΩ|及其明确参考尺度度量，电子电流检查以独立端面和其余边界的净流出衡量。既有1%通量门槛和0.1%四/五阶热源积分差门槛是本任务的数值验收标准，不是期刊规定或全场误差界。

解析验证优先复用与当前代码兼容、几何及边界有效的既有记录；只有核到实际误差及运行日志后才写为已完成。对无法复核原场、结果存在几何失效或仅有旧版输出的项目，保留阻塞记录。本轮的小型输运基准及既有解析验证实际结果见配套证据表；没有闭合的验证不由示意图或“PASS”字样替代。

本轮小型复核实际通过V01、V02、G02、F1M和G03入口，原始输出与退出码保留。V02解析验证采用1、2、4、8、16个单元的序列，电势相对L2误差依次为6.317758×10⁻³、8.008457×10⁻⁴、1.004809×10⁻⁴、1.257206×10⁻⁵和1.571883×10⁻⁶，最细层H1误差为9.012461×10⁻⁵；相应观测阶约为3.000和2.000。这证明该解析验证条件下的离散误差按预期下降，不能替代含孔ZnO生产算例的空间验收。G02、F1M与G03的材料和单元一致切线检查通过，也不证明候选材料已经实验定标。

**本轮V02验证数值表（与图2历史A1/A2不同，不混用算例身份）**

| V02单元数 | 电势相对L2误差 | 电势相对H1误差 |
|---:|---:|---:|
| 1 | 6.317758×10⁻³ | 2.265987×10⁻² |
| 2 | 8.008457×10⁻⁴ | 5.741188×10⁻³ |
| 4 | 1.004809×10⁻⁴ | 1.440378×10⁻³ |
| 8 | 1.257206×10⁻⁵ | 3.604174×10⁻⁴ |
| 16 | 1.571883×10⁻⁶ | 9.012461×10⁻⁵ |

来源为本次validation/runtime_retry/V02.txt，表中有效位数仅为排版舍入，完整值保留在源文件。

![图2 历史A1耦合合成条带和A2绝缘圆孔合成解析误差；非3 V生产场验证](C:/Users/Administrator/Desktop/codex-workspaces/Console3/solver_v2/paper_sprint/20261008_2110/figures/fig02_analytic_errors.png)

**图2** 历史A1耦合合成条带与A2绝缘圆孔合成算例的实际解析误差。Fig. 2 Recorded analytical errors of the historical coupled-bar and insulated-hole synthetic tests. A1/A2来自F1_20260906保存日志，本轮重新提取曲线并核对源文件；它们不是上表V02重放，也不是3 V含孔ZnO生产状态。A1最细层φ相对L2误差为4.943250×10⁻⁶；A2的L2层φ及E相对L2误差分别为2.782470×10⁻⁴、5.408310×10⁻³。这些误差限定在各自合成条件内，不认证主工况全场。图源为data/A1_error_source.csv与A2_error_source.csv，原日志身份及哈希见audit/postprocess/analytic_sources.json。

### 4.2 一维漂移—扩散核对式

为隔离输运离散问题，考虑常E、μ和D的一维稳态电子方程J=qμEn+qDdn/dx=常数。取x∈[0,L]、n(0)=nL、n(L)=nR，记b=μE/D。当b≠0时，解析浓度及电流为

$$
n(x)=n_L+(n_R-n_L)\frac{1-e^{-bx}}{1-e^{-bL}},\qquad
J=q\mu E\frac{n_R-n_Le^{-bL}}{1-e^{-bL}}. \tag{14}
$$

在E→0时，浓度退化为线性分布，J=qD(nR−nL)/L。较大正bL时，应核对边界层位置和符号；不将电子粒子通量与常规电流的方向混用。解析式用于程序路径的测试规格，本身并不证明三维程序已经通过该测试。

本轮复用并实际重放了现有R22生产HEX20单元核输运条带，结果与其既有标准输出逐字节一致。在域Pe=1000、A=bL=2000、128个单元时，电子采样浓度和电流最大归一偏差分别为7.595833×10⁻⁵和2.518331×10⁻³，空穴对应偏差为4.047920×10⁻⁵和1.662881×10⁻³。这些误差来自每单元11点采样；浓度用2归一，电流以max(1,|精确电流|)归一，不是独立相对L2误差。该条带使用生产材料点与单元残差/Jacobian，但采用专用缩约装配与稠密求解，未覆盖完整全局CSR路径，也未完成同种电子的反向场测试。

因此，小型基准仍需补充物理浓度n的独立相对L2积分、积分细化稳定性、完整内部剖面、单独端面通量及自由载流子行残差。粗网格在高Péclet条件下出现较大误差是可报告的离散局限，不要求所有粗网格通过。三维Pe诊断还必须写明采用的长度与热电压定义，不能将单元边电势跨度自动等同于唯一的流向Pe。式(14)—(16)的额外代数核对只验证解析关系，不被冒充三维程序测试。

### 4.3 常系数电热参考关系

设独立参考问题满足J₀=σ₀E、σ₀和k均为常数、∇·J₀=0，且左右端面保持φ=0、−V*及θ=0，其余边界电绝缘、绝热。定义ψ=−φ∈[0,V*]，则Δψ=0，而热方程为−kΔθ=σ₀|∇ψ|²。由链式法则，函数

$$
\theta_{\rm ref}(\psi)=\frac{\sigma_0}{2k}\psi(V_*-\psi) \tag{15}
$$

满足热方程与同样的边界条件。对连接两端电极的连续材料域，ψ沿连接路径必经V*/2，因此该参考问题的最高温升为σ₀V*²/(8k)。这一结论依赖常系数、两端等温及匹配的绝缘绝热侧壁条件；不直接推广到张量系数不匹配、载流子非均匀或非等温接触问题。

式(15)表明，在上述参考模型内，圆孔可改变总电流、总功率和温度空间分布，而不必改变峰温。因此，孔径变化造成峰温变化较弱，本身不足以证明特殊耦合机制。本文将实际θ与按σ₀=q(n₀μn₀+p₀μp₀)构造的参考量比较；偏离反映全部模型差别及离散误差的合成，不能直接命名为纯热释电、纯压电或纯迁移率贡献。LQ仍包含电—力—载流子反馈和分阶段冻结，亦不等同于这个严格欧姆参考问题。

### 4.4 功率恒等式及局部抵消

以Jn和Jp的常规电流约定，定义

$$
P_{\rm term}=-\int_{\partial\Omega}\phi(\mathbf J_n+\mathbf J_p)\!\cdot\boldsymbol\nu\,d\Gamma,
\quad P_d=\int_\Omega Q_d\,d\Omega,
\quad P_{\rm diff}=\int_\Omega q(\mathbf D_n\nabla n-\mathbf D_p\nabla p)\!\cdot\mathbf E\,d\Omega. \tag{16}
$$

对稳态守恒电流，由分部积分得Pterm=∫(Jn+Jp)·EdΩ=Pd+Pdiff。当前热方程则要求Qout=∮qh·νdΓ=Pd。两个关系检验不同对象：端面电功率与漂移热相差时，可能来自扩散功及离散积分误差，不能立即把它当成热守恒失败；Qout与Pd的差也不能通过更换端面功率分母掩盖。

局部诊断同时保留Jn,d、Jn,diff以及空穴对应分项的方向和大小。在两项均非零的区域，可用|Jd+Jdiff|/(|Jd|+|Jdiff|)辅助描述抵消，分母接近零的点单独处理。若扩散贡献不可忽略，本文仅报告给定漂移热模型下的温度结果，暂不将其外推为完整器件能量转化或热耗散预测。

## 5 结果与讨论

### 5.1 主工况与可用证据

主比较限定为3 V的H0/H1及现有LQ/N0/NT模型。已定位的历史F1_FAST_HOLE_L0_20260907记录包含六个L0状态，其原有标签为PASS_NOT_MESH_ACCEPTANCE，表示相应历史求解与通量门槛通过，并非网格验收。当前工作不更改这一标签，也不将最新程序版本名称追溯贴到旧结果上。本轮已核对现存input、mesh、exe及54个归档源文件身份；但早期H0三例缺少当时同时记录的逐case运行provenance，环境信息属于后补回溯，来源状态仍为PENDING。H1有当时运行记录。现存档案哈希匹配不能补造历史进程来源，也不是重新编译后的二进制可复现证明。

**表2 主比较的证据层级（随本轮审计更新）**  
Table 2 Evidence levels of the principal comparisons.

| 比较 | 已有证据 | 尚缺证据与允许表述 |
|---|---|---|
| H0/H1、LQ/N0/NT，3 V，历史L0 | 六例输入、网格、exe与源快照身份已核；保留历史标签 | H0缺当时逐case运行来源记录；六例空间验收未完成，仅探索性比较 |
| H1、3 V、后续定向网格 | 有独立网格/积分审计记录 | 不同离散形式、半域及约束须分开；不据此认证旧六例全部场 |
| H1 NT−N0 | 同条件指标可建立对照 | 同系列网格上的模型差值稳定性待核验 |
| H0/H1机械响应 | 有位移及积分点应力输出入口 | 相对位移、厚度变化和应力分位数的空间证据待核验 |

### 5.2 电热响应与热源解释

在已核到的H1—NT历史L0状态中，漂移热总功率为1.594472711 mW，端面电功率为1.513271195 mW，外边界热流为1.584440099 mW，采样峰温升为3.906877953 K，原有独立热流相对缺口为0.6292119%。**这些数值均为历史粗网格探索性结果，未完成空间验收，不能作为已验证的温升或功率精度结论。** 原始指标来自该状态的metrics记录；本轮重构的功率分项已与这些记录交叉核对。

本轮从保存的HEX20场重构五阶体积分，并与原metrics中已有的独立边界积分对照；本轮没有重新计算端面通量。得到H1—NT的Pdiff=−0.08118642 mW，体积积分的总电流功为1.513286294 mW，与原记录端面电功率的相对缺口约9.97794×10⁻⁶。H0—NT的Pd=1.850211132 mW，Pdiff=−0.107199672 mW。上述分项说明，在这些**未通过空间验收的探索性L0场**中，扩散功的量级已不能与漂移功混同；这并不是对完整载流子能量模型的验证。

**表3 同3 V主工况指标：全部为历史L0探索性结果，未完成空间验收**  
Table 3 Principal metrics at 3 V; all entries are exploratory L0 results without spatial acceptance.

| 几何/模型 | 采样峰温升/K | Pd/mW | Pdiff/mW | Pterm/mW | Qout/mW |
|---|---:|---:|---:|---:|---:|
| H0/LQ | 3.69319 | 1.77361 | 约0 | 1.77361 | 1.77361 |
| H0/N0 | 4.02515 | 1.88627 | −0.112659 | 1.77360 | 1.87222 |
| H0/NT | 3.96088 | 1.85021 | −0.107200 | 1.74300 | 1.83695 |
| H1/LQ | 3.69463 | 1.53951 | 约0 | 1.53951 | 1.53951 |
| H1/N0 | 3.96721 | 1.62465 | −0.0851411 | 1.53950 | 1.61405 |
| H1/NT | 3.90688 | 1.59447 | −0.0811864 | 1.51327 | 1.58444 |

表3的功率均为全域量；Pdiff保留符号，“约0”对应原始积分在数值舍入量级的值而非删除分项。完整精度、原场位置与积分定义见data/metrics_summary.csv及后处理契约。LQ数据仅表示三阶段顺序模型结果。H0三例另有第5.1节所述同时运行来源记录缺口；本表与下列H0/H1图仍限于现存归档场的探索性比较。表3不为任何一列提供网格误差界。

![图4 同3 V条件下H0/H1 NT电场、电流及温升，中面单元内侧有限元值；历史L0探索性，未完成空间验收](C:/Users/Administrator/Desktop/codex-workspaces/Console3/solver_v2/paper_sprint/20261008_2110/figures/fig04_midplane_fields.png)

**图4** 同3 V条件下H0/H1 NT的电场幅值、电流幅值和温升，中面z=0⁺的原始单元内有限元采样；同列采用共同色标。Fig. 4 Midplane electric field, current and temperature rise of H0 and H1 at 3 V, with common column scales. **历史L0探索性，未完成空间验收。** 图中孔周电场与电流绕流结构以及温度变化可从原始场重构；峰值精度仍受网格限制。

### 5.3 温度相关输运的比较边界

对于表1中的βn=βp=2，NT中迁移率随T⁻²变化，扩散系数随T⁻¹变化。因而在给定温升下，漂移和扩散系数均减小，但最终电流或温升变化还受到电场和载流子重分布的反馈，不能仅凭幂律导数判断完整耦合解的变化幅度。N0提供冻结两个输运系数的对照，其余几何、边界、材料和数值设置应保持一致。

NT−N0的峰温、端电流、局部浓度及应力差值需要在同网格序列上检查。旧文档中的约1.5%等差值尚不能替代这一检验，故本稿不将这些百分比写入摘要和结论。若模型差异与相邻网格敏感性处于同一量级，应报告“当前分辨率下未充分分辨”，而非通过增加电压或改变指数人为放大信号。

表3给出H1—N0与H1—NT的实际L0温升3.96721 K和3.90688 K，可作为后续同系列网格差值检查的起点。本文暂不将其转换为精确反馈百分比，也不单独制作暗示差值已分辨的“机制贡献图”。图5保留为待补的差值敏感性图位，所需证据为相同离散形式、边界和网格序列下的N0/NT对照。

### 5.4 局部载流子层与空间分辨率

现有审计提示，整体指标通过与局部载流子指标通过并不总是一致。在**3 V、R23 log-weighted、y半域**的L1/L2对照中，局部电子和空穴扰动差分别约3.00654%和3.00977%，超过该局部检查预设的2%门槛。该失败应保留，其身份不同于表3的F1全域普通Galerkin六例；它不自动否定所有整体量，也不能反过来作为旧六例已经完成网格验证的证据。

厚度剖面需在孔外材料内固定平面位置，再比较从下表面到上表面的实际场值。中面值接近参考态并不能排除表面载流子响应，插值曲线平滑也不能排除未解析的薄层。当前参考屏蔽尺度与表面首层厚度均需结合具体网格登记，重点核查表层、内部粗区和沿场方向是否同时足以分辨相关观测量。仅改变某个方向的网格时，结论限定为该方向敏感性。

![图3 节点电势跨度Pe诊断；并非单元内部极值或唯一流向Pe](C:/Users/Administrator/Desktop/codex-workspaces/Console3/solver_v2/paper_sprint/20261008_2110/figures/fig03_peclet_regions.png)

**图3** 已有3 V状态的节点电势跨度Pe区域诊断。Fig. 3 Region-wise potential-span Péclet diagnostics at 3 V. 此Pe使用单元节点电势范围、NT单元节点平均温度，N0/LQ用T₀；区域依次按接触、上下表层、孔周及其余体内划分，体积比例采用detJ权重。它不是单元内部电势极值搜索，不能单独认证局部场精度。

![图3b 3 V R23 log-weighted半域局部载流子场比较失败证据](C:/Users/Administrator/Desktop/codex-workspaces/Console3/solver_v2/paper_sprint/20261008_2110/figures/fig03b_spatial_failure.png)

**图3b** 3 V R23 log-weighted半域L1→L2的局部场映射比较，电子和空穴扰动超过2%冻结门槛。Fig. 3b Local-field comparison retaining the failed carrier criteria. 数据源为field_shape_r23_L1_L2_V030.json的实际比较记录。

![图6 孔外材料内部厚度剖面；L0探索性](C:/Users/Administrator/Desktop/codex-workspaces/Console3/solver_v2/paper_sprint/20261008_2110/figures/fig06_thickness_profiles.png)

**图6** (x,y)=(0,3.25 μm)处电子扰动、温升、Ez与电热参考偏离的厚度剖面，整条取样线位于孔外材料内。Fig. 6 Through-thickness profiles along a material line outside the hole. 曲线由原始HEX20插值得到，保留单元分段和场变化；**历史L0探索性，未完成空间验收**。此处上下表面存在载流子响应，而中面接近参考态，说明仅中面图不足以覆盖厚度响应；精确峰值仍待局部分辨率核验。

### 5.5 机械响应与几何效应

温升通过热应力项直接作用于力学方程，并经热释电、电场及逆压电效应间接影响位移和应力。圆孔还改变应力传递路径。然而，在未完成关闭对照的情况下，这些项之间的作用不能由一幅位移云图分离。本文首先报告总耦合响应，只有主结论确实依赖某一机制时才考虑单一且物理定义明确的关闭对照。

本轮位移评价在(x,y)=(0,3.25 μm)厚度线上减去中面位移，以保留真实厚度变化。H0—NT与H1—NT的上下表面厚度差分别为1.0701113×10⁻¹¹ m和1.0604187×10⁻¹¹ m。圆周应力取r=3.25 μm、z=0的单元内侧迹，不做节点平滑；H1—NT的该圆周采样最大|σθθ|为5.338699×10⁵ Pa。**这些均为探索性L0采样量，未完成空间验收，不是收敛的孔壁最大应力。** 未计算可靠公共体积环带增强比，故不报告此类比值。

![图7 固定材料圆周应力与厚度线相对位移，未作节点平滑，L0探索性](C:/Users/Administrator/Desktop/codex-workspaces/Console3/solver_v2/paper_sprint/20261008_2110/figures/fig07_mechanics.png)

**图7** r=3.25 μm、z=0的周向应力与厚度线相对中点的位移幅值；LQ为顺序三阶段参考。Fig. 7 Circumferential stress on a material circle and displacement relative to the thickness-line midpoint. **历史L0探索性；机械数值的空间敏感性待补。**

### 5.6 普通电热参照与物理外推限制

对照式(15)时，首先比较θ与电势变量ψ之间的关系，再比较空间分布。若无孔与有孔的峰温接近，而总功率存在差别，这可与普通电热参考相容，并不自动意味着圆孔对热响应无作用。若实际解明显偏离参考，可能来源包括载流子非均匀、温度相关输运、电—力反馈、不同热源定义及离散误差。本文将这些来源保留为竞争解释，不从相关性直接建立单项机制。

![图8 功率分项与实际温升对常系数电热参考的偏离，L0探索性](C:/Users/Administrator/Desktop/codex-workspaces/Console3/solver_v2/paper_sprint/20261008_2110/figures/fig08_power_and_reference.png)

**图8** 漂移热、带符号/绝对值扩散功、端面电功率，以及温升相对式(15)的体积加权均方根偏离。Fig. 8 Power components and volume-weighted RMS departure from the electrothermal reference. 所有模型使用相同σ₀、k和3 V；**历史L0探索性，未完成空间验收**。H1—NT的参考偏离均方根为0.1332175 K，H0—NT为0.1822167 K；它们是合成模型偏离指标，不称为纯机制贡献。

本研究的物理边界还受到材料定标、热源和接触条件约束。表1参数不是同一1 μm厚ZnO样品的统一测量集，空穴参数也未针对双极输运标定；接触假定同时固定电势、浓度与温度，孔壁采用理想自然边界。现有模型可以用于研究这些约定下的耦合响应及数值可靠性，但尚不足以给出真实器件寿命、击穿、效率或工艺优化结论。

## 6 结论

本文建立了无孔与圆孔ZnO板在固定电压和等温接触下的三维热—电—力及双载流子有限元表述，明确区分LQ顺序参考、N0常输运非线性模型及NT温度相关输运模型。模型热源为漂移电导项，NT−N0包含迁移率和扩散系数的共同温度反馈，LQ最终组合状态具有分阶段求解的适用限制。

常系数电热关系表明，在电绝缘、绝热侧壁和两端等温的特定参考问题中，几何改变总功率并不必然改变峰温。端面电功率、扩散功、漂移热和外边界热流的分别核对，是解释当前约化热源结果的必要条件。这些分析关系提供结果解释的参照，而不是有限元全场已验证的替代证据。

现有主工况的归档文件身份已核，H0三例仍缺当时逐case运行来源记录，局部载流子及模型差值的空间证据也未全部闭合。当前工作稿不宣称完整场网格无关，不给出已验证的圆孔增强幅度或独立温度机制贡献，也不声称完成实体实验验证。待同条件几何比较、NT−N0差值及机械观测量的证据闭合后，再据实际数据修订摘要和定量结论。

## 参考文献

[1] Zhao M, Yan X, Wang B, Zhang Q. Finite element formulation for piezoelectric semiconductor plates[J]. Materials Today Communications, 2022, 30: 103098. DOI: [10.1016/j.mtcomm.2021.103098](https://doi.org/10.1016/j.mtcomm.2021.103098). 出版方可检索的[摘要与正文片段](https://www.sciencedirect.com/science/article/pii/S2352492821010825)已核对；直接全文访问受限，本文仅采用已可读范围支持引言比较。

[2] Kong D, He Z, Liu C, Zhang C. Analysis of multi-field coupling behaviors of sandwich piezoelectric semiconductor beams under thermal loadings[J]. Applied Mathematics and Mechanics (English Edition), 2025, 46(8): 1571–1590. DOI: [10.1007/s10483-025-3284-7](https://doi.org/10.1007/s10483-025-3284-7). [出版方页面](https://link.springer.com/article/10.1007/s10483-025-3284-7)已核对作者、题名、年份、卷页及摘要支持范围；未据订阅全文添加未经核读的公式或结果。

## 附录A 写作与证据说明（投稿前移入补充材料或内部审查记录）

本稿的一句话论证是：在两端等温、电压驱动的无孔/圆孔ZnO板中，用三维耦合有限元及独立电热参照建立可核验的比较，同时把未解析的局部场和小模型差值保留为未闭合结论。正文依次回答研究对象、模型、可信性、比较与适用边界；当前证据不足以将论文定位为已发现强增强规律。

术语统一为“漂移电导型热源”“温度相关输运”“顺序三阶段LQ参考”“空间敏感性”“探索性结果”。π只表示热释电系数，p只表示空穴浓度，κ表示介电张量，k表示热导率，D表示电位移，Dₙ/Dₚ表示扩散张量。T与θ分别是绝对温度和相对温升。

正文数值以本次证据登记与新后处理文件为索引。F1历史L0数值在第5节紧邻标明未完成空间验收，不进入摘要和定量结论；后续定向网格结果与历史L0的离散形式、全/半域差异单列。未提供作者、机构、基金、实验或未核实参考文献。引用[1]、[2]支持引言中的模型比较，不支持本程序已获COMSOL或实验验证的表述。

本轮已将八幅实际图件、六例全域L0主指标、V02解析误差和小型条带重放结果写入正文，其中图1为几何示意，图2为历史合成验证，另六图为真实场后处理。图源为data/metrics_summary.csv、data/peclet_regions.csv、data/*midplane.csv、data/A1_error_source.csv、data/A2_error_source.csv及data/profiles/中的实际文件；大体积原始场保留原位置。当前未完成事项包括：H0早期三例缺同时记录的逐case运行provenance，且本轮未重新编译证明归档二进制可复现；NT−N0同系列网格差值及图5；局部载流子空间失败项；机械观测量空间敏感性；小型输运全规格以及中文期刊当前模板。详细状态以audit/open_items.md及证据登记为准，避免以正文概括覆盖逐项失败。

《固体力学学报》官方入口已查阅，当前可访问页面主要给出英文刊链接。学会网站可访问的旧征稿简则支持采用中英文题名/摘要、双语图表题、SI单位与顺序编码参考文献，但其年代较早，不能据此声称已满足现行字数或格式硬性要求。本稿暂按清晰A4工作稿组织，当前中文模板核验保留为PENDING。

