# G05 M0 解析解验证冻结规范 v1

## 1. 目的与边界

G05 只验证理论 v2.1 第 15 节冻结的一阶小扰动模型、HEX20 插值与积分、全局 CSR、边界条件、缩放和 PARDISO 求解链。它不把小扰动方程冒充 M1/M2 生产方程，也不证明非线性 Newton 已经正确。

M0 使用与主程序相同的七自由度顺序

\[
[u_1,u_2,u_3,\phi,\eta_n,\eta_p,\theta],
\]

但在一阶参考态附近解释为

\[
\eta_n=\Delta n/n_0,\qquad \eta_p=\Delta p/p_0.
\]

因此 M0 电流严格删除二阶项 \(\Delta nE\) 和 \(\Delta pE\)。所有算例均为一维场嵌入三维直六面体，场只随 \(x\) 变化；横截面只承担三维 HEX20 积分，不引入横向物理。

## 2. 冻结参数（SI）

\[
\begin{aligned}
L&=10^{-6}\ \mathrm m, & C&=10^{11}\ \mathrm{Pa}, & e&=1\ \mathrm{C/m^2},\\
\kappa&=10^{-8}\ \mathrm{F/m}, & n_0&=10^{21}\ \mathrm{m^{-3}}, & p_0&=8\times10^{20}\ \mathrm{m^{-3}},\\
\mu_n&=0.05\ \mathrm{m^2/(Vs)}, & D_n&=1.2926\times10^{-3}\ \mathrm{m^2/s},\\
\mu_p&=0.02\ \mathrm{m^2/(Vs)}, & D_p&=5.1704\times10^{-4}\ \mathrm{m^2/s},\\
k&=100\ \mathrm{W/(mK)}, & \alpha&=10^{-7}\ \mathrm{K^{-1}}, & \pi&=10^{-8}\ \mathrm{C/(m^2K)}.
\end{aligned}
\]

电子与空穴扩散系数取相同的 \(\mu/D\) 比例。电压幅值取 \(V=10^{-3}\ \mathrm V\)，热耦合温差取 \(\Delta T=1\ \mathrm K\)，确保 \(|\eta_n|,|\eta_p|\ll1\)。

## 3. 算例定义

### V02：力—电屏蔽耦合

以零电子电流关系消去 \(\Delta n\)：

\[
\Delta n=A_n\phi,\qquad A_n=n_0\mu_n/D_n.
\]

有限元弱式保持原始符号：

\[
R_u=\int v'(Cu'+e\phi')\,\mathrm d\Omega,
\]

\[
R_\phi=\int \psi'(eu'-\kappa\phi')\,\mathrm d\Omega
-\int \psi qA_n\phi\,\mathrm d\Omega.
\]

边界为 \(u(0)=0,\phi(0)=0,\phi(L)=V\)，右端力学边界为自然零应力。精确解为

\[
\phi=V\frac{\sinh(\gamma_nx)}{\sinh(\gamma_nL)},\qquad
u=-\frac eC\phi,
\]

\[
\gamma_n^2=\frac{qn_0\mu_n}{D_n(\kappa+e^2/C)}.
\]

### V03：电子漂移—扩散符号

规定线性电势 \(\phi=Vx/L\)，求解

\[
R_n=\int w'qn_0(-\mu_n\phi'+D_n\eta_n')\,\mathrm d\Omega=0,
\]

并取 \(\eta_n(0)=0\)、右端零电子电流。精确关系为

\[
\eta_n=(\mu_n/D_n)\phi>0,
\qquad J_n=0.
\]

### V04：空穴漂移—扩散符号

同样规定 \(\phi=Vx/L\)，求解

\[
R_p=\int w'qp_0(-\mu_p\phi'-D_p\eta_p')\,\mathrm d\Omega=0,
\]

并取 \(\eta_p(0)=0\)、右端零空穴电流。精确关系为

\[
\eta_p=-(\mu_p/D_p)\phi<0,
\qquad J_p=0.
\]

### V05：一维导热

包含两个子算例：

1. 无热源、两端温差：\(\theta=\Delta T x/L\)；
2. 均匀热源、两端零温升：\(\theta=Q_0x(L-x)/(2k)\)。

第二个子算例同时比较左右端数值热流之和与 \(Q_0AL\)。

### V06：热—压电—电子联合解

联合弱式为

\[
\begin{aligned}
R_u&=\int v'(Cu'+e\phi'-C\alpha\theta),\\
R_\phi&=\int\psi'(eu'-\kappa\phi'+\pi\theta)-\int\psi qn_0\eta_n,\\
R_n&=\int w'qn_0(-\mu_n\phi'+D_n\eta_n'),\\
R_T&=\int\zeta'k\theta'.
\end{aligned}
\]

边界为 \(u(0)=0\)，\(\phi(0)=\phi(L)=0\)，\(\eta_n(0)=0\)，\(\theta(0)=0\)，\(\theta(L)=\Delta T\)。精确解为

\[
\theta=\Delta T x/L,
\]

\[
\phi=\frac{s}{\gamma_n^2}
\left[\frac{\cosh\{\gamma_n(x-L/2)\}}{\cosh(\gamma_nL/2)}-1\right],
\quad
s=\frac{(\pi+e\alpha)\Delta T}{(\kappa+e^2/C)L},
\]

\[
\eta_n=(\mu_n/D_n)\phi,qquad
u=\frac{\alpha\Delta T}{2L}x^2-\frac eC\phi.
\]

## 4. 预先冻结的通过条件

- V02 和 V06 使用 \(1,2,4,8,16\) 个轴向 HEX20 单元；末级相对 \(L_2\) 误差小于 \(2\times10^{-5}\)，相对 \(H^1\) 半范误差小于 \(2\times10^{-3}\)。
- V02、V06 末两级观测阶：\(L_2>2.7\)，\(H^1>1.8\)。
- V03、V04 的浓度变量相对 \(L_2\) 误差与归一化最大电流均小于 \(10^{-10}\)，并分别满足末端符号正、负。
- V05 两个多项式解的相对误差小于 \(10^{-11}\)；均匀热源能量不平衡小于 \(10^{-11}\)。
- 所有算例必须连续两次完全相同，并由命令行构建和 Visual Studio 构建交叉复算。

阈值在首次运行前冻结；若未通过，必须检查方程、边界、缩放或实现，不能通过放宽阈值获得 PASS。

## 5. 执行审计记录

- 初始计划使用 \(1,2,4,8\) 级网格。首次用与求解积分相同的 3×3×3 Gauss 计算误差时，V06 最粗网格出现非单调假象，确认属于误差积分混叠。
- 误差测量随后改为独立 5×5×5 Gauss；方程和数值解不变。独立积分下，V06 的 8 单元 \(L_2\) 误差为 \(2.18\times10^{-5}\)，高于冻结阈值 \(2.0\times10^{-5}\)，因此该次运行判为失败。
- 不放宽阈值、不改变材料参数；只把收敛序列延伸到 16 单元，以确认预期渐近阶并达到原阈值。此修改提高验证分辨率，不删除上述失败记录。
