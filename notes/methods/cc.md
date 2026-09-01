# $\alpha$ 衰变的耦合道描述

### 物理体系与通道基

考虑母核 $P$ 衰变为子核 $D$ 与 $\alpha$ 粒子的过程。$P$ 的总角动量、投影和宇称为 $J$、$M$ 和 $\pi_P$，$D$ 的相应量为 $I$、$M_I$ 和 $\pi_D$。$\alpha$ 粒子视为无内部自由度的 $0^+$ 粒子：

$$
P(J^{\pi_P})\rightarrow D(I^{\pi_D})+\alpha(0^+;l).
$$

取子核质心为实验室系原点，$\mathbf r_\alpha$ 为 $\alpha$ 粒子的坐标：

$$
r_\alpha=|\mathbf r_\alpha|,\qquad \hat{\mathbf r}_\alpha=\frac{\mathbf r_\alpha}{r_\alpha}.
$$

总角动量 $J$ 及其投影 $M$ 守恒。$(l,m_l)$ 为 $\alpha$ 粒子的轨道角动量及其投影，$\Phi_{IM_I}(D)$ 为子核态。每组可耦合到 $J$ 的 $\{l,I\}$ 构成一个通道 $c$：

$$
c=\{l,I\},\qquad |I-l|\leq J\leq I+l,\qquad \pi_P=\pi_D(-1)^l.
$$

记 $\Phi_{lm_l}^{(c)}(\mathbf r_\alpha)$ 为通道 $c$ 中的 $\alpha$ 粒子波函数。通道波函数为

$$
\Psi^{(c)}_{JM}(\mathbf r_\alpha)\equiv\big[\Phi_l^{(c)}(\mathbf r_\alpha)\otimes\Phi_I(D)\big]_{JM}=\sum_{m_l,M_I}\langle lm_l;IM_I|JM\rangle\Phi_{lm_l}^{(c)}(\mathbf r_\alpha)\Phi_{IM_I}(D).
$$

其中 $u_c$ 为约化径向波函数，$Y_{lm_l}$ 为球谐函数：

$$
\Phi_{lm_l}^{(c)}(\mathbf r_\alpha)=\frac{u_c(r_\alpha)}{r_\alpha}Y_{lm_l}(\hat{\mathbf r}_\alpha).
$$

对所有允许通道求和：

$$
\boxed{\Psi_{JM}(\mathbf r_\alpha)=\sum_{c}a_c\Psi_{JM}^{(c)}(\mathbf r_\alpha)=\sum_{c}a_c\frac{u_c(r_\alpha)}{r_\alpha}\big[Y_l(\hat{\mathbf r}_\alpha)\otimes\Phi_I(D)\big]_{JM}}.
$$

### 耦合道薛定谔方程

以 $\hat H_D$、$\hat T_\alpha$ 和 $\hat V(\mathbf r_\alpha)$ 分别表示子核内部哈密顿量、$\alpha$ 粒子动能和相对运动势场。系统哈密顿量为

$$
\hat H=\hat H_D+\hat T_\alpha(\mathbf r_\alpha)+\hat V(\mathbf r_\alpha).
$$

子核态的本征能量为 $E_I$：

$$
\hat H_D\Phi_{IM_I}(D)=E_I\Phi_{IM_I}(D).
$$

以 $\mu$ 表示约化质量，$\hat{\mathbf L}$ 表示 $\alpha$ 粒子的轨道角动量算符。动能算符为

$$
\hat T_\alpha=-\frac{\hbar^2}{2\mu}\left[\frac{1}{r_\alpha^2}\frac{d}{dr_\alpha}\left(r_\alpha^2\frac{d}{dr_\alpha}\right)-\frac{\hat{\mathbf L}^2}{\hbar^2r_\alpha^2}\right].
$$

对于通道 $c=\{l,I\}$ 和 $c'=\{l',I'\}$，定义势场矩阵元

$$
V_{cc'}^{JM}(r_\alpha)\equiv\left\langle\big[Y_l(\hat{\mathbf r}_\alpha)\otimes\Phi_I(D)\big]_{JM}\right|\hat V(\mathbf r_\alpha)\left|\big[Y_{l'}(\hat{\mathbf r}_\alpha)\otimes\Phi_{I'}(D)\big]_{JM}\right\rangle.
$$

以 $E$ 表示系统总能量，通道能量为

$$
E_c\equiv E-E_I.
$$

以 $\delta_{cc'}$ 表示 Kronecker 符号。将薛定谔方程投影到通道 $c$，得到

$$
\boxed{\sum_{c'}\left[-\frac{\hbar^2}{2\mu}\left(\frac{d^2}{dr_\alpha^2}-\frac{l(l+1)}{r_\alpha^2}\right)\delta_{cc'}+V_{cc'}^{JM}(r_\alpha)\right]u_{c'}(r_\alpha)=E_cu_c(r_\alpha)}.
$$

### 多极势场与耦合矩阵元

以 $\lambda$ 和 $\mu$ 表示多极阶数及其分量，$\hat V_{\lambda\mu}(r_\alpha)$ 为作用于子核态空间的球张量分量。势场展开为

$$
\hat V(\mathbf r_\alpha)=\sum_{\lambda,\mu}\hat V_{\lambda\mu}(r_\alpha)Y_{\lambda\mu}^*(\hat{\mathbf r}_\alpha)=\sum_{\lambda,\mu}(-1)^\mu\hat V_{\lambda\mu}(r_\alpha)Y_{\lambda,-\mu}(\hat{\mathbf r}_\alpha).
$$

对于通道 $c=\{l,I\}$ 和 $c'=\{l',I'\}$，势场矩阵元为

$$
V_{cc'}^{JM}(r_\alpha)=\left\langle\big[Y_l(\hat{\mathbf r}_\alpha)\otimes\Phi_I(D)\big]_{JM}\right|\sum_{\lambda,\mu}(-1)^\mu\hat V_{\lambda\mu}(r_\alpha)Y_{\lambda,-\mu}(\hat{\mathbf r}_\alpha)\left|\big[Y_{l'}(\hat{\mathbf r}_\alpha)\otimes\Phi_{I'}(D)\big]_{JM}\right\rangle.
$$

展开角动量耦合：

$$
V_{cc'}^{JM}(r_\alpha)=\sum_{\lambda,\mu}\sum_{m_l,M_I}\sum_{m_l',M_I'}(-1)^\mu\langle JM|lm_l;IM_I\rangle\langle l'm_l';I'M_I'|JM\rangle\langle lm_l|Y_{\lambda,-\mu}|l'm_l'\rangle\left\langle\Phi_{IM_I}(D)\right|\hat V_{\lambda\mu}(r_\alpha)\left|\Phi_{I'M_I'}(D)\right\rangle.
$$

由 Wigner–Eckart 定理：

$$
\left\langle\Phi_{IM_I}(D)\right|\hat V_{\lambda\mu}(r_\alpha)\left|\Phi_{I'M_I'}(D)\right\rangle=(-1)^{I-M_I}\begin{pmatrix}I&\lambda&I'\\-M_I&\mu&M_I'\end{pmatrix}\left\langle\Phi_I(D)\right\|\hat V_\lambda(r_\alpha)\left\|\Phi_{I'}(D)\right\rangle.
$$

球谐函数矩阵元为

$$
\langle lm_l|Y_{\lambda,-\mu}|l'm_l'\rangle=(-1)^{l-m_l}\begin{pmatrix}l&\lambda&l'\\-m_l&-\mu&m_l'\end{pmatrix}\langle l\|Y_\lambda\|l'\rangle,
$$

其中

$$
\langle l\|Y_\lambda\|l'\rangle=(-1)^l\sqrt{\frac{(2l+1)(2\lambda+1)(2l'+1)}{4\pi}}\begin{pmatrix}l&\lambda&l'\\0&0&0\end{pmatrix}.
$$

完成磁量子数求和：

$$
\boxed{V_{cc'}^{JM}(r_\alpha)=\sum_\lambda(-1)^{l'+I+J}\begin{Bmatrix}l&I&J\\I'&l'&\lambda\end{Bmatrix}\langle l\|Y_\lambda\|l'\rangle\left\langle\Phi_I(D)\right\|\hat V_\lambda(r_\alpha)\left\|\Phi_{I'}(D)\right\rangle}.
$$

非零耦合满足

$$
|l-l'|\leq\lambda\leq l+l',\qquad |I-I'|\leq\lambda\leq I+I',\qquad l+l'+\lambda\in2\mathbb Z.
$$

相关约定见 [Wigner Symbols](../mathematical_physics/wigner_symbols.md)。

### 子核的轴对称转动态

子核采用轴对称强耦合转子模型：

- 内禀对称轴为 $z'$；
- 强耦合：子核角动量在 $z'$ 轴上的投影 $K$ 守恒；
- 集体转动与内禀运动满足绝热分离。

以 $\Omega_D=(\varphi_D,\theta_D,\psi_D)$ 表示子核取向的 Euler 角。内禀态记为 $\chi_K(D)$，其时间反演伙伴为 $\chi_{-K}(D)$。

只考虑一个固定 $K$ 和固定宇称 $\pi_D$ 的转动带。第一节中的子核态记号表示

$$
\Phi_{IM_I}(D)\equiv\Phi_{IM_IK}(D).
$$

对于 $K>0$，归一化子核态为

$$
\boxed{\Phi_{IM_IK}(D)=\sqrt{\frac{2I+1}{16\pi^2}}\left[D_{M_IK}^{I*}(\Omega_D)\chi_K(D)+(-1)^{I+K}D_{M_I,-K}^{I*}(\Omega_D)\chi_{-K}(D)\right]}.
$$

对于 $K=0$，归一化子核态为

$$
\boxed{\Phi_{IM_I0}(D)=\sqrt{\frac{2I+1}{8\pi^2}}D_{M_I0}^{I*}(\Omega_D)\chi_0(D)}.
$$

允许的角动量及其投影满足

$$
I=|K|,|K|+1,\ldots,\qquad M_I=-I,-I+1,\ldots,I.
$$

对于偶偶子核的 $K=0$ 基态转动带，

$$
I=0,2,4,\ldots.
$$

Wigner-$D$ 函数的定义与约定见 [SO(3) and SU(2)](../group_theory/so3_su2.md#wigner-d-functions)。

### 内禀系到实验室系的几何变换

$\mathbf r_\alpha$ 与 $\mathbf r_\alpha'$ 表示同一 $\alpha$ 粒子坐标在实验室系与子核内禀系中的分量，两套坐标系均以子核质心为原点。坐标变换保持径向长度不变：

$$
r_\alpha=|\mathbf r_\alpha|=|\mathbf r_\alpha'|.
$$

以 $\hat{\mathbf n}_D$ 表示实验室系中的子核对称轴，以 $\hat{\mathbf z}'$ 表示内禀对称轴。内禀系中的夹角 $\theta'$ 满足

$$
\boxed{\cos\theta'=\hat{\mathbf r}_\alpha'\cdot\hat{\mathbf z}'=\hat{\mathbf r}_\alpha\cdot\hat{\mathbf n}_D}.
$$

以 $\hat V(\mathbf r_\alpha)$ 和 $\hat V'(\mathbf r_\alpha')$ 分别表示同一势场在实验室系与内禀系中的形式。轴对称性使内禀系展开只保留 $\mu=0$：

$$
\hat V'(\mathbf r_\alpha')=\sum_\lambda\hat V'_{\lambda0}(r_\alpha)Y_{\lambda0}(\hat{\mathbf r}_\alpha')=\sum_\lambda\hat V'_{\lambda0}(r_\alpha)\sqrt{\frac{2\lambda+1}{4\pi}}P_\lambda(\cos\theta').
$$

其中 $P_\lambda$ 为 Legendre 多项式。球谐函数加法定理给出

$$
P_\lambda(\hat{\mathbf r}_\alpha\cdot\hat{\mathbf n}_D)=\frac{4\pi}{2\lambda+1}\sum_\mu Y_{\lambda\mu}^*(\hat{\mathbf r}_\alpha)Y_{\lambda\mu}(\hat{\mathbf n}_D).
$$

代入内禀系展开：

$$
\hat V'(\mathbf r_\alpha')=\sum_{\lambda,\mu}\left[\sqrt{\frac{4\pi}{2\lambda+1}}\hat V'_{\lambda0}(r_\alpha)Y_{\lambda\mu}(\hat{\mathbf n}_D)\right]Y_{\lambda\mu}^*(\hat{\mathbf r}_\alpha).
$$

同一势场的两种表示满足

$$
\boxed{\hat V(\mathbf r_\alpha)=\hat V'(\mathbf r_\alpha')=\sum_{\lambda,\mu}\hat V_{\lambda\mu}(r_\alpha)Y_{\lambda\mu}^*(\hat{\mathbf r}_\alpha)=\sum_\lambda\hat V'_{\lambda0}(r_\alpha)Y_{\lambda0}(\hat{\mathbf r}_\alpha')}.
$$

比较实验室系中的球谐函数系数：

$$
\hat V_{\lambda\mu}(r_\alpha)=\sqrt{\frac{4\pi}{2\lambda+1}}\hat V'_{\lambda0}(r_\alpha)Y_{\lambda\mu}(\hat{\mathbf n}_D).
$$

利用

$$
Y_{\lambda\mu}(\hat{\mathbf n}_D)=\sqrt{\frac{2\lambda+1}{4\pi}}D_{\mu0}^{\lambda *}(\Omega_D),
$$

得到

$$
\boxed{\hat V_{\lambda\mu}(r_\alpha)=\hat V'_{\lambda0}(r_\alpha)D_{\mu0}^{\lambda *}(\Omega_D)}.
$$

### 转动态的势场矩阵元

记

$$
\hat I\equiv\sqrt{2I+1},\qquad \hat I'\equiv\sqrt{2I'+1},
$$

并定义内禀矩阵元

$$
V'_{\lambda0;K,K'}(r_\alpha)\equiv\left\langle\chi_K(D)\right|\hat V'_{\lambda0}(r_\alpha)\left|\chi_{K'}(D)\right\rangle.
$$

对于 $K,K'>0$，将转动态与

$$
\hat V_{\lambda\mu}(r_\alpha)=\hat V'_{\lambda0}(r_\alpha)D_{\mu0}^{\lambda *}(\Omega_D)
$$

代入矩阵元，得到

$$
\begin{aligned}
&\left\langle\Phi_{IM_IK}(D)\right|\hat V_{\lambda\mu}(r_\alpha)\left|\Phi_{I'M_I'K'}(D)\right\rangle\\
={}&\frac{\hat I\hat I'}{16\pi^2}\left[\int d\Omega_D\,D^I_{M_IK}D^{\lambda *}_{\mu0}D^{I'*}_{M_I'K'}\right]V'_{\lambda0;K,K'}(r_\alpha)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I'+K'}\left[\int d\Omega_D\,D^I_{M_IK}D^{\lambda *}_{\mu0}D^{I'*}_{M_I',-K'}\right]V'_{\lambda0;K,-K'}(r_\alpha)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I+K}\left[\int d\Omega_D\,D^I_{M_I,-K}D^{\lambda *}_{\mu0}D^{I'*}_{M_I'K'}\right]V'_{\lambda0;-K,K'}(r_\alpha)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I+K+I'+K'}\left[\int d\Omega_D\,D^I_{M_I,-K}D^{\lambda *}_{\mu0}D^{I'*}_{M_I',-K'}\right]V'_{\lambda0;-K,-K'}(r_\alpha).
\end{aligned}
$$

利用 Wigner-$D$ 函数的复共轭关系

$$
D_{mk}^{j*}(\Omega_D)=(-1)^{m-k}D_{-m,-k}^{j}(\Omega_D),
$$

上式化为

$$
\begin{aligned}
&\left\langle\Phi_{IM_IK}(D)\right|\hat V_{\lambda\mu}(r_\alpha)\left|\Phi_{I'M_I'K'}(D)\right\rangle\\
={}&\frac{\hat I\hat I'}{16\pi^2}(-1)^{\mu+M_I'-K'}\left[\int d\Omega_D\,D^I_{M_IK}D^\lambda_{-\mu,0}D^{I'}_{-M_I',-K'}\right]V'_{\lambda0;K,K'}(r_\alpha)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I'+K'}(-1)^{\mu+M_I'+K'}\left[\int d\Omega_D\,D^I_{M_IK}D^\lambda_{-\mu,0}D^{I'}_{-M_I',K'}\right]V'_{\lambda0;K,-K'}(r_\alpha)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I+K}(-1)^{\mu+M_I'-K'}\left[\int d\Omega_D\,D^I_{M_I,-K}D^\lambda_{-\mu,0}D^{I'}_{-M_I',-K'}\right]V'_{\lambda0;-K,K'}(r_\alpha)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I+K+I'+K'}(-1)^{\mu+M_I'+K'}\left[\int d\Omega_D\,D^I_{M_I,-K}D^\lambda_{-\mu,0}D^{I'}_{-M_I',K'}\right]V'_{\lambda0;-K,-K'}(r_\alpha).
\end{aligned}
$$

三重 Wigner-$D$ 函数积分满足

$$
\int d\Omega_D\,D^{(j_1)}_{\mu_1\nu_1}(\Omega_D)D^{(j_2)}_{\mu_2\nu_2}(\Omega_D)D^{(j_3)}_{\mu_3\nu_3}(\Omega_D)=8\pi^2\begin{pmatrix}j_1&j_2&j_3\\\mu_1&\mu_2&\mu_3\end{pmatrix}\begin{pmatrix}j_1&j_2&j_3\\\nu_1&\nu_2&\nu_3\end{pmatrix}.
$$

取

$$
(j_1,j_2,j_3)=(I,\lambda,I'),
$$

并将各项对应的 $\mu_i$ 与 $\nu_i$ 分别代入，得到

$$
\begin{aligned}
&\left\langle\Phi_{IM_IK}(D)\right|\hat V_{\lambda\mu}(r_\alpha)\left|\Phi_{I'M_I'K'}(D)\right\rangle\\
={}&\frac{\hat I\hat I'}{2}\begin{pmatrix}I&\lambda&I'\\M_I&-\mu&-M_I'\end{pmatrix}\Bigg\{(-1)^{\mu+M_I'-K'}\begin{pmatrix}I&\lambda&I'\\K&0&-K'\end{pmatrix}V'_{\lambda0;K,K'}(r_\alpha)\\
&\qquad+(-1)^{I'+K'}(-1)^{\mu+M_I'+K'}\begin{pmatrix}I&\lambda&I'\\K&0&K'\end{pmatrix}V'_{\lambda0;K,-K'}(r_\alpha)\\
&\qquad+(-1)^{I+K}(-1)^{\mu+M_I'-K'}\begin{pmatrix}I&\lambda&I'\\-K&0&-K'\end{pmatrix}V'_{\lambda0;-K,K'}(r_\alpha)\\
&\qquad+(-1)^{I+K+I'+K'}(-1)^{\mu+M_I'+K'}\begin{pmatrix}I&\lambda&I'\\-K&0&K'\end{pmatrix}V'_{\lambda0;-K,-K'}(r_\alpha)\Bigg\}.
\end{aligned}
$$

四个体轴 Wigner-$3j$ 符号的下排之和必须为零，因此四项分别要求

$$
K=K',\qquad K=-K',\qquad K=-K',\qquad K=K'.
$$

由于 $K,K'>0$，第二项和第三项为零，并有 $K'=K$。矩阵元化为

$$
\begin{aligned}
&\left\langle\Phi_{IM_IK}(D)\right|\hat V_{\lambda\mu}(r_\alpha)\left|\Phi_{I'M_I'K}(D)\right\rangle=\frac{\hat I\hat I'}{2}\begin{pmatrix}I&\lambda&I'\\M_I&-\mu&-M_I'\end{pmatrix}\Bigg[(-1)^{\mu+M_I'-K}\begin{pmatrix}I&\lambda&I'\\K&0&-K\end{pmatrix}V'_{\lambda0;K,K}(r_\alpha)\\
&\qquad+(-1)^{I+I'+2K}(-1)^{\mu+M_I'+K}\begin{pmatrix}I&\lambda&I'\\-K&0&K\end{pmatrix}V'_{\lambda0;-K,-K}(r_\alpha)\Bigg].
\end{aligned}
$$

第一个 Wigner-$3j$ 符号要求 $M_I-\mu-M_I'=0$。再利用反转性质

$$
\begin{pmatrix}I&\lambda&I'\\-K&0&K\end{pmatrix}=(-1)^{I+\lambda+I'}\begin{pmatrix}I&\lambda&I'\\K&0&-K\end{pmatrix},
$$

得到

$$
\left\langle\Phi_{IM_IK}(D)\right|\hat V_{\lambda\mu}(r_\alpha)\left|\Phi_{I'M_I'K}(D)\right\rangle
=\frac{\hat I\hat I'}{2}(-1)^{M_I-K}\begin{pmatrix}I&\lambda&I'\\M_I&-\mu&-M_I'\end{pmatrix}\begin{pmatrix}I&\lambda&I'\\K&0&-K\end{pmatrix}\left[V'_{\lambda0;K,K}(r_\alpha)+(-1)^\lambda V'_{\lambda0;-K,-K}(r_\alpha)\right].
$$

若势场具有时间反演对称性、内禀矩阵元取实，并且仅保留偶数 $\lambda$，则

$$
V'_{\lambda0;-K,-K}(r_\alpha)=V'_{\lambda0;K,K}(r_\alpha),
$$

从而

$$
\left\langle\Phi_{IM_IK}(D)\right|\hat V_{\lambda\mu}(r_\alpha)\left|\Phi_{I'M_I'K}(D)\right\rangle=\hat I\hat I'(-1)^{M_I-K}\begin{pmatrix}I&\lambda&I'\\M_I&-\mu&-M_I'\end{pmatrix}\begin{pmatrix}I&\lambda&I'\\K&0&-K\end{pmatrix}V'_{\lambda0;K,K}(r_\alpha).
$$

等价地，

$$
=\hat I\hat I'(-1)^{I-M_I}\begin{pmatrix}I&\lambda&I'\\-M_I&\mu&M_I'\end{pmatrix}(-1)^{I-K}\begin{pmatrix}I&\lambda&I'\\-K&0&K\end{pmatrix}V'_{\lambda0;K,K}(r_\alpha).
$$

因此，约化矩阵元为

$$
\left\langle\Phi_I(D)\left\|\hat V_\lambda(r_\alpha)\right\|\Phi_{I'}(D)\right\rangle=\hat I\hat I'(-1)^{I-K}\begin{pmatrix}I&\lambda&I'\\-K&0&K\end{pmatrix}V'_{\lambda0;K,K}(r_\alpha).
$$

对于 $K=0$，约化矩阵元化为

$$
\left\langle\Phi_I(D)\left\|\hat V_\lambda(r_\alpha)\right\|\Phi_{I'}(D)\right\rangle=\hat I\hat I'(-1)^I\begin{pmatrix}I&\lambda&I'\\0&0&0\end{pmatrix}V'_{\lambda0;0,0}(r_\alpha).
$$

### 库伦耦合势

以 $\mathbf r_D'=r_D\hat{\mathbf r}_D'$ 表示子核内相对于其质心的内禀坐标。子核电荷密度 $\rho_D$ 归一化为

$$
\int d^3r_D'\,\rho_D(\mathbf r_D')=1.
$$

内禀系中的库伦势为

$$
V^{\prime C}(\mathbf r_\alpha')=\frac{e^2Z_\alpha Z_D}{4\pi\epsilon_0}\int d^3r_D'\,\frac{\rho_D(\mathbf r_D')}{|\mathbf r_\alpha'-\mathbf r_D'|}.
$$

利用

$$
\frac{1}{|\mathbf r_\alpha'-\mathbf r_D'|}=\sum_{\lambda,\mu}\frac{4\pi}{2\lambda+1}\frac{r_{<}^\lambda}{r_{>}^{\lambda+1}}Y_{\lambda\mu}^*(\hat{\mathbf r}_\alpha')Y_{\lambda\mu}(\hat{\mathbf r}_D'),\qquad r_{<}=\min(r_\alpha,r_D),\quad r_{>}=\max(r_\alpha,r_D),
$$

对于积分域内处处满足 $r_\alpha\geq r_D$ 的子核外部区域，定义

$$
V^{\prime C}(\mathbf r_\alpha')=\frac{eZ_\alpha}{4\pi\epsilon_0}\sum_{\lambda,\mu}\frac{4\pi}{2\lambda+1}\frac{Q_{\lambda\mu}^D}{r_\alpha^{\lambda+1}}Y_{\lambda\mu}^*(\hat{\mathbf r}_\alpha'),\qquad Q_{\lambda\mu}^D=eZ_D\int d^3r_D'\,\rho_D(\mathbf r_D')r_D^\lambda Y_{\lambda\mu}(\hat{\mathbf r}_D').
$$

因此

$$
V_{\lambda0}^{\prime C}(r_\alpha)=\frac{eZ_\alpha}{4\pi\epsilon_0}\frac{4\pi}{2\lambda+1}\frac{Q_{\lambda0}^D}{r_\alpha^{\lambda+1}}.
$$

考虑均匀带电球体，其表面允许轴对称形变：

$$
\rho_D(\mathbf r_D')=\rho_0\Theta\!\left(R_D(\hat{\mathbf r}_D')-r_D\right),\qquad 1=\rho_0\frac{4\pi}{3}(R_0^D)^3,
$$

$$
R_D(\hat{\mathbf r}_D')=R_0^D\left[1+\sum_{\lambda>0}\beta_{\lambda0}^D Y_{\lambda0}(\hat{\mathbf r}_D')\right].
$$

完成径向积分并保留形变参数的一阶项，

$$
Q_{\lambda0}^D\simeq eZ_D\frac{3}{4\pi}(R_0^D)^\lambda\left[\frac{\sqrt{4\pi}}{3}\delta_{\lambda0}+\beta_{\lambda0}^D\right].
$$

对于 $\lambda>0$，

$$
V_{\lambda0;K,K'}^{\prime C}(r_\alpha)=\delta_{KK'}\frac{e^2Z_DZ_\alpha}{4\pi\epsilon_0}\frac{3}{2\lambda+1}\frac{(R_0^D)^\lambda}{r_\alpha^{\lambda+1}}\beta_{\lambda0}^D.
$$

单极项为

$$
V_{00;K,K'}^{\prime C}(r_\alpha)=\delta_{KK'}\frac{e^2Z_DZ_\alpha}{4\pi\epsilon_0}\frac{\sqrt{4\pi}}{r_\alpha}.
$$

将 $V_{\lambda0;K,K'}^{\prime C}(r_\alpha)$ 代入前述 $\left\langle\Phi_I(D)\left\|\hat V_\lambda^C(r_\alpha)\right\|\Phi_{I'}(D)\right\rangle$，即可得到库伦耦合矩阵元。

### 核势耦合

采用 Woods–Saxon 形式描述 $\alpha$ 粒子与子核之间的核势：

$$
\hat V^{\prime N}(\mathbf r_\alpha')=\frac{V_0}{1+\exp\!\left([r_\alpha-R_0^\alpha-\hat R(D)]/a\right)},
$$

其中 $V_0$ 为势深，$a$ 为表面弥散参数，$R_0^\alpha$ 为 $\alpha$ 粒子半径。全部角向依赖由子核半径算符给出：

$$
\hat R(D)=R_0^D\left[1+\sum_{\lambda,\mu}\hat a_{\lambda\mu}Y_{\lambda\mu}^*(\hat{\mathbf r}_\alpha')\right].
$$

仿照前述转动态约化矩阵元，

$$
\left\langle\Phi_{IK}(D)\left\|\hat R_{D,\lambda}\right\|\Phi_{I'K'}(D)\right\rangle=\hat I\hat I'(-1)^{I-K}\begin{pmatrix}I&\lambda&I'\\-K&0&K'\end{pmatrix}\left\langle\chi_K(D)\right|\hat R'_{D,\lambda0}\left|\chi_{K'}(D)\right\rangle.
$$

内禀半径算符的矩阵元为

$$
\left\langle\chi_K(D)\right|\hat R'_{D,\lambda0}\left|\chi_{K'}(D)\right\rangle=\delta_{KK'}R_0^D\left[\sqrt{4\pi}\,\delta_{\lambda0}+\beta_{\lambda0}^D\right].
$$

定义球形 Woods–Saxon 势

$$
V_0^N(r_\alpha)\equiv\frac{V_0}{1+\exp\!\left([r_\alpha-R_0^\alpha-R_0^D]/a\right)}.
$$

对于较大形变，可先构造并对角化半径矩阵，再计算 Woods–Saxon 势，而不作一阶展开。
