# $\alpha$ 衰变的耦合道描述

### 物理体系与通道基

母核态 $P$ 的总角动量、投影和宇称为 $(J,M,\pi_P)$；子核态 $D$ 的相应量为 $(I,M_I,\pi_I)$；$\alpha$ 粒子处于 $0^+$ 基态。相对运动的轨道角动量量子数为 $(l,m)$。

$$
P(J^{\pi_P}M)\rightarrow D(I^{\pi_I}M_I)+\alpha(0^+;lm).
$$

设 $\mathbf r_\alpha$ 和 $\mathbf r_D$ 为两个碎片的质心坐标，则

$$
\mathbf r=\mathbf r_\alpha-\mathbf r_D,\qquad r=|\mathbf r|,\qquad \hat{\mathbf r}=\mathbf r/r.
$$

子核与 $\alpha$ 粒子的内部坐标分别为 $\xi_D$ 和 $\xi_\alpha$，归一化内部态分别为 $\Phi_{I M_I}(\xi_D)$ 和 $\Phi_\alpha(\xi_\alpha)$。

固定 $(J,M,\pi_P)$，并假设每个 $I$ 只对应一个纳入计算的子核态。定义通道

$$
c=(I,l).
$$

允许通道集合 $\mathcal C_{J\pi_P}$ 满足

$$
|I-l|\leq J\leq I+l,\qquad \pi_P=\pi_I\pi_\alpha(-1)^l,\qquad \pi_\alpha=+1.
$$

以 $Y_{lm}$ 表示球谐函数，以 $\langle lm;IM_I|JM\rangle$ 表示 Clebsch–Gordan 系数，则

$$
\mathcal Y_{JM}^{(c)}=\big[Y_l\otimes\Phi_I\big]_{JM}\Phi_\alpha=\sum_{m,M_I}\langle lm;IM_I|JM\rangle Y_{lm}(\hat{\mathbf r})\Phi_{I M_I}(\xi_D)\Phi_\alpha(\xi_\alpha).
$$

以 $u_c(r)$ 表示通道 $c$ 的约化径向波函数，并定义相应的通道分量

$$
\Psi_{JM}^{(c)}(\mathbf r,\xi_D,\xi_\alpha)\equiv\frac{u_c(r)}{r}\mathcal Y_{JM}^{(c)}(\hat{\mathbf r},\xi_D,\xi_\alpha).
$$

总波函数为

$$
\boxed{\Psi_{JM}=\sum_{c\in\mathcal C_{J\pi_P}}\Psi_{JM}^{(c)}}.
$$

径向坐标保持显式，不包含在通道标记 $c$ 中。

### 耦合道薛定谔方程

以 $\hat H_D$、$\hat H_\alpha$、$\hat T_r$ 和 $\hat V$ 分别表示子核内部哈密顿量、$\alpha$ 粒子内部哈密顿量、相对运动动能和碎片间相互作用，则

$$
\hat H=\hat H_D(\xi_D)+\hat H_\alpha(\xi_\alpha)+\hat T_r(\mathbf r)+\hat V(\mathbf r,\xi_D,\xi_\alpha).
$$

子核态与 $\alpha$ 粒子态的本征能量分别为 $E_I$ 和 $E_\alpha$：

$$
\hat H_D\Phi_{IM_I}=E_I\Phi_{IM_I},\qquad \hat H_\alpha\Phi_\alpha=E_\alpha\Phi_\alpha.
$$

设 $m_D$ 和 $m_\alpha$ 为碎片质量，$\mu=m_Dm_\alpha/(m_D+m_\alpha)$ 为约化质量，$\hat{\mathbf L}$ 为相对轨道角动量算符，$\hbar$ 为约化普朗克常数，则

$$
\hat T_r=-\frac{\hbar^2}{2\mu}\left[\frac{1}{r^2}\frac{d}{dr}\left(r^2\frac{d}{dr}\right)-\frac{\hat{\mathbf L}^2}{\hbar^2r^2}\right].
$$

以 $l_c$ 表示通道 $c$ 的轨道角动量，则

$$
\hat T_r\Psi_{JM}^{(c)}=-\frac{\hbar^2}{2\mu r}\left[\frac{d^2}{dr^2}-\frac{l_c(l_c+1)}{r^2}\right]u_c(r)\mathcal Y_{JM}^{(c)}.
$$

总能量 $E$ 满足

$$
(\hat H-E)\Psi_{JM}=0.
$$

对于通道 $c=(I,l)$ 和 $c'=(I',l')$，定义耦合势矩阵元

$$
V^{JM}_{cc'}(r)\equiv\big\langle \mathcal Y_{JM}^{(c)}\big|\hat V\big|\mathcal Y_{JM}^{(c')}\big\rangle_{\hat{\mathbf r},\xi_D,\xi_\alpha}.
$$

定义通道能量

$$
E_c\equiv E-E_I-E_\alpha.
$$

以 $\delta_{cc'}$ 表示 Kronecker 符号。将薛定谔方程投影到通道 $c$，得到

$$
\boxed{\sum_{c'}\left[-\frac{\hbar^2}{2\mu}\left(\frac{d^2}{dr^2}-\frac{l_c(l_c+1)}{r^2}\right)\delta_{cc'}+V^{JM}_{cc'}(r)\right]u_{c'}(r)=E_cu_c(r)}.
$$

### 多极相互作用与耦合矩阵元

以 $\lambda$ 和 $\mu$ 表示多极阶数及其分量，$\hat V_{\lambda\mu}(r,\xi_D,\xi_\alpha)$ 表示作用于碎片内部空间的球张量分量。相互作用展开为

$$
\hat V=\sum_{\lambda,\mu}\hat V_{\lambda\mu}(r)Y_{\lambda\mu}^*(\hat{\mathbf r})=\sum_{\lambda,\mu}(-1)^\mu\hat V_{\lambda\mu}(r)Y_{\lambda,-\mu}(\hat{\mathbf r}).
$$

对于 $c=(I,l)$ 和 $c'=(I',l')$，耦合势矩阵元为

$$
V^{JM}_{cc'}(r)=\left\langle\big[Y_l\otimes\Phi_I\big]_{JM}\Phi_\alpha\right|\sum_{\lambda,\mu}(-1)^\mu\hat V_{\lambda\mu}(r)Y_{\lambda,-\mu}\left|\big[Y_{l'}\otimes\Phi_{I'}\big]_{JM}\Phi_\alpha\right\rangle.
$$

展开角动量耦合，

$$
V^{JM}_{cc'}(r)=\sum_{\lambda,\mu}\sum_{m,M_I}\sum_{m',M_I'}(-1)^\mu\langle JM|lm;IM_I\rangle\langle l'm';I'M_I'|JM\rangle\langle lm|Y_{\lambda,-\mu}|l'm'\rangle\left\langle\Phi_{IM_I}\Phi_\alpha\right|\hat V_{\lambda\mu}(r)\left|\Phi_{I'M_I'}\Phi_\alpha\right\rangle.
$$

定义内部约化矩阵元

$$
\mathcal V^\lambda_{II'}(r)\equiv\left\langle\Phi_I\Phi_\alpha\right\|\hat V_\lambda(r)\left\|\Phi_{I'}\Phi_\alpha\right\rangle.
$$

由 Wigner–Eckart 定理，

$$
\left\langle\Phi_{IM_I}\Phi_\alpha\right|\hat V_{\lambda\mu}(r)\left|\Phi_{I'M_I'}\Phi_\alpha\right\rangle=(-1)^{I-M_I}\begin{pmatrix}I&\lambda&I'\\-M_I&\mu&M_I'\end{pmatrix}\mathcal V^\lambda_{II'}(r).
$$

球谐函数矩阵元为

$$
\langle lm|Y_{\lambda,-\mu}|l'm'\rangle=(-1)^{l-m}\begin{pmatrix}l&\lambda&l'\\-m&-\mu&m'\end{pmatrix}\langle l\|Y_\lambda\|l'\rangle.
$$

其中

$$
\langle l\|Y_\lambda\|l'\rangle=(-1)^l\sqrt{\frac{(2l+1)(2\lambda+1)(2l'+1)}{4\pi}}\begin{pmatrix}l&\lambda&l'\\0&0&0\end{pmatrix}.
$$

代回并完成磁量子数求和，

$$
\boxed{V^{JM}_{cc'}(r)=\sum_\lambda(-1)^{l'+I+J}\begin{Bmatrix}l&I&J\\I'&l'&\lambda\end{Bmatrix}\langle l\|Y_\lambda\|l'\rangle\mathcal V^\lambda_{II'}(r)}.
$$

非零耦合满足

$$
|l-l'|\leq\lambda\leq l+l',\qquad |I-I'|\leq\lambda\leq I+I',\qquad l+l'+\lambda\in2\mathbb Z.
$$

相关约定见 [Wigner Symbols](../mathematical_physics/wigner_symbols.md)。

### 子核的轴对称转动态

采用以下假设：

- 轴对称：内禀对称轴为 $z'$；
- 强耦合：子核角动量在 $z'$ 轴上的投影 $K$ 守恒；
- 绝热近似：集体转动与内禀运动可分离。

以 $\Omega_D=(\varphi_D,\theta_D,\psi_D)$ 表示子核的取向，以 $\eta_D$ 表示内禀坐标，则 $\xi_D=(\Omega_D,\eta_D)$。内禀态记为 $\chi_K(\eta_D)$，其时间反演伙伴为 $\chi_{-K}(\eta_D)$。

对于 $K>0$，归一化子核态为

$$
\boxed{\Phi_{IM_IK}(\Omega_D,\eta_D)=\sqrt{\frac{2I+1}{16\pi^2}}\left[D_{M_IK}^{I*}(\Omega_D)\chi_K(\eta_D)+(-1)^{I+K}D_{M_I,-K}^{I*}(\Omega_D)\chi_{-K}(\eta_D)\right]}.
$$

对于 $K=0$，

$$
\boxed{\Phi_{IM_I0}(\Omega_D,\eta_D)=\sqrt{\frac{2I+1}{8\pi^2}}D_{M_I0}^{I*}(\Omega_D)\chi_0(\eta_D)}.
$$

允许的角动量满足

$$
I=|K|,|K|+1,\ldots,\qquad M_I=-I,-I+1,\ldots,I.
$$

常见的 $K$ 取值为：

- 偶偶核基态转动带：$K=0$；
- 奇质量核：$K=\Omega$，由未配对核子决定；
- 奇奇核：$K=|\Omega_p\pm\Omega_n|$。

其中 $\Omega$、$\Omega_p$ 和 $\Omega_n$ 分别表示未配对核子、质子和中子的内禀角动量投影。Wigner-$D$ 函数的定义与约定见 [SO(3) and SU(2)](../group_theory/so3_su2.md#wigner-d-functions)。

### 内禀系到实验室系的几何变换

$\mathbf r$ 与 $\mathbf r'$ 表示同一相对坐标在实验室系与内禀系中的分量。

设 $\hat{\mathbf n}_D$ 为实验室系中的子核对称轴，$\theta'$ 为内禀系中 $\hat{\mathbf r}'$ 与对称轴 $\hat{\mathbf z}'$ 的夹角。坐标变换不改变物理夹角，因此

$$
\boxed{\cos\theta'=\hat{\mathbf r}'\cdot\hat{\mathbf z}'=\hat{\mathbf r}\cdot\hat{\mathbf n}_D}.
$$

轴对称性使内禀系展开只保留 $\mu=0$：

$$
\hat V'(\mathbf r')=\sum_\lambda\hat V'_{\lambda0}(r)Y_{\lambda0}(\hat{\mathbf r}')=\sum_\lambda\hat V'_{\lambda0}(r)\sqrt{\frac{2\lambda+1}{4\pi}}P_\lambda(\cos\theta').
$$

利用球谐函数加法定理，

$$
P_\lambda(\hat{\mathbf r}\cdot\hat{\mathbf n}_D)=\frac{4\pi}{2\lambda+1}\sum_\mu Y_{\lambda\mu}^*(\hat{\mathbf r})Y_{\lambda\mu}(\hat{\mathbf n}_D).
$$

将夹角关系与加法定理代入内禀系展开，

$$
\hat V'(\mathbf r')=\sum_{\lambda,\mu}\left[\sqrt{\frac{4\pi}{2\lambda+1}}\hat V'_{\lambda0}(r)Y_{\lambda\mu}(\hat{\mathbf n}_D)\right]Y_{\lambda\mu}^*(\hat{\mathbf r}).
$$

同一个物理相互作用在两套坐标系中的多极展开满足

$$
\boxed{\hat V(\mathbf r)=\hat V'(\mathbf r')=\sum_{\lambda,\mu}\hat V_{\lambda\mu}(r)Y_{\lambda\mu}^*(\hat{\mathbf r})=\sum_\lambda\hat V'_{\lambda0}(r)Y_{\lambda0}(\hat{\mathbf r}')}.
$$

比较实验室系展开与上式中的球谐函数系数，得到

$$
\boxed{\hat V_{\lambda\mu}(r)=\sqrt{\frac{4\pi}{2\lambda+1}}\hat V'_{\lambda0}(r)Y_{\lambda\mu}(\hat{\mathbf n}_D)}.
$$

最后利用

$$
Y_{\lambda\mu}(\hat{\mathbf n}_D)=\sqrt{\frac{2\lambda+1}{4\pi}}D_{\mu0}^{\lambda *}(\Omega_D),
$$

得到

$$
\hat V_{\lambda\mu}(r)=\hat V'_{\lambda0}(r)D_{\mu0}^{\lambda *}(\Omega_D).
$$

Wigner-$D$ 形式由夹角不变与球谐函数加法定理导出，而不是作为变换起点。

### 转动态的相互作用矩阵元

记

$$
\hat I\equiv\sqrt{2I+1},\qquad \hat I'\equiv\sqrt{2I'+1},
$$

并定义内禀矩阵元

$$
V'_{\lambda0;\kappa\kappa'}(r)\equiv\left\langle\chi_\kappa\Phi_\alpha\right|\hat V'_{\lambda0}(r)\left|\chi_{\kappa'}\Phi_\alpha\right\rangle.
$$

对于 $K,K'>0$，将转动态与

$$
\hat V_{\lambda\mu}(r)=\hat V'_{\lambda0}(r)D_{\mu0}^{\lambda *}(\Omega_D)
$$

代入矩阵元，得到

$$
\begin{aligned}
&\left\langle\Phi_{IM_IK}\Phi_\alpha\right|\hat V_{\lambda\mu}(r)\left|\Phi_{I'M_I'K'}\Phi_\alpha\right\rangle\\
={}&\frac{\hat I\hat I'}{16\pi^2}\left[\int d\Omega_D\,D^I_{M_IK}D^{\lambda *}_{\mu0}D^{I'*}_{M_I'K'}\right]V'_{\lambda0;K,K'}(r)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I'+K'}\left[\int d\Omega_D\,D^I_{M_IK}D^{\lambda *}_{\mu0}D^{I'*}_{M_I',-K'}\right]V'_{\lambda0;K,-K'}(r)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I+K}\left[\int d\Omega_D\,D^I_{M_I,-K}D^{\lambda *}_{\mu0}D^{I'*}_{M_I'K'}\right]V'_{\lambda0;-K,K'}(r)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I+K+I'+K'}\left[\int d\Omega_D\,D^I_{M_I,-K}D^{\lambda *}_{\mu0}D^{I'*}_{M_I',-K'}\right]V'_{\lambda0;-K,-K'}(r).
\end{aligned}
$$

利用 Wigner-$D$ 函数的复共轭关系

$$
D_{mk}^{j*}(\Omega_D)=(-1)^{m-k}D_{-m,-k}^{j}(\Omega_D),
$$

上式化为

$$
\begin{aligned}
&\left\langle\Phi_{IM_IK}\Phi_\alpha\right|\hat V_{\lambda\mu}(r)\left|\Phi_{I'M_I'K'}\Phi_\alpha\right\rangle\\
={}&\frac{\hat I\hat I'}{16\pi^2}(-1)^{\mu+M_I'-K'}\left[\int d\Omega_D\,D^I_{M_IK}D^\lambda_{-\mu,0}D^{I'}_{-M_I',-K'}\right]V'_{\lambda0;K,K'}(r)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I'+K'}(-1)^{\mu+M_I'+K'}\left[\int d\Omega_D\,D^I_{M_IK}D^\lambda_{-\mu,0}D^{I'}_{-M_I',K'}\right]V'_{\lambda0;K,-K'}(r)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I+K}(-1)^{\mu+M_I'-K'}\left[\int d\Omega_D\,D^I_{M_I,-K}D^\lambda_{-\mu,0}D^{I'}_{-M_I',-K'}\right]V'_{\lambda0;-K,K'}(r)\\
&+\frac{\hat I\hat I'}{16\pi^2}(-1)^{I+K+I'+K'}(-1)^{\mu+M_I'+K'}\left[\int d\Omega_D\,D^I_{M_I,-K}D^\lambda_{-\mu,0}D^{I'}_{-M_I',K'}\right]V'_{\lambda0;-K,-K'}(r).
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
&\left\langle\Phi_{IM_IK}\Phi_\alpha\right|\hat V_{\lambda\mu}(r)\left|\Phi_{I'M_I'K'}\Phi_\alpha\right\rangle\\
={}&\frac{\hat I\hat I'}{2}\begin{pmatrix}I&\lambda&I'\\M_I&-\mu&-M_I'\end{pmatrix}\Bigg\{(-1)^{\mu+M_I'-K'}\begin{pmatrix}I&\lambda&I'\\K&0&-K'\end{pmatrix}V'_{\lambda0;K,K'}(r)\\
&\qquad+(-1)^{I'+K'}(-1)^{\mu+M_I'+K'}\begin{pmatrix}I&\lambda&I'\\K&0&K'\end{pmatrix}V'_{\lambda0;K,-K'}(r)\\
&\qquad+(-1)^{I+K}(-1)^{\mu+M_I'-K'}\begin{pmatrix}I&\lambda&I'\\-K&0&-K'\end{pmatrix}V'_{\lambda0;-K,K'}(r)\\
&\qquad+(-1)^{I+K+I'+K'}(-1)^{\mu+M_I'+K'}\begin{pmatrix}I&\lambda&I'\\-K&0&K'\end{pmatrix}V'_{\lambda0;-K,-K'}(r)\Bigg\}.
\end{aligned}
$$

四个体轴 Wigner-$3j$ 符号的下排之和必须为零，因此四项分别要求

$$
K=K',\qquad K=-K',\qquad K=-K',\qquad K=K'.
$$

由于 $K,K'>0$，第二项和第三项为零，并有 $K'=K$。矩阵元化为

$$
\begin{aligned}
&\left\langle\Phi_{IM_IK}\Phi_\alpha\right|\hat V_{\lambda\mu}(r)\left|\Phi_{I'M_I'K}\Phi_\alpha\right\rangle=\frac{\hat I\hat I'}{2}\begin{pmatrix}I&\lambda&I'\\M_I&-\mu&-M_I'\end{pmatrix}\Bigg[(-1)^{\mu+M_I'-K}\begin{pmatrix}I&\lambda&I'\\K&0&-K\end{pmatrix}V'_{\lambda0;K,K}(r)\\
&\qquad+(-1)^{I+I'+2K}(-1)^{\mu+M_I'+K}\begin{pmatrix}I&\lambda&I'\\-K&0&K\end{pmatrix}V'_{\lambda0;-K,-K}(r)\Bigg].
\end{aligned}
$$

第一个 Wigner-$3j$ 符号要求 $M_I-\mu-M_I'=0$。再利用反转性质

$$
\begin{pmatrix}I&\lambda&I'\\-K&0&K\end{pmatrix}=(-1)^{I+\lambda+I'}\begin{pmatrix}I&\lambda&I'\\K&0&-K\end{pmatrix},
$$

得到

$$
\left\langle\Phi_{IM_IK}\Phi_\alpha\right|\hat V_{\lambda\mu}(r)\left|\Phi_{I'M_I'K}\Phi_\alpha\right\rangle
=\frac{\hat I\hat I'}{2}(-1)^{M_I-K}\begin{pmatrix}I&\lambda&I'\\M_I&-\mu&-M_I'\end{pmatrix}\begin{pmatrix}I&\lambda&I'\\K&0&-K\end{pmatrix}\left[V'_{\lambda0;K,K}(r)+(-1)^\lambda V'_{\lambda0;-K,-K}(r)\right].
$$

若相互作用具有时间反演对称性、内禀矩阵元取实，并且仅保留偶数 $\lambda$，则

$$
V'_{\lambda0;-K,-K}(r)=V'_{\lambda0;K,K}(r),
$$

从而

$$
\left\langle\Phi_{IM_IK}\Phi_\alpha\right|\hat V_{\lambda\mu}(r)\left|\Phi_{I'M_I'K}\Phi_\alpha\right\rangle=\hat I\hat I'(-1)^{M_I-K}\begin{pmatrix}I&\lambda&I'\\M_I&-\mu&-M_I'\end{pmatrix}\begin{pmatrix}I&\lambda&I'\\K&0&-K\end{pmatrix}V'_{\lambda0;K,K}(r).
$$

等价地，

$$
=\hat I\hat I'(-1)^{I-M_I}\begin{pmatrix}I&\lambda&I'\\-M_I&\mu&M_I'\end{pmatrix}(-1)^{I-K}\begin{pmatrix}I&\lambda&I'\\-K&0&K\end{pmatrix}V'_{\lambda0;K,K}(r).
$$

因此，约化矩阵元为

$$
\mathcal V^\lambda_{II'}(r)\equiv\left\langle IK\left\|\hat V_\lambda(r)\right\|I'K\right\rangle=\hat I\hat I'(-1)^{I-K}\begin{pmatrix}I&\lambda&I'\\-K&0&K\end{pmatrix}V'_{\lambda0;K,K}(r).
$$

对于 $K=0$，约化矩阵元化为

$$
\mathcal V^\lambda_{II'}(r)\equiv\left\langle I0\left\|\hat V_\lambda(r)\right\|I'0\right\rangle=\hat I\hat I'(-1)^I\begin{pmatrix}I&\lambda&I'\\0&0&0\end{pmatrix}V'_{\lambda0;0,0}(r).
$$

### 库伦耦合势

设子核电荷密度 $\rho_D$ 归一化为

$$
\int d^3r_D'\,\rho_D(\mathbf r_D')=1.
$$

内禀系中的库伦势为

$$
V^{\prime C}(\mathbf r')=\frac{e^2Z_\alpha Z_D}{4\pi\epsilon_0}\int d^3r_D'\,\frac{\rho_D(\mathbf r_D')}{|\mathbf r'-\mathbf r_D'|}.
$$

利用

$$
\frac{1}{|\mathbf r'-\mathbf r_D'|}=\sum_{\lambda,\mu}\frac{4\pi}{2\lambda+1}\frac{r_<^\lambda}{r_>^{\lambda+1}}Y_{\lambda\mu}^*(\hat{\mathbf r}')Y_{\lambda\mu}(\hat{\mathbf r}_D'),
$$

对于 $r\geq r_D$，定义

$$
V^{\prime C}(\mathbf r')=\frac{eZ_\alpha}{4\pi\epsilon_0}\sum_{\lambda,\mu}\frac{4\pi}{2\lambda+1}\frac{Q_{\lambda\mu}^D}{r^{\lambda+1}}Y_{\lambda\mu}^*(\hat{\mathbf r}')\qquad Q_{\lambda\mu}^D=eZ_D\int d^3r_D'\,\rho_D(\mathbf r_D')r_D^\lambda Y_{\lambda\mu}(\hat{\mathbf r}_D').
$$

因此

$$
V_{\lambda0}^{\prime C}(r)=\frac{eZ_\alpha}{4\pi\epsilon_0}\frac{4\pi}{2\lambda+1}\frac{Q_{\lambda0}^D}{r^{\lambda+1}}.
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
V_{\lambda0;K,K'}^{\prime C}(r)=\delta_{KK'}\frac{e^2Z_DZ_\alpha}{4\pi\epsilon_0}\frac{3}{2\lambda+1}\frac{(R_0^D)^\lambda}{r^{\lambda+1}}\beta_{\lambda0}^D.
$$

单极项为

$$
V_{00;K,K'}^{\prime C}(r)=\delta_{KK'}\frac{e^2Z_DZ_\alpha}{4\pi\epsilon_0}\frac{\sqrt{4\pi}}{r}.
$$

将 $V_{\lambda0;K,K'}^{\prime C}(r)$ 代入前述 $\mathcal V_{II'}^\lambda(r)$，即可得到库伦耦合矩阵元。

### 核势耦合

采用 Woods–Saxon 形式描述 $\alpha$ 粒子与子核之间的核势：

$$
\hat V^{\prime N}(\mathbf r')=\frac{V_0}{1+\exp\!\left([r-R_0^\alpha-\hat R_D(\hat{\mathbf r}',\xi_D)]/a\right)},
$$

其中 $V_0$ 为势深，$a$ 为表面弥散参数，$R_0^\alpha$ 为 $\alpha$ 粒子半径。全部角向依赖由子核半径算符给出：

$$
\hat R_D(\hat{\mathbf r}',\xi_D)=R_0^D\left[1+\sum_{\lambda,\mu}\hat a_{\lambda\mu}(\xi_D)Y_{\lambda\mu}^*(\hat{\mathbf r}')\right].
$$

仿照前述转动态约化矩阵元，

$$
\left\langle\Phi_{IK}\left\|\hat R_{D,\lambda}\right\|\Phi_{I'K'}\right\rangle=\hat I\hat I'(-1)^{I-K}\begin{pmatrix}I&\lambda&I'\\-K&0&K'\end{pmatrix}\left\langle\chi_K\right|\hat R'_{D,\lambda0}\left|\chi_{K'}\right\rangle.
$$

内禀半径算符的矩阵元为

$$
\left\langle\chi_K\right|\hat R'_{D,\lambda0}\left|\chi_{K'}\right\rangle=\delta_{KK'}R_0^D\left[\sqrt{4\pi}\,\delta_{\lambda0}+\beta_{\lambda0}^D\right].
$$

定义球形 Woods–Saxon 势

$$
V_0^N(r)\equiv\frac{V_0}{1+\exp\!\left([r-R_0^\alpha-R_0^D]/a\right)}.
$$

在内禀系中以 $\beta_{\lambda0}^D$ 代替形变算符，并将核势展开到一阶：

$$
V^{\prime N}(\mathbf r')\simeq V_0^N(r)-R_0^D\frac{dV_0^N(r)}{dr}\sum_{\lambda>0}\beta_{\lambda0}^D Y_{\lambda0}(\hat{\mathbf r}').
$$

因此，核势的内禀多极分量为

$$
V_{\lambda0;K,K'}^{\prime N}(r)=\delta_{KK'}\begin{cases}\sqrt{4\pi}\,V_0^N(r),&\lambda=0,\\-R_0^D\beta_{\lambda0}^D\dfrac{dV_0^N(r)}{dr},&\lambda>0.\end{cases}
$$

将 $V_{\lambda0;K,K'}^{\prime N}(r)$ 代入前述 $\mathcal V_{II'}^\lambda(r)$，即可得到核势耦合矩阵元。对于较大形变，可先构造并对角化半径矩阵，再计算 Woods–Saxon 势，而不作一阶展开。
