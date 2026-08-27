# 耦合道方法研究 \(\alpha\)-衰变

### 1. 背景与研究问题

##### 1. 重核和超重核

> - 合成新核素，探索原子核质量存在极限是超重核研究的重要科学问题。
> - 其中的关键在于是否存在相对稳定的或长寿命的超重核素，这主要受两个因素影响：
> - 一方面：在超重核区，库伦排斥会显著增强其裂变趋势，从而降低超重核稳定性；
> - 另一方面：壳结构效应又能够提高原子核的稳定性，使其具有相对较长的寿命。
> - 目前，发现的超重核大多是不稳定的，而alpha衰变是其重要衰变模式之一。

![超重核壳修正能分布](assets/10.1088_0034-4885_78_3_036301.png)

![已知最重核素及其衰变模式](assets/10.1088_1402-4896_aa53c1.png)

##### 2. \(\alpha\)-衰变

> - 在理论上，\(\alpha\)-衰变通常描述为：α 团簇首先在母核内部预形成，随后通过量子隧穿穿过 α–子核势垒并发射到核外。
> - 在实验上，通常通过从未知核素逐级衰变到已知核素的级联 \(\alpha\)-衰变链鉴别新合成核素；因此，\(\alpha\)-衰变链构成探索未知重核区的重要实验观测入口。

![α-衰变势垒与量子隧穿示意图](assets/10.1016_j.ppnp.2018.11.003.png)

![超重核半衰期与衰变模式分布](assets/10.1103_RevModPhys.72.3.png)

##### 3. 现有研究方法与耦合道动机

> - 目前，\(\alpha\)-衰变理论研究大体可分为宏观系统学研究与微观研究两类。
> - 宏观系统学方法，如 VSF 和 UDL，将已有实验规律参数化为经验或半经验公式，通常具有较好的整体预测精度；但核结构信息往往被平均化或吸收进参数，难以揭示内部结构与衰变机制。
> - \(\alpha\)-衰变的微观研究主要涉及单粒子高激发态、配对关联和四核子团簇化等物理机制，但目前仍难以精确计算 α 团簇形成振幅。
> - 任中洲、倪冬冬等采用耦合道方法研究 \(\alpha\)-衰变精细结构与分支比；我们的工作首先复现该方法，再在此基础上进行改进。

![α 团簇形成振幅与配对能隙](assets/10.1103_PhysRevLett.110.242502.png)

![α-衰变宏观系统学规律](assets/10.1103_PhysRevC.80.044326.png)

![耦合道方法计算的 α-衰变分支比](assets/10.1103_PhysRevC.86.054608.png)

### 2. 耦合道方法的形式理论

##### 1. 耦合道薛定谔方程

> 1. 母核 \(P\) 衰变为子核 \(D\) 与 α 粒子：
>
> \[
> P(J^{\pi_P})\rightarrow D(I^{\pi_I})+\alpha(0^+;l)\qquad c=\{l,I\}.
> \]
>
> 2. 定义通道 \(c=(I,l)\) 及通道基：
>
> \[
> \mathcal Y_{JM}^{(c)}
> =\big[Y_l\otimes\Phi_I\big]_{JM}\Phi_\alpha
> =\sum_{m,M_I}\langle lm;IM_I|JM\rangle
> Y_{lm}(\hat{\mathbf r})\Phi_{IM_I}(\xi_D)\Phi_\alpha(\xi_\alpha).
> \]
>
> 3. 以 \(u_c(r)\) 表示约化径向波函数，定义通道分量：
>
> \[
> \Psi_{JM}^{(c)}(\mathbf r,\xi_D,\xi_\alpha)
> =\frac{u_c(r)}{r}\mathcal Y_{JM}^{(c)}(\hat{\mathbf r},\xi_D,\xi_\alpha).
> \]
>
> 4. 总波函数为所有允许通道分量之和：
>
> \[
> \boxed{\Psi_{JM}=\sum_{c\in\mathcal C_{J\pi_P}}\Psi_{JM}^{(c)}}.
> \]
>
> 5. 体系哈密顿量包括子核、α 粒子、相对运动和相互作用：
>
> \[
> \hat H=\hat H_D+\hat H_\alpha+\hat T_r+\hat V.
> \]
>
> 6. 定义通道间的耦合势矩阵元：
>
> \[
> V_{cc'}^{JM}(r)
> =\left\langle\mathcal Y_{JM}^{(c)}
> \middle|\hat V\middle|
> \mathcal Y_{JM}^{(c')}\right\rangle_{\hat{\mathbf r},\xi_D,\xi_\alpha}.
> \]
>
> 7. 通道能量为：
>
> \[
> E_c\equiv E-E_I-E_\alpha.
> \]
>
> 8. 将薛定谔方程投影到通道 \(c\)，得到：
>
> \[
> \boxed{
> \sum_{c'}
> \left[
> -\frac{\hbar^2}{2\mu}
> \left(
> \frac{d^2}{dr^2}-\frac{l_c(l_c+1)}{r^2}
> \right)\delta_{cc'}
> +V_{cc'}^{JM}(r)
> \right]u_{c'}(r)
> =E_cu_c(r)
> }.
> \]

##### 2. 耦合矩阵元的一般形式

> 1. 以 \(\lambda\) 和 \(\mu\) 表示多极阶数及其分量，相互作用展开为：
>
> \[
> \hat V
> =\sum_{\lambda,\mu}\hat V_{\lambda\mu}(r)Y_{\lambda\mu}^{*}(\hat{\mathbf r})
> =\sum_{\lambda,\mu}(-1)^\mu\hat V_{\lambda\mu}(r)Y_{\lambda,-\mu}(\hat{\mathbf r}).
> \]
>
> 2. 对于 \(c=(I,l)\) 和 \(c'=(I',l')\)，耦合势矩阵元为：
>
> \[
> V_{cc'}^{JM}(r)
> =\left\langle
> \big[Y_l\otimes\Phi_I\big]_{JM}\Phi_\alpha
> \right|
> \sum_{\lambda,\mu}(-1)^\mu\hat V_{\lambda\mu}(r)Y_{\lambda,-\mu}
> \left|
> \big[Y_{l'}\otimes\Phi_{I'}\big]_{JM}\Phi_\alpha
> \right\rangle.
> \]
>
> 3. 展开角动量耦合：
>
> \[
> V_{cc'}^{JM}(r)
> =\sum_{\lambda,\mu}
> \sum_{m,M_I}
> \sum_{m',M_I'}
> (-1)^\mu
> \langle JM|lm;IM_I\rangle
> \langle l'm';I'M_I'|JM\rangle
> \langle lm|Y_{\lambda,-\mu}|l'm'\rangle
> \left\langle
> \Phi_{IM_I}\Phi_\alpha
> \right|
> \hat V_{\lambda\mu}(r)
> \left|
> \Phi_{I'M_I'}\Phi_\alpha
> \right\rangle.
> \]
>
> 4. 定义内部约化矩阵元，并应用 Wigner–Eckart 定理：
>
> \[
> \mathcal V_{II'}^\lambda(r)
> \equiv
> \left\langle
> \Phi_I\Phi_\alpha
> \right\|
> \hat V_\lambda(r)
> \left\|
> \Phi_{I'}\Phi_\alpha
> \right\rangle,
> \]
>
> \[
> \left\langle
> \Phi_{IM_I}\Phi_\alpha
> \right|
> \hat V_{\lambda\mu}(r)
> \left|
> \Phi_{I'M_I'}\Phi_\alpha
> \right\rangle
> =
> (-1)^{I-M_I}
> \begin{pmatrix}
> I&\lambda&I'\\
> -M_I&\mu&M_I'
> \end{pmatrix}
> \mathcal V_{II'}^\lambda(r).
> \]
>
> 5. 球谐函数矩阵元为：
>
> \[
> \langle lm|Y_{\lambda,-\mu}|l'm'\rangle
> =
> (-1)^{l-m}
> \begin{pmatrix}
> l&\lambda&l'\\
> -m&-\mu&m'
> \end{pmatrix}
> \langle l\|Y_\lambda\|l'\rangle,
> \]
>
> \[
> \langle l\|Y_\lambda\|l'\rangle
> =
> (-1)^l
> \sqrt{\frac{(2l+1)(2\lambda+1)(2l'+1)}{4\pi}}
> \begin{pmatrix}
> l&\lambda&l'\\
> 0&0&0
> \end{pmatrix}.
> \]
>
> 6. 代回并完成磁量子数求和，得到：
>
> \[
> \boxed{
> V_{cc'}^{JM}(r)
> =
> \sum_\lambda
> (-1)^{l'+I+J}
> \begin{Bmatrix}
> l&I&J\\
> I'&l'&\lambda
> \end{Bmatrix}
> \langle l\|Y_\lambda\|l'\rangle
> \mathcal V_{II'}^\lambda(r)
> }.
> \]

##### 3. 从内禀系到实验室系的变换

> 1. \(\mathbf r\) 与 \(\mathbf r'\) 表示同一相对坐标在实验室系和内禀系中的分量。坐标变换保持物理夹角不变：
>
> \[
> \boxed{
> \cos\theta'
> =\hat{\mathbf r}'\cdot\hat{\mathbf z}'
> =\hat{\mathbf r}\cdot\hat{\mathbf n}_D
> }.
> \]
>
> 2. 轴对称性使内禀系的多极展开仅保留 \(\mu=0\)：
>
> \[
> \hat V'(\mathbf r')
> =\sum_\lambda\hat V'_{\lambda0}(r)Y_{\lambda0}(\hat{\mathbf r}')
> =\sum_\lambda\hat V'_{\lambda0}(r)
> \sqrt{\frac{2\lambda+1}{4\pi}}
> P_\lambda(\cos\theta').
> \]
>
> 3. 球谐函数加法定理给出：
>
> \[
> P_\lambda(\hat{\mathbf r}\cdot\hat{\mathbf n}_D)
> =\frac{4\pi}{2\lambda+1}
> \sum_\mu
> Y_{\lambda\mu}^{*}(\hat{\mathbf r})
> Y_{\lambda\mu}(\hat{\mathbf n}_D).
> \]
>
> 4. 将夹角关系与加法定理代入内禀系展开：
>
> \[
> \hat V'(\mathbf r')
> =\sum_{\lambda,\mu}
> \left[
> \sqrt{\frac{4\pi}{2\lambda+1}}
> \hat V'_{\lambda0}(r)
> Y_{\lambda\mu}(\hat{\mathbf n}_D)
> \right]
> Y_{\lambda\mu}^{*}(\hat{\mathbf r}).
> \]
>
> 5. 同一相互作用在两套坐标系中的多极展开满足：
>
> \[
> \boxed{
> \hat V(\mathbf r)
> =\hat V'(\mathbf r')
> =\sum_{\lambda,\mu}
> \hat V_{\lambda\mu}(r)
> Y_{\lambda\mu}^{*}(\hat{\mathbf r})
> =\sum_\lambda
> \hat V'_{\lambda0}(r)
> Y_{\lambda0}(\hat{\mathbf r}')
> }.
> \]
>
> 6. 比较球谐函数系数，并利用球谐函数与 Wigner-\(D\) 函数的关系：
>
> \[
> \hat V_{\lambda\mu}(r)
> =\sqrt{\frac{4\pi}{2\lambda+1}}
> \hat V'_{\lambda0}(r)
> Y_{\lambda\mu}(\hat{\mathbf n}_D),
> \]
>
> \[
> Y_{\lambda\mu}(\hat{\mathbf n}_D)
> =\sqrt{\frac{2\lambda+1}{4\pi}}
> D_{\mu0}^{\lambda *}(\Omega_D),
> \]
>
> \[
> \boxed{
> \hat V_{\lambda\mu}(r)
> =\hat V'_{\lambda0}(r)
> D_{\mu0}^{\lambda *}(\Omega_D)
> }.
> \]

##### 4. 内禀系约化矩阵元

> 1. 对于 \(K>0\)，子核转动态为：
>
> \[
> \boxed{
> \Phi_{IM_IK}(\Omega_D,\eta_D)
> =
> \sqrt{\frac{2I+1}{16\pi^2}}
> \left[
> D_{M_IK}^{I*}(\Omega_D)\chi_K(\eta_D)
> +(-1)^{I+K}
> D_{M_I,-K}^{I*}(\Omega_D)\chi_{-K}(\eta_D)
> \right]
> }.
> \]
>
> 2. 记 \(\hat I=\sqrt{2I+1}\)、\(\hat I'=\sqrt{2I'+1}\)，转动态相互作用矩阵元可写为：
>
> \[
> \left\langle
> \Phi_{IM_IK}\Phi_\alpha
> \right|
> \hat V_{\lambda\mu}(r)
> \left|
> \Phi_{I'M_I'K}\Phi_\alpha
> \right\rangle
> =
> \frac{\hat I\hat I'}{2}
> (-1)^{M_I-K}
> \begin{pmatrix}
> I&\lambda&I'\\
> M_I&-\mu&-M_I'
> \end{pmatrix}
> \begin{pmatrix}
> I&\lambda&I'\\
> K&0&-K
> \end{pmatrix}
> \left[
> V'_{\lambda0;K,K}(r)
> +(-1)^\lambda V'_{\lambda0;-K,-K}(r)
> \right].
> \]
>
> 3. 若相互作用具有时间反演对称性、内禀矩阵元为实数，并且仅保留偶数 \(\lambda\)，则：
>
> \[
> V'_{\lambda0;-K,-K}(r)
> =
> V'_{\lambda0;K,K}(r),
> \]
>
> \[
> \left\langle
> \Phi_{IM_IK}\Phi_\alpha
> \right|
> \hat V_{\lambda\mu}(r)
> \left|
> \Phi_{I'M_I'K}\Phi_\alpha
> \right\rangle
> =
> \hat I\hat I'
> (-1)^{I-M_I}
> \begin{pmatrix}
> I&\lambda&I'\\
> -M_I&\mu&M_I'
> \end{pmatrix}
> (-1)^{I-K}
> \begin{pmatrix}
> I&\lambda&I'\\
> -K&0&K
> \end{pmatrix}
> V'_{\lambda0;K,K}(r).
> \]
>
> 4. 与 Wigner–Eckart 定理比较，分离出约化矩阵元：
>
> \[
> \boxed{
> \mathcal V_{II'}^\lambda(r)
> \equiv
> \left\langle
> IK
> \left\|
> \hat V_\lambda(r)
> \right\|
> I'K
> \right\rangle
> =
> \hat I\hat I'
> (-1)^{I-K}
> \begin{pmatrix}
> I&\lambda&I'\\
> -K&0&K
> \end{pmatrix}
> V'_{\lambda0;K,K}(r)
> }.
> \]

### 3. 耦合模型与阶段性结果

##### 1. 库仑耦合

> 1. 设子核电荷密度 \(\rho_D\) 归一化为 \(1\)，内禀系中的库仑势为：
>
> \[
> \int d^3r_D'\,\rho_D(\mathbf r_D')=1,
> \]
>
> \[
> V^{\prime C}(\mathbf r')
> =
> \frac{e^2Z_\alpha Z_D}{4\pi\epsilon_0}
> \int d^3r_D'\,
> \frac{\rho_D(\mathbf r_D')}
> {|\mathbf r'-\mathbf r_D'|}.
> \]
>
> 2. 库仑核的多极展开为：
>
> \[
> \frac{1}{|\mathbf r'-\mathbf r_D'|}
> =
> \sum_{\lambda,\mu}
> \frac{4\pi}{2\lambda+1}
> \frac{r_<^\lambda}{r_>^{\lambda+1}}
> Y_{\lambda\mu}^{*}(\hat{\mathbf r}')
> Y_{\lambda\mu}(\hat{\mathbf r}_D').
> \]
>
> 3. 对于核外区域，定义子核的电多极矩：
>
> \[
> Q_{\lambda\mu}^D
> =
> eZ_D
> \int d^3r_D'\,
> \rho_D(\mathbf r_D')
> r_D^\lambda
> Y_{\lambda\mu}(\hat{\mathbf r}_D').
> \]
>
> 4. 将多极展开与电多极矩代入库仑势：
>
> \[
> V^{\prime C}(\mathbf r')
> =
> \frac{eZ_\alpha}{4\pi\epsilon_0}
> \sum_{\lambda,\mu}
> \frac{4\pi}{2\lambda+1}
> \frac{Q_{\lambda\mu}^D}{r^{\lambda+1}}
> Y_{\lambda\mu}^{*}(\hat{\mathbf r}').
> \]
>
> 5. 对于轴对称子核，仅保留内禀系中的 \(\mu=0\) 分量：
>
> \[
> \boxed{
> V_{\lambda0}^{\prime C}(r)
> =
> \frac{eZ_\alpha}{4\pi\epsilon_0}
> \frac{4\pi}{2\lambda+1}
> \frac{Q_{\lambda0}^D}{r^{\lambda+1}}
> }.
> \]
>
> 6. 对于均匀带电形变球体，保留形变参数的一阶项：
>
> \[
> \boxed{
> Q_{\lambda0}^D
> \simeq
> eZ_D\frac{3}{4\pi}(R_0^D)^\lambda
> \left[
> \frac{\sqrt{4\pi}}{3}\delta_{\lambda0}
> +\beta_{\lambda0}^D
> \right]
> }.
> \]

##### 2. 核势耦合

> 1. 子核表面形变由半径算符描述：
>
> \[
> \hat R_D(\hat{\mathbf r}',\xi_D)
> =
> R_0^D
> \left[
> 1+
> \sum_{\lambda,\mu}
> \hat a_{\lambda\mu}(\xi_D)
> Y_{\lambda\mu}^{*}(\hat{\mathbf r}')
> \right].
> \]
>
> 2. 转动态之间的半径算符约化矩阵元为：
>
> \[
> \left\langle
> \Phi_{IK}
> \left\|
> \hat R_{D,\lambda}
> \right\|
> \Phi_{I'K'}
> \right\rangle
> =
> \hat I\hat I'
> (-1)^{I-K}
> \begin{pmatrix}
> I&\lambda&I'\\
> -K&0&K'
> \end{pmatrix}
> \left\langle
> \chi_K
> \right|
> \hat R'_{D,\lambda0}
> \left|
> \chi_{K'}
> \right\rangle.
> \]
>
> 3. 内禀半径算符的矩阵元为：
>
> \[
> \left\langle
> \chi_K
> \right|
> \hat R'_{D,\lambda0}
> \left|
> \chi_{K'}
> \right\rangle
> =
> \delta_{KK'}R_0^D
> \left[
> \sqrt{4\pi}\,\delta_{\lambda0}
> +\beta_{\lambda0}^D
> \right].
> \]
>
> 4. 将半径算符代入 Woods–Saxon 形式，得到 α 粒子与子核之间的核势：
>
> \[
> \boxed{
> \hat V^{\prime N}(\mathbf r')
> =
> \frac{V_0}
> {1+\exp\!\left(
> [r-R_0^\alpha-\hat R_D(\hat{\mathbf r}',\xi_D)]/a
> \right)}
> }.
> \]
>
> 其中，\(V_0\) 为势深，\(a\) 为表面弥散参数，\(R_0^\alpha\) 为 α 粒子半径。

##### 3. 数值计算与阶段性结果

> 1. 将库仑耦合与核势耦合矩阵元代入多通道薛定谔方程，数值求解各通道径向波函数。
>
> 2. 采用打靶法匹配内、外区波函数的对数导数，并通过奇异值分解（SVD）搜索匹配矩阵的最小奇异值：
>
> \[
> \sigma_{\min}\!\left[\mathbf M(V_0)\right]\rightarrow 0.
> \]
>
> 由此确定核势深度 \(V_0\)，并由对应的右奇异向量构造匹配波函数 \(\mathbf u_{\mathrm{match}}\)。
>
> 3. 对 \(\mathbf u_{\mathrm{match}}\) 进行归一化，随后计算各通道的衰变宽度与分支比：
>
> \[
> \Gamma_{\mathrm{ch}}
> = P_\alpha \hbar j_{\rm ch}(r_\max)
> =
> P_\alpha\hbar v_{\mathrm{ch}}
> \left|
> u_{\mathrm{ch}}(r_{\max})
> \right|^2,
> \]
>
> \[
> \mathrm{BR}_{\mathrm{ch}}
> =
> \frac{\Gamma_{\mathrm{ch}}}
> {\displaystyle\sum_{\mathrm{ch}'}\Gamma_{\mathrm{ch}'}}.
> \]
>
> 其中，\(P_\alpha\) 为 α 粒子预形成概率，\(v_{\mathrm{ch}}\) 为对应通道的渐近速度。

### 4. 总结与展望

##### 1. 阶段性总结

> - 学习 \(\alpha\)-衰变的耦合道形式理论与从内禀系到实验室系的变换。
> - 构造库仑耦合与核势耦合矩阵元，并代入多通道径向薛定谔方程。
> - 初步复现数值求解流程，获得 \(\alpha\)-衰变分支比与半衰期计算结果。

##### 2. 下一步工作

> - 结合投影壳模型，引入角动量投影与 \(K\) 混合，改进当前耦合道方法中固定 \(M\)、\(K\) 的近似处理；
> - 将现有理论框架推广至奇 \(A\) 核和奇奇核，发展相应的耦合道计算方法。
