# 代码文件 `hfbtho_gauss.f90` 说明

本文档对模块 `HFBTHO_Gauss`（1113 行）中的全部子程序与函数进行说明。该模块为 HFBTHO 提供完整的高斯求积（Gauss quadrature）功能，包括柱坐标谐振子（HO）基函数及其导数的计算、三种正交多项式积分网格的构建、以及 Hermite 多项式乘积展开系数。

---

### `gaupol`

##### 函数定义

- 函数位置: `193-352 lines`

```fortran
Subroutine gaupol(lpr)
```

##### 物理公式

- 物理含义: 在柱坐标 $(\rho,\theta,z)$ 下计算 HO 波函数 $\psi_{n_r\Lambda n_z}(\rho,\theta,z)$ 及其一阶导数在 Gauss 积分节点上的值，并检验正交归一性精度。

完整 HO 波函数为
$$
\psi_{n_r\Lambda n_z}(\rho,\theta,z)
= \psi_{n_r}^{|\Lambda|}(\eta)\,\psi_{n_z}(\xi)\,
\frac{e^{i\Lambda\theta}}{\sqrt{2\pi}}\,\chi_{\Sigma}(\sigma),
$$
其中无量纲坐标 $\eta=\beta_{\perp}^{2}\rho^{2}$，$\xi=\beta_{z}z$。

**z 方向**
$$
\psi_{n_z}(z)=\mathcal{N}_{n_z}\,\beta_z^{1/2}\,e^{-\xi^{2}/2}\,H_{n_z}(\xi),
\qquad
\mathcal{N}_{n_z}=\left(\frac{1}{\sqrt{\pi}\,2^{n_z}n_z!}\right)^{1/2},
$$
$$
\frac{\partial\psi_{n_z}}{\partial z}
=\beta_z\left[\sqrt{2n_z}\,\psi_{n_z-1}(\xi)-\xi\,\psi_{n_z}(\xi)\right].
$$

**$\rho$ 方向**
$$
\psi_{n_r}^{\Lambda}(\eta)
=\mathcal{N}_{n_r}^{\Lambda}\,\beta_{\perp}\sqrt{2}\,
\eta^{|\Lambda|/2}e^{-\eta/2}L_{n_r}^{|\Lambda|}(\eta),
\qquad
\mathcal{N}_{n_r}^{\Lambda}
=\left(\frac{n_r!}{(n_r+|\Lambda|)!}\right)^{1/2},
$$
$$
\frac{\partial\psi_{n_r}^{\Lambda}}{\partial\rho}
=\frac{\beta_{\perp}}{\sqrt{\eta}}
\Bigl[(2n_r+|\Lambda|-\eta)\psi_{n_r}^{\Lambda}(\eta)
-2\sqrt{n_r(n_r+|\Lambda|)}\,\psi_{n_r-1}^{\Lambda}(\eta)\Bigr].
$$

- 实现步骤:
  1. 利用 `xh(ih)`、`wh(ih)` 和递推关系计算 Hermite 多项式相关的 `qh(n,ih)`（$\psi_{n_z}/\beta_z^{1/2}$）和 `qh1(n,ih)`（$\partial_z\psi_{n_z}/\beta_z^{3/2}$）。
  2. 利用 `xl(il)`、`wl(il)` 和递推关系计算 Laguerre 多项式相关的 `ql(n,l,il)`（$\psi_{n_r}^{\Lambda}/(\beta_{\perp}\sqrt{2})$）和 `ql1(n,l,il)`（$\sqrt{\eta}\,\partial_\rho\psi_{n_r}^{\Lambda}/(\beta_{\perp}^{2}\sqrt{2})$）。
  3. 对 `qh` 和 `ql` 分别做正交归一性数值检验，输出最大非正交偏差与归一化偏差。

##### 变量对应

- 输入变量1: `lpr (Logical)`，控制是否打印网格信息与精度测试结果
- 输入变量2: `xh(ngh)`，Gauss-Hermite 无量纲节点 $\xi_h$
- 输入变量3: `wh(ngh)`，已含权函数的高斯-厄米权重 $w_h$
- 输入变量4: `xl(ngl)`，Gauss-Laguerre 无量纲节点 $\eta_l$
- 输入变量5: `wl(ngl)`，已含权函数的高斯-拉盖尔权重 $w_l$
- 输出变量1: `qh(0:nzm,1:ngh)`，z 方向归一化 HO 波函数值（无量纲）
- 输出变量2: `qh1(0:nzm,1:ngh)`，z 方向 HO 波函数一阶导数（无量纲）
- 输出变量3: `ql(0:nrm,0:nlm,1:ngl)`，$\rho$ 方向归一化 HO 波函数值（无量纲）
- 输出变量4: `ql1(0:nrm,0:nlm,1:ngl)`，$\rho$ 方向 HO 波函数一阶导数（无量纲）

##### 函数数据流
- 调用情况: `preparer <- gaupol`
- 被 `preparer`（`hfbtho_solver.f90`）在基组构建完成后调用，用于填充全局数组 `qh`、`qh1`、`ql`、`ql1`

---

### `coordinateLST_HO`

##### 函数定义

- 函数位置: `357-375 lines`

```fortran
Subroutine coordinateLST_HO()
```

##### 物理公式

- 物理含义: 将无量纲 Gauss 节点转换为有量纲柱坐标 $(z_h, \rho_l)$，并计算二维积分体积元权重 $\mathrm{d}v = \pi\,b_z\,b_{\perp}^{2}\,w_h\,w_l$。

网格点映射关系：
$$
z_h = b_z\,\xi_h, \qquad \rho_l = b_{\perp}\sqrt{\eta_l},
$$
体积元权重：
$$
w_{hl} = \pi\,w_h\,w_l\,b_z\,b_{\perp}^{2}.
$$

##### 变量对应

- 输入变量1: `xh(ngh)`，Hermite 无量纲节点 $\xi_h$
- 输入变量2: `xl(ngl)`，Laguerre 无量纲节点 $\eta_l$
- 输入变量3: `wh(ngh)`，Hermite 权重 $w_h$
- 输入变量4: `wl(ngl)`，Laguerre 权重 $w_l$
- 输入变量5: `bz (Real(pr))`，z 方向谐振子长度 $b_z$
- 输入变量6: `bp (Real(pr))$，$\rho$ 方向谐振子长度 $b_{\perp}$
- 输出变量1: `fh(nghl)`，有量纲 z 节点 $z_h$
- 输出变量2: `fl(nghl)`，有量纲 $\rho$ 节点 $\rho_l$
- 输出变量3: `wdcor(nghl)`，体积元权重 $w_{hl}$
- 输出变量4: `wdcori(nghl)`，体积元权重倒数 $1/w_{hl}$

##### 函数数据流
- 调用情况: `preparer <- coordinateLST_HO <- optHFBTHO_HO`
- 被 `preparer`（`hfbtho_solver.f90`）调用；内部自动调用 `optHFBTHO_HO` 预计算优化数组

---

### `optHFBTHO_HO`

##### 函数定义

- 函数位置: `382-423 lines`

```fortran
Subroutine optHFBTHO_HO()
```

##### 物理公式

- 物理含义: 预计算 HO 基函数值 $f_i$、径向导数 $\nabla_r f_i$、z 向导数 $\nabla_z f_i$ 以及 Laplacian $\Delta f_i$ 在所有 Gauss 网格点上的值，供密度与场计算时直接查表，避免重复计算。

对于第 $i$ 个基函数（量子数 $n_z,n_r,\Lambda,\Sigma$），在网格点 $(\xi_h,\eta_l)$ 上有：
$$
f_i = \psi_{n_z}(\xi_h)\,\psi_{n_r}^{\Lambda}(\eta_l),
\qquad
\nabla_r f_i = \frac{2\sqrt{\eta_l}}{b_{\perp}}\,\psi_{n_z}\,\partial_\eta\psi_{n_r}^{\Lambda},
\qquad
\nabla_z f_i = \frac{1}{b_z}\,\partial_\xi\psi_{n_z}\,\psi_{n_r}^{\Lambda},
$$
$$
\Delta f_i
=\Bigl[(\xi_h^{2}-2n_z-1)\,b_z^{-2}
+4\bigl(\tfrac{1}{4}-(2n_r+|\Lambda|+1)\eta_l^{-1}+\Lambda^{2}\eta_l^{-2}\bigr)\eta_l\,b_{\perp}^{-2}\Bigr]f_i.
$$

- 实现步骤:
  1. 分配优化数组 `QHLA_opt`、`FI1R_opt`、`FI1Z_opt`、`FI2D_opt`、`y_opt`。
  2. 遍历所有 $K^{\pi}$ 分块与基函数，在每个 Gauss 网格点上利用已计算的 `qh`、`qh1`、`ql`、`ql1` 组合得到 $f_i$、$\nabla_r f_i$、$\nabla_z f_i$、$\Delta f_i$。

##### 变量对应

- 输入变量1: `qh`、`qh1`、`ql`、`ql1`，由 `gaupol` 预先计算
- 输入变量2: `bp`、`bz`，HO 长度参数
- 输入变量3: `ntx`、`NB`、`ID(ib)`、`ia(ib)`，基组分块信息
- 输入变量4: `nl(ja)`、`nr(ja)`、`nz(ja)`、`ns(ja)`，各基函数的量子数
- 输出变量1: `QHLA_opt(ntx,nghl)`，基函数值 $f_i$
- 输出变量2: `FI1R_opt(ntx,nghl)`，径向导数 $\nabla_r f_i$
- 输出变量3: `FI1Z_opt(ntx,nghl)`，z 向导数 $\nabla_z f_i$
- 输出变量4: `FI2D_opt(ntx,nghl)`，Laplacian $\Delta f_i$
- 输出变量5: `y_opt(nghl)`，$1/\rho_l$（避免除零的预计算量）

##### 函数数据流
- 调用情况: `coordinateLST_HO <- optHFBTHO_HO`
- 仅由 `coordinateLST_HO` 内部调用

---

### `gausspoints`

##### 函数定义

- 函数位置: `440-498 lines`

```fortran
Subroutine gausspoints()
```

##### 物理公式

- 物理含义: 分别计算 Gauss-Hermite、Gauss-Laguerre 和 Gauss-Legendre 积分的节点与权重。对于含权函数的积分，权重已并入 $w_h$、$w_l$，使得后续密度计算无需再显式写出权函数。

Gauss-Hermite（全空间积分，权函数 $e^{-\xi^{2}}$）：
$$
\int_{-\infty}^{+\infty} f(\xi)\,g(\xi)\,\mathrm{d}\xi
=\sum_{h} w_h\,f(\xi_h)\,g(\xi_h),
\qquad w_h = \omega_h\,e^{\xi_h^{2}}.
$$

Gauss-Laguerre（半无限空间，权函数 $e^{-\eta}\eta^{\alpha}$，此处 $\alpha=0$）：
$$
\int_{0}^{+\infty} f(\eta)\,g(\eta)\,\mathrm{d}\eta
=\sum_{l} w_l\,f(\eta_l)\,g(\eta_l),
\qquad w_l = \omega_l\,e^{\eta_l}.
$$

- 实现步骤:
  1. 若 `Parity=.True.`，Hermite 网格点数加倍（$N=2\,\text{ngh}$），仅保留正节点并加倍权重；否则直接取 $N=\text{ngh}$。
  2. 调用 `Gaussq` 计算 Hermite、Laguerre、Legendre 的节点 $t$ 和原始权重 $w$。
  3. 将权函数指数因子并入最终权重 `wh`、`wl`。

##### 变量对应

- 输入变量1: `ngh`、`ngl`、`nleg`，各类积分点数
- 输入变量2: `Parity (Logical)`，是否守恒宇称
- 输出变量1: `xh(ngh)`，Hermite 无量纲节点 $\xi_h$
- 输出变量2: `wh(ngh)`，含权 Hermite 权重 $w_h$
- 输出变量3: `xl(ngl)`，Laguerre 无量纲节点 $\eta_l$
- 输出变量4: `wl(ngl)`，含权 Laguerre 权重 $w_l$
- 输出变量5: `sxl(ngl)`，$\sqrt{\eta_l}$
- 输出变量6: `xleg(nleg)`，Legendre 节点 $\mu_{\text{leg}}$
- 输出变量7: `wleg(nleg)`，Legendre 权重 $w_{\text{leg}}$

##### 函数数据流
- 调用情况: `preparer <- gausspoints <- Gaussq`
- 被 `preparer`（`hfbtho_solver.f90`）调用；内部调用 `Gaussq` 完成实际计算

---

### `Gaussq`

##### 函数定义

- 函数位置: `550-584 lines`

```fortran
Subroutine Gaussq(kindi,N,alpha,beta,kpts,endpts,b,t,w)
```

##### 物理公式

- 物理含义: 通用 Gauss 求积节点与权重计算。基于 Golub–Welsch 算法：将正交多项式的递推系数组织成对称三对角 Jacobi 矩阵，其特征值即为节点，特征向量模方乘以零阶矩 $\mu_0$ 即得权重。

递推关系（正交归一多项式）：
$$
b_j\,p_j(x)=(x-a_j)\,p_{j-1}(x)-b_{j-1}\,p_{j-2}(x),
\qquad
\mu_0=\int w(x)\,\mathrm{d}x.
$$

- 实现步骤:
  1. 调用 `Class` 获取指定正交多项式族的递推系数 $a_j$、$b_j$ 和零阶矩 $\mu_0$。
  2. 若 `kpts=0`（普通 Gauss），直接对三对角矩阵调用 `GBTQL2` 求特征值/向量。
  3. 若 `kpts=1` 或 `2`（Gauss-Radau / Gauss-Lobatto），调用 `GBSLVE` 修正边界节点后，再调用 `GBTQL2`。
  4. 权重由 $w_j=\mu_0\,[v_j^{(1)}]^2$ 得到（$v^{(1)}$ 为首分量归一化特征向量）。

##### 变量对应

- 输入变量1: `kindi (Integer)`，多项式类型（1=Legendre，2=Chebyshev I，3=Chebyshev II，4=Hermite，5=Jacobi，6=Laguerre）
- 输入变量2: `N (Integer)`，求积阶数
- 输入变量3: `alpha`、`beta (Real(pr))`，Jacobi/Laguerre 参数
- 输入变量4: `kpts (Integer)`，固定端点个数（0/1/2）
- 输入变量5: `endpts(2) (Real(pr))`，固定端点值
- 输入/输出变量6: `b(N) (Real(pr))`，工作数组，也存放递推系数 $b_j$
- 输出变量1: `t(N) (Real(pr))`，求积节点 $t_j$
- 输出变量2: `w(N) (Real(pr))`，求积权重 $w_j$

##### 函数数据流
- 调用情况: `gausspoints <- Gaussq <- Class / GBSLVE / GBTQL2`
- 被 `gausspoints` 调用；内部调用 `Class`、`GBSLVE`、`GBTQL2`

---

### `GBSLVE`

##### 函数定义

- 函数位置: `588-603 lines`

```fortran
Real(pr) Function GBSLVE(SHIFT,N,A,B)
```

##### 物理公式

- 物理含义: 在 Golub–Welsch 算法中，给定移位量 `SHIFT`（固定端点），计算修正三对角矩阵所需的比例因子。通过对递推关系的连分式求值实现：
$$
\alpha_1 = a_1 - \text{SHIFT},
\qquad
\alpha_i = a_i - \text{SHIFT} - \frac{b_{i-1}^{2}}{\alpha_{i-1}}
\;(i=2,\dots,N-1),
$$
返回 $1/\alpha_{N-1}$。

##### 变量对应

- 输入变量1: `SHIFT (Real(pr))`，待固定的端点值
- 输入变量2: `N (Integer)`，矩阵维度
- 输入变量3: `A(N) (Real(pr))`，三对角矩阵对角元 $a_i$
- 输入变量4: `B(N) (Real(pr))`，三对角矩阵次对角元 $b_i$
- 输出变量1: 函数返回值，$1/\alpha_{N-1}$

##### 函数数据流
- 调用情况: `Gaussq <- GBSLVE`
- 仅由 `Gaussq` 在 Gauss-Radau / Gauss-Lobatto 情况下调用

---

### `Class`

##### 函数定义

- 函数位置: `630-699 lines`

```fortran
Subroutine Class(kindi,N,ALPHA,BETA,B,A,MUZERO)
```

##### 物理公式

- 物理含义: 为六种经典正交多项式族提供递推系数 $a_j$、$b_j$ 和零阶矩 $\mu_0$。各项结果如下表：

| `kindi` | 多项式 | $\mu_0$ | $a_j$ | $b_j$ |
|---------|--------|---------|-------|-------|
| 1 | Legendre | $2$ | $0$ | $j/\sqrt{4j^{2}-1}$ |
| 2 | Chebyshev I | $\pi$ | $0$ | $1/2$（$b_1=1/\sqrt{2}$） |
| 3 | Chebyshev II | $\pi/2$ | $0$ | $1/2$ |
| 4 | Hermite | $\sqrt{\pi}$ | $0$ | $\sqrt{j/2}$ |
| 5 | Jacobi | 见 $\Gamma$ 函数表达式 | $\frac{\beta^{2}-\alpha^{2}}{(2j+\alpha+\beta)(2j+\alpha+\beta-2)}$ | 复杂有理式 |
| 6 | Laguerre | $\Gamma(\alpha+1)$ | $2j-1+\alpha$ | $\sqrt{j(j+\alpha)}$ |

- 实现步骤:
  1. 根据 `kindi` 进入对应分支。
  2. 对 Legendre、Chebyshev、Hermite 使用解析公式；对 Jacobi 和 Laguerre 需调用 `pr_gamma` 计算 $\Gamma$ 函数。

##### 变量对应

- 输入变量1: `kindi (Integer)`，多项式类型（1–6）
- 输入变量2: `N (Integer)`，阶数
- 输入变量3: `ALPHA`、`BETA (Real(pr))`，Jacobi/Laguerre 参数
- 输出变量1: `A(N) (Real(pr))`，递推系数 $a_j$
- 输出变量2: `B(N) (Real(pr))`，递推系数 $b_j$
- 输出变量3: `MUZERO (Real(pr))`，零阶矩 $\mu_0$

##### 函数数据流
- 调用情况: `Gaussq <- Class <- pr_gamma`
- 仅由 `Gaussq` 调用；在 Jacobi/Laguerre 分支中调用 `pr_gamma`

---

### `GBTQL2`

##### 函数定义

- 函数位置: `703-786 lines`

```fortran
Subroutine GBTQL2(N,D,E,Z,IERR)
```

##### 物理公式

- 物理含义: 隐式 QL 算法（Implicit QL Algorithm）求解实对称三对角矩阵的全部特征值与特征向量。Jacobi 矩阵为
$$
J = \begin{pmatrix}
d_1 & e_1 & & \\
e_1 & d_2 & e_2 & \\
& \ddots & \ddots & \ddots \\
& & e_{N-1} & d_N
\end{pmatrix},
$$
算法通过一系列正交相似变换将 $J$ 对角化，得到特征值 $d_j$（即求积节点）和首分量归一化特征向量 $z_j$（用于计算权重）。

- 实现步骤:
  1. 对次对角元进行零判断，确定当前需处理的未收敛子块。
  2. 选取 Wilkinson 位移，执行隐式 QL 迭代。
  3. 迭代收敛后按特征值从小到大排序，并同步交换特征向量分量。

##### 变量对应

- 输入变量1: `N (Integer)`，矩阵维度
- 输入/输出变量2: `D(N) (Real(pr))`，输入为对角元，输出为排序后的特征值
- 输入/输出变量3: `E(N) (Real(pr))`，输入为次对角元（$E(N)=0$），迭代过程中被覆盖
- 输入/输出变量4: `Z(N) (Real(pr))`，输入初值为 $Z(1)=1$、其余为 $0$，输出为特征向量首分量
- 输出变量5: `IERR (Integer)`，错误标志（0=成功，非零=第 `IERR` 个特征值未收敛）

##### 函数数据流
- 调用情况: `Gaussq <- GBTQL2`
- 仅由 `Gaussq` 调用

---

### `gauleg`

##### 函数定义

- 函数位置: `801-841 lines`

```fortran
Subroutine gauleg(x1,x2,x,w,n)
```

##### 物理公式

- 物理含义: 直接利用 Newton 迭代法计算 $N$ 阶 Gauss-Legendre 求积节点与权重，积分区间为 $[x_1, x_2]$。

变换关系：
$$
x_i = x_m - x_l\,\mu_i,
\qquad
x_{N+1-i} = x_m + x_l\,\mu_i,
\qquad
w_i = w_{N+1-i} = \frac{2x_l}{(1-\mu_i^{2})\,[P_N'(\mu_i)]^{2}},
$$
其中 $x_m=(x_1+x_2)/2$，$x_l=(x_2-x_1)/2$，$\mu_i$ 为 $P_N(\mu)$ 的正零点。

- 实现步骤:
  1. 以 $\cos\bigl[\pi(i-1/4)/(N+1/2)\bigr]$ 为初值。
  2. 用 Legendre 多项式递推计算 $P_N(\mu)$ 及其导数 $P_N'(\mu)$。
  3. Newton 迭代修正零点直至收敛。

##### 变量对应

- 输入变量1: `x1`、`x2 (Real(pr))`，积分区间端点
- 输入变量2: `n (Integer)`，求积阶数
- 输出变量1: `x(n) (Real(pr))`，节点
- 输出变量2: `w(n) (Real(pr))`，权重

##### 函数数据流
- 调用情况: `coulom / recompute_coulomb_expansion <- gauleg`
- 被 `coulom` 和 `recompute_coulomb_expansion`（`hfbtho_solver.f90`）调用，用于库仑势的 Gauss-Legendre 积分

---

### `D_HERM`

##### 函数定义

- 函数位置: `846-871 lines`

```fortran
Subroutine D_HERM(X,N,HER,DHER,NDIM)
```

##### 物理公式

- 物理含义: 计算物理 Hermite 多项式 $H_n(x)$（满足 $H_0=1, H_1=2x$）及其导数 $H_n'(x)=2nH_{n-1}(x)$ 在单点 $x$ 处的值，直到 $N$ 阶。

递推关系：
$$
H_{n+1}(x)=2xH_n(x)-2nH_{n-1}(x),
\qquad
H_n'(x)=2nH_{n-1}(x).
$$

##### 变量对应

- 输入变量1: `X (Real(pr))`，自变量 $x$
- 输入变量2: `N (Integer)`，最高阶数
- 输入变量3: `NDIM (Integer)`，数组维度
- 输出变量1: `HER(1:NDIM) (Real(pr))`，$H_n(x)$
- 输出变量2: `DHER(1:NDIM) (Real(pr))`，$H_n'(x)$

##### 函数数据流
- 调用情况: `DEVHER <- D_HERM`
- 仅由 `DEVHER` 调用

---

### `DEVHER`

##### 函数定义

- 函数位置: `879-945 lines`

```fortran
Subroutine DEVHER(NOSACT)
```

##### 物理公式

- 物理含义: 将两个归一化 Hermite 多项式的乘积 $H_n(x)H_m(x)$ 展开为 Hermite 多项式的线性组合，计算展开系数；同时预计算归一化因子。该展开用于 Gogny 有限程力中 Gauss-Hermite 卷积型积分的解析处理。

展开式：
$$
H_n(x)\,H_m(x)=\sum_{k=0}^{n+m} c_{k}^{(n,m)}\,H_k(x),
$$
系数由正交性得到：
$$
c_{k}^{(n,m)}=\frac{1}{\sqrt{\pi}\,2^{k}k!}
\int_{-\infty}^{+\infty}e^{-x^{2}}H_n(x)H_m(x)H_k(x)\,\mathrm{d}x.
$$
在代码中利用 Gauss-Hermite 节点直接数值计算：
$$
c_{k}^{(n,m)}=\sum_{h} w_h\,e^{-x_h^{2}}\,
\tilde{H}_n(x_h)\,\tilde{H}_m(x_h)\,\tilde{H}_k(x_h),
$$
其中 $\tilde{H}_n(x)=H_n(x)/\sqrt{\sqrt{\pi}\,2^{n}n!}$ 为归一化多项式。

- 实现步骤:
  1. 计算归一化因子 `HERFAC(k)`。
  2. 在所有 Gauss-Hermite 节点上调用 `D_HERM`，得到归一化多项式值 `HERPLN`。
  3. 利用 Gauss-Hermite 数值积分计算三重积分的展开系数 `COEF00(k,m,n)`。
  4. 仅当 $k$ 与 $n+m$ 同奇偶时系数非零，因此循环步长为 $-2$。

##### 变量对应

- 输入变量1: `NOSACT (Integer)`，展开最高阶数
- 输入变量2: `ngh`、`xh(ngh)`、`wh(ngh)`，Gauss-Hermite 网格
- 输出变量1: `COEF00(0:2*NOSACT,0:NOSACT,0:NOSACT) (Real(pr))`，展开系数 $c_{k}^{(n,m)}$
- 输出变量2: `HERFAC(0:2*NOSACT) (Real(pr))`，归一化因子 $\sqrt{\sqrt{\pi}\,2^{k}k!}$

##### 函数数据流
- 调用情况: `preparer <- DEVHER <- D_HERM`
- 被 `preparer`（`hfbtho_solver.f90`）在 Gogny 力计算分支中调用

---

### `pr_gamma`

##### 函数定义

- 函数位置: `982-1109 lines`

```fortran
Real(pr) Function pr_gamma(x)
```

##### 物理公式

- 物理含义: 实数自变量 $\Gamma(x)$ 的双精度计算。采用 Cody & Stoltz 的 rational minimax 逼近：
  - $0<x<1$：利用 $\Gamma(x)=\Gamma(x+1)/x$ 递推到 $(1,2)$ 区间。
  - $1\le x<12$：直接对 $(1,2)$ 区间做有理逼近（8 阶分子/分母）。
  - $x\ge 12$：利用 Stirling 型渐近展开：
$$
\ln\Gamma(x)\approx\Bigl(x-\tfrac{1}{2}\Bigr)\ln x-x+\tfrac{1}{2}\ln(2\pi)+\frac{1}{x}\sum_{i}\frac{c_i}{x^{2i}}.
$$
  - $x<0$：利用反射公式 $\Gamma(x)=\pi/[\sin(\pi x)\,\Gamma(1-x)]$。

##### 变量对应

- 输入变量1: `x (Real(pr))`，自变量
- 输出变量1: 函数返回值，$\Gamma(x)$

##### 函数数据流
- 调用情况: `Class <- pr_gamma`
- 仅由 `Class` 在 Jacobi 和 Laguerre 分支中调用

---

## 模块全局变量

| 变量名 | 类型 | 说明 |
|--------|------|------|
| `ngh` | `Integer(ipr)` | Gauss-Hermite 积分点数 |
| `ngl` | `Integer(ipr)` | Gauss-Laguerre 积分点数 |
| `nghl` | `Integer(ipr)` | 二维网格总点数 `ngh*ngl` |
| `nleg` | `Integer(ipr)` | Gauss-Legendre 积分点数 |
| `xh(ngh)` | `Real(pr)` | Hermite 无量纲节点 $\xi_h$ |
| `wh(ngh)` | `Real(pr)` | 含权 Hermite 权重 $w_h$ |
| `xl(ngl)` | `Real(pr)` | Laguerre 无量纲节点 $\eta_l$ |
| `sxl(ngl)` | `Real(pr)` | $\sqrt{\eta_l}$ |
| `wl(ngl)` | `Real(pr)` | 含权 Laguerre 权重 $w_l$ |
| `xleg(nleg)` | `Real(pr)` | Legendre 节点 $\mu_{\text{leg}}$ |
| `wleg(nleg)` | `Real(pr)` | Legendre 权重 $w_{\text{leg}}$ |
| `fh(nghl)` | `Real(pr)` | 有量纲 z 节点 $z_h$ |
| `fl(nghl)` | `Real(pr)` | 有量纲 $\rho$ 节点 $\rho_l$ |
| `wdcor(nghl)` | `Real(pr)` | 体积元权重 $w_{hl}$ |
| `wdcori(nghl)` | `Real(pr)` | $1/w_{hl}$ |
| `COEF00` | `Real(pr)` | Hermite 乘积展开系数 |
| `HERFAC` | `Real(pr)` | Hermite 归一化因子 |
| `qh`、`qh1` | `Real(pr)` | z 方向 HO 波函数与导数 |
| `ql`、`ql1` | `Real(pr)` | $\rho$ 方向 HO 波函数与导数 |
| `QHLA_opt` 等 | `Real(pr)` | 优化查表数组 |

---

## 模块数据流总览

```
preparer
  ├── gausspoints
  │     └── Gaussq
  │           ├── Class ── pr_gamma
  │           ├── GBSLVE
  │           └── GBTQL2
  ├── gaupol
  ├── coordinateLST_HO
  │     └── optHFBTHO_HO
  └── DEVHER
        └── D_HERM

coulom / recompute_coulomb_expansion
  └── gauleg
```
