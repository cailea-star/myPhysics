# 代码文件 `hfbtho_fission.f90` 说明

本文档对模块 `HFBTHO_fission_fragments`（1317 行）中的全部子程序与函数进行说明。该模块负责裂变碎片相关物理量的计算，包括颈部位置确定、高斯颈部算符、碎片质心、碎片多极矩以及碎片内角动量期望值的评估。

---

### `neck_computeField`

##### 函数定义

- 函数位置: `131-187 lines`

```fortran
Subroutine neck_computeField(ib, gaussian_neck)
```

##### 物理公式

- 物理含义: 在当前 $K^{\pi}$ 分块内，计算高斯颈部算符 $\hat{Q}_{\rm N}$ 的矩阵元。算符形式为
$$
\hat{Q}_{\rm N}(z)=\exp\!\left[-\left(\frac{z-Z_{\rm neck}\,b_{z}}{A_{N}}\right)^{2}\right],
$$
其中 $Z_{\rm neck}$ 为无量纲颈部位置，`AN_VAL` 为高斯宽度（单位 fm），$b_{z}$ 为 z 方向谐振子长度。

矩阵元计算采用 HO 基函数的 Gauss 网格积分：
$$
Q_{\rm N,ij}=\sum_{h,l}w_{hl}\,\phi_{i}(\mathbf{r}_{hl})\,\phi_{j}(\mathbf{r}_{hl})\,
\hat{Q}_{\rm N}(z_{h}),
$$
仅当两态自旋投影满足 $m_{s,i}+m_{s,j}=0$ 时（即同为自旋上或同为自旋下），矩阵元非零。

- 实现步骤:
  1. 在所有 `nghl` 个 Gauss 网格点上计算高斯颈部势值 `Vmom`。
  2. 遍历分块内所有基函数对 $(n_1,n_2)$，利用预计算的 `QHLA_opt` 查表得到 HO 波函数值。
  3. 按自旋投影条件累加得到上三角存储的矩阵元数组 `gaussian_neck`。

##### 变量对应

- 输入变量1: `ib (Integer(ipr))`，当前 $K^{\pi}$ 分块索引
- 输入变量2: `fh(nghl)`，有量纲 z 坐标节点
- 输入变量3: `Z_NECK`、`AN_VAL`、`bz`，颈部参数
- 输入变量4: `QHLA_opt(ntx,nghl)`，预计算 HO 基函数值
- 输出变量1: `gaussian_neck(:) (Real(pr))`，当前分块内 $\hat{Q}_{\rm N}$ 的上三角矩阵元（按行优先压缩存储）

##### 函数数据流
- 调用情况: `neck_expectation <- neck_computeField`
- 被 `neck_expectation` 调用，也可由 `gamdel`（`hfbtho_solver.f90`）在约束计算中通过场构建间接涉及

---

### `QNFIND`

##### 函数定义

- 函数位置: `195-313 lines`

```fortran
Subroutine QNFIND()
```

##### 物理公式

- 物理含义: 确定裂变颈部位置 $Z_{\rm neck}$，使其对应于密度沿 z 轴的极小值点。在极小点处，高斯颈部期望值 $\langle\hat{Q}_{\rm N}\rangle$ 对 $z_{N}$ 的导数应为零：
$$
\frac{\partial}{\partial z_{N}}\langle\hat{Q}_{\rm N}\rangle=0.
$$

- 实现步骤:
  1. 检验系统是否具有反射对称性（`Parity` 标志 + 奇多极矩是否为零）。
  2. 若对称性破缺，在区间 $[-5,+5]$（无量纲）内用 100 点扫描 $\langle\hat{Q}_{\rm N}\rangle$，由数值导数符号变化定位极小值区间。
  3. 在粗定位区间内调用 `ZBRENT`（割线/二分混合法）精确求根。
  4. 若激活 `neck_constraints`，对 $Z_{\rm neck}$ 做线性混合（`xmix`）。
  5. 若对称性守恒，直接取 $Z_{\rm neck}=0$。
  6. 计算最终颈部期望值 `Q_NECK = QMOM_K(Z_NECK, 0)`。

##### 变量对应

- 输入变量1: `iiter (Integer)`，当前迭代次数
- 输入变量2: `Parity (Logical)`，是否守恒反射对称性
- 输入变量3: `neck_constraints (Logical)`，是否对颈部做约束混合
- 输入变量4: `xmix`、`neckLag`，混合参数与 Lagrange 乘子
- 输出变量1: `Z_NECK (Real(pr))`，无量纲颈部位置
- 输出变量2: `Q_NECK (Real(pr))`，高斯颈部算符期望值
- 输出变量3: `neckValue (Real(pr))`，输出到全局的颈部值

##### 函数数据流
- 调用情况: `expect <- QNFIND <- DERIVE / QMOM_K <- ZBRENT`
- 被 `expect`（`hfbtho_solver.f90`）在每次 HFB 迭代后调用

---

### `DERIVE`

##### 函数定义

- 函数位置: `318-332 lines`

```fortran
Real(pr) Function DERIVE(Z_POSI)
```

##### 物理公式

- 物理含义: 计算高斯颈部期望值对颈部位置 $z_{N}$ 的解析导数：
$$
\frac{\partial}{\partial z_{N}}\langle\hat{Q}_{\rm N}\rangle
=\frac{2b_{z}^{2}}{A_{N}^{2}}
\Bigl[I_{1}(z_{N})-z_{N}\,I_{0}(z_{N})\Bigr],
$$
其中 $I_{k}(z_{N})$ 为 `QMOM_K` 计算的 k 阶矩。

##### 变量对应

- 输入变量1: `Z_POSI (Real(pr))`，当前颈部位置（无量纲）
- 输入变量2: `bz`、`AN_VAL`，HO 长度与高斯宽度
- 输出变量1: 函数返回值，导数值

##### 函数数据流
- 调用情况: `QNFIND <- DERIVE <- QMOM_K`
- 仅由 `QNFIND` 调用，作为 `ZBRENT` 的目标函数

---

### `QMOM_K`

##### 函数定义

- 函数位置: `340-367 lines`

```fortran
Real(pr) Function QMOM_K(Z_POSI, KORDER)
```

##### 物理公式

- 物理含义: 计算高斯颈部算符的 $k$ 阶矩：
$$
I_{k}(z_{N})=\int\mathrm{d}^{3}r\,\xi^{k}\,\rho(\mathbf{r})\,
\exp\!\left[-\left(\frac{(\xi-z_{N})b_{z}}{A_{N}}\right)^{2}\right],
\qquad \xi=\frac{z}{b_{z}}.
$$
在代码中利用二维 Gauss 网格数值积分：
$$
I_{k}=\sum_{h,l}w_{hl}\,\bigl[\rho_{n}(\mathbf{r}_{hl})+\rho_{p}(\mathbf{r}_{hl})\bigr]\,
\exp(\cdots)\,\xi_{h}^{k}.
$$

##### 变量对应

- 输入变量1: `Z_POSI (Real(pr))`，颈部位置（无量纲）
- 输入变量2: `KORDER (Integer(ipr))`，矩的阶数 $k$
- 输入变量3: `ro(nghl,2)`，中子/质子密度
- 输入变量4: `wdcor(nghl)`、`xh(ngh)`、`bz`、`AN_VAL`
- 输出变量1: 函数返回值，$I_{k}(z_{N})$

##### 函数数据流
- 调用情况: `QNFIND / DERIVE <- QMOM_K`
- 被 `QNFIND` 和 `DERIVE` 调用

---

### `test_density`

##### 函数定义

- 函数位置: `371-388 lines`

```fortran
Real(pr) Function test_density()
```

##### 物理公式

- 物理含义: 数值检验总密度的归一化：
$$
N_{\text{tot}}=\int\mathrm{d}^{3}r\,\bigl[\rho_{n}(\mathbf{r})+\rho_{p}(\mathbf{r})\bigr]
=\sum_{h,l}w_{hl}\,\bigl[\rho_{n}(\mathbf{r}_{hl})+\rho_{p}(\mathbf{r}_{hl})\bigr].
$$

##### 变量对应

- 输入变量1: `ro(nghl,2)`，中子/质子密度
- 输入变量2: `wdcor(nghl)`，体积元权重
- 输出变量1: 函数返回值，总粒子数

##### 函数数据流
- 调用情况: 独立工具函数，用于调试

---

### `neck_expectation`

##### 函数定义

- 函数位置: `396-430 lines`

```fortran
Subroutine neck_expectation(ib, qval, rho, dd)
```

##### 物理公式

- 物理含义: 在当前分块内计算高斯颈部算符的期望值：
$$
\langle\hat{Q}_{\rm N}\rangle_{K}=\sum_{i,j=1}^{n_{K}}
\hat{Q}_{\rm N,ij}\,\rho_{ji}^{K},
$$
其中 $\rho^{K}$ 为当前分块的密度矩阵，通过准粒子 $V$ 系数构造：
$
\rho^{K}=V^{K}(V^{K})^{T}.
$

- 实现步骤:
  1. 调用 `neck_computeField` 构建该分块的 `gaussian_neck` 矩阵元。
  2. 将压缩存储展开为完整对称矩阵 `dblmul`。
  3. 用 `dgemm` 计算 $\rho\cdot\hat{Q}_{\rm N}$，取对角迹得到期望值。

##### 变量对应

- 输入变量1: `ib (Integer(ipr))`，分块索引
- 输入变量2: `rho(nd,nd) (Real(pr))`，密度矩阵
- 输出变量1: `qval (Real(pr))`，累积的颈部期望值（输入输出，累加模式）
- 输出变量2: `dd (Real(pr))`，累积的密度矩阵迹（输入输出，累加模式）

##### 函数数据流
- 调用情况: `neck <- neck_expectation <- neck_computeField`
- 被 `neck` 调用

---

### `neck`

##### 函数定义

- 函数位置: `435-494 lines`

```fortran
Subroutine neck()
```

##### 物理公式

- 物理含义: 在全空间组态空间中计算高斯颈部算符的期望值，并验证密度矩阵的迹等于总粒子数：
$$
\langle\hat{Q}_{\rm N}\rangle=\text{Tr}\bigl(\hat{Q}_{\rm N}\,\rho\bigr),
\qquad
\text{Tr}(\rho)=N.
$$
分别对中子（`it=1`）和质子（`it=2`）遍历所有 $K^{\pi}$ 分块，利用准粒子 $U$、$V$ 矩阵构造密度矩阵后调用 `neck_expectation`。

##### 变量对应

- 输入变量1: `nb`、`id(ib)`、`kd(ib,it)`、`ka(ib,it)`，分块结构信息
- 输入变量2: `RVqpN/UqpN`、`RVqpP/UqpP`，中子/质子准粒子波函数
- 输出变量1: 屏幕输出 `qvala`（2倍，含中子+质子）与 `dd2`（总粒子数检验）

##### 函数数据流
- 调用情况: `expect <- neck <- neck_expectation <- neck_computeField`
- 被 `expect`（`hfbtho_solver.f90`）在调试模式下调用

---

### `wave_localization`

##### 函数定义

- 函数位置: `498-534 lines`

```fortran
Subroutine wave_localization(ZPOINT)
```

##### 物理公式

- 物理含义: 计算单准粒子波函数在颈部某一侧的局域化概率。利用 `SFACTO(nz1,nz2)` 矩阵（由 `DEFMAS` 计算）表示 z 方向 HO 波函数在颈部一侧的重叠积分，结合径向量子数相同的条件，统计各准粒子态在左/右碎片的分布权重。

##### 变量对应

- 输入变量1: `ZPOINT (Real(pr))`，颈部位置（无量纲）
- 输入变量2: `SFACTO(0:nzx,0:nzx)`，z 方向截断重叠矩阵
- 输出变量: 无直接返回值，中间量用于后续碎片分析

##### 函数数据流
- 调用情况: `wave_localization <- DEVHER / DEFMAS`
- 被外部分析流程调用

---

### `DEFMAS`

##### 函数定义

- 函数位置: `544-593 lines`

```fortran
Subroutine DEFMAS(NZMAXX, ZPOINT, SFACTO)
```

##### 物理公式

- 物理含义: 计算 z 方向 HO 波函数乘积在颈部一侧的截断积分矩阵：
$$
S_{m_{z}n_{z}}=\sum_{k=0}^{n_{z}+m_{z}}C_{n_{z}m_{z}}^{k}\,
\int_{z_{N}}^{+\infty}\!\mathrm{d}\xi\,H_{k}(\xi)\,e^{-\xi^{2}},
$$
其中 $C_{n_{z}m_{z}}^{k}$ 为 `COEF00`（`DEVHER` 预计算），$H_{k}$ 为 Hermite 多项式。

积分利用误差函数关系解析处理：
$$
\int_{z_{N}}^{+\infty}H_{k}(\xi)e^{-\xi^{2}}\mathrm{d}\xi
=\begin{cases}
\dfrac{H_{k}(z_{N})e^{-z_{N}^{2}}}{\text{HERFAC}(k)}, & k\text{ 奇}\\[6pt]
\dfrac{H_{k}(z_{N})e^{-z_{N}^{2}}-H_{k}(0)}{\text{HERFAC}(k)}, & k\text{ 偶}
\end{cases}
$$
（符号随 $z_{N}$ 正负调整）。

##### 变量对应

- 输入变量1: `NZMAXX (Integer(ipr))`，z 方向最大量子数
- 输入变量2: `ZPOINT (Real(pr))`，颈部位置（无量纲）
- 输入变量3: `COEF00`、`HERFAC`、`D_HERM`
- 输出变量1: `SFACTO(0:NZMAXX,0:NZMAXX) (Real(pr))`，截断积分矩阵 $S_{m_{z}n_{z}}$

##### 函数数据流
- 调用情况: `wave_localization <- DEFMAS <- D_HERM / PROINT`
- 被 `wave_localization` 调用

---

### `center_of_mass`

##### 函数定义

- 函数位置: `598-739 lines`

```fortran
Subroutine center_of_mass(Z_POSI, CENLEF, CENRIG)
```

##### 物理公式

- 物理含义: 分别计算左、右裂变碎片的质心位置（单位 fm）：
$$
Z_{\text{CoM}}^{L}=\frac{\displaystyle\int_{-\infty}^{z_{N}}z\,\rho(z)\,\mathrm{d}z}
{\displaystyle\int_{-\infty}^{z_{N}}\rho(z)\,\mathrm{d}z},
\qquad
Z_{\text{CoM}}^{R}=\frac{\displaystyle\int_{z_{N}}^{+\infty}z\,\rho(z)\,\mathrm{d}z}
{\displaystyle\int_{z_{N}}^{+\infty}\rho(z)\,\mathrm{d}z},
$$
其中纵向密度 $\rho(z)$ 通过对横向坐标 $(x,y)$ 的 Gauss-Hermite 积分得到。

- 实现步骤:
  1. 对每个 `ngl` 的径向层，沿 z 方向用 `SPLINE`+`SPLINT` 对密度插值到密网格（`NPOINT=500`）。
  2. 对左、右片段分别用 Simpson 3/8 规则积分 $\rho(z)$ 和 $z\rho(z)$。
  3. 比值得到质心位置。

##### 变量对应

- 输入变量1: `Z_POSI (Real(pr))`，颈部位置（无量纲）
- 输入变量2: `ro(nghl,2)`，密度数组
- 输入变量3: `xh(ngh)`、`wl(ngl)`、`bz`、`bp`
- 输出变量1: `CENLEF (Real(pr))`，左碎片质心（fm）
- 输出变量2: `CENRIG (Real(pr))`，右碎片质心（fm）

##### 函数数据流
- 调用情况: `expect <- center_of_mass <- SPLINE / SPLINT / integrate_simpson`
- 被 `expect`（`hfbtho_solver.f90`）在裂变模式激活时调用

---

### `QLMFRA`

##### 函数定义

- 函数位置: `745-948 lines`

```fortran
Subroutine QLMFRA(Z_POSI, LAMACT, QLMLEF, QLMRIG, CENLEF, CENRIG, I_TYPE)
```

##### 物理公式

- 物理含义: 在左、右裂变碎片的质心系中，计算多极矩算符 $\hat{Q}_{\ell m}$ 的期望值：
$$
Q_{\ell m}^{L/R}=\int_{V_{L/R}}\rho(\mathbf{r})\,r^{\ell}\,Y_{\ell m}(\theta,\phi)\,\mathrm{d}^{3}r,
$$
其中积分区域以颈部 $z_{N}$ 为界，坐标原点平移到各碎片的质心。在轴对称情况下仅计算 $m=0$ 分量，此时 $Y_{\ell 0}\propto P_{\ell}(\cos\theta)$。

- 实现步骤:
  1. 与 `center_of_mass` 类似，先对密度沿 z 方向做样条插值。
  2. 对每个插值 z 点和径向层计算球坐标 $r=\sqrt{\rho_{l}^{2}+(z-Z_{\text{CoM}})^{2}}$ 和 $\cos\theta=(z-Z_{\text{CoM}})/r$。
  3. 调用 `DEFSPH` 计算 $P_{\ell}(\cos\theta)$，组合得到被积函数后用 Simpson 3/8 规则积分。
  4. `I_TYPE=1` 时用总密度，`I_TYPE=2` 时用质子（电荷）密度。

##### 变量对应

- 输入变量1: `Z_POSI (Real(pr))`，颈部位置
- 输入变量2: `LAMACT (Integer(ipr))`，多极阶数 $\ell$
- 输入变量3: `I_TYPE (Integer(ipr))`，密度类型（1=总密度，2=质子密度）
- 输入变量4: `CENLEF`、`CENRIG (Real(pr))`，左/右碎片质心
- 输入变量5: `ro(nghl,2)`、`xl(ngl)`、`wl(ngl)`、`bz`、`bp`
- 输出变量1: `QLMLEF (Real(pr))`，左碎片多极矩
- 输出变量2: `QLMRIG (Real(pr))`，右碎片多极矩

##### 函数数据流
- 调用情况: `expect <- QLMFRA <- DEFSPH <- SPLINE / SPLINT / integrate_simpson`
- 被 `expect`（`hfbtho_solver.f90`）调用

---

### `DEFSPH`

##### 函数定义

- 函数位置: `953-985 lines`

```fortran
Subroutine DEFSPH(LAMACT, COSTHE, SPHHAR)
```

##### 物理公式

- 物理含义: 计算连带 Legendre 函数 $P_{\ell}^{m}(\cos\theta)$ 在 $m=0$ 时的值，并结合多极矩单位化因子得到球谐函数值：
$$
Y_{\ell 0}(\theta)=\sqrt{\frac{2\ell+1}{4\pi}}\,P_{\ell}(\cos\theta)\,q_{\text{units}}(\ell).
$$
其中 $q_{\text{units}}(\ell)$ 来自 `HFBTHO_multipole_moments` 模块的单位转换。

##### 变量对应

- 输入变量1: `LAMACT (Integer(ipr))`，阶数 $\ell$
- 输入变量2: `COSTHE (Real(pr))`，$\cos\theta$
- 输入变量3: `DEFLEG(LAMACT,0,COSTHE)`，连带 Legendre 函数值
- 输出变量1: `SPHHAR(1:2) (Real(pr))`，$Y_{\ell 0}(\theta)$ 及占位零分量

##### 函数数据流
- 调用情况: `QLMFRA <- DEFSPH <- DEFLEG`
- 仅由 `QLMFRA` 调用

---

### `print_moments`

##### 函数定义

- 函数位置: `989-1014 lines`

```fortran
Subroutine print_moments(iw)
```

##### 物理公式

- 物理含义: 格式化输出裂变碎片的多极矩、电荷数、质量数、角动量等物理量。输出量包括：
  - 高斯颈部算符值 `neckValue` 与颈部位置 `Z_NECK`
  - 左/右碎片质心 `CENLEF` / `CENRIG`
  - 电荷数 `QLMPRO(0,:)` 与质量数 `QLMTOT(0,:)`
  - 多极矩 $q_{10}$–$q_{80}$
  - 角动量平方期望值 $\langle J_{x}^{2}\rangle$（来自 $\rho$ 与 $\kappa$ 的贡献）

##### 变量对应

- 输入变量1: `iw (Integer(ipr))`，Fortran 输出单元号
- 输入变量2: `QLMTOT(0:8,0:1)`、`QLMPRO(0:8,0:1)`，碎片多极矩
- 输入变量3: `jxsq_rho_left/right`、`jxsq_kappa_left/right`，角动量分量
- 输出变量: 无返回值，向单元 `iw` 写入格式化文本

##### 函数数据流
- 调用情况: `expect <- print_moments`
- 被 `expect`（`hfbtho_solver.f90`）在裂变计算完成后调用

---

### `calc_qh`

##### 函数定义

- 函数位置: `1018-1040 lines`

```fortran
Subroutine calc_qh(zpoint)
```

##### 物理公式

- 物理含义: 在任意指定 z 坐标点（无量纲）计算归一化 Hermite 波函数值 $\psi_{n_z}(z)/\sqrt{b_z}$，用于碎片内角动量计算中 z 方向波函数的插值。递推公式：
$$
\psi_{0}=\pi^{-1/4}e^{-z^{2}/2}/\sqrt{b_z},
\qquad
\psi_{1}=\sqrt{2}\,z\,\psi_{0},
\qquad
\psi_{n}=\frac{1}{\sqrt{n}}\bigl(\sqrt{2}\,z\,\psi_{n-1}-\sqrt{n-1}\,\psi_{n-2}\bigr).
$$

##### 变量对应

- 输入变量1: `zpoint (Real(pr))`，无量纲 z 坐标
- 输入变量2: `bz`、`nzm`
- 输出变量1: `qh_zneck(0:nzm) (Real(pr))`，Hermite 波函数值数组

##### 函数数据流
- 调用情况: `average_J_in_fragments <- calc_qh`
- 被 `average_J_in_fragments` 调用

---

### `calc_lx`

##### 函数定义

- 函数位置: `1045-1131 lines`

```fortran
Subroutine calc_lx(nz1, nr1, nl1, nz2, nr2, nl2, zc_right, rlx_right, zc_left, rlx_left)
```

##### 物理公式

- 物理含义: 在 HO 基下计算轨道角动量 $x$ 分量 $\hat{\ell}_x$ 的矩阵元。利用柱坐标分解：
$$
\hat{\ell}_x=-i\Bigl(y\frac{\partial}{\partial z}-z\frac{\partial}{\partial y}\Bigr),
$$
由于轴对称性，$\hat{\ell}_x$ 仅当 $|\Lambda_1-\Lambda_2|=1$ 时矩阵元非零。选择定则和升降算符关系导致矩阵元可分解为径向积分和 z 方向积分的乘积。

对于右碎片（颈部右侧 $z>z_N$），z 方向积分利用 `SFACTO` 截断矩阵；对于左碎片（$z<z_N$），需用全空间结果减去右碎片贡献。

##### 变量对应

- 输入变量1: `nz1,nr1,nl1`，bra 态量子数
- 输入变量2: `nz2,nr2,nl2`，ket 态量子数
- 输入变量3: `zc_right`、`zc_left`，右/左碎片质心（fm）
- 输入变量4: `SFACTO`、`qh_zneck`、`ql`、`ql1`、`bz`、`bp`
- 输出变量1: `rlx_right (Real(pr))`，右碎片 $\langle\hat{\ell}_x\rangle$ 矩阵元
- 输出变量2: `rlx_left (Real(pr))`，左碎片 $\langle\hat{\ell}_x\rangle$ 矩阵元

##### 函数数据流
- 调用情况: `average_J_in_fragments <- calc_lx`
- 被 `average_J_in_fragments` 调用

---

### `average_J_in_fragments`

##### 函数定义

- 函数位置: `1135-1313 lines`

```fortran
Subroutine average_J_in_fragments()
```

##### 物理公式

- 物理含义: 计算裂变碎片中总角动量平方的期望值 $\langle\hat{J}^{2}\rangle=2\langle\hat{J}_{x}^{2}\rangle$（轴对称下 $J_y=J_z=0$）。其中
$$
\hat{J}_x=\hat{\ell}_x+\hat{s}_x,
$$
期望值分解为密度矩阵贡献和配对张量贡献：
$$
\langle J_x^{2}\rangle_{\rho}=2\,\text{Tr}\bigl(J_x\,(1-\rho)\,J_x\,\rho\bigr),
\qquad
\langle J_x^{2}\rangle_{\kappa}=-2\,\text{Tr}\bigl(J_x\,\kappa\,J_x\,\kappa\bigr).
$$
分别对中子和质子、左碎片和右碎片进行计算。

- 实现步骤:
  1. 调用 `calc_qh` 在颈部位置预计算 z 方向波函数。
  2. 遍历全部 HO 基函数对，计算 $\hat{\ell}_x$（`calc_lx`）和 $\hat{s}_x$ 矩阵元，组装 $J_x$ 矩阵。
  3. 由准粒子 $U$、$V$ 系数构建全空间的密度矩阵 $\rho$ 和配对张量 $\kappa$。
  4. 对 $\rho$ 部分计算 $(1-\rho)$，通过 `dgemm` 链式乘法得到 $J_x(1-\rho)J_x\rho$ 的迹。
  5. 对 $\kappa$ 部分计算 $J_x\kappa J_x\kappa$ 的迹。
  6. 累加中子+质子结果，存入全局变量 `jxsq_rho_left/right` 和 `jxsq_kappa_left/right`。

##### 变量对应

- 输入变量1: `RVqpN/UqpN`、`RVqpP/UqpP`，准粒子波函数
- 输入变量2: `ka`、`kd`、`id`、`ia`、`nr`、`nl`、`nz`、`ns`，分块与量子数信息
- 输入变量3: `SFACTO`、`qh_zneck`、`CENLEF`、`CENRIG`
- 输出变量1: `jxsq_rho_left/right (Real(pr))`，密度矩阵贡献的 $\langle J_x^{2}\rangle$
- 输出变量2: `jxsq_kappa_left/right (Real(pr))`，配对张量贡献的 $\langle J_x^{2}\rangle$

##### 函数数据流
- 调用情况: `expect <- average_J_in_fragments <- calc_lx / calc_qh / dgemm`
- 被 `expect`（`hfbtho_solver.f90`）在裂变计算分支中调用

---

## 模块全局变量

| 变量名 | 类型 | 说明 |
|--------|------|------|
| `AN_VAL` | `Real(pr)` | 高斯颈部算符宽度（fm） |
| `Q_NECK` | `Real(pr)` | 高斯颈部算符期望值 |
| `Z_NECK` | `Real(pr)` | 颈部无量纲位置 |
| `CENLEF` | `Real(pr)` | 左碎片质心（fm） |
| `CENRIG` | `Real(pr)` | 右碎片质心（fm） |
| `QLMTOT(0:8,0:1)` | `Real(pr)` | 总密度碎片多极矩（左/右） |
| `QLMPRO(0:8,0:1)` | `Real(pr)` | 质子密度碎片多极矩（左/右） |
| `SFACTO(0:nzx,0:nzx)` | `Real(pr)` | z 方向截断重叠积分矩阵 |
| `qh_zneck(0:nzm)` | `Real(pr)` | 颈部处 Hermite 波函数值 |
| `qp_occupation` | `Logical` | 是否使用 qp 占据数 |

---

## 模块数据流总览

```
expect (hfbtho_solver)
  ├── QNFIND
  │     ├── DERIVE ── QMOM_K
  │     └── ZBRENT
  ├── center_of_mass
  │     ├── SPLINE / SPLINT
  │     └── integrate_simpson
  ├── QLMFRA
  │     ├── DEFSPH ── DEFLEG
  │     └── SPLINE / SPLINT / integrate_simpson
  ├── average_J_in_fragments
  │     ├── calc_qh
  │     └── calc_lx ── SFACTO / ql / ql1
  ├── print_moments
  └── neck (debug)
        └── neck_expectation
              └── neck_computeField
```
