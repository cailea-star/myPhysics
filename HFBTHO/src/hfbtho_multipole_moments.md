# 代码文件 `hfbtho_multipole_moments.f90` 说明

本模块定义并计算轴对称多极矩算符 $\hat{Q}_{\lambda 0}$ 在坐标空间和组态空间的矩阵元及期望值。多极矩在柱坐标 $(\rho, z)$ 下表示，归一化常数与 HFODD 程序一致，最终单位为 $\text{barn}^\lambda$。同时包含角动量算符 $\hat{J}_x^2$ 期望值的计算。

---

### `moments_setUnits`

##### 函数定义

- 函数位置: `133-150 lines`

```fortran
Subroutine moments_setUnits()
```

##### 物理公式

- 物理含义: 定义多极矩的标准归一化单位 $q\_units(\lambda)$，使得所有多极矩以 $\text{barn}^\lambda$ 为单位输出。

先设定前置系数：

$
\mathcal{N}_0 = \sqrt{4\pi}, \quad
\mathcal{N}_1 = \sqrt{\frac{4\pi}{3}}, \quad
\mathcal{N}_2 = 2\sqrt{\frac{4\pi}{5}} = \sqrt{\frac{16\pi}{5}}
$

对任意 $\lambda \in [0, \lambda_{\max}]$：

$
q\_units(\lambda) = \frac{\mathcal{N}_\lambda}{10^\lambda}
$

其中当 $\lambda > 2$ 时，初始值 $\mathcal{N}_\lambda = 1$。

##### 变量对应

- 输出变量: `q_units(0:lambdaMax) (Real(pr))`, 公式对应 $q\_units(\lambda)$

##### 函数数据流

- 调用情况: `hfbtho_solver.f90 <- moments_setUnits`

---

### `moments_computeValue`

##### 函数定义

- 函数位置: `154-184 lines`

```fortran
Subroutine moments_computeValue()
```

##### 物理公式

- 物理含义: 在 Gauss 积分网格上计算多极矩的坐标空间期望值。

$
Q_\lambda^{(q)} = \int d^3r \, \rho_q(\mathbf{r}) \, \hat{Q}_\lambda(\mathbf{r})
              = \sum_{i=1}^{nghl} wdcor_i \, \rho_q(i) \, Q_\lambda(z_i, \rho_i^2)
$

其中 $q=1$ (中子), $q=2$ (质子), $q=3$ (总和)。

- 实现步骤:
  1. 遍历 Gauss 积分网格点 `ihli`
  2. 提取权重 $w_i$、密度 $\rho_q$、坐标 $z_i$ 和径向坐标平方 $\rho_i^2 = \text{fl}^2$
  3. 对每个网格点调用 `moments_valueMesh(z, rrr, Qval)` 计算该点处 $Q_\lambda$
  4. 用权重和密度加权累加得到 `qmoment(lambda, q)`，最后求和得到总矩

##### 变量对应

- 输入变量 (模块变量):
  - `ro(ihli, q) (Real(pr))`, 公式对应 $\rho_q(i)$ — 网格点上的中子/质子密度
  - `wdcor(ihli) (Real(pr))`, 公式对应 $w_i$ — 积分权重（已含体积元）
  - `fh(ihli) (Real(pr))`, 公式对应 $z_i$ — 网格点 $z$ 坐标
  - `fl(ihli) (Real(pr))`, 公式对应 $\rho_i$ — 网格点径向坐标（传入时 `rrr = fl**2` 即 $\rho_i^2$）
- 输出变量: `qmoment(0:lambdaMax, 1:3) (Real(pr))`, 公式对应 $Q_\lambda^{(q)}$

##### 函数数据流

- 调用情况: `hfbtho_solver.f90 / hfbtho_pnp.f90 <- moments_computeValue <- moments_valueMesh`

---

### `moments_valueMesh`

##### 函数定义

- 函数位置: `192-227 lines`

```fortran
Subroutine moments_valueMesh(z, rrr, Qval)
```

##### 物理公式

- 物理含义: 在单个网格点 $(\rho, z)$ 上计算多极矩算符 $\hat{Q}_{\lambda}(\rho, z)$ 的值。输入 `rrr` 实际为 $\rho^2$。最终返回的 `Qval` 已乘以 $q\_units(\lambda)$，单位为 $\text{barn}^\lambda$。

低阶矩的显式表达式：

$
Q_0 = 1
$

$
Q_1 = \frac{z}{10}
$

$
Q_2 = \frac{1}{100} \cdot \frac{1}{2}\left(2z^2 - \rho^2\right)
$

$
Q_3 = \frac{1}{1000} \cdot \frac{1}{2}\left(2z^3 - 3z\rho^2\right)
$

一般形式基于球谐函数与 Legendre 多项式：

$
\hat{Q}_{\lambda}(\rho, z)
= \mathcal{N}_\lambda \sqrt{\frac{2\lambda+1}{4\pi}} \, r^\lambda P_\lambda(\cos\theta)
$

其中 $r^2 = \rho^2 + z^2$，$\cos\theta = z/r$。更高阶矩（$\lambda = 4,\dots,8$）使用对应的 $P_\lambda(\cos\theta)$ 展开。若 `Parity = .True.`，则奇数 $\lambda$ 被置零。

##### 变量对应

- 输入变量:
  - `z (Real(pr))`, 公式对应 $z$ — 轴向坐标
  - `rrr (Real(pr))`, 公式对应 $\rho^2$ — 径向坐标平方
- 输出变量: `Qval(0:lambdaMax) (Real(pr))`, 公式对应 $\hat{Q}_{\lambda}(\rho, z)$

##### 函数数据流

- 调用情况: `moments_computeValue / moments_computeField / hfbtho_solver.f90 <- moments_valueMesh`

---

### `moments_computeField`

##### 函数定义

- 函数位置: `235-298 lines`

```fortran
Subroutine moments_computeField(lambda, ib, multMatElems, debug)
```

##### 物理公式

- 物理含义: 计算多极矩算符在谐振子（HO）基下的矩阵元。

$
(Q_\lambda)_{n_1 n_2} = \int d^3r \, \phi_{n_1}(\mathbf{r}) \, \hat{Q}_\lambda(\mathbf{r}) \, \phi_{n_2}(\mathbf{r})
                      = \sum_{ihli} w_{ihli} \, \hat{Q}_\lambda(ihli) \, \phi_{n_1}(ihli) \phi_{n_2}(ihli)
$

- 实现步骤:
  1. 遍历网格点，调用 `moments_valueMesh` 提取指定 $\lambda$ 的 $Q_\lambda(ihli)$，存入 `Vmom`
  2. 对 HO 基函数 `QHLA_opt` 按自旋投影分解为 `OMPFIU`（自旋向上）和 `OMPFID`（自旋向下）
  3. 双重循环遍历基态 $n_1, n_2$，计算 $vh = \phi_{n_1} \phi_{n_2}$（仅当两态自旋投影相同且非零时计入）
  4. 累加 `multMatElems(i) += vh * vnhl`，使用 packed 格式存储对称矩阵

##### 变量对应

- 输入变量:
  - `lambda (Integer(ipr))`, 公式对应 $\lambda$ — 多极阶数
  - `ib (Integer(ipr))`, K 块（block）指标
  - `debug (Logical, Optional)`, 可选调试标志；若传入则只累加 $vh$（不含 $Q_\lambda$）
- 输入/输出变量:
  - `multMatElems(:) (Real(pr), Allocatable)`, 公式对应 $(Q_\lambda)_{n_1 n_2}$ — packed 存储的矩阵元
- 内部变量:
  - `Vmom(ihli) (Real(pr))`, 公式对应 $\hat{Q}_\lambda(ihli)$
  - `vnhl (Real(pr))`, 当前网格点的 $\hat{Q}_\lambda$ 值
  - `vh (Real(pr))`, 公式对应 $\phi_{n_1} \phi_{n_2}$（已匹配自旋）

##### 函数数据流

- 调用情况: `moments_expectation / energy_moments (hfbtho_collective.f90) / hfbtho_solver.f90 <- moments_computeField <- moments_valueMesh`

---

### `moments_expectation`

##### 函数定义

- 函数位置: `303-337 lines`

```fortran
Subroutine moments_expectation(lambda, it, ib, qval, rho, dd, multMatElems)
```

##### 物理公式

- 物理含义: 在组态空间中计算多极矩算符的期望值 $\langle \hat{Q}_\lambda \rangle$，并检查粒子数守恒。

$
dd = \text{Tr}(\rho) = \sum_{n=1}^{nd} \rho_{nn}
$

$
\langle \hat{Q}_\lambda \rangle = \text{Tr}(\rho \, Q_\lambda)
$

其中矩阵元 $Q_\lambda$ 由 `moments_computeField` 计算，通过 `dgemm` 完成矩阵乘法后取迹。

##### 变量对应

- 输入变量:
  - `lambda (Integer(ipr))`, 公式对应 $\lambda$
  - `it (Integer(ipr))`, 粒子种类（1=中子, 2=质子）
  - `ib (Integer(ipr))`, 块指标
  - `rho(:,:) (Real(pr), Allocatable)`, 公式对应 $\rho$ — 密度矩阵
- 输入/输出变量:
  - `qval (Real(pr))`, 公式对应 $\langle \hat{Q}_\lambda \rangle$ — 期望值（累加）
  - `dd (Real(pr))`, 公式对应 $dd$ — 粒子数（累加）
  - `multMatElems(:) (Real(pr), Allocatable)`, 公式对应 $(Q_\lambda)_{n_1 n_2}$

##### 函数数据流

- 调用情况: `hfbtho_collective.f90 <- moments_expectation <- moments_computeField`

---

### `J_computeValue`

##### 函数定义

- 函数位置: `342-461 lines`

```fortran
Subroutine J_computeValue()
```

##### 物理公式

- 物理含义: 计算总角动量平方 $\hat{J}_x^2$ 的期望值，包含正常密度 $\rho$ 和配对张量 $\kappa$ 两部分的贡献。

总角动量算符：

$
\hat{J}_x = \hat{L}_x + \hat{S}_x
$

正常部分（Wick 定理）：

$
\langle \hat{J}_x^2 \rangle_\rho = 2 \, \text{Tr}\big[ j_x \, (1-\rho) \, j_x \, \rho \big]
$

配对部分：

$
\langle \hat{J}_x^2 \rangle_\kappa = 2 \, \text{Tr}\big[ j_x \, \kappa \, j_x \, \kappa \big]
$

- 实现步骤:
  1. 构建 $j_x$ 矩阵元：
     - 轨道部分 `rlx_tot` 由 `calc_lx_tot` 给出，要求 $s_1 = s_2$
     - 自旋部分 `rsx_tot = 0.5$ 当 $n_{r1}=n_{r2},\ n_{z1}=n_{z2},\ n_{l1}=n_{l2},\ s_1 = -s_2$
  2. 从 HFB 系数 $U, V$ 构建密度矩阵 $\rho$ 和配对张量 $\kappa$（中子和质子分别计算）
  3. 计算 $\rho^m = 1 - \rho$
  4. 通过 `dgemm` 链式矩阵乘法计算 $\rho$ 和 $\kappa$ 对 $\langle J_x^2 \rangle$ 的贡献，结果分别存入：
     - `jxsq_rhoN`, `jxsq_rhoP` — 正常部分
     - `jxsq_kappaN`, `jxsq_kappaP` — 配对部分

##### 变量对应

- 输入变量 (模块变量):
  - `RVqpN`, `RUqpN`, `RVqpP`, `RUqpP (Real(pr))`, HFB 变换系数 $V, U$
  - `KpwiN`, `KpwiP`, `ka (Integer(ipr))`, 准粒子索引映射
- 输出变量:
  - `jxsq_rhoN`, `jxsq_rhoP (Real(pr))`, 公式对应 $2 \, \text{Tr}[j_x (1-\rho) j_x \rho]$
  - `jxsq_kappaN`, `jxsq_kappaP (Real(pr))`, 公式对应 $2 \, \text{Tr}[j_x \kappa j_x \kappa]$

##### 函数数据流

- 调用情况: `hfbtho_solver.f90 <- J_computeValue <- calc_lx_tot`

---

### `calc_lx_tot`

##### 函数定义

- 函数位置: `466-518 lines`

```fortran
Subroutine calc_lx_tot(nz1, nr1, nl1, nz2, nr2, nl2, rlx_tot)
```

##### 物理公式

- 物理含义: 在柱坐标谐振子基下计算轨道角动量算符 $\hat{l}_x$ 的矩阵元。

选择定则：$n_{l1} = n_{l2} \pm 1$ 且 $n_{z1} = n_{z2} \pm 1$，否则矩阵元为零。

$
\langle n_{z1} n_{r1} n_{l1} | \hat{l}_x | n_{z2} n_{r2} n_{l2} \rangle = \text{fac1} + \text{fac2} + \text{fac3}
$

其中：

$
\text{fac1} = -\frac{1}{2}(\delta_{n_{l1},n_{l2}+1} - \delta_{n_{l1},n_{l2}-1}) \cdot \text{facr1} \cdot \text{facz1}
$

$
\text{fac2} = -\frac{1}{2}(\delta_{n_{l1},n_{l2}+1} + \delta_{n_{l1},n_{l2}-1}) \cdot n_{l2} \cdot \text{facr2} \cdot \text{facz2}
$

$
\text{fac3} = +\frac{1}{2}(\delta_{n_{l1},n_{l2}+1} - \delta_{n_{l1},n_{l2}-1}) \cdot \text{facr3} \cdot \text{facz2}
$

各项物理来源：
- `facr1`：径向积分 $\sim \int d\rho \, \phi_{n_{r1}n_{l1}} \cdot \rho \cdot \phi_{n_{r2}n_{l2}}$
- `facr2`：径向积分 $\sim \int d\rho \, \phi_{n_{r1}n_{l1}} \cdot \rho^{-1} \cdot \phi_{n_{r2}n_{l2}}$
- `facr3`：径向积分 $\sim \int d\rho \, \phi_{n_{r1}n_{l1}} \cdot \partial_\rho \phi_{n_{r2}n_{l2}}$
- `facz1`：$z$ 方向导数矩阵元（来自 $\partial_z$ 项）
- `facz2`：$z$ 方向坐标矩阵元（来自 $z$ 项）

##### 变量对应

- 输入变量:
  - `nz1, nr1, nl1 (Integer(ipr))`, 左态（bra）量子数 $n_z, n_r, n_l$
  - `nz2, nr2, nl2 (Integer(ipr))`, 右态（ket）量子数 $n_z, n_r, n_l$
- 输出变量:
  - `rlx_tot (Real(pr))`, 公式对应 $\langle n_{z1}n_{r1}n_{l1} | \hat{l}_x | n_{z2}n_{r2}n_{l2} \rangle$

##### 函数数据流

- 调用情况: `J_computeValue <- calc_lx_tot`
