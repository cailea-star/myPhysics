# 代码文件 `hfbtho_solver.f90` 说明

本模块 `HFBTHO_solver` 是 HFBTHO 代码的 **核心 DFT 求解器**，实现轴对称形变的约束/无约束 Hartree-Fock-Bogoliubov（HFB）自洽计算。支持 Skyrme/UNEDF/Gogny 能量密度泛函、HO/THO 基函数、反射对称/非对称、Lipkin-Nogami 修正、温度效应、阻塞（blocking）、多极约束、颈部约束、粒子数投影（PNP）等功能。模块包含约 30 个子程序，按功能可分为：顶层驱动、初始化、基函数构造、自洽迭代循环、HFB 矩阵对角化、密度计算、平均场与配对场构建、Broyden 混合、期望值计算、约束处理及阻塞处理。

---

### `HFBTHO_DFT_SOLVER`

##### 函数定义

- 函数位置: `143-430 lines`

```fortran
Subroutine HFBTHO_DFT_SOLVER(filename_binary)
```

##### 物理公式

- 物理含义: HFBTHO 求解器的顶层入口，负责根据用户输入的 `*_INI` 参数初始化求解器，处理约束、阻塞、投影等设置，并分阶段调用 HFB+HO 和 HFB+THO 计算。计算流程遵循：初始化 → 约束判断 → 阻塞循环 → HFB+HO 自洽迭代 →（可选）HFB+THO 自洽迭代 → 结果输出。

- 实现步骤:
  1. 调用 `initialize_HFBTHO_SOLVER()` 从 `*_INI` 值初始化所有运行参数。
  2. 调用 `Constraint_or_not()` 判断是否存在多极矩约束。
  3. 进入 **阻塞主循环**（`Do While(irestart /= 0)`）：处理奇核子时，先计算偶偶核芯，再逐个尝试阻塞候选准粒子态。
  4. 设置粒子数（阻塞时调整中子/质子数 ±1）。
  5. **HFB+HO 阶段**（`iLST <= 0`）:
     - 调用 `preparer(.True.)` 构建基函数和网格。
     - 尝试从二进制文件读取初始波函数（`inout(1,...)`）。
     - 若读取失败（`iexit > 0`），调用 `start()` 从 Woods-Saxon 势产生初始波函数。
     - 处理阻塞（`handle_blocking`）。
     - 若有 kickoff 约束，先用低精度约束迭代（`iter(.True.,.True.)`），再恢复常规精度。
     - 调用 `iter(.True.,.False.)` 进行常规 HFB+HO 自洽迭代。
     - 调用 `resu()` 输出结果并写入二进制文件。
  6. **HFB+THO 阶段**（`iLST < 0`，即从 HO 结果出发做 THO）:
     - 设置 `iLST1=1`，调用 `coordinateLST()` 构建 THO 基。
     - 调用 `densit()`、`field()` 重新计算密度和场。
     - 调用 `iter(.True.,.False.)` 进行 HFB+THO 自洽迭代。
     - 调用 `resu()` 输出结果。
  7. **纯 THO 阶段**（`iLST > 0`）: 与 HO 阶段类似，但使用 THO 基。

##### 变量对应

- 输入变量: `filename_binary (Character(Len=256), Optional)`，二进制重启/输出文件名
- 输出: 通过全局变量和文件输出计算结果，包括 `ehfb`, `etot`, `rms`, `q2`, `def` 等

##### 函数数据流

- 调用情况: `Main_Program <- HFBTHO_DFT_SOLVER <- initialize_HFBTHO_SOLVER / Constraint_or_not / preparer / inout / start / handle_blocking / iter / resu / coordinateLST / densit / field`

---

### `heading`

##### 函数定义

- 函数位置: `434-503 lines`

```fortran
Subroutine heading()
```

##### 物理公式

- 物理含义: 向输出文件打印程序标题、版本信息（git hash）、核名称、反射对称性设置、MPI/OpenMP 并行信息及时间戳。

##### 变量对应

- 输入: 全局变量 `nucname`, `npr`, `Parity`, `do_print`, `HFB_cores`
- 输出: 写入 `lout`（屏幕）和 `lfile`（`thoout.dat`）

##### 函数数据流

- 调用情况: `preparer <- heading`

---

### `thoalloc`

##### 函数定义

- 函数位置: `508-643 lines`

```fortran
Subroutine thoalloc()
```

##### 物理公式

- 物理含义: 根据当前壳层数 `n00` 和 Gauss 积分点数 `ngh`, `ngl` 分配全部动态数组。计算关键维度参数：

$$
\begin{aligned}
N_{oscil} &= \frac{(n_{00}+1)(n_{00}+2)(n_{00}+3)}{6} \\
nghl &= ngh \cdot ngl \\
ndx &= \frac{(n_{00}+2)^2}{4}, \quad ndx2 = 2 \cdot ndx \\
nhfbx &= ndx + ndx, \quad nqx = ndx \cdot ndx \\
nb2x &= nbx + nbx
\end{aligned}
$$

- 实现步骤:
  1. 根据反射对称性设置 `ngh`（`Parity=.False.` 时 `ngh=2*ngh_INI`）。
  2. 计算各维度上限 `nzx, nrx, nlx, ndx` 等。
  3. 分配 Gauss 网格数组：`xh`, `wh`, `xl`, `wl`, `fl`, `fh`, `wdcor`, `wdcori` 等。
  4. 分配密度数组：`ro`, `tau`, `dro`, `dj`, `aka`, `NABLAR`, `NABLAZ`, `SZFI`, `SFIZ`, `SRFI`, `SFIR`。
  5. 分配约束场 `qfield`、有效质量 `MEFFn/p`。
  6. 分配基函数数组 `qh`, `ql`, `qh1`, `ql1` 及量子数数组 `nz`, `nr`, `nl`, `ns`, `npar`。
  7. 分配 HFB 矩阵与特征值数组 `hfb`, `zhfb`, `evvk`, `erhfb`, `drhfb`。
  8. 分配准粒子 U/V 数组 `RUqpN/P`, `RVqpN/P`, `REqpN/P` 及索引数组 `KqpN/P`, `KpwiN/P`。
  9. 分配 LAPACK 工作空间 `alwork`, `lwork`。
  10. 若为有限程力，调用 `allocate_fr()`。

##### 变量对应

- 输入: 全局变量 `ngh_INI`, `ngl_INI`, `nleg_INI`, `n00`, `Parity`, `nbx`, `ntx`, `nzx`, `nrx`, `nlx`, `ndx`, `nqp`, `nuv`, `lambdaMax`
- 输出: 分配全部模块级动态数组并置零

##### 函数数据流

- 调用情况: `preparer <- thoalloc <- allocate_fr (条件)`

---

### `realloc`

##### 函数定义

- 函数位置: `649-734 lines`

```fortran
Subroutine realloc()
```

##### 物理公式

- 物理含义: 与 `thoalloc` 类似，但用于从二进制文件读取后重新分配数组（当 `functional='READ'` 时）。此时 `n00`, `b0`, `bp`, `bz` 等可能已从文件中读取，与输入不同，需重新计算基空间并分配数组。

##### 变量对应

- 输入: 全局变量（与 `thoalloc` 相同）
- 输出: 重新分配数组

##### 函数数据流

- 调用情况: `HFBTHO_DFT_SOLVER <- realloc`

---

### `preparer`

##### 函数定义

- 函数位置: `740-1048 lines`

```fortran
Subroutine preparer(lpr)
```

##### 物理公式

- 物理含义: 计算前准备子程序。确定核名称、 oscillator 长度、基函数变形、库仑参数、投影网格点、配对窗口、THO 参数等，然后依次调用基函数构建、内存分配、Gauss 点计算、库仑展开、HO 量子数设置、基函数值计算、坐标 LST 变换及 Gogny 矩阵元计算。

关键公式：
- 基函数变形与长度关系：

$$
q = \exp\left(3\sqrt{\frac{5}{16\pi}} \beta_0\right), \quad b_\perp = b_0 q^{-1/6}, \quad b_z = b_0 q^{1/3}
$$

- 若 `b0 <= 0`，使用经验公式：

$$
\hbar\omega = 41 \cdot A^{-1/3} \cdot r_0 \, [\mathrm{MeV}], \quad b_0 = \sqrt{\frac{2\hbar^2}{m\hbar\omega}}
$$

- 库仑交换能系数：$coex = -e^2(3/\pi)^{1/3}$，$cex = -0.75 \cdot coex$

- 质心系修正（若 `use_cm_cor=.True.`）：$hb0 = hbzero \cdot (1 - 1/A)$

- 实现步骤:
  1. 调用 `nucleus()` 确定元素符号。
  2. 若 `lpr=.True.` 且 `do_print=1`，调用 `heading()` 和 `print_functional_parameters()`，输出所有 namelist 内容。
  3. 计算粒子数 `tz`, `amas`，处理阻塞时的偶偶核芯质量数。
  4. 计算库仑常数 `chargee2`, `coex`, `cex`。
  5. 计算 `hb0`, `hb0n`, `hb0p`（含质心系修正）。
  6. 计算基函数变形参数 `q`, `bp`, `bz`, `bpp`。
  7. 设置 PNP 网格点数 `keypjn`, `keypjp`。
  8. 初始化 THO 渐近参数 `ass=0`, `iasswrong=0`。
  9. 初始化迭代变量 `etot`, `varmas`, `rms`, `del`, `alast`, `siold` 等。
  10. 若 `lpr=.True.`，输出运行特征信息（壳层数、HO 长度、基变形、最大迭代次数、初始波函数来源、泛函类型、库仑选项、配对截断、温度、Broyden/线性混合等）。
  11. 依次调用 `gfv()`, `base0(lpr)`, `set_all_gaussians(icou)`, `thoalloc()`, `gausspoints()`, `recompute_coulomb_expansion()`, `base(lpr)`, `gaupol(lpr)`, `coordinateLST()`, `gogny_matrix_elements()`。

##### 变量对应

- 输入变量: `lpr (Logical)`，控制是否打印运行信息
- 输出: 全局变量 `b0`, `bp`, `bz`, `q`, `hb0`, `hb0n`, `hb0p`, `chargee2`, `coex`, `cex`, `ass`, 以及分配的所有数组

##### 函数数据流

- 调用情况: `HFBTHO_DFT_SOLVER / start <- preparer <- heading / nucleus / gfv / base0 / set_all_gaussians / thoalloc / gausspoints / recompute_coulomb_expansion / base / gaupol / coordinateLST / gogny_matrix_elements`

---

### `initialize_HFBTHO_NAMELIST`

##### 函数定义

- 函数位置: `2124-2224 lines`

```fortran
Subroutine initialize_HFBTHO_NAMELIST()
```

##### 物理公式

- 物理含义: 初始化全部 13 组 namelist 的默认值。包括：
  - `HFBTHO_GENERAL`: 壳层数 `n00=10`, 长度 `b0=-1`, 变形 `q=0`, 质子数 `Z=24`, 中子数 `N=26`, 计算类型 `kindhfb=1`
  - `HFBTHO_INITIAL`: 初始变形 `beta2=0`, `beta3=0`, `beta4=0`
  - `HFBTHO_ITERATIONS`: 最大迭代 `maxi=100`, 精度 `epsi=1e-5`
  - `HFBTHO_FUNCTIONAL`: 泛函 `SLY4`, 库仑类型 `icou=2`
  - `HFBTHO_PAIRING`: 配对强度 `V0n=V0p=-300`, 截断 `pwi=60`
  - `HFBTHO_CONSTRAINTS`: 约束多极矩（全零）
  - `HFBTHO_BLOCKING`: 阻塞信息（全零）
  - `HFBTHO_PROJECTION`: 投影参数（关闭）
  - `HFBTHO_TEMPERATURE`: 温度（关闭）
  - `HFBTHO_FEATURES`: 各种高级特性（关闭）
  - `HFBTHO_NECK`: 颈部约束（关闭）
  - `HFBTHO_DEBUG`: 调试参数
  - `HFBTHO_RESTORATION`: 对称性恢复参数

##### 变量对应

- 输入: 无
- 输出: 全局 namelist 变量的默认值

##### 函数数据流

- 调用情况: `read_HFBTHO_NAMELIST / HFBTHO_DFT_SOLVER (间接) <- initialize_HFBTHO_NAMELIST`

---

### `read_HFBTHO_NAMELIST`

##### 函数定义

- 函数位置: `2229-2366 lines`

```fortran
Subroutine read_HFBTHO_NAMELIST(filename_nml)
```

##### 物理公式

- 物理含义: 从输入文件读取全部 13 组 namelist。默认文件名为 `hfbtho_NAMELIST.dat`，也可通过参数指定。读取失败时设置错误标志。

##### 变量对应

- 输入变量: `filename_nml (Character(Len=256), Optional)`，输入文件名
- 输出: 填充全局 namelist 变量

##### 函数数据流

- 调用情况: `Main_Program / HFBTHO_DFT_SOLVER <- read_HFBTHO_NAMELIST`

---

### `check_consistency`

##### 函数定义

- 函数位置: `2370-2564 lines`

```fortran
Subroutine check_consistency()
```

##### 物理公式

- 物理含义: 对输入参数进行一致性检查，确保各参数在有效范围内。检查项包括：
  - 壳层数 `n00_INI` 在 `[1, 50]`
  - 质子/中子数大于 1
  - 计算类型 `kindhfb_INI` 为 `±1`
  - 精度 `epsi_INI > 0`
  - 重启标记 `inin_INI` 在 `{1,2,3}` 中
  - 泛函名称在 45 个预设列表中
  - 配对截断 `pwi_INI >= 0`
  - 表面-体积参数 `CpV1_INI` 在 `[0,1]`
  - 库仑选项 `icou_INI` 在 `[-5, 2]`
  - THO 选项 `iLST_INI` 在 `[-1, 1]`
  - 投影网格点数合法
  - 温度非负
  - Gauss/Laguerre/Legendre 点数合法

##### 变量对应

- 输入: 全局 `*_INI` 变量
- 输出: `ierror_flag`, `ierror_info`

##### 函数数据流

- 调用情况: `Main_Program <- check_consistency`

---

### `initialize_HFBTHO_SOLVER`

##### 函数定义

- 函数位置: `2569-2709 lines`

```fortran
Subroutine initialize_HFBTHO_SOLVER()
```

##### 物理公式

- 物理含义: 从 `*_INI` 变量向运行变量赋值，初始化阻塞状态、混合参数、物理常数、Woods-Saxon 参数、角动量标签等。设置循环范围 `itmin`, `itmax`（处理纯中子滴情况）。最后调用 `set_functional_parameters()` 设置 EDF 参数和 `moments_setUnits()` 设置多极矩单位。

关键初始化值：
- 混合参数：`xmix0=0.1`, `xmix=0.3`, `xmax=0.9`
- Woods-Saxon 参数：`v0ws=-71.28`, `akv=0.4616`, `r0v=1.2334`, `av=0.6150`, `vso=11.1175`
- 物理常数：`amn=938.90590`, `amu=931.4940130`, `hqc=197.328910`, `alphi=137.036020`

##### 变量对应

- 输入: 全局 `*_INI` 变量
- 输出: 全局运行变量 `epsi`, `icou`, `Parity`, `xmix`, `ala`, `pwi` 等

##### 函数数据流

- 调用情况: `HFBTHO_DFT_SOLVER / start <- initialize_HFBTHO_SOLVER <- set_functional_parameters / moments_setUnits`

---

### `adjust_basis`

##### 函数定义

- 函数位置: `2716-2753 lines`

```fortran
Subroutine adjust_basis(q2val,initial_deformations,zz,nn)
```

##### 物理公式

- 物理含义: 基于裂变路径的经验拟合自动调整基函数变形和 oscillator 频率（目前主要针对 $^{240}\mathrm{Pu}$）。

$$
\omega_0 = \begin{cases}
0.1 \exp(-0.02 Q_2) Q_2 + \omega_{large}, & Q_2 > 30 \\
\omega_{small}, & Q_2 \le 30
\end{cases}
$$

基函数长度和变形：

$$
b_0 = \frac{\hbar c}{\sqrt{\omega_0 m_n}}, \quad q_{def} = 0.05 \sqrt{|Q_2|}
$$

##### 变量对应

- 输入变量:
  - `q2val (Real(pr))`: 请求的 $Q_2$ 约束值
  - `initial_deformations (Logical)`: 是否同时调整初始 Woods-Saxon 变形
  - `zz (Integer)`: 质子数
  - `nn (Integer)`: 中子数
- 输出: 全局变量 `oscillator_length`, `basis_deformation`, `beta2_deformation`, `beta4_deformation`

##### 函数数据流

- 调用情况: `Main_Program <- adjust_basis`

---

### `base0`

##### 函数定义

- 函数位置: `2757-2914 lines`

```fortran
Subroutine base0(lpr)
```

##### 物理公式

- 物理含义: 计数并排序柱坐标 HO 基态，确定基空间截断能量 `EBASECUT` 和各维度上限。HO 能量为：

$$
E_{n_z,n_r,m} = \frac{\hbar^2}{m b_z^2}\left(n_z + \frac{1}{2}\right) + \frac{\hbar^2}{m b_\perp^2}\left(2n_r + m + 1\right)
$$

- 实现步骤:
  1. 对 `n00max` 内的全部量子数组合计算能量并排序。
  2. 取第 `NOSCIL=(n00+1)(n00+2)(n00+3)/6` 个能量为截断 `EBASECUT`。
  3. 仅保留能量低于 `EBASECUT` 的态，计数得到 `nbx`（块数）、`ntx`（总态数）、`ndx`（最大块维度）、`nqp`, `nuv`。

##### 变量对应

- 输入变量: `lpr (Logical)`，是否打印基函数统计
- 输出: `nbx`, `ntx`, `ndx`, `nqp`, `nuv`, `nzx`, `nrx`, `nlx`, `nnx`, `EBASECUT`

##### 函数数据流

- 调用情况: `preparer <- base0 <- ord`

---

### `base`

##### 函数定义

- 函数位置: `2918-3114 lines`

```fortran
Subroutine base(lpr)
```

##### 物理公式

- 物理含义: 为每个保留的 HO 态分配量子数 `(nz, nr, nl, ns, npar)`，按 $K^\pi$ 分块（`nb` 块，`nt` 总态数），并计算 Broyden/线性混合所需的存储维度。

块结构：
- 反射对称时（`Parity=.True.`）：每块对应固定 $K = m + 1/2$ 和宇称 $\pi$
- 反射非对称时（`Parity=.False.`）：每块仅对应固定 $K$

Broyden 存储维度：

$$
nhhdim = \sum_{ib=1}^{NB} \frac{nd(ib)(nd(ib)+1)}{2}
$$

##### 变量对应

- 输入变量: `lpr (Logical)`
- 输出: `nb`, `nt`, `ia(ib)`, `id(ib)`, `ikb(ib)`, `ipb(ib)`, `nz(i)`, `nr(i)`, `nl(i)`, `ns(i)`, `npar(i)`, `tb(i)`, `txb(ib)`, `nhhdim`, `nhhdim2`, `nhhdim3`, `nhhdim4`, `brin`, `brout`

##### 函数数据流

- 调用情况: `preparer / HFBTHO_DFT_SOLVER <- base`

---

### `nucleus`

##### 函数定义

- 函数位置: `3375-3413 lines`

```fortran
Subroutine nucleus(is,npr2,te)
```

##### 物理公式

- 物理含义: 简单的查表子程序，实现质子数与元素符号之间的双向查找。支持 `Z=0`（中子）到 `Z=133`（超重新元素）。

##### 变量对应

- 输入变量:
  - `is (Integer(ipr))`: `1`=由质子数查符号，`2`=由符号查质子数
  - `npr2 (Integer(ipr))`: `is=1` 时为输入质子数，`is=2` 时为输出质子数
  - `te (Character(Len=2))`: `is=1` 时为输出元素符号，`is=2` 时为输入元素符号

##### 函数数据流

- 调用情况: `preparer <- nucleus`

---

### `start`

##### 函数定义

- 函数位置: `3183-3371 lines`

```fortran
Subroutine start()
```

##### 物理公式

- 物理含义: 当无法从二进制文件读取初始波函数时，通过求解形变 Woods-Saxon（WS）势产生初始密度和 HFB 矩阵元。WS 势参数来自 von Koepf 和 Ring（1991）：

$$
\begin{aligned}
V_{WS}(r) &= \frac{V_0}{1 + \exp((r - R_{WS})/a)} \\
V_{LS}(r) &= -\frac{\hbar^2}{2m^2c^2}\frac{V_0 \kappa_{so}}{r}\frac{d}{dr}\frac{1}{1 + \exp((r - R_{so})/a_{so})}
\end{aligned}
$$

其中 $V_0 = -71.28$ MeV，$\kappa = 0.4616$，$r_0 = 1.2334$ fm，$a = 0.615$ fm。

核表面采用多极展开（至 $\lambda=4$）：

$$
R(\theta) = R_0 c(\beta)\left[1 + \sum_{\lambda=1}^{4} \beta_\lambda Y_{\lambda 0}(\theta)\right]
$$

体积守恒系数：

$$
c(\beta) = \left[\frac{4\pi}{\int f^3(\theta) \sin\theta d\theta}\right]^{1/3}
$$

- 实现步骤:
  1. 重新初始化求解器，调用 `initialize_HFBTHO_SOLVER()` 和 `preparer(.False.)`。
  2. 限制初始变形参数范围：$|\beta_2| \le 3$, $|\beta_3| \le 1.2$, $|\beta_4| \le 1$。
  3. 计算形变 WS 势在 Gauss 网格上的值，包括中心势 $u$ 和自旋-轨道势 $w$。
  4. 计算初始密度 `ro` 和配对密度 `aka`。
  5. 计算库仑势（均匀带电球近似）。
  6. 设置初始密度矩阵元 `rk=0.1`, `ak=0.1`，调用 `gamdel(.true.,.true.)` 构建初始 HFB 矩阵。

##### 变量对应

- 输入: 全局变量 `npr_INI`, `b0_INI`, `q_INI`, `b2_0`, `b3_0`, `b4_0`, `v0ws`, `akv`, `r0v`, `av`, `vso`, `rso`, `aso`
- 输出: 初始场 `vn`, `vp`, `vsn`, `vsp`, `vhbn`, `vhbp`, `dvn`, `dvp`, `cou`, 初始密度 `ro`, `aka`, 初始 HFB 矩阵元 `brin`

##### 函数数据流

- 调用情况: `HFBTHO_DFT_SOLVER <- start <- initialize_HFBTHO_SOLVER / preparer / gamdel`

---

### `coordinateLST`

##### 函数定义

- 函数位置: `1054-1065 lines`

```fortran
Subroutine coordinateLST()
```

##### 物理公式

- 物理含义: HO/THO 坐标 LST 变换的调度器。根据 `iLST1` 选择调用 HO 版本 `coordinateLST_HO()` 或 THO 版本 `coordinateLST_THO()`。

##### 变量对应

- 输入: `iLST1`
- 输出: 调用对应子程序

##### 函数数据流

- 调用情况: `preparer <- coordinateLST <- coordinateLST_HO / coordinateLST_THO`

---

### `coulom1`

##### 函数定义

- 函数位置: `3418-3480 lines`

```fortran
Subroutine coulom1()
```

##### 物理公式

- 物理含义: 基于 Vautherin  prescription（椭圆积分法）计算轴对称库仑势的直接部分。库仑势核为：

$$
v_{C}(\boldsymbol{r}, \boldsymbol{r}') = \frac{e^2}{2\pi} \int_0^{2\pi} \frac{d\phi}{\sqrt{(z-z')^2 + \rho^2 + \rho'^2 - 2\rho\rho'\cos\phi}}
$$

利用椭圆积分表示：

$$
v_{ik} = \frac{e^2}{2\pi} \sqrt{d_{ik}} \cdot K_2(x_{ik}), \quad d_{ik} = (\rho_i + \rho_k)^2 + (z_i - z_k)^2, \quad x_{ik} = \frac{4\rho_i\rho_k}{d_{ik}}
$$

其中 $K_2$ 为第二类完全椭圆积分 `CompleteEllipticFunction_2nd`。

对于无反射对称的情况，需额外计算 $z < 0$ 部分的贡献。

##### 变量对应

- 输入: `fl(i)`, `fh(i)`, `wdcor(i)`, `chargee2`, `Parity`
- 输出: `vc(i,k)`（库仑势核矩阵），`cou(i)`（质子密度上的库仑场）

##### 函数数据流

- 调用情况: `densit <- coulom1 <- CompleteEllipticFunction_2nd`

---

### `coulom`

##### 函数定义

- 函数位置: `3485-3564 lines`

```fortran
Subroutine coulom()
```

##### 物理公式

- 物理含义: 基于 Gogny prescription 计算库仑势，将 $1/|\boldsymbol{r}-\boldsymbol{r}'|$ 展开为 Gauss-Legendre 积分：

$$
\frac{1}{|\boldsymbol{r}-\boldsymbol{r}'|} = \frac{e^2}{\sqrt{\pi}} \int_0^{\infty} \frac{d\mu}{\mu^2} e^{-(\boldsymbol{r}-\boldsymbol{r}')^2/\mu^2}
$$

通过变量替换 $\mu = b\sqrt{(1-\xi^2)}/\xi$（$b=50$ fm），将积分转化为 Gauss-Legendre 求积：

$$
v_{ik} = \frac{e^2}{\sqrt{\pi}} \frac{1}{b} \sum_{\ell=1}^{nleg} w_\ell \frac{e^{u_\ell^2 \cdot 2\rho\rho'} \left(e^{-u_\ell^2 \cdot d_1} + e^{-u_\ell^2 \cdot d_2}\right)}{(1-\xi_\ell^2)^{3/2}}
$$

其中 $u_\ell = x_\ell / (b\sqrt{1-x_\ell^\beta})^\alpha$，$d_1 = (\rho-\rho')^2 + (z-z')^2$，$d_2 = (\rho-\rho')^2 + (z+z')^2$。

使用修正贝塞尔函数 $e^{-x}I_0(x)$ 计算 `besei0`。

##### 变量对应

- 输入: `fl`, `fh`, `wdcor`, `chargee2`, `Parity`, `nleg`, `xleg`, `wleg`
- 输出: `vc(i,k)`, `cou(i)`

##### 函数数据流

- 调用情况: `densit <- coulom <- besei0`

---

### `coulom_test`

##### 函数定义

- 函数位置: `3568-3628 lines`

```fortran
Subroutine coulom_test()
```

##### 物理公式

- 物理含义: 测试库仑势计算的数值精度，通过固定 Gauss-Legendre 节点计算库仑能量贡献并输出。

##### 变量对应

- 输入: 与 `coulom` 相同
- 输出: 屏幕打印各节点的能量贡献

##### 函数数据流

- 调用情况: （调试调用）

---

### `recompute_coulomb_expansion`

##### 函数定义

- 函数位置: `4742-4761 lines`

```fortran
Subroutine recompute_coulomb_expansion()
```

##### 物理公式

- 物理含义: 为库仑高斯展开计算 Gauss-Legendre 节点和权重。变换关系：

$$
\mu_j = b_{max} \frac{\sqrt{1-t_j^2}}{t_j}, \quad V_j = \frac{2}{\sqrt{\pi}} \frac{w_j}{b_{max}(1-t_j^2)^{3/2}}
$$

其中 $b_{max} = \max(b_\perp, b_z)$，`gauleg` 在 `[0,1]` 上产生节点 $t_j$ 和权重 $w_j$。

##### 变量对应

- 输入: `bp`, `bz`, `n_g_coul`
- 输出: `mu_g_coul(j)`, `V_g_coul(j)`

##### 函数数据流

- 调用情况: `preparer <- recompute_coulomb_expansion <- gauleg`

---

### `iter`

##### 函数定义

- 函数位置: `1069-1271 lines`

```fortran
Subroutine iter(lpr,history)
```

##### 物理公式

- 物理含义: **HFB 自洽迭代主循环**。每次迭代执行：
  1. `hfbdiag(it,0)`：HFB 矩阵对角化（最小正则基）。
  2. `expect(.False.)`：计算密度和期望值。
  3. `getLagrange(ite)`：更新约束的 Lagrange 乘子。
  4. `field()`：计算新平均场和配对场。
  5. `gamdel(.false., do_broyden)`：构建新 HFB 矩阵并进行 Broyden/线性混合。

收敛判据：
- 主要判据：$si = \max|v_{out} - v_{in}| < \epsilon$（`epsi`）
- 额外判据：迭代次数超过 1000 且化学势为正时中断
- 达到 `maxi` 时中断
- 硬中断：若化学势搜索连续失败次数 `iter_fail > iter_fail_max`，强制中断并返回错误

混合参数自适应调整：
- 若 `si` 减小：`xmix = \min(xmax, xmix \times 1.10)`
- 若 `si` 增大：`xmix = 0.2`

- 实现步骤:
  1. 打印迭代表头。
  2. 主循环 `ite = mini, maxi+2`：
     - 打印当前迭代信息（`iiter`, `si`, `xmix`, `bet`, `etot`, `varmas`, `rms`, `ept`, `del`, `alast` 等）。
     - 检查收敛条件或中断条件。
     - 调用 `hfbdiag` → `expect` → `getLagrange` → `field` → `gamdel`。
     - 调整 `xmix`。
  3. 收敛/中断后打印总结信息。

##### 变量对应

- 输入变量:
  - `lpr (Logical)`: 是否打印迭代信息
  - `history (Logical)`: `.True.` 时保留历史（用于 kickoff 约束阶段），`.False.` 时最终迭代不做混合
- 输出: `iiter`（实际迭代次数），`si`（收敛残差），`etot`, `varmas`, `rms`, `ept`, `del`, `alast`, `iError_in_HO/THO`

##### 函数数据流

- 调用情况: `HFBTHO_DFT_SOLVER <- iter <- hfbdiag / expect / getLagrange / field / gamdel`

---

### `hfbdiag`

##### 函数定义

- 函数位置: `1276-1828 lines`

```fortran
Subroutine hfbdiag(it,icanon)
```

##### 物理公式

- 物理含义: **HFB 矩阵块对角化**。对每个 $K^\pi$ 块构建 HFB 矩阵并调用 LAPACK（`DSYEVD`/`DSYEVX`）对角化：

$$
\mathcal{H}_{HFB} = \begin{pmatrix} h - \lambda & \Delta \\ \Delta^\dagger & -(h - \lambda) \end{pmatrix}
$$

其中 $h$ 为单粒子哈密顿量（来自 `brin` 的上三角存储），$\Delta$ 为配对场，$\lambda$ 为化学势。

对角化后得到准粒子能量 $E_k$ 和波函数 $(U_k, V_k)$。密度矩阵和配对张量：

$$
\rho_{ab} = \sum_k V_{ak} V_{bk}^*, \quad \kappa_{ab} = \sum_k U_{ak} V_{bk}^*
$$

阻塞处理：对奇核子，通过最大化新旧阻塞态的重叠来确定实际阻塞的准粒子态。

温度效应：Fermi-Dirac 占据数

$$
f_T(E_k) = \frac{1}{2}\left(1 - \tanh\frac{E_k}{2T}\right)
$$

粒子数归一化：

$$
N = 2\sum_k v_k^2 + 2\sum_k (1-2v_k^2)f_T(E_k)
$$

- 实现步骤:
  1. 初始化指针指向中子/质子数组。
  2. 进入化学势搜索循环（`norm_to_improve`）：
     - 对每个块 `ib` 构建 HFB 矩阵（OpenMP 并行）。
     - 调用 `DSYEVD` 对角化。
     - 处理阻塞（若 `iparenti(it) /= 0`）。
     - 计算每个准粒子态的占据概率 $v_k^2$ 和准粒子能量。
     - 应用配对窗口截断（`pwi` + Fermi 截止函数）。
     - 存储准粒子信息到全局数组。
     - 计算密度矩阵 `rk` 和配对张量 `ak`（含温度修正和阻塞修正）。
     - 调用 `Canonical()` 进行正则基变换（若 `icanon=1`）。
  3. 调用 `ALambda()` 搜索化学势 $\lambda$ 使粒子数守恒。
  4. 计算渐近衰减常数：`ass(it) = 2\sqrt{|(E_{min} - \lambda_{eff})/hb0_{n/p}|}`，其中 $hb0_{n/p}$ 为含质心系修正的有效质量参数；若启用 Lipkin-Nogami（`kindhfb < 0`），化学势需额外包含 `ala2` 修正项。

##### 变量对应

- 输入变量:
  - `it (Integer(ipr))`: `1`=中子，`2`=质子
  - `icanon (Integer(ipr))`: `0`=最小正则基，`1`=完整正则基（用于结果输出）
- 输出: `REqpN/P`, `RUqpN/P`, `RVqpN/P`, `KqpN/P`, `KpwiN/P`, `rk`, `ak`, `ala`, `alast`, `ass`, `erhfb`, `drhfb`, `uk`, `ek`, `dk`, `vk`, `numax`, `ddc` 等

##### 函数数据流

- 调用情况: `iter / resu <- hfbdiag <- ALambda / Canonical / requested_blocked_level`

---

### `ALambda`

##### 函数定义

- 函数位置: `1833-1936 lines`

```fortran
Subroutine ALambda(al,it,kl)
```

##### 物理公式

- 物理含义: **化学势（Fermi 能）搜索**。通过牛顿迭代法求解粒子数守恒方程：

$$
N(\lambda) = \sum_{k=1}^{kl} 2v_k^2(\lambda) = N_{target}
$$

其中 $v_k^2 = \frac{1}{2}\left(1 - \frac{e_k - \lambda}{\sqrt{(e_k - \lambda)^2 + \Delta_k^2}}\right)$，$e_k$ 为参考单粒子能，$\Delta_k$ 为配对隙。

温度修正：

$$
v_k^2 \to v_k^2 + (1 - 2v_k^2)f_T(E_k), \quad f_T = \frac{1}{2}(1 - \tanh(E_k/2T))
$$

导数：

$$\frac{\partial N}{\partial \lambda} = \sum_k 2(1-2f_T)\frac{\partial v_k^2}{\partial \lambda} + 2(1-2v_k^2)\frac{\partial f_T}{\partial \lambda}
$$

- 实现步骤:
  1. 若无配对且泛函不为 D1S（`CpV0=0 .And. trim(skyrme) /= 'D1S'`），取第 $N/2$ 与 $N/2+1$ 个单粒子能的中点为化学势。
  2. 有配对时，以牛顿法迭代求解（最多 500 次）：
     - 计算当前 $\lambda$ 下的粒子数 $sn$ 及其导数 `dez`, `dfz`。
     - 调整上下界 `xinf`, `xsup`。
     - 更新 $\lambda = \lambda - (sn - N)/(dez + dfz)$。
     - 若超出边界，取中点。
     - 收敛判据：$|sn - N|/N < 10^{-10}$。

##### 变量对应

- 输入变量:
  - `it (Integer(ipr))`: 同位旋
  - `kl (Integer(ipr))`: 活跃准粒子态数
- 输出变量: `al (Real(pr))`: 化学势 $\lambda$

##### 函数数据流

- 调用情况: `hfbdiag <- ALambda`

---

### `densit`

##### 函数定义

- 函数位置: `3632-4058 lines`

```fortran
Subroutine densit()
```

##### 物理公式

- 物理含义: **在 Gauss 积分网格上计算所有 Skyrme 密度**。利用准粒子波函数 $U_k, V_k$ 和预计算的 HO/THO 基函数值 `QHLA_opt`, `FI1R_opt`, `FI1Z_opt`, `FI2D_opt` 构建坐标空间密度。

对于每个块 `ib` 和同位旋 `it`，准粒子波函数在坐标空间为：

$$
\phi_k^{(V)}(\boldsymbol{r}) = \sum_a V_{ak} \varphi_a(\boldsymbol{r}), \quad \phi_k^{(U)}(\boldsymbol{r}) = \sum_a U_{ak} \varphi_a(\boldsymbol{r})
$$

各密度分量的计算公式（显式含温度权重 $f_T$，$T=0$ 时 $f_T=0$ 退化为仅 $V$ 贡献）：

$$
\begin{aligned}
\rho(\boldsymbol{r}) &= \sum_k \left[(1-f_T)|\phi_k^{(V)}|^2 + f_T|\phi_k^{(U)}|^2\right] \\
\kappa(\boldsymbol{r}) &= \sum_k \phi_k^{(U)*} \phi_k^{(V)} \\
\tau(\boldsymbol{r}) &= \sum_k \left[(1-f_T)|\nabla\phi_k^{(V)}|^2 + f_T|\nabla\phi_k^{(U)}|^2\right] \\
\Delta\rho(\boldsymbol{r}) &= \sum_k \left[(1-f_T)\left(\phi_k^{(V)*}\Delta\phi_k^{(V)} + |\nabla\phi_k^{(V)}|^2\right) + f_T\left(\phi_k^{(U)*}\Delta\phi_k^{(U)} + |\nabla\phi_k^{(U)}|^2\right)\right] \\
\boldsymbol{J}(\boldsymbol{r}) &= \sum_k \left[(1-f_T)\phi_k^{(V)*}\hat{\boldsymbol{j}}\phi_k^{(V)} + f_T\phi_k^{(U)*}\hat{\boldsymbol{j}}\phi_k^{(U)}\right]
\end{aligned}
$$

其中 $f_T(E_k) = \frac{1}{2}\left(1 - \tanh\frac{E_k}{2T}\right)$ 为准粒子热占据数。

阻塞修正：从总密度中减去被阻塞态的贡献，并加上单粒子波函数的修正项。

- 实现步骤:
  1. 初始化所有密度数组为零。
  2. OpenMP 并行循环遍历块 `ib` 和同位旋 `it`：
     - 读取准粒子波函数 `RVqpN/P`, `RUqpN/P` 到局部数组 `OMPAN`, `OMPANK`。
     - 读取温度占据数 `fn_T`, `fp_T`。
     - 对每个 Gauss 点 `ihil`，通过 `DAXPY` 累加基函数贡献得到坐标空间波函数。
     - 计算各类密度分量（`TROIHIL`, `TTAUIHIL`, `TAKAIHIL`, `TDROIHIL`, `TDJIHIL` 等）。
     - 应用阻塞修正。
  3. OpenMP 归约：将各块密度累加到全局数组。
  4. 去除积分权重并乘以 Jacobi：`Tro = Tro * WGT` 等。
  5. 计算密度归一化因子 `DNFactor(it) = tz(it) / (2*Sum(tro))`。
  6. 调用 `coulom1()` 或 `coulom()` 计算库仑场。

##### 变量对应

- 输入: `RVqpN/P`, `RUqpN/P`, `KpwiN/P`, `fn_T`, `fp_T`, `QHLA_opt`, `FI1R_opt`, `FI1Z_opt`, `FI2D_opt`, `y_opt`, `wdcori`
- 输出: `ro`, `tau`, `dro`, `dj`, `aka`, `SZFI`, `SFIZ`, `SRFI`, `SFIR`, `NABLAR`, `NABLAZ`, `varmas`, `Sumnz`, `DNFactor`, `cou`

##### 函数数据流

- 调用情况: `expect / iter <- densit <- coulom1 / coulom`

---

### `field`

##### 函数定义

- 函数位置: `4103-4438 lines`

```fortran
Subroutine field()
```

##### 物理公式

- 物理含义: **计算 Skyrme 平均场和配对场**。基于当前密度调用 `calculate_U_parameters()` 计算 $U$ 振幅，然后构建各分量场。

平均场的标量部分（以中子为例）：

$$
\begin{aligned}
U_n &= U^{\rho\rho} + U^{\rho\tau}\tau + U^{\rho\Delta\rho}\Delta\rho + U^{\rho\nabla J}\nabla\cdot J + U^{JJ}J^2 + U^{J\nabla\rho}J\cdot\nabla\rho \\
&\quad + V_{Coul} + V_{ext} + V_{constrain}
\end{aligned}
$$

有效质量项：

$$
U^{\tau}_n = \frac{\hbar^2}{2m_n^*} = U^{\rho\tau}(0,0)\rho_0 + U^{\rho\tau}(3,0)\rho_1 + \frac{\hbar^2}{2m_n}(1 - 1/A)
$$

自旋-轨道场：

$$
\boldsymbol{W}_n = U^{\rho\nabla J}\nabla\rho + U^{J\nabla\rho}\boldsymbol{J}
$$

配对场（delta 配对，无正规化）：

$$
\Delta_n(\boldsymbol{r}) = V_0^{(n)}(1 - \eta_n \rho/\rho_c) \kappa_n(\boldsymbol{r})
$$

配对正规化（`pairing_regularization=.True.`）：

$$
\Delta_n = \frac{\kappa_n}{g_{eff}^{-1}}, \quad g_{eff}^{-1} = \frac{1}{g} - \frac{1}{4\pi^2}\left(\frac{k_c}{m^*} - \frac{k_F}{2m^*}\ln\frac{k_c+k_F}{k_c-k_F}\right)
$$

其中 $k_c = \sqrt{2m^*(\lambda + E_c - U)}$，$k_F = \sqrt{2m^*(\lambda - U)}$。

**注意**：上述公式仅适用于 $pUr + ec - \lambda < 0$ 且 $\lambda > pUr$ 的参数区间。代码中还包含 $\lambda \le pUr$ 等其他分支（对应反正切项或取 $g_{eff}^{-1}=1/g$），具体见 `field()` 中多分支判断逻辑。

约束势：

$$V_{constrain} = -\sum_\lambda \lambda_\lambda Q_\lambda(\boldsymbol{r}) - \lambda_{neck} e^{-(z-Z_{neck})^2/\sigma^2}$$

- 实现步骤:
  1. 对每个 Gauss 点 `ihli`，构建同位旋标量/矢量密度组合 `RHO_0/1`, `TAU_0/1`, `DRHO_0/1`, `DJ_0/1` 等。
  2. 调用 `calculate_U_parameters()` 计算全部 $U$ 振幅。
  3. 按密度依赖的 $U$ 振幅公式计算各偏导数场 `tUr`, `tUt`, `tUDr`, `tUNr`, `tUNz`, `tUDj`, `tUFIZ`, `tUZFI`, `tUFIR`, `tURFI`。
  4. 转换回质子-中子表象 `pUr`, `pUt`, `pUDr` 等。
  5. 加入约束势、库仑势（直接+交换 Slater）、外部势。
  6. 若启用配对正规化，计算有效配对强度 `geff_inv` 和正规化配对场 `dvn`, `dvp`；否则使用标准 delta 配对。
  7. 存储所有场分量到全局数组 `vn`, `vp`, `vhbn`, `vhbp`, `vrn`, `vrp`, `vzn`, `vzp`, `vdn`, `vdp`, `vsn`, `vsp`, `vSFIZn/p`, `vSZFIn/p`, `vSFIRn/p`, `vSRFIn/p`。

##### 变量对应

- 输入: `ro`, `tau`, `dro`, `dj`, `aka`, `NABLAR`, `NABLAZ`, `SZFI`, `SFIZ`, `SRFI`, `SFIR`, `cou`, `vDHartree`, `multLag`, `neckLag`
- 输出: 全部势场数组 `vn`, `vp`, `vhbn`, `vhbp`, `vrn`, `vrp`, `vzn`, `vzp`, `vdn`, `vdp`, `vsn`, `vsp`, `dvn`, `dvp` 及张量场

##### 函数数据流

- 调用情况: `iter <- field <- calculate_U_parameters / moments_valueMesh`

---

### `constraining_field`

##### 函数定义

- 函数位置: `4442-4472 lines`

```fortran
Subroutine constraining_field()
```

##### 物理公式

- 物理含义: 预计算约束多极场和颈部高斯场在 Gauss 网格上的值，存入 `qfield(ihli,lambda)`。

##### 变量对应

- 输入: `multLambda`, `numberCons`
- 输出: `qfield(ihli,lambda)`

##### 函数数据流

- 调用情况: （独立调用）

---

### `gamdel`

##### 函数定义

- 函数位置: `4477-4707 lines`

```fortran
Subroutine gamdel(WoodsSaxon,DoMixing)
```

##### 物理公式

- 物理含义: **在组态空间（HO/THO 基）中构建 HFB 平均场和配对场矩阵元**。通过 Gauss 积分将坐标空间场投影到 HO 基上：

$$
\begin{aligned}
h_{ab} &= \int d^3r \, \varphi_a^*(\boldsymbol{r}) \left[-\frac{\hbar^2}{2m^*}\nabla^2 + U(\boldsymbol{r})\right] \varphi_b(\boldsymbol{r}) \\
\Delta_{ab} &= \int d^3r \, \varphi_a^*(\boldsymbol{r}) \Delta(\boldsymbol{r}) \varphi_b(\boldsymbol{r})
\end{aligned}
$$

利用预计算的基函数导数 `FI1R_opt`, `FI1Z_opt`, `FI2D_opt`，矩阵元通过解析积分公式直接计算。

Lipkin-Nogami 修正：若 `kindhfb < 0`，在哈密顿量对角元上添加 $2\lambda_2 \rho_{aa}$。

- 实现步骤:
  1. 若启用颈部约束，调用 `QNFIND()` 确定颈部位置。
  2. OpenMP 并行遍历各块 `ib` 和 Gauss 点 `ihli`：
     - 读取坐标空间场值到局部变量。
     - 对每个基态对 `(N1,N2)` 计算重叠积分和导数矩阵元。
     - 根据自旋组合（UpUp, DoDo, UpDo, DoUp）分别计算不同密度通道的贡献。
     - 累加到场矩阵元 `brout`（中子：1–`nhhdim`，质子：`nhhdim+1`–`nhhdim2`，配对：`nhhdim2+1`–`nhhdim3`/`nhhdim4`）。
     - 添加 LN 修正（若 `kindhfb < 0`）。
  3. 若为有限程力，调用 `gamdel_gogny()` 添加 Gogny 矩阵元。
  4. 将约束 Lagrange 乘子存入 `brout(nhhdim4+lambda)`。
  5. 若 `DoMixing=.True.`，调用 `broyden_min()` 进行 Broyden/线性混合，更新 `brin`；否则直接令 `brin=brout`。

##### 变量对应

- 输入变量:
  - `WoodsSaxon (Logical)`: `.True.` 时跳过 Gogny 矩阵元（用于初始 WS 势）
  - `DoMixing (Logical)`: 是否启用 Broyden/线性混合
- 输入: `vn`, `vp`, `vhbn`, `vhbp`, `vrn`, `vrp`, `vzn`, `vzp`, `vdn`, `vdp`, `vsn`, `vsp`, `dvn`, `dvp` 及张量场
- 输出: `brin`（混合后的新场矩阵元），`si`（残差），`multLag`, `neckLag`

##### 函数数据流

- 调用情况: `iter / start <- gamdel <- gamdel_gogny (条件) / broyden_min`

---

### `broyden_min`

##### 函数定义

- 函数位置: `4776-4868 lines`

```fortran
Subroutine broyden_min(N,vout,vin,alpha,si,iter,M,bbroyden)
```

##### 物理公式

- 物理含义: **改进的 Broyden 混合方法**（Johnson, 1988）。用于加速 HFB 自洽迭代的收敛。定义残差 $f^{(n)} = v_{out}^{(n)} - v_{in}^{(n)}$，Broyden 更新为：

$$
v_{in}^{(n+1)} = v_{in}^{(n)} + G^{(n)} f^{(n)}
$$

其中 $G^{(n)}$ 为拟牛顿逆 Jacobian 的近似。通过存储前 $M$ 步的 $df^{(m)} = f^{(m+1)} - f^{(m)}$ 和 $dv^{(m)} = v_{in}^{(m+1)} - v_{in}^{(m)}$ 构造低秩更新。

归一化：$\Delta f^{(m)} = df^{(m)} / ||df^{(m)}||$。

矩阵 $\beta_{ij} = (\Delta f^{(i)}, \Delta f^{(j)}) + w_0^2 \delta_{ij}$，$w_0 = 0.01$。

曲率判据：若 $(f, \text{curv}) > -1$，接受 Broyden 步；否则回退到线性混合。

- 实现步骤:
  1. 计算残差 `vout = vout - vin`，$si = \max|vout|$。
  2. `M <= 0`：线性混合 `vin = vin + xmix * vout`（使用全局混合参数 `xmix`）；`iter == 0`：线性混合 `vin = vin + alpha * vout`。
  3. Broyden 步：
     - 存储前 $M$ 步的归一化 $df$, $dv$。
     - 构建并分解 $\beta$ 矩阵（`DSYTRF`/`DPOTRF`）。
     - 计算曲率向量 `curv = alpha*vout - sum(gamma_j*(dv_j + alpha*df_j))`。
     - 若曲率判据通过，更新 `vin = vin + curv`；否则线性混合 `vin = vin + 0.5*alpha*vout`。

##### 变量对应

- 输入变量:
  - `N (Integer(ipr))`: 向量维度
  - `vout(N) (Real(pr))`: 当前输出场（**入出**，先被减去 `vin` 变为残差）
  - `vin(N) (Real(pr))`: 当前输入场（**入出**，被更新为新输入）
  - `alpha (Real(pr))`: 混合因子
  - `iter (Integer(ipr))`: 当前迭代数
  - `M (Integer(ipr))`: Broyden 历史长度
- 输出变量:
  - `si (Real(pr))`: 最大残差 $\max|vout-vin|$
  - `bbroyden (Character(1))`: `'B'`=Broyden 混合，`'L'`=线性混合

##### 函数数据流

- 调用情况: `gamdel <- broyden_min`

---

### `expect`

##### 函数定义

- 函数位置: `4872-5562 lines`

```fortran
Subroutine expect(lpr,type_basis)
```

##### 物理公式

- 物理含义: **计算可观测量的期望值**，包括能量逐项分解、半径、多极矩、变形、角动量、熵等。总能量：

$$
\begin{aligned}
E_{tot} &= E_{kin} + E_{vol} + E_{surf} + E_{so} + E_{Coul} + E_{pair} + E_{tens} + E_{extra} + E_{ext} \\
E_{kin} &= \frac{\hbar^2}{2m_n}\tau_n + \frac{\hbar^2}{2m_p}\tau_p \\
E_{vol} &= \int d^3r \left[U^{\rho\tau}\rho\tau + U^{\rho\rho}\rho^2 + \frac{1}{2}V_{Hartree}^{dir}\rho\right] \\
E_{surf} &= \int d^3r \left[U^{\rho\Delta\rho}\rho\Delta\rho + U^{\nabla\rho\nabla\rho}(\nabla\rho)^2\right] \\
E_{so} &= \int d^3r \left[U^{\rho\nabla J}\rho\nabla\cdot J + U^{J\nabla\rho}\boldsymbol{J}\cdot\nabla\rho\right] \\
E_{tens} &= \int d^3r \left[U^{JJ}J^2 + U^{J_{ab}J_{ba}}J_{ab}J_{ba}\right]
\end{aligned}
$$

配对能（delta 配对）：

$$E_{pair} = \int d^3r \, \Delta(\boldsymbol{r})\kappa(\boldsymbol{r}) = \sum_{it} \int \Delta_{it} \kappa_{it} \, d^3r$$

库仑能：直接部分 $E_{Coul}^{dir} = \frac{1}{2}\int \rho_p V_{Coul} d^3r$，交换部分采用 Slater 近似 $E_{Coul}^{ex} = -C_{ex}\int \rho_p^{4/3} d^3r$。

rms 半径：$R_{rms} = \sqrt{\langle r^2 \rangle}$，电荷半径：$R_c = \sqrt{R_p^2 + 0.769 + (N/Z)(-0.1161) + 0.033}$ fm。

多极矩：$Q_\lambda = \int r^\lambda P_\lambda(\cos\theta) \rho d^3r$。

变形参数：$\beta_2 = \sqrt{\pi/5} \, Q_2 / (A R_{rms}^2)$。

- 实现步骤:
  1. 调用 `densit()` 计算密度（若 `basis=1`）。
  2. 在 Gauss 网格上逐项积分能量密度，计算各能量分量、rms、多极矩。
  3. 计算有限程力贡献（`trace_product`）。
  4. 计算库仑高斯展开贡献（`trace_product_2`）。
  5. 汇总总能量 `etot = ekt + evol + esurf + eso + ecoul + ept + etens + ...`
  6. 若启用 LN，调用 `tracesln()` 并修正能量和密度。
  7. 调用 `moments_computeValue()` 计算多极矩，`J_computeValue()` 计算角动量。
  8. 若启用集体惯性且 `lpr=.True.`，调用 `calculate_collective_mass()`。
  9. 若启用裂变碎片分析，调用相关子程序计算碎片性质。
  10. 若 `lpr=.True.`，向输出文件打印全部结果。
  11. 若启用 PNP，调用 `expectpj()` 进行粒子数投影后的期望值计算。

##### 变量对应

- 输入变量:
  - `lpr (Logical)`: `.True.` 时打印并存储最终结果，`.False.` 时仅计算期望值供迭代使用
  - `type_basis (Integer(ipr), Optional)`: 基类型选择（`1`=Bogoliubov 基）
- 输出: `etot`, `ehfb`, `rms`, `q2`, `q4`, `def`, `bet`, `ept`, `del`, `ecoul`, `eso`, `esurf`, `evol`, `ekt`, `entropy`, `eresu`, `eresl` 等

##### 函数数据流

- 调用情况: `iter / resu <- expect <- densit / moments_computeValue / J_computeValue / calculate_collective_mass / tracesln / expectpj / trace_product / trace_product_2`

---

### `resu`

##### 函数定义

- 函数位置: `1940-2120 lines`

```fortran
Subroutine resu(filename_binary)
```

##### 物理公式

- 物理含义: 收敛后输出最终结果。首先写入二进制文件（若 `write_hel=.True.`），然后调用 `hfbdiag(it,1)` 进行完整正则基对角化，再次写入二进制文件，最后打印准粒子能谱、正则单粒子能谱，并调用 `expect(.True.,1)` 计算并打印全部期望值。

- 实现步骤:
  1. 判断是否需要写入二进制文件。
  2. 调用 `inout(2,...)` 写入预对角化数据。
  3. 调用 `hfbdiag(it,1)` 对每个同位旋进行完整正则基对角化。
  4. 调用 `inout(3,...)` 写入后对角化数据。
  5. 打印准粒子能谱：$E_k^{qp}$, $e_k^{ref}$, $v_k^2$, $\Delta_k$。
  6. 打印正则单粒子能谱：$\varepsilon_k^{canon}$, $v_k^2$, $\Delta_k^{canon}$。
  7. 计算 Lipkin-Nogami 相关求和量 `ssln`。
  8. 调用 `expect(.True.,1)` 打印全部期望值。

##### 变量对应

- 输入变量: `filename_binary (Character(Len=256), Optional)`
- 输出: 屏幕/文件输出，二进制文件

##### 函数数据流

- 调用情况: `HFBTHO_DFT_SOLVER <- resu <- inout / hfbdiag / expect`

---

### `Constraint_or_not`

##### 函数定义

- 函数位置: `5566-5585 lines`

```fortran
Subroutine Constraint_or_not(inin_INI0,inin0,icstr0)
```

##### 物理公式

- 物理含义: 判断是否存在激活的多极矩约束。遍历 `lambda_active(l)`，若有正值则 `icstr0=1`，否则为 0。

##### 变量对应

- 输入变量: `inin_INI0 (Integer(ipr))`
- 输出变量: `inin0 (Integer(ipr))`, `icstr0 (Integer(ipr))`

##### 函数数据流

- 调用情况: `HFBTHO_DFT_SOLVER / start <- Constraint_or_not`

---

### `getLagrange`

##### 函数定义

- 函数位置: `5592-6086 lines`

```fortran
Subroutine getLagrange(ite)
```

##### 物理公式

- 物理含义: **基于 QRPA 微扰 Cranking 近似更新约束的 Lagrange 乘子**。约束算符 $\hat{F}_i$ 在准粒子基下的矩阵元为：

$$
\tilde{F}_i^{12} = U^\dagger F_i V - V^\dagger F_i U, \quad \tilde{F}_i^{11} = U^\dagger F_i U - V^\dagger F_i V
$$

约束关联矩阵（极化率）：

$$
M_{ij} = \sum_{kl} \frac{\tilde{F}_i^{12}(k,l) \tilde{F}_j^{12}(k,l)}{E_k + E_l} \quad (T=0)
$$

温度修正：

$$
M_{ij} = \sum_{kl} \tilde{F}_i^{12}(k,l) \tilde{F}_j^{12}(k,l) \frac{1 + f_k + f_l}{E_k + E_l} + \frac{1}{2}\sum_{kl} \tilde{F}_i^{11}(k,l) \tilde{F}_j^{11}(k,l) \frac{f_k - f_l}{E_k - E_l}
$$

Lagrange 乘子更新：

$$
\lambda_i^{new} = \lambda_i^{old} + \sum_j (M^{-1})_{ij} (Q_j^{requested} - Q_j^{current})
$$

- 实现步骤:
  1. 计算当前多极矩与目标值的偏差 `cnsvec`。
  2. OpenMP 并行遍历块 `ib` 和同位旋 `it`：
     - 对每个约束算符，调用 `moments_computeField()` 或 `neck_computeField()` 计算其在 HO 基下的矩阵元。
     - 通过 `dgemm` 变换到准粒子基得到 $\tilde{F}^{12}$ 和 $\tilde{F}^{11}$。
     - 计算约束关联矩阵 `cnsmat`（含温度修正）。
  3. 对所有块求和得到总关联矩阵 `cnsorg`。
  4. 调用 `dgetrf` 和 `dgetri` 求逆。
  5. 计算 Lagrange 乘子修正并更新 `multLag` 和 `neckLag`。

##### 变量对应

- 输入变量: `ite (Integer(ipr))`，当前迭代数
- 输入: `qmoment`, `multRequested`, `neckRequested`, `brin`, `RUqpN/P`, `RVqpN/P`, `REqpN/P`
- 输出: `multLag(lambda)`, `neckLag`, `brout(nhhdim4+...)`

##### 函数数据流

- 调用情况: `iter <- getLagrange <- moments_computeField / neck_computeField`

---

### `requested_blocked_level`

##### 函数定义

- 函数位置: `6090-6133 lines`

```fortran
Subroutine requested_blocked_level(ib,it)
```

##### 物理公式

- 物理含义: 根据用户指定的量子数（$\Omega$, $\pi$, $N=n_z+2n_r+n_l$, $n_z$, $n_l$）在 HFB 对角化后的波函数中搜索匹配的阻塞候选态。通过最大化 $|U|$ 或 $|V|$ 分量确定主分量对应的量子数。

##### 变量对应

- 输入变量:
  - `ib (Integer(ipr))`: 当前 $K^\pi$ 块编号
  - `it (Integer(ipr))`: 同位旋
- 输入: `nkblo(it,1:5)`（用户指定的阻塞量子数）
- 输出: `keyblo(it)`, `bloblo(keyblo,it)`, `blo123(keyblo,it)`

##### 函数数据流

- 调用情况: `hfbdiag <- requested_blocked_level`

---

### `handle_blocking`

##### 函数定义

- 函数位置: `6137-6230 lines`

```fortran
Subroutine handle_blocking(irestart,iexit)
```

##### 物理公式

- 物理含义: 管理阻塞计算的循环逻辑。对于奇核子，先计算偶偶核芯（`irestart=0`），然后逐个尝试阻塞候选态（`irestart=1,2,...`）。支持特定阻塞（用户指定量子数）和自动阻塞（在能量窗口内遍历所有候选）。

- 实现步骤:
  1. 检查阻塞模式一致性（不允许同时自动和特定阻塞）。
  2. 若从文件读取且为首次迭代，设置 `irestart=1` 触发阻塞循环。
  3. `irestart=1`：初始化自动阻塞循环变量 `nloops`。
  4. `irestart >= 2`：递增阻塞候选索引 `keyblo(prim)`，处理奇-奇核时同时递增第二个同位旋索引。
  5. 当所有候选遍历完毕后设置 `irestart=0` 退出循环。

##### 变量对应

- 输入变量:
  - `irestart (Integer(ipr))`: 入出，控制阻塞循环状态
  - `iexit (Integer(ipr))`: 文件读取状态
- 输出: `irestart`, `keyblo`, `blocking_never_done`

##### 函数数据流

- 调用情况: `HFBTHO_DFT_SOLVER <- handle_blocking`
