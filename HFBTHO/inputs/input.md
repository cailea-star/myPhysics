# HFBTHO 输入文件说明

HFBTHO 的输入文件采用 Fortran **Namelist** 格式，文件名默认为 `hfbtho_NAMELIST.dat`（可通过命令行参数覆盖）。文件由若干以 `&HFBTHO_` 开头的 namelist 组构成，每组包含若干 `变量名 = 值` 的键值对，以 `/` 结束。

> **注意**：所有 namelist 组均为可选。若输入文件中未出现某 namelist 组，程序将使用内部硬编码的默认值（见下表“默认值”列）。

---

## 文件格式示例

```fortran
&HFBTHO_GENERAL
 number_of_shells = 10, oscillator_length = -1.0, basis_deformation = 0.0,
 proton_number = 20, neutron_number = 28, type_of_calculation = 1 /
&HFBTHO_INITIAL
 beta2_deformation = 0.0, beta3_deformation = 0.0, beta4_deformation = 0.0 /
&HFBTHO_ITERATIONS
 number_iterations = 50, accuracy = 1.E-5, restart_file = -1 /
&HFBTHO_FUNCTIONAL
 functional = 'SKM*', add_initial_pairing = F, type_of_coulomb = 2 /
&HFBTHO_PAIRING
 user_pairing = T, vpair_n = -280.0, vpair_p = -280.0,
 pairing_cutoff = 60.0, pairing_feature = 0.5 /
&HFBTHO_CONSTRAINTS
 lambda_values = 1, 2, 3, 4, 5, 6, 7, 8,
 lambda_active = 0, 0, 0, 0, 0, 0, 0, 0,
 expectation_values = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /
&HFBTHO_BLOCKING
 proton_blocking = 0, 0, 0, 0, 0, neutron_blocking = -5, 0, 2, 0, 2 /
&HFBTHO_PROJECTION
 switch_to_THO = 0, projection_is_on = 0,
 gauge_points = 1, delta_Z = 0, delta_N = 0 /
&HFBTHO_TEMPERATURE
 set_temperature = F, temperature = 0.0 /
&HFBTHO_FEATURES
 collective_inertia = F, fission_fragments = F, pairing_regularization = F,
 automatic_basis = F, localization_functions = F /
&HFBTHO_NECK
 set_neck_constrain = F, neck_value = 13.00 /
&HFBTHO_DEBUG
 number_Gauss = 40, number_Laguerre = 40, number_Legendre = 80,
 compatibility_HFODD = F, number_states = 500, force_parity = T,
 print_time = 0 /
&HFBTHO_RESTORATION
 PNP_is_on = 0, number_of_gauge_points = 9, delta_neutrons = 6, delta_protons = 6,
 AMP_is_on = 0, number_of_rotational_angles = 16, maximal_angular_momentum = 14 /
```

---

## Namelist 变量详解

### `&HFBTHO_GENERAL` — 全局计算参数

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `number_of_shells` | `Integer` | 谐振子基截断主壳层数 $N_{\max}$ | `10` | [1, 50] | `check_consistency`, `execute_HFBTHO`, `preparer`, `thoalloc`, `base0` |
| `oscillator_length` | `Real` | 谐振子长度 $b_0$（fm） | `-1.0` | `-1.0`=自动计算；`>0`=用户指定 | `adjust_basis`, `execute_HFBTHO`, `preparer` |
| `basis_deformation` | `Real` | 基函数的初始四极形变 $q_{20}^{(\text{ini})}$ | `0.0` | 任意实数 | `adjust_basis`, `compute_PES`, `compute_driplines`, `preparer` |
| `proton_number` | `Integer` | 质子数 $Z$ | `24` | >1 | `Main_Program`, `compute_PES`, `compute_driplines`, `check_consistency` |
| `neutron_number` | `Integer` | 中子数 $N$ | `26` | >1 | `Main_Program`, `compute_PES`, `compute_driplines`, `check_consistency` |
| `type_of_calculation` | `Integer` | 计算类型 | `1` | `1`=HFB；`-1`=HFB+LN | `check_consistency`, `execute_HFBTHO`, `HFBTHO_DFT_SOLVER` |

### `&HFBTHO_INITIAL` — 初始 Woods-Saxon 势形变

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `beta2_deformation` | `Real` | 初始四极形变 $\beta_2$ | `0.0` | 任意实数 | `adjust_basis`, `compute_PES`, `execute_HFBTHO`, `start` |
| `beta3_deformation` | `Real` | 初始八极形变 $\beta_3$ | `0.0` | 任意实数 | `compute_PES`, `execute_HFBTHO`, `start` |
| `beta4_deformation` | `Real` | 初始十六极形变 $\beta_4$ | `0.0` | 任意实数 | `adjust_basis`, `compute_PES`, `execute_HFBTHO`, `start` |

### `&HFBTHO_ITERATIONS` — 迭代控制

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `number_iterations` | `Integer` | 最大自洽迭代次数 | `100` | >0 | `execute_HFBTHO`, `iter` |
| `accuracy` | `Real` | 能量收敛判据 $\epsilon_{\text{conv}}$ | `1.D-5` | >0 | `ALambda`, `HFBTHO_DFT_SOLVER`, `check_consistency`, `iter` |
| `restart_file` | `Integer` | 初始波函数来源 | `-1` | `1`=球形；`2`=长球形；`3`=扁球形；负值=二进制重启文件 | `check_consistency`, `execute_HFBTHO`, `inout`, `start` |

### `&HFBTHO_FUNCTIONAL` — 能量密度泛函与库仑

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `functional` | `Character(30)` | 能量泛函 / Skyrme / Gogny 名称 | `'SLY4'` | 见下方预设泛函列表 | `HFBTHO_DFT_SOLVER`, `check_consistency`, `set_functional_parameters`, `preparer` |
| `add_initial_pairing` | `Logical` | 从重启文件读取时是否叠加初始配对场 | `.False.` | `T` / `F` | `execute_HFBTHO`, `initialize_HFBTHO_SOLVER` |
| `type_of_coulomb` | `Integer` | 库仑相互作用处理方式 | `2` | 见下方库仑选项 | `check_consistency`, `execute_HFBTHO`, `coulom`, `field`, `preparer` |
| `include_3N_force` | `Logical` | 是否包含三体力（特定 DME 泛函） | `.False.` | `T` / `F` | `execute_HFBTHO`, `set_functional_parameters` |

### `&HFBTHO_PAIRING` — 配对参数

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `user_pairing` | `Logical` | 是否使用用户自定义配对强度 | `.False.` | `T` / `F` | `execute_HFBTHO`, `initialize_HFBTHO_SOLVER` |
| `vpair_n` | `Real` | 中子配对强度 $V_0^{(n)}$（MeV） | `-300.0` | 任意实数 | `execute_HFBTHO`, `initialize_HFBTHO_SOLVER` |
| `vpair_p` | `Real` | 质子配对强度 $V_0^{(p)}$（MeV） | `-300.0` | 任意实数 | `execute_HFBTHO`, `initialize_HFBTHO_SOLVER` |
| `pairing_cutoff` | `Real` | 单粒子能截断（MeV） | `60.0` | >=0 | `check_consistency`, `execute_HFBTHO`, `hfbdiag`, `field` |
| `pairing_feature` | `Real` | 体积-表面配对混合参数 | `0.5` | [0.0, 1.0] | `check_consistency`, `execute_HFBTHO`, `initialize_HFBTHO_SOLVER` |

### `&HFBTHO_CONSTRAINTS` — 多极约束

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `lambda_values(1:8)` | `Integer(8)` | 约束多极阶数 $\lambda$ | `0×8` | 1~8 | `HFBTHO_DFT_SOLVER`, `Constraint_or_not`, `preparer` |
| `lambda_active(1:8)` | `Integer(8)` | 各阶约束激活标志 | `0×8` | `0`=关闭；`>0`=激活 | `Constraint_or_not`, `HFBTHO_DFT_SOLVER`, `Main_Program`, `getLagrange` |
| `expectation_values(1:8)` | `Real(8)` | 约束期望值 $\langle Q_\lambda \rangle$ | `0.0×8` | 任意实数 | `HFBTHO_DFT_SOLVER`, `Main_Program`, `compute_PES`, `getLagrange`, `expect` |

### `&HFBTHO_BLOCKING` — 奇核子阻塞

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `proton_blocking(1:5)` | `Integer(5)` | 质子阻塞态量子数 | `0×5` | $(n_r, n_z, m_l, \pm, k)$ | `execute_HFBTHO`, `handle_blocking`, `requested_blocked_level` |
| `neutron_blocking(1:5)` | `Integer(5)` | 中子阻塞态量子数 | `0×5` | $(n_r, n_z, m_l, \pm, k)$ | `execute_HFBTHO`, `handle_blocking`, `requested_blocked_level` |

### `&HFBTHO_PROJECTION` — 投影到不同核（旧版接口）

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `switch_to_THO` | `Integer` | THO 基模式 | `0` | `0`=HO；`-1`=HO→THO；`1`=THO | `check_consistency`, `execute_HFBTHO`, `HFBTHO_DFT_SOLVER`, `coordinateLST` |
| `projection_is_on` | `Integer` | 是否投影到不同核子数的核 | `0` | `0`=关闭；非零=开启 | `execute_HFBTHO`, `check_consistency`, `preparer` |
| `gauge_points` | `Integer` | PNP 规范角积分点数 | `1` | >=1 | `check_consistency`, `execute_HFBTHO`, `preparer`, `HFBTHO_PNP` |
| `delta_Z` | `Integer` | 投影目标质子数偏差 $\Delta Z$ | `0` | 整数（$Z+\Delta Z \ge 1$） | `check_consistency`, `execute_HFBTHO`, `preparer` |
| `delta_N` | `Integer` | 投影目标中子数偏差 $\Delta N$ | `0` | 整数（$N+\Delta N \ge 1$） | `check_consistency`, `execute_HFBTHO`, `preparer` |

### `&HFBTHO_TEMPERATURE` — 有限温度

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `set_temperature` | `Logical` | 是否启用有限温度计算 | `.False.` | `T` / `F` | `execute_HFBTHO`, `iter`, `densit` |
| `temperature` | `Real` | 温度 $T$（MeV） | `0.0` | >=0 | `check_consistency`, `densit`, `execute_HFBTHO`, `iter`, `hfbdiag` |

### `&HFBTHO_FEATURES` — 计算特性开关

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `collective_inertia` | `Logical` | 计算集体惯性质量与零点能 | `.False.` | `T` / `F` | `execute_HFBTHO`, `expect`, `preparer` |
| `fission_fragments` | `Logical` | 计算裂变碎片特征 | `.False.` | `T` / `F` | `execute_HFBTHO`, `expect`, `preparer` |
| `pairing_regularization` | `Logical` | 激活配对正则化 | `.False.` | `T` / `F` | `execute_HFBTHO`, `expect`, `field`, `preparer` |
| `automatic_basis` | `Logical` | 自动优化基参数 | `.False.` | `T` / `F` | `Main_Program`, `compute_PES`, `execute_HFBTHO`, `adjust_basis` |
| `localization_functions` | `Logical` | 计算局域化函数 | `.False.` | `T` / `F` | `execute_HFBTHO`, `localization` module |
| `canonical_wavefunctions` | `Logical` | 输出正则波函数 | `.False.` | `T` / `F` | `resu`, `canonical` module |

### `&HFBTHO_NECK` — 颈部约束

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `set_neck_constrain` | `Logical` | 是否激活颈部约束 | `.False.` | `T` / `F` | `execute_HFBTHO`, `Constraint_or_not`, `read_data` |
| `neck_value` | `Real` | 颈部约束期望值 $Q_{\text{neck}}$ | `0.5` | 任意实数 | `execute_HFBTHO`, `getLagrange`, `expect` |

### `&HFBTHO_DEBUG` — 调试与数值参数

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `number_Gauss` | `Integer` | $z$ 方向 Gauss-Hermite 积分点数 | `40` | [1, 100] | `check_consistency`, `execute_HFBTHO`, `thoalloc`, `gausspoints` |
| `number_Laguerre` | `Integer` | $\rho$ 方向 Gauss-Laguerre 积分点数 | `40` | [1, 100] | `check_consistency`, `execute_HFBTHO`, `thoalloc`, `gausspoints` |
| `number_Legendre` | `Integer` | 库仑/PNP Gauss-Legendre 积分点数 | `80` | $(-\infty, 100]$ | `check_consistency`, `execute_HFBTHO`, `preparer` |
| `compatibility_HFODD` | `Logical` | 强制使用与 HFODD 相同的基函数排序 | `.False.` | `T` / `F` | `execute_HFBTHO`, `base0`, `base` |
| `number_states` | `Integer` | 基函数中输出/保存的最大状态数 | `500` | >0 | `check_consistency`, `execute_HFBTHO`, `base0` |
| `force_parity` | `Logical` | 强制反射对称性 | `.True.` | `T` / `F` | `execute_HFBTHO`, `thoalloc`, `preparer` |
| `write_hel` | `Logical` | 是否写入二进制重启文件 | `.True.` | `T` / `F` | `resu`, `inout`, `write_data` |
| `print_time` | `Integer` | 屏幕输出详细程度 / 计时级别 | `0` | `0`=最少；更大则更详细 | `execute_HFBTHO`, `preparer`, `heading` |

### `&HFBTHO_RESTORATION` — 对称性恢复（PNP / AMP）

| 变量名 | 类型 | 物理含义 | 默认值 | 取值选项 / 范围 | 关键调用函数 |
|--------|------|----------|--------|-----------------|--------------|
| `PNP_is_on` | `Integer` | 粒子数投影开关 | `0` | `0`=关闭；`1`=开启 | `calculate_densities`, `calculate_energies`, `calculate_overlaps`, `execute_HFBTHO` |
| `number_of_gauge_points` | `Integer` | PNP 规范角积分点数 | `1` | >=1 | `initialize_projections`, `preparer` |
| `delta_neutrons` | `Integer` | PNP 中子数偏差 $\Delta N$ | `0` | 任意整数 | `initialize_projections`, `preparer` |
| `delta_protons` | `Integer` | PNP 质子数偏差 $\Delta Z$ | `0` | 任意整数 | `initialize_projections`, `preparer` |
| `AMP_is_on` | `Integer` | 角动量投影开关 | `0` | `0`=关闭；`1`=开启 | `calculate_overlaps`, `execute_HFBTHO`, `initialize_projections`, `preparer` |
| `number_of_rotational_angles` | `Integer` | AMP 欧拉角积分点数 | `1` | >=1 | `initialize_projections`, `preparer` |
| `maximal_angular_momentum` | `Integer` | AMP 投影最大角动量 $J_{\max}$ | `0` | >=0 | `initialize_projections`, `preparer` |

#### 预设泛函列表（`functional`）

| 类别 | 可选值 |
|------|--------|
| **Skyrme** | `SIII`, `SKM*`, `SKP`, `SLY4`, `SLY5`, `SLY6`, `SLY7`, `SKI3`, `SKO`, `SKX`, `SKM*mod`, `SKOP`, `SLY4mod` |
| **UNEDF** | `UNE0`, `UNE1`, `UNE2` |
| **DME** | `N0LO`, `N1LO`, `N2LO`, `DME_LO`, `DME_NLO`, `DME_N2LO`, `DME_NLOD`, `DME_N2LOD` |
| **Reg.优化** | `FITS`, `REG_LO`, `REG_NLO`, `REG_N2LO`, `REG_NLOD`, `REG_N2LOD` |
| **Gogny** | `D1`, `D1p`, `D1S`, `D1N`, `T0X0` |
| **NEDF** | `NEDF`, `SeaLL1`, `NEDF1`, `NEDF2`, `NEDF3`, `NEDF4`, `NEDF5` |
| **其他** | `HFB1`, `READ`（从二进制文件读取泛函参数）, `RAND` |

#### 库仑选项（`type_of_coulomb`）

| 值 | 含义 |
|----|------|
| `0` | 完全关闭库仑力 |
| `1` | 直接库仑用 substitution 方法 |
| `2` | 直接库仑用 substitution 方法，交换库仑用 Slater 近似 |
| `-1` | 直接库仑力仅由高斯求和计算 |
| `-2` | 直接库仑用 substitution 方法，交换库仑由高斯求和计算 |
| `-3` | 直接库仑由高斯求和计算，交换库仑用 Slater 近似 |
| `-4` | 直接+交换库仑均由高斯求和计算 |
| `-5` | 直接+交换+配对库仑均由高斯求和计算 |

> **注意**：`HFBTHO_RESTORATION` 是现代 HFBTHO 版本中推荐使用的 PNP/AMP 控制接口；`HFBTHO_PROJECTION` 保留了向后兼容性，用于投影到不同核子数（`delta_Z`, `delta_N`）。两者可同时存在于输入文件中。

---

## 变量快速参考索引

| Namelist | 核心物理量 |
|----------|-----------|
| `HFBTHO_GENERAL` | 核子数、基空间截断、计算类型 |
| `HFBTHO_INITIAL` | 初始 Woods-Saxon 形变 |
| `HFBTHO_ITERATIONS` | 收敛标准、初始波函数来源 |
| `HFBTHO_FUNCTIONAL` | EDF/Skyrme/Gogny 泛函、库仑处理 |
| `HFBTHO_PAIRING` | 配对强度、截断、类型 |
| `HFBTHO_CONSTRAINTS` | 多极约束阶数与期望值 |
| `HFBTHO_BLOCKING` | 奇核子阻塞量子数 |
| `HFBTHO_PROJECTION` | THO 模式、投影到不同核 |
| `HFBTHO_TEMPERATURE` | 有限温度开关与温度值 |
| `HFBTHO_FEATURES` | 集体惯性、裂变碎片、正则化等 |
| `HFBTHO_NECK` | 颈部约束 |
| `HFBTHO_DEBUG` | 积分网格、反射对称、输出控制 |
| `HFBTHO_RESTORATION` | PNP/AMP 对称性恢复参数 |
