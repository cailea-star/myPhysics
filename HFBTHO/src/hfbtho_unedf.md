# 代码文件 `hfbtho_unedf.f90` 说明

本模块定义 Skyrme/UNEDF 能量密度泛函（EDF）的参数、耦合常数以及不同参数化形式之间的转换。同时计算 $U$ 振幅（能量密度对密度双线性组合的偏导数），用于在求解器中构建 HFB 平均场。

---

### `calculate_U_parameters`

##### 函数定义

- 函数位置: `198-442 lines`

```fortran
Subroutine calculate_U_parameters(rho0_in,rho1_in,tau0_in,tau1_in,laprho0,laprho1,nablarho0s,nablarho1s)
```

##### 物理公式

- 物理含义: 计算 $U$ 振幅——能量密度 $\mathcal{H}$ 对各类局部密度双线性组合（$\rho\rho$、$\rho\tau$、$\rho\Delta\rho$、$JJ$ 等）的偏导数。结果存储在公共数组 `Urhorho`、`Urhotau`、`UrhoDrho`、`UJJ`、`UJabJba`、`UJnablarho`、`UrhonablaJ`、`Unablarho` 以及配对振幅数组 `Urhorhopr` 中。

**2N（类 Skyrme）项**: 对于同位标量 ($t=0$, $\phi_t=+1$) 和同位矢量 ($t=1$, $\phi_t=-1$),

$$
U^{\rho\rho}_{t,0} = C^{\rho\rho}_t + C^{\rho\rho\sigma}_t \rho_0^{\sigma} + \frac{1}{2}(a_{\rho\rho} + \phi_t b_{\rho\rho})\cdot\text{mevfm}
$$

$$
U^{\rho\tau}_{t,0} = C^{\rho\tau}_t + \frac{1}{2}(a_{\rho\tau} + \phi_t b_{\rho\tau})\cdot\text{mevfm}
$$

$$
U^{\rho\Delta\rho}_{t,0} = C^{\rho\Delta\rho}_t + 4(\alpha^2-\alpha+\tfrac12)\cdot\frac{1}{2}(a_{\rho\Delta\rho} + \phi_t b_{\rho\Delta\rho})\cdot\text{mevfm}
$$

$$
U^{JJ}_{t,0} = C^{JJ}_t + \frac{1}{2}(a_{JJ} + \phi_t b_{JJ})\cdot\text{mevfm}
$$

$$
U^{\nabla\rho}_{t,0} = C^{\nabla\rho}_t + \frac{1}{2}(a_{\nabla\rho} + \phi_t b_{\nabla\rho})\cdot\text{mevfm}
$$

**密度导数** ($j=1\dots7$ 覆盖 $\partial\rho_0$、$\partial\rho_1$、$\partial^2\rho_0^2$、$\partial^2\rho_1^2$、$\partial^2\rho_0\rho_1$、$\partial\tau_0$、$\partial\Delta\rho_0$):

$$
U^{\rho\rho}_{t,1} = \frac{\sigma C^{\rho\rho\sigma}_t \rho_0^{\sigma}}{\rho_0+\epsilon} + \frac{1}{2}(\partial a_{\rho\rho} + \phi_t \partial b_{\rho\rho})\partial u \cdot\text{mevfm}
$$

**配对 (pp) 振幅:**

$$
U^{\rho\rho,\text{pp}}_{0,0} = \frac{1}{16}\bigl[C^{pV0}_0(1-C^{pV1}_0\rho_0/0.16) + C^{pV0}_1(1-C^{pV1}_1\rho_0/0.16)\bigr]
$$

**3N (DME) 项:** 若 `use_DME3N_terms=.True.`，累加如下贡献

$$
U^{\rho\rho}_{0,0} \mathrel{+}= h_{\rho_0\rho_0}\,\rho_0, \qquad
U^{\rho\tau}_{0,0} \mathrel{+}= h_{\rho_0\tau_0}\,\rho_0, \qquad
U^{JJ}_{0,0} \mathrel{+}= h_{J_0J_0}\,\rho_0
$$

以及相应的一阶/二阶密度导数。

**NEDF 分支:** 若 `is_nedf=.True.`，$U^{\rho\rho}$ 振幅被覆盖为

$$
U^{\rho\rho}_{t,0} = \frac{a_t}{\rho_0^{1/3}} + b_t + c_t\,\rho_0^{1/3}, \qquad t=0,1
$$

$$
U^{\rho\rho}_{2,0} = \Bigl(\frac{a_2}{\rho_0^{10/3}} + \frac{b_2}{\rho_0^{3}} + \frac{c_2}{\rho_0^{8/3}}\Bigr)\rho_1^{3}
$$

并附带对 $\rho_0$ 和 $\rho_1$ 的解析导数。

- 实现步骤:
  1. 将输入密度存入模块级数组 `rho(0:1)`、`tau(0:1)` 等。
  2. 将所有 $a$、$b$、$h$ 辅助系数及所有 $U$ 数组清零。
  3. 在 $t=0,1$ 循环中计算 2N $U$ 振幅及其导数。
  4. 计算配对振幅 `Urhorhopr`。
  5. 若 `use_DME3N_terms`，向 $U$ 及导数添加 3N/DME 贡献。
  6. 若 `is_nedf`，用 NEDF 解析公式覆盖 $U^{\rho\rho}$。
  7. 若 `.Not.use_j2terms`，将 $U^{JJ}$ 和 $U^{J\bar{J}}$ 置零。

##### 变量对应

- 输入变量 `rho0_in (real(pr))`: 公式对应 $\rho_0$（同位标量密度）
- 输入变量 `rho1_in (real(pr))`: 公式对应 $\rho_1$（同位矢量密度）
- 输入变量 `tau0_in (real(pr))`: 公式对应 $\tau_0$（同位标量动能密度）
- 输入变量 `tau1_in (real(pr))`: 公式对应 $\tau_1$（同位矢量动能密度）
- 输入变量 `laprho0 (real(pr), optional)`: 公式对应 $\Delta\rho_0$
- 输入变量 `laprho1 (real(pr), optional)`: 公式对应 $\Delta\rho_1$
- 输入变量 `nablarho0s (real(pr), optional)`: 公式对应 $(\nabla\rho_0)^2$
- 输入变量 `nablarho1s (real(pr), optional)`: 公式对应 $(\nabla\rho_1)^2$
- 输出数组 `Urhorho(0:3,0:7)`, `Urhotau(0:3,0:7)`, `UrhoDrho(0:3,0:7)`, `UJJ(0:3,0:7)` 等: 公式对应 $U^{\cdot}_{t,t'}$ 及其各阶导数

##### 函数数据流

- 调用情况: `hfbtho_solver::field`, `hfbtho_solver::expect`, `calculate_C_from_NM`, `find_NM_RHOC`, `calculate_NM_properties` `<- calculate_U_parameters`

---

### `default_UNEDF_NAMELIST`

##### 函数定义

- 函数位置: `446-509 lines`

```fortran
Subroutine default_UNEDF_NAMELIST()
```

##### 物理公式

- 物理含义: 设置 UNEDF namelist 的所有默认参数值（默认采用 UNEDF1 参数集）。

##### 变量对应

- 无输入参数
- 输出: 设置模块级公共变量默认值，包括 `Crho`、`Cdrho`、`Ctau`、`CrDr`、`CrdJ`、`CJ`、`CpV0`、`CpV1`、`sigma`、`hbzero`、`e2charg`、`E_NM`、`K_NM`、`RHO_NM` 等

##### 函数数据流

- 调用情况: `read_UNEDF_NAMELIST <- default_UNEDF_NAMELIST`

---

### `read_UNEDF_NAMELIST`

##### 函数定义

- 函数位置: `513-1189 lines`

```fortran
Subroutine read_UNEDF_NAMELIST(fname,noForces,filename_nml)
```

##### 物理公式

- 物理含义: 读取功能名称并设置耦合常数。支持大量预定义泛函：UNE0、UNE1、UNE2、HFB1、N0LO、N1LO、N2LO、NEDF、SeaLL1、NEDF1–5、FITS、READ、RAND 等；也可通过 `UNEDF_NAMELIST.dat` 文件自定义读取。若 `fname` 匹配传统 Skyrme 力名称，则调用 `skforce`。

- 实现步骤:
  1. 调用 `default_UNEDF_NAMELIST` 重置所有参数。
  2. 调用 `skforce` 尝试匹配传统 Skyrme/Gogny 力名称。
  3. 若未匹配，通过 `Select Case` 根据功能名称硬编码各泛函参数。
  4. 若 `use_Namelist=.True.`，从 `UNEDF_NAMELIST.dat` 读取 namelist 并校验。

##### 变量对应

- 输入变量 `fname (character(len=30), inout)`: 公式对应功能名称（如 'UNE1'、'SLY4'）
- 输出变量 `noForces (integer(ipr))`: 若成功设置力参数则为 1，否则为 0
- 输入变量 `filename_nml (character(len=256), optional)`: namelist 文件路径，默认 'UNEDF_NAMELIST.dat'

##### 函数数据流

- 调用情况: `hfbtho_solver::Main_Program <- read_UNEDF_NAMELIST <- default_UNEDF_NAMELIST, skforce`

---

### `gogny_force`

##### 函数定义

- 函数位置: `1193-1262 lines`

```fortran
Subroutine gogny_force(fname)
```

##### 物理公式

- 物理含义: 加载 Gogny 有限程核力的参数集（D1、D1S、D1p、D1N）。每组参数包含两个高斯项的力程 $\mu_g$ 及交换系数 $W_g$、$B_g$、$H_g$、$M_g$。

##### 变量对应

- 输入变量 `fname (character(len=30), inout)`: Gogny 力名称
- 输出: 分配并设置模块级数组 `mu_g`、`W_g`、`B_g`、`H_g`、`M_g`

##### 函数数据流

- 调用情况: `skforce <- gogny_force`

---

### `set_all_gaussians`

##### 函数定义

- 函数位置: `1266-1282 lines`

```fortran
Subroutine set_all_gaussians(icoul)
```

##### 物理公式

- 物理含义: 将 Gogny 核力的高斯力程参数与库仑相互作用的高斯展开力程合并到统一数组 `mu_g_all` 中，供后续有限程计算使用。

##### 变量对应

- 输入变量 `icoul (integer(ipr))`: 库仑处理方式标记；`icoul < 0` 时启用高斯库仑展开
- 输出: 分配并设置 `mu_g_all(1:n_g_all)`

##### 函数数据流

- 调用情况: `hfbtho_solver::Main_Program <- set_all_gaussians`

---

### `skforce`

##### 函数定义

- 函数位置: `1287-1679 lines`

```fortran
Subroutine skforce(fname,noForces)
```

##### 物理公式

- 物理含义: 为传统 Skyrme 核力（SIII、SKM*、SKP、SLy4–7、SKI3、SKO、SKO'、SKX、T0X0 等）和 Gogny 核力（D1、D1S、D1p、D1N）设置参数。对 Skyrme 力设置 $(t_0,t_1,t_2,t_3,x_0,x_1,x_2,x_3,\sigma)$ 等，并调用 `C_from_t()` 转换为 $C$ 耦合常数。同时预计算频繁出现的组合量 `tv1`–`tv6`、`ts1`、`ts2`、`t0s`、`t0a`、`drs`、`dra` 等。

- 实现步骤:
  1. 根据 `fname` 匹配具体力参数，通过 `Select Case` 硬编码各 Skyrme/Gogny 参数。
  2. 对 Gogny 力调用 `gogny_force`。
  3. 调用 `C_from_t` 将 $(t,x)$ 参数转换为 $C$ 耦合常数。
  4. 计算能量和势中频繁出现的组合量。

##### 变量对应

- 输入变量 `fname (character(len=30), inout)`: 力名称
- 输出变量 `noForces (integer(ipr))`: 若成功则为 1，否则为 0
- 输出: 设置模块级私有变量 `t0`–`t3`、`x0`–`x3`、`sigma`、`hbzero` 及组合量 `tv1`–`tv6`、`ts1`、`ts2` 等

##### 函数数据流

- 调用情况: `read_UNEDF_NAMELIST <- skforce <- gogny_force, C_from_t`

---

### `set_functional_parameters`

##### 函数定义

- 函数位置: `1683-1710 lines`

```fortran
Subroutine set_functional_parameters(fname,npr)
```

##### 物理公式

- 物理含义: EDF 参数主设置入口。若 `use_INM=.True.` 则调用 `calculate_C_from_NM` 从核物质性质反解体积耦合常数；否则直接使用已读取的耦合常数。随后计算核物质性质和自然单位。

- 实现步骤:
  1. 若 `use_INM`，调用 `calculate_C_from_NM`。
  2. 调用 `calculate_NM_properties` 校验核物质性质。
  3. 调用 `calculate_natural_units` 计算自然单位。

##### 变量对应

- 输入变量 `fname (character(len=30), inout)`: 功能名称
- 输入变量 `npr (integer(ipr), dimension(3))`: 粒子数数组（1=中子，2=质子，3=总和）

##### 函数数据流

- 调用情况: `hfbtho_solver::Main_Program <- set_functional_parameters <- calculate_C_from_NM, calculate_NM_properties, calculate_natural_units`

---

### `print_functional_parameters`

##### 函数定义

- 函数位置: `1714-1844 lines`

```fortran
Subroutine print_functional_parameters(fname)
```

##### 物理公式

- 物理含义: 将所有功能参数格式化输出到日志文件。输出内容包括 $C$ 耦合常数、$(t,x)$ 参数（如适用）、Gogny 参数、DME 参数、自然单位制下的耦合常数、以及核物质性质（$E_{\text{NM}}$、$K$、$\rho_{\text{NM}}$、$a_{\text{sym}}$、$L$ 等）。

##### 变量对应

- 输入变量 `fname (character(len=30), optional)`: 功能名称，若未提供则使用 `FunctionalName`

##### 函数数据流

- 调用情况: `hfbtho_solver::Main_Program <- print_functional_parameters <- t_from_C`

---

### `calculate_natural_units`

##### 函数定义

- 函数位置: `1848-1859 lines`

```fortran
Subroutine calculate_natural_units()
```

##### 物理公式

- 物理含义: 将耦合常数转换为自然单位，以 $f_\pi$ 和 $\Lambda$ 为基准进行量纲缩放。

$$
\nu C^{\rho\rho} = C^{\rho\rho} \frac{f_\pi^2}{\text{mevfm}^3}, \quad
\nu C^{\rho\tau} = C^{\rho\tau} \frac{(f_\pi\Lambda)^2}{\text{mevfm}^5}
$$

##### 变量对应

- 无输入参数
- 输出: 设置模块级私有数组 `nuCrho`、`nuCtau`、`nuCrDr`、`nuCrdJ`、`nuCJ`、`nuCpV0`、`nuCpV1` 等

##### 函数数据流

- 调用情况: `set_functional_parameters <- calculate_natural_units`

---

### `calculate_C_from_NM`

##### 函数定义

- 函数位置: `1873-1954 lines`

```fortran
Subroutine calculate_C_from_NM(E,K,SMASS,RHO,ASS,LASS,VMASS,npr,sigma_NM)
```

##### 物理公式

- 物理含义: 从无限核物质（INM）性质反解体积耦合常数 $C^{\rho\rho}_t$、$C^{\rho\rho\sigma}_t$、$C^{\rho\tau}_t$ 以及密度依赖指数 $\sigma$。输入的核物质性质包括能量、不可压缩系数、标量有效质量、饱和密度、对称能、对称能斜率、矢量有效质量。若存在有限程贡献（Gogny），通过 `calculate_finite_range_NM` 进行扣除。

**同位标量耦合常数:**

$$
C^{\rho\rho}_0 = \frac{1}{3\sigma\rho}\bigl[\tau_c\hbar^2(-3+(2-3\sigma)(m_s^*-1)) + 3(1+\sigma)E + 3P/\rho - 3\sigma\rho a_{\rho_0\rho_0} + \dots\bigr]
$$

$$
C^{\rho\rho\sigma}_0 = \frac{\rho^{-1-\sigma}}{3\sigma}\bigl[\tau_c\hbar^2(3-2(m_s^*-1)) - 3(E+P/\rho) - 3\rho^2 h_{\rho_0\rho_0} - 3\tau_c\rho^2 h_{\rho_0\tau_0} - \dots\bigr]
$$

$$
C^{\rho\tau}_0 = \frac{\hbar^2(m_s^*-1) - \rho(a_{\rho_0\tau_0}+\rho h_{\rho_0\tau_0})}{\rho}
$$

**同位矢量耦合常数**（类似长公式，代码中显式给出）。

- 实现步骤:
  1. 在核物质密度下调用 `calculate_U_parameters` 获取 $a$、$b$、$h$ 系数。
  2. 若 `finite_range`，调用 `calculate_finite_range_NM` 扣除有限程贡献。
  3. 计算 $a$、$b$、$h$ 系数的同位标量/矢量组合。
  4. 计算密度依赖指数 $\sigma$（若未显式提供）。
  5. 计算同位标量和同位矢量体积耦合常数。

##### 变量对应

- 输入变量 `E (real(pr))`: 公式对应 $E_{\text{NM}}$（每核子能量）
- 输入变量 `K (real(pr))`: 公式对应 $K_{\text{NM}}$（不可压缩系数）
- 输入变量 `SMASS (real(pr))`: 公式对应 $1/m_s^*$（标量有效质量倒数）
- 输入变量 `RHO (real(pr))`: 公式对应 $\rho_{\text{NM}}$（饱和密度）
- 输入变量 `ASS (real(pr))`: 公式对应 $a_{\text{sym}}$（对称能）
- 输入变量 `LASS (real(pr))`: 公式对应 $L_{\text{sym}}$（对称能斜率）
- 输入变量 `VMASS (real(pr))`: 公式对应 $1/m_v^*$（矢量有效质量倒数）
- 输入变量 `npr (integer(ipr), dimension(3))`: 粒子数数组
- 输入变量 `sigma_NM (real(pr), optional)`: 可选的密度依赖指数
- 输出: 设置模块级公共数组 `Crho`、`Cdrho`、`Ctau` 以及 `sigma`

##### 函数数据流

- 调用情况: `set_functional_parameters <- calculate_C_from_NM <- calculate_U_parameters, calculate_finite_range_NM`

---

### `calculate_finite_range_NM`

##### 函数定义

- 函数位置: `2024-2076 lines`

```fortran
Subroutine calculate_finite_range_NM(rho,E,P,K,SMASS,ASS,LASS,KASS,VMASS)
```

##### 物理公式

- 物理含义: 计算 Gogny 有限程力对核物质性质的贡献。采用 Sellahewa & Bender (2014) 的解析公式，对每个高斯项 $i$ 求和。

$$
E^{\text{NM}} = \frac{1}{2}\sum_i\bigl[A_0^i\rho + B_0^i g_0(\mu_i k_F)\bigr], \quad
P^{\text{NM}} = \sum_i\bigl[\tfrac{1}{2}A_0^i\rho^2 + B_0^i p_0(\mu_i k_F)\rho\bigr]
$$

$$
K^{\text{NM}} = -3\sum_i B_0^i k_0(\mu_i k_F), \quad
M_s = \frac{2m_N}{\hbar^2}\frac{1}{2}\sum_i B_0^i m_0(\mu_i k_F)\frac{3}{5k_F^2}
$$

其中 $A_0^i = \frac{\pi^{3/2}\mu_i^3}{4}(4W_i+2B_i-2H_i-M_i)$，$B_0^i = B_{nn}^i + B_{np}^i$，$k_F = (3\pi^2\rho/2)^{1/3}$。函数 $g_0$、$p_0$、$k_0$、$m_0$、$s_1$、$s_2$、$l_1$、$l_2$ 等为解析定义的辅助函数。

##### 变量对应

- 输入变量 `rho (real(pr))`: 公式对应 $\rho$（核物质密度）
- 输出变量 `E (real(pr))`: 公式对应 $E^{\text{NM}}$（每核子能量）
- 输出变量 `P (real(pr))`: 公式对应 $P^{\text{NM}}$（压强）
- 输出变量 `K (real(pr))`: 公式对应 $K^{\text{NM}}$（不可压缩系数）
- 输出变量 `SMASS (real(pr))`: 公式对应 $1/m_s^*$
- 输出变量 `ASS (real(pr))`: 公式对应 $a_{\text{sym}}$
- 输出变量 `LASS (real(pr))`: 公式对应 $L_{\text{sym}}$
- 输出变量 `KASS (real(pr))`: 公式对应 $\Delta K^{\text{NM}}$
- 输出变量 `VMASS (real(pr))`: 公式对应 $1/m_v^*$

##### 函数数据流

- 调用情况: `calculate_C_from_NM`, `calculate_NM_properties <- calculate_finite_range_NM`

---

### `P_SNM_FR`

##### 函数定义

- 函数位置: `2083-2109 lines`

```fortran
Function P_SNM_FR(rho) result(P)
```

##### 物理公式

- 物理含义: 计算有限程力对对称核物质压强的贡献。

$$
P^{\text{FR}}(\rho) = \sum_{i=1}^{n_g}\Bigl\{\frac{1}{2}A_0^i\rho^2 + B_0^i p_0(\mu_i k_F)\rho\Bigr\}
$$

##### 变量对应

- 输入变量 `rho (real(pr))`: 公式对应 $\rho$（核物质密度）
- 输出变量 `P (real(pr))`: 公式对应 $P^{\text{FR}}$（有限程压强贡献）

##### 函数数据流

- 调用情况: `find_NM_RHOC <- P_SNM_FR`

---

### `calculate_NM_properties`

##### 函数定义

- 函数位置: `2113-2179 lines`

```fortran
Subroutine calculate_NM_properties()
```

##### 物理公式

- 物理含义: 基于当前耦合常数计算核物质性质，包括对称核物质（SNM）的能量 $E_{\text{NM}}$、压强 $P_{\text{NM}}$、不可压缩系数 $K_{\text{NM}}$、标量有效质量 $m_s^*$，以及非对称核物质（ANM）的对称能 $a_{\text{sym}}$、斜率 $L_{\text{sym}}$、曲率 $K_a$、矢量有效质量 $m_v^*$。首先通过 `find_NM_RHOC` 确定饱和密度。

**对称核物质:**

$$
E_{\text{NM}} = \tau_c\hbar^2 + \rho(a_{\rho_0\rho_0}+\rho h_{\rho_0\rho_0}+C^{\rho\rho}_0+\rho^{\sigma}C^{\rho\rho\sigma}_0) + \tau_c\rho(a_{\rho_0\tau_0}+\rho h_{\rho_0\tau_0}+C^{\rho\tau}_0) + E_{\text{FR}}
$$

##### 变量对应

- 无输入参数（使用模块级耦合常数）
- 输出: 设置模块级公共变量 `E_NM`、`P_NM`、`K_NM`、`SMASS_NM`、`ASS_NM`、`LASS_NM`、`KA_NM`、`VMASS_NM`

##### 函数数据流

- 调用情况: `set_functional_parameters <- calculate_NM_properties <- find_NM_RHOC, calculate_finite_range_NM, calculate_U_parameters`

---

### `find_NM_RHOC`

##### 函数定义

- 函数位置: `2183-2232 lines`

```fortran
Real(pr) Function find_NM_RHOC()
```

##### 物理公式

- 物理含义: 用割线法（secant method）求解压强为零的密度，即饱和密度 $\rho_{\text{NM}}$。

$$
P(\rho) = \frac{1}{3}\rho^2\Bigl[\frac{2\tau_c\hbar^2}{\rho} + 3a_{\rho_0\rho_0} + 5\tau_c a_{\rho_0\tau_0} + 6\rho h_{\rho_0\rho_0} + 8\tau_c\rho h_{\rho_0\tau_0} + 3C^{\rho\rho}_0 + 3(1+\sigma)\rho^{\sigma}C^{\rho\rho\sigma}_0 + \dots\Bigr] + P_{\text{FR}} = 0
$$

- 实现步骤:
  1. 以 $\rho=0.170$ 为初始猜测，调用 `calculate_U_parameters` 计算压强。
  2. 用割线法迭代：$\rho_{n+1} = \rho_n - P_n(\rho_n-\rho_{n-1})/(P_n-P_{n-1})$。
  3. 当步长 $|{\rm step}| < 100\epsilon$ 时收敛。

##### 变量对应

- 输出: 函数返回值即为 $\rho_{\text{NM}}$（饱和密度）

##### 函数数据流

- 调用情况: `calculate_NM_properties <- find_NM_RHOC <- calculate_U_parameters, P_SNM_FR`

---

### `C_from_t`

##### 函数定义

- 函数位置: `2275-2311 lines`

```fortran
Subroutine C_from_t()
```

##### 物理公式

- 物理含义: 将传统 Skyrme $(t,x)$ 参数化转换为能量密度泛函的 $C$ 耦合常数。

**时-even 耦合常数:**

$$
C^{\rho\rho}_0 = \frac{3}{8}t_0, \quad C^{\rho\rho}_1 = -\frac{1}{4}t_0\bigl(\tfrac{1}{2}+x_0\bigr)
$$

$$
C^{\rho\rho\sigma}_0 = \frac{1}{16}t_3, \quad C^{\rho\rho\sigma}_1 = -\frac{1}{24}t_3\bigl(\tfrac{1}{2}+x_3\bigr)
$$

$$
C^{\rho\tau}_0 = \frac{3}{16}t_1 + \frac{1}{4}t_2\bigl(\tfrac{5}{4}+x_2\bigr), \quad
C^{\rho\Delta\rho}_0 = \frac{1}{16}t_2\bigl(\tfrac{5}{4}+x_2\bigr) - \frac{9}{64}t_1
$$

$$
C^{\rho\nabla J}_0 = -b_4 - \tfrac{1}{2}b_4', \quad C^{JJ}_0 = -\frac{1}{16}\bigl[t_1(2x_1-1) + t_2(2x_2+1) - 5t_e - 15t_o\bigr]
$$

**时-odd 耦合常数**由规范不变性确定，如 $C^{ss}_t$、$C^{s\Delta s}_t$、$C^{sj}_t = -C^{\rho\tau}_t$ 等。

##### 变量对应

- 无输入参数（使用模块级私有变量 `t0`–`t3`、`x0`–`x3`、`b4`、`b4p`、`te`、`to`）
- 输出: 设置模块级公共数组 `Crho`、`Cdrho`、`Ctau`、`CrDr`、`CrdJ`、`CJ`、`CJbar`、`Cs`、`Cds`、`Csj`、`CsDs`、`CsdJ`、`cT` 等

##### 函数数据流

- 调用情况: `skforce <- C_from_t`

---

### `t_from_C`

##### 函数定义

- 函数位置: `2316-2340 lines`

```fortran
Subroutine t_from_C()
```

##### 物理公式

- 物理含义: $C$ 耦合常数到 $(t,x)$ 参数化的逆转换。

$$
t_0 = \frac{8}{3}C^{\rho\rho}_0, \quad
x_0 = -\frac{1}{2}\Bigl(3\frac{C^{\rho\rho}_1}{C^{\rho\rho}_0} + 1\Bigr)
$$

$$
t_1 = \frac{4}{3}\bigl(C^{\rho\tau}_0 - 4C^{\rho\Delta\rho}_0\bigr), \quad
t_3 = 16C^{\rho\rho\sigma}_0
$$

$$
b_4 = C^{\rho\nabla J}_1 - C^{\rho\nabla J}_0, \quad
b_4' = -2C^{\rho\nabla J}_1
$$

$$
t_e = \frac{4}{15}\bigl(3C^{JJ}_0 - 9C^{JJ}_1 - 4C^{\rho\Delta\rho}_0 + 12C^{\rho\Delta\rho}_1 - 2C^{\rho\tau}_0 + 6C^{\rho\tau}_1\bigr)
$$

##### 变量对应

- 无输入参数（使用模块级公共数组 `Crho`、`Cdrho`、`Ctau`、`CrDr`、`CrdJ`、`CJ`）
- 输出: 设置模块级私有变量 `t0`–`t3`、`x0`–`x3`、`b4`、`b4p`、`te`、`to` 以及公共输出变量 `t0_pub`–`to_pub`

##### 函数数据流

- 调用情况: `print_functional_parameters`, `hfbtho_solver <- t_from_C`

---

### `Vexternal`

##### 函数定义

- 函数位置: `2351-2358 lines`

```fortran
Elemental Function Vexternal(t,x,y,z)
```

##### 物理公式

- 物理含义: 定义外部谐振子势场（当前实现返回 0，未启用）。若启用，形式为 $V_{\text{ext}}(\mathbf{r}) = E_0 + \frac{1}{2}\hbar\Omega_z z^2 + \hbar\Omega_\perp \rho^2$。

##### 变量对应

- 输入变量 `t (integer(ipr))`: 同位旋指标（0=同位标量，1=同位矢量）
- 输入变量 `x, y, z (real(pr))`: 笛卡尔坐标
- 输出: 函数返回值 $V_{\text{ext}}$（当前恒为 0）

##### 函数数据流

- 调用情况: `hfbtho_solver::field`, `hfbtho_solver::expect <- Vexternal`
