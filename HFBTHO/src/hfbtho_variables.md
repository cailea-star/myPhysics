# 代码文件 ```hfbtho_variables.f90``` 说明

### `Module HFBTHO`

##### 模块定义

- 模块位置: `80-396 lines`

```f90
Module HFBTHO
```

##### 物理含义

本模块为 HFBTHO 求解器的**全局变量中心**，声明了程序运行所需的全部全局标量、数组、可分配数组及 Namelist 输入接口。所有物理量（密度、势场、能量、形变等）均在此定义，供 `hfbtho_solver`、`hfbtho_unedf`、`hfbtho_projections` 等模块共享。

---

##### 关键变量组

###### 1. HFB 迭代输入变量（`_INI` 后缀）

- 变量位置: `86-94 lines`

这些变量通过 Namelist 读入，用于初始化单次 HFB 计算的基与粒子数配置：

- `n00_INI`: 谐振子主壳层数截断 $N_{\max}$
- `npr_INI(3)`: 质子数 $Z$、中子数 $N$、总核子数 $A$
- `b0_INI, bz_INI, bp_INI`: 谐振子长度参数 $b_0, b_z, b_\perp$
- `q_INI`: 初始四极约束期望值 $q_{20}^{(\text{ini})}$
- `skyrme_INI`: Skyrme/泛函名称标识
- `pwi_INI, V0n_INI, V0p_INI`: 配对窗口（MeV）与初始配对强度 $V_0^{(n)}, V_0^{(p)}$
- `basis_HFODD_INI, Add_Pairing_INI`: 基类型与初始配对开关

###### 2. 回归优化输出变量（能量密度泛函分解）

- 变量位置: `96-98 lines`

用于核质量表拟合时，将总能量 $E_{\text{tot}}$ 按 EDF 分量分解输出：

$$
E_{\text{tot}} = \underbrace{\mathcal{E}_{\rho\rho}}_{\texttt{efit\_rhorho}} + \underbrace{\mathcal{E}_{\rho\tau}}_{\texttt{efit\_rhotau}} + \underbrace{\mathcal{E}_{\rho\Delta\rho}}_{\texttt{efit\_rhoDrho}} + \underbrace{\mathcal{E}_{\rho\nabla J}}_{\texttt{efit\_rhonablaJ}} + \underbrace{\mathcal{E}_{JJ}}_{\texttt{efit\_JJ}} + \underbrace{\mathcal{E}_{\text{pair}}}_{\texttt{efitV0}} + \cdots
$$

###### 3. 形变转换常数

- 变量位置: `108 lines`

在 `hfbtho_solver.f90` 中初始化，用于将多极矩 $Q_{\lambda\mu}$ 转换为无量纲形变参数 $\beta_\lambda$：

$$
\beta_2 = \frac{\text{ffdef}_6 \cdot Q_2}{A \, r_0^2}, \quad \text{ffdef}_6 = \frac{\sqrt{5\pi}}{3}
$$

$$
\beta_4 = \frac{\text{ffdef}_7 \cdot Q_4}{A \, r_0^4}, \quad \text{ffdef}_7 = \frac{\sqrt{\pi}}{4}
$$

其中 $r_0 = 1.2\,\text{fm}$ 为核半径参数。

###### 4. 密度与势场数组（坐标空间网格）

- 变量位置: `172-184 lines`

这些可分配数组存储在圆柱坐标 $(r, z)$ 网格上的微观物理量，是 HFB 能量密度泛函计算的核心数据结构：

| 变量名 | 物理意义 | 公式对应 |
|--------|----------|----------|
| `ro(:,:)` | 局域粒子密度 $\rho_q(\boldsymbol{r})$ | $\rho_q = \sum_k |V_{qk}(\boldsymbol{r})|^2$ |
| `tau(:,:)` | 局域动能密度 $\tau_q(\boldsymbol{r})$ | $\tau_q = \sum_k |\nabla V_{qk}(\boldsymbol{r})|^2$ |
| `dro(:,:)` | 密度梯度相关的量 | — |
| `dj(:,:)` | 自旋-轨道密度 $\boldsymbol{J}_q(\boldsymbol{r})$ | — |
| `aka(:,:)` | 配对密度（异常密度）$\tilde{\rho}_q(\boldsymbol{r})$ | $\tilde{\rho}_q = \sum_k U_{qk}^*(\boldsymbol{r}) V_{qk}(\boldsymbol{r})$ |
| `SZFI(:,:), SFIZ(:,:), SRFI(:,:), SFIR(:,:)` | 自旋密度张量分量 $J_{\mu\nu}$ | — |
| `NABLAR(:,:), NABLAZ(:,:)` | 密度梯度分量 $\nabla_r\rho, \nabla_z\rho$ | — |
| `vhbn(:), vrn(:), vzn(:), vdn(:), vsn(:)` | 中子平均场分量：$V_{\text{HB}}, V_r, V_z, V_{\text{div}}, V_{\text{so}}$ | — |
| `vhbp(:), vrp(:), vzp(:), vdp(:), vsp(:)` | 质子平均场分量（同上）| — |
| `qfield(:,:)` | 约束拉格朗日场 $V_{\text{constr}}^{(\lambda)}(\boldsymbol{r})$ | — |

下标 $q \in \{n, p\}$ 由数组第二维或独立数组（后缀 `n`/`p`）区分。

###### 5. HFB 对角化相关数组

- 变量位置: `194-195 lines`

- `hfb(:,:)`: HFB 准粒子哈密顿矩阵 $\mathcal{H}_{\text{HFB}}$
- `zhfb(:)`: 工作数组（本征矢量）
- `evvk(:)`: 准粒子能量 $E_k$
- `erhfb(:), drhfb(:)`: HFB 能量期望值及其迭代差分

HFB 矩阵形式为：

$$
\mathcal{H}_{\text{HFB}} = \begin{pmatrix} h - \lambda & \Delta \\ -\Delta^* & -(h - \lambda)^* \end{pmatrix}
$$

其中 $h$ 为单粒子哈密顿量，$\Delta$ 为配对场，$\lambda$ 为化学势。

###### 6. Broyden 混合参数

- 变量位置: `214-220 lines`

用于自洽迭代加速：

- `alphamix=0.70_pr`: 线性混合参数 $\alpha_{\text{mix}}$
- `nbroyden=7`: Broyden 历史长度 $m_{\text{Broyd}}$
- `brout(:), brin(:)`: 输出/输入密度向量，满足迭代映射 $\boldsymbol{\rho}^{(n+1)} = \boldsymbol{\rho}^{(n)} + \alpha \, \delta\boldsymbol{\rho}$

###### 7. 温度相关变量

- 变量位置: `167-170 lines`

- `switch_on_temperature`: 有限温度计算开关
- `temper`: 温度 $T$（单位：MeV）
- `entropy(3)`: 熵 $S_n, S_p, S_{\text{tot}}$
- `fn_T(:), fp_T(:)`: 热占据因子 $f_k = \bigl[1 + \exp(E_k / T)\bigr]^{-1}$

###### 8. 错误处理

- 变量位置: `243-244 lines`

- `ierror_flag=0`: 全局错误标志（`0` 表示正常）
- `ierror_info(0:11)`: 错误信息字符串数组

---

##### Namelist 输入接口

以下为程序运行时从输入文件 `hfbtho_NAMELIST.dat` 中读取的全部 Namelist 组，各组变量与上述全局 `_INI` 变量直接映射。

###### `HFBTHO_GENERAL`

- 变量位置: `250-253 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `number_of_shells` | `Integer(ipr)` | $N_{\max}$（谐振子壳层截断）|
| `oscillator_length` | `Real(pr)` | $b_0$（振荡器长度，fm）|
| `basis_deformation` | `Real(pr)` | $q_{20}^{(\text{ini})}$（初始四极形变）|
| `proton_number` | `Integer(ipr)` | $Z$ |
| `neutron_number` | `Integer(ipr)` | $N$ |
| `type_of_calculation` | `Integer(ipr)` | 计算类型（0=HFB, 1=HF, ...）|

###### `HFBTHO_INITIAL`

- 变量位置: `255-256 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `beta2_deformation` | `Real(pr)` | $\beta_2$（初始四极形变）|
| `beta3_deformation` | `Real(pr)` | $\beta_3$（初始八极形变）|
| `beta4_deformation` | `Real(pr)` | $\beta_4$（初始十六极形变）|

###### `HFBTHO_ITERATIONS`

- 变量位置: `258-260 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `number_iterations` | `Integer(ipr)` | 最大自洽迭代次数 |
| `accuracy` | `Real(pr)` | 收敛判据 $\epsilon_{\text{conv}}$ |
| `restart_file` | `Integer(ipr)` | 重启文件编号 |

###### `HFBTHO_FUNCTIONAL`

- 变量位置: `262-265 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `functional` | `Character(30)` | EDF/Skyrme/Gogny 泛函名称 |
| `add_initial_pairing` | `Logical` | 初始配对开关 |
| `type_of_coulomb` | `Integer(ipr)` | 库仑处理类型 |
| `include_3N_force` | `Logical` | 三体力开关 |

###### `HFBTHO_PAIRING`

- 变量位置: `267-269 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `user_pairing` | `Logical` | 用户自定义配对开关 |
| `vpair_n, vpair_p` | `Real(pr)` | 中子/质子配对强度 $V_0^{(n)}, V_0^{(p)}$ |
| `pairing_cutoff` | `Real(pr)` | 单粒子能截断（MeV）|
| `pairing_feature` | `Real(pr)` | 配对特征参数 |

###### `HFBTHO_CONSTRAINTS`

- 变量位置: `271-273 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `lambda_values(1:8)` | `Integer(ipr)` | 约束多极阶数 $\lambda$ |
| `lambda_active(1:8)` | `Integer(ipr)` | 各阶约束激活标志 |
| `expectation_values(1:8)` | `Real(pr)` | 约束期望值 $\langle Q_\lambda \rangle$ |

###### `HFBTHO_BLOCKING`

- 变量位置: `275-276 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `proton_blocking(1:5)` | `Integer(ipr)` | 质子阻塞量子数 $(n_r, n_z, m_l, \pm, k)$ |
| `neutron_blocking(1:5)` | `Integer(ipr)` | 中子阻塞量子数（同上）|

###### `HFBTHO_PROJECTION`

- 变量位置: `278-279 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `switch_to_THO` | `Integer(ipr)` | THO 基开关 |
| `projection_is_on` | `Integer(ipr)` | 投影开关 |
| `gauge_points` | `Integer(ipr)` | 规范角积分点数 $N_\varphi$ |
| `delta_Z, delta_N` | `Integer(ipr)` | 投影粒子数偏差 $\Delta Z, \Delta N$ |

###### `HFBTHO_TEMPERATURE`

- 变量位置: `281-283 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `set_temperature` | `Logical` | 有限温度开关 |
| `temperature` | `Real(pr)` | 温度 $T$（MeV）|

###### `HFBTHO_FEATURES`

- 变量位置: `285-288 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `collective_inertia` | `Logical` | 集体惯性计算开关 |
| `fission_fragments` | `Logical` | 裂变碎片计算开关 |
| `pairing_regularization` | `Logical` | 配对正则化开关 |
| `automatic_basis` | `Logical` | 自动优化基开关 |
| `localization_functions` | `Logical` | 局域化函数开关 |
| `canonical_wavefunctions` | `Logical` | 正则波函数输出开关 |

###### `HFBTHO_FISSION`

- 变量位置: `290-292 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `fragment_properties` | `Logical` | 碎片性质计算开关 |
| `AMP_fragments, PNP_fragments` | `Logical` | 碎片 AMP/PNP 投影开关 |
| `real_Z, real_N` | `Real(pr)` | 碎片真实质子/中子数 |
| `filter` | `Logical` | 碎片滤波开关 |

###### `HFBTHO_EXPERT`

- 变量位置: `294-295 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `density_constraint` | `Logical` | 密度约束开关 |
| `restart_from_canonical` | `Logical` | 从正则基重启开关 |
| `record_HFB_matrix` | `Logical` | 记录 HFB 矩阵开关 |

###### `HFBTHO_LATTICE`

- 变量位置: `297-300 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `spinors_on_grid` | `Logical` | 格点波函数输出开关 |
| `Nx_lattice, Ny_lattice, Nz_lattice` | `Integer(ipr)` | 笛卡尔格点数 |
| `lattice_step` | `Real(pr)` | 格点步长 $dx$（fm）|

###### `HFBTHO_NECK`

- 变量位置: `302-304 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `set_neck_constrain` | `Logical` | 颈部约束开关 |
| `neck_value` | `Real(pr)` | 颈部约束期望值 $Q_{\text{neck}}$ |

###### `HFBTHO_DEBUG`

- 变量位置: `306-310 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `number_Gauss, number_Laguerre, number_Legendre` | `Integer(ipr)` | 高斯/拉盖尔/勒让德积分点数 $N_g, N_l, N_{\text{leg}}$ |
| `compatibility_HFODD` | `Logical` | HFODD 兼容模式 |
| `number_states` | `Integer(ipr)` | 输出状态数 |
| `force_parity` | `Logical` | 强制宇称开关 |
| `write_hel` | `Logical` | 二进制重启文件写入开关 |
| `print_time` | `Integer(ipr)` | 时间打印级别 |

###### `HFBTHO_RESTORATION`

- 变量位置: `312-315 lines`

| 变量 | 类型 | 物理含义 |
|------|------|----------|
| `PNP_is_on` | `Integer(ipr)` | 粒子数投影开关 |
| `number_of_gauge_points` | `Integer(ipr)` | 规范角点数 $N_\varphi$ |
| `delta_neutrons, delta_protons` | `Integer(ipr)` | 中子/质子数偏差 $\Delta N, \Delta Z$ |
| `AMP_is_on` | `Integer(ipr)` | 角动量投影开关 |
| `number_of_rotational_angles` | `Integer(ipr)` | 欧拉角点数 |
| `maximal_angular_momentum` | `Integer(ipr)` | 投影最大角动量 $J_{\max}$ |

---

##### 大规模计算条件编译变量

- 变量位置: `317-373 lines`

根据编译宏 `DO_MASSTABLE=1`、`DO_PES=1`、`DRIP_LINES=1` 定义不同的大规模并行变量：

| 宏定义 | 功能 | 关键变量 |
|--------|------|----------|
| `DO_MASSTABLE` | 质量表计算 | `Z_masstable(:), N_masstable(:), Q20_masstable(:), beta_masstable(:)` |
| `DO_PES` | 势能面计算 | `Z_PES(:), N_PES(:), bet2_PES(:), bet3_PES(:), bet4_PES(:)` |
| `DRIP_LINES` | 滴线计算 | `Z_stable_line(:), N_stable_line(:), direction_sl(:), beta_step` |

---

##### 变量与公式的全局对应关系

| 变量组 | 核心公式 |
|--------|----------|
| 密度场 | $\rho_q(\boldsymbol{r}) = \sum_k V_{qk}^*(\boldsymbol{r}) V_{qk}(\boldsymbol{r})$ |
| 配对场 | $\tilde{\rho}_q(\boldsymbol{r}) = \sum_k U_{qk}^*(\boldsymbol{r}) V_{qk}(\boldsymbol{r})$ |
| 四极矩 | $Q_{20} = \int d^3r \, (2z^2 - r_\perp^2) \, \rho(\boldsymbol{r})$ |
| 形变参数 | $\beta_2 = \dfrac{\sqrt{5\pi}}{3} \dfrac{Q_{20}}{A r_0^2}$ |
| 总能量 | $E = \mathcal{E}[\rho,\tau,\boldsymbol{J},\tilde{\rho}] + E_{\text{Coul}} + E_{\text{pair}} + E_{\text{cm}}$ |
| 热占据 | $f_k = \dfrac{1}{1 + \exp(E_k / T)}$ |

---

##### 模块数据流

- 依赖关系: `HFBTHO_utilities <- HFBTHO`
- 被使用: `hfbthoprog / hfbtho_solver / hfbtho_unedf / hfbtho_gogny / hfbtho_projections / hfbtho_canonical / hfbtho_pnp / hfbtho_collective / hfbtho_fission / hfbtho_io / hfbtho_large_scale / hfbtho_tho / ... <- HFBTHO`
