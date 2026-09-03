# 代码文件 ```hfbtho_gogny.f90``` 说明

### `allocate_fr`

##### 函数定义

- 函数位置: `176-189 lines`

```f90
Subroutine allocate_fr()
```

##### 物理含义

分配 Gogny 有限程相互作用及高斯型库仑相互作用所需的全部工作数组，包括直接/交换平均场矩阵 `gamma_g_dir`/`gamma_g_exc`、配对场矩阵 `delta_g_dir`、以及库仑场矩阵 `coulf_g_dir`/`coulf_g_exc`/`coulf_d_dir`。

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `gamma_g_dir(ndx²,2nbx)` | 隐式输出 | 输出 | Gogny 直接平均场 $\Gamma^{(\text{dir})}_{nm}$ |
| `gamma_g_exc(ndx²,2nbx)` | 隐式输出 | 输出 | Gogny 交换平均场 $\Gamma^{(\text{exc})}_{nm}$ |
| `delta_g_dir(ndx²,2nbx)` | 隐式输出 | 输出 | Gogny 配对场 $\Delta_{nm}$ |
| `coulf_g_dir(ndx²,2nbx)` | 隐式输出 | 输出 | 库仑直接平均场 |
| `coulf_g_exc(ndx²,2nbx)` | 隐式输出 | 输出 | 库仑交换平均场 |
| `coulf_d_dir(ndx²,2nbx)` | 隐式输出 | 输出 | 库仑配对场 |

##### 函数数据流

- 调用情况: `gamdel / field <- allocate_fr`

---

### `gogny_matrix_elements`

##### 函数定义

- 函数位置: `895-929 lines`

```f90
Subroutine gogny_matrix_elements()
```

##### 物理公式

- 物理含义: 计算 Gogny 有限程高斯势的两体矩阵元，并将其分离为**轴向分量** $V^z$ 与**径向分量** $V^r$ 的乘积：

$$
\langle n_{z_i} n_{r_i} \Lambda_i | V_i(r) | n_{z_j} n_{r_j} \Lambda_j \rangle
= V^z_{n_{z_i} n_{z_j} n_{z_k} n_{z_l}} \cdot V^r_{n_{r_i} \Lambda_i n_{r_j} \Lambda_j n_{r_k} \Lambda_k n_{r_l} \Lambda_l}
$$

其中高斯势的形式为 $V_i(r) = e^{-(\mathbf{r}_1 - \mathbf{r}_2)^2 / \mu_i^2}$，$\mu_i$ 为高斯展开宽度。

- 实现步骤:
  1. 调用 `calculateTz()` 计算轴向导数相关量
  2. 根据编译宏 `GOGNY_HYPER` 选择计算方法：
     - `GOGNY_HYPER=1`：使用超几何函数 $_2F_1$ 展开（大基下精度更高）
     - 其他：使用直接展开（小基下更快，大基下精度下降）
  3. 调用 `CalculateVzGogny()` 计算并存储所有轴向矩阵元 $V^z$
  4. 调用 `CalculateVrGogny()` 计算并存储所有径向矩阵元 $V^r$

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `VzGogny(:,:)` | 隐式输出 | 输出 | 轴向矩阵元 $V^z$（`GOGNY_SYMMETRIES=1` 时压缩存储）|
| `VrGogny(:,:)` | 隐式输出 | 输出 | 径向矩阵元 $V^r$ |
| `matrix_elements_calculated` | 隐式输出 | 输出 | 计算完成标志（避免重复计算）|

##### 函数数据流

- 调用情况: `field <- gogny_matrix_elements <- calculateTz / calculateTr / calculateCpolar2cartesian / calculateME1D / CalculateVzGogny / CalculateVrGogny`

---

### `gamdel_gogny`

##### 函数定义

- 函数位置: `207-874 lines`

```f90
Subroutine gamdel_gogny()
```

##### 物理公式

- 物理含义: 基于预先计算的 Gogny 两体矩阵元 $V^z$ 和 $V^r$，构建 HFB 平均场 $\Gamma$ 和配对场 $\Delta$。Gogny 相互作用的一般形式为：

$$
V(\mathbf{r}_1, \mathbf{r}_2) = \sum_{i=1}^{N_g} (W_i + B_i \hat{P}_\sigma - H_i \hat{P}_\tau - M_i \hat{P}_\sigma \hat{P}_\tau) e^{-(\mathbf{r}_1 - \mathbf{r}_2)^2 / \mu_i^2}
$$

其中 $W_i, B_i, H_i, M_i$ 为 Gogny 力参数，$\hat{P}_\sigma$ 和 $\hat{P}_\tau$ 分别为自旋和同位旋交换算符。

平均场和配对场的直接/交换项为：

$$
\Gamma^{(\text{dir})}_{ac} = \sum_{bd} \bar{V}_{abcd} \, \rho_{db}, \quad
\Gamma^{(\text{exc})}_{ac} = \sum_{bd} \bar{V}_{abdc} \, \rho_{db}
$$

$$
\Delta^{(\text{dir})}_{ac} = \sum_{bd} \bar{V}_{acbd} \, \kappa_{db}
$$

其中 $\bar{V}$ 包含自旋-同位旋通道的组合系数 $(W \pm B \mp H \mp M)$ 等。

- 实现步骤:
  1. **轴向收缩**（`GOGNY_SYMMETRIES=1`）：先将轴向矩阵元 $V^z$ 与密度/配对矩阵 $\rho, \kappa$ 收缩，得到中间量 $Z_\rho, Z_\kappa$
  2. **径向收缩**：再将径向矩阵元 $V^r$ 与 $Z_\rho, Z_\kappa$ 收缩，得到 $\Gamma^{(\text{dir})}, \Gamma^{(\text{exc})}, \Delta^{(\text{dir})}$
  3. **自旋通道处理**：根据 $s_a + s_c$ 和 $s_b + s_d$ 的取值选择不同的自旋组合（up-up/down-down 或 up-down/down-up）
  4. **库仑场**：若 `icou` 标志激活，同时计算高斯型库仑场的直接/交换项
  5. **DME 近似**：若 `force_is_dme=.True.`，仅保留 Hartree（直接）部分
  6. 将结果累加到全局输出数组 `brout`（用于 Broyden 混合迭代）

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `rk(:,:)` | 隐式输入 | 输入 | 密度矩阵 $\rho_{db}$ |
| `ak(:,:)` | 隐式输入 | 输入 | 配对矩阵 $\kappa_{db}$ |
| `W_g, B_g, H_g, M_g` | 隐式输入 | 输入 | Gogny 力参数 $W_i, B_i, H_i, M_i$ |
| `VzGogny(:,:)` | 隐式输入 | 输入 | 轴向矩阵元 $V^z$ |
| `VrGogny(:,:)` | 隐式输入 | 输入 | 径向矩阵元 $V^r$ |
| `gamma_g_dir(:,:)` | 隐式输出 | 输出 | Gogny 直接平均场 $\Gamma^{(\text{dir})}$ |
| `gamma_g_exc(:,:)` | 隐式输出 | 输出 | Gogny 交换平均场 $\Gamma^{(\text{exc})}$ |
| `delta_g_dir(:,:)` | 隐式输出 | 输出 | Gogny 配对场 $\Delta^{(\text{dir})}$ |
| `coulf_g_dir/exc/d_dir` | 隐式输出 | 输出 | 库仑场（直接/交换/配对）|
| `brout(:)` | 隐式输出 | 输出 | 混合向量（供 Broyden 迭代）|

##### 函数数据流

- 调用情况: `gamdel <- gamdel_gogny <- gogny_matrix_elements / zindex / rindex`

---

### `MatrixElement_z`

##### 函数定义

- 函数位置: `1659-1721 lines`

```f90
Function MatrixElement_z(ni, nj, nk, nl, mu, b) result(Vz)
```

##### 物理公式

- 物理含义: 计算高斯势在轴向上的两体矩阵元：

$$
V^z_{n_i n_j n_k n_l} = \langle n_{z_i} n_{z_j} | e^{-(z_1 - z_2)^2 / \mu^2} | n_{z_k} n_{z_l} \rangle
$$

**标准展开**（`GOGNY_HYPER=0`）：

$$
V^z = \sqrt{\frac{G_z - 1}{G_z + 1}} \sum_{m_z} \sum_{n_z} T_z(n_i, n_k, m_z) \, T_z(n_j, n_l, n_z) \, \bar{I}(m_z, n_z, G_z)
$$

其中 $G_z = 1 + (\mu / b_z)^2$，$T_z$ 为 Talmi-Moshinsky 系数，$\bar{I}$ 为辅助积分。

**超几何展开**（`GOGNY_HYPER=1`）：

$$
V^z = \frac{\mu}{\sqrt{2\pi^3} \, b_z} \sum_{n_z} T_z(n_j, n_l, n_z) \, \bar{F}(n_i, n_k, n_z, z)
$$

其中 $\bar{F}$ 包含 Gamma 函数与 $_2F_1$ 超几何函数。

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `ni, nj, nk, nl` | `Integer(ipr), In` | 输入 | 轴向 HO 量子数 $n_{z_i}, n_{z_j}, n_{z_k}, n_{z_l}$ |
| `mu` | `Real(pr), In` | 输入 | 高斯宽度 $\mu$ |
| `b` | `Real(pr), In` | 输入 | 轴向振荡器长度 $b_z$ |
| `Vz` | `Real(pr)` | 输出 | 轴向矩阵元 $V^z$ |

##### 函数数据流

- 调用情况: `CalculateVzGogny <- MatrixElement_z <- Ibarz / HyperGeom2F1 / T_z`

---

### `trace_product`

##### 函数定义

- 函数位置: `2901-2931 lines`

```f90
Function trace_product(A, B) result(tr)
```

##### 物理公式

- 物理含义: 计算两个块对角对称矩阵乘积的迹，用于求 Gogny 相互作用的能量贡献：

$$
\text{tr} = \sum_{q=n,p} \sum_{\text{block } ib} \sum_{n_1 \ge n_2} f_{n_1 n_2} \, A_{n_1 n_2}^{(ib,q)} \, B_{n_1 n_2}^{(ib,q)}
$$

其中权重因子 $f_{n_1 n_2} = 1$（对角元）或 $2$（非对角元），以利用矩阵对称性。

- 实现步骤:
  1. 遍历所有 $\Omega$ 块 `ib`（中子和质子块通过 `ibx = ib + nbx` 区分）
  2. 对块内下三角元素求和
  3. 非对角元乘以因子 `2`

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `A(:,:)` | `Real(pr), Allocatable, In` | 输入 | 第一个块对角对称矩阵 |
| `B(:,:)` | `Real(pr), Allocatable, In` | 输入 | 第二个块对角对称矩阵 |
| `tr` | `Real(pr)` | 输出 | 乘积的迹 $\text{Tr}(\mathbf{A}\mathbf{B})$ |

##### 函数数据流

- 调用情况: `expect / RESU <- trace_product`
