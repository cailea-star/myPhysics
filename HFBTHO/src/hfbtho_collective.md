# 代码文件 ```hfbtho_collective.f90``` 说明

### `calculate_collective_mass`

##### 函数定义

- 函数位置: `135-276 lines`

```f90
Subroutine calculate_collective_mass()
```

##### 物理公式

- 物理含义: 基于微扰推转近似（perturbative cranking），计算 ATDHFB 与 GCM 两种框架下的**集体惯性质量张量** $M_{ab}$ 及**零点能振动修正** $\epsilon_0$。

核心公式依赖于能量矩 $\mathbf{M}^{(K)}$（由 `energy_moments` 计算）：

**GCM 度规（metric）**

$$
\mathbf{G} = \frac{1}{2} \bigl[\mathbf{M}^{(1)}\bigr]^{-1} \mathbf{M}^{(2)} \bigl[\mathbf{M}^{(1)}\bigr]^{-1}
$$

**ATDHFB 质量张量**

$$
\mathbf{M}_{\mathrm{ATD}} = \hbar^{2} \bigl[\mathbf{M}^{(1)}\bigr]^{-1} \mathbf{M}^{(3)} \bigl[\mathbf{M}^{(1)}\bigr]^{-1}
$$

**ATDHFB 零点能修正**

$$
\epsilon_{\mathrm{ATD}} = -\frac{\hbar^2}{2} \mathrm{Tr}\bigl(\mathbf{B}_{\mathrm{ATD}} \, \mathbf{G}\bigr), \quad \mathbf{B}_{\mathrm{ATD}} = \mathbf{M}_{\mathrm{ATD}}^{-1}
$$

**GCM 惯性张量**

$$
\mathbf{B}_{\mathrm{GCM}} = \frac{1}{4} \mathbf{G}^{-1} \bigl[\mathbf{M}^{(1)}\bigr]^{-1} \mathbf{G}^{-1}
$$

**GCM 质量张量**

$$
\mathbf{M}_{\mathrm{GCM}} = \mathbf{B}_{\mathrm{GCM}}^{-1}
$$

**GCM 零点能修正**

$$
\epsilon_{\mathrm{GCM}} = \frac{\hbar^2}{2} \mathrm{Tr}\bigl(\mathbf{G} \, \mathbf{B}_{\mathrm{GCM}}\bigr)
$$

- 实现步骤:
  1. 调用 `define_collective_space()` 确定集体变量数 `NOMULT` 和张量独立元数 `NOCOMP`
  2. 对每个同位旋 `it`（中子、质子）调用 `energy_moments(it)` 计算 $M^{(1)}, M^{(2)}, M^{(3)}$
  3. 由 `INVMUL` 计算 $[M^{(1)}]^{-1} M^{(K)} [M^{(1)}]^{-1}$ 形式的矩阵
  4. 通过 `INVERT` 求逆得到 $B = M^{-1}$
  5. 通过 `TRACEM` 计算迹，得到零点能
  6. 对中子+质子求和得到总质量张量与总零点能

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 公式对应 |
|--------|------|-----------|----------|
| `SK1(:,it)` | 隐式输出 | 输出 | 一阶能量矩 $M^{(1)}$ |
| `SK2(:,it)` | 隐式输出 | 输出 | 二阶能量矩 $M^{(2)}$ |
| `SK3(:,it)` | 隐式输出 | 输出 | 三阶能量矩 $M^{(3)}$ |
| `ATDMAS(:,it)` | 隐式输出 | 输出 | ATDHFB 质量张量 $M_{\mathrm{ATD}}$ |
| `GCMMAS(:,it)` | 隐式输出 | 输出 | GCM 质量张量 $M_{\mathrm{GCM}}$ |
| `E0_ATD(it)` | 隐式输出 | 输出 | ATDHFB 零点能修正 $\epsilon_{\mathrm{ATD}}$ |
| `E0_GCM(it)` | 隐式输出 | 输出 | GCM 零点能修正 $\epsilon_{\mathrm{GCM}}$ |
| `COVMET(:,:,it)` | 隐式输出 | 输出 | GCM 度规 $\mathbf{G}$ |

##### 函数数据流

- 调用情况: `expect / RESU <- calculate_collective_mass <- define_collective_space / energy_moments / INVMUL / INVERT / TRACEM`

---

### `energy_moments`

##### 函数定义

- 函数位置: `345-496 lines`

```f90
Subroutine energy_moments(it)
```

##### 物理公式

- 物理含义: 计算能量矩 $M_{ab}^{(K)}$（$K=1,2,3$），它是集体惯性张量的核心输入量。矩的定义为两准粒子激发能量的负幂次求和：

$$
M_{ab}^{(K)} = 2 \times \mathfrak{Re} \sum_{\mu\nu}
\frac{\langle \mu\nu | \hat{Q}_{a} | 0\rangle \langle 0 | \hat{Q}_{b} | \mu\nu\rangle}{(E_{\mu}+E_{\nu})^{K}}
$$

其中因子 $2$ 来源于只遍历正单态（positive simplex）的准粒子；$E_\mu$ 为准粒子能量；$\hat{Q}_a$ 为约束多极算符。

- 实现步骤:
  1. 遍历所有集体变量对 $(a, b)$，由 `INDMAS` 映射到多极阶数 $(\lambda, \lambda')$
  2. 对每个 $\Omega$ 块，提取准粒子能量 `Eqp` 和 HFB 系数 `Umatr`, `Vmatr`
  3. 调用 `moments_computeField` 计算多极算符在 HO 基下的矩阵元 `multMatElems`
  4. 调用 `qp_basis_F20` 将算符转换到准粒子基，得到矩阵 $F^{(a)}_{\mu\nu}$ 和 $F^{(b)}_{\mu\nu}$
  5. 对块内所有准粒子对 $(\mu, \nu)$ 求和：
     $\text{aux}(K) \mathrel{+}= F^{(a)}_{\mu\nu} F^{(b)}_{\mu\nu} / (E_\mu + E_\nu)^K$
  6. 乘以因子 $2$ 并存储到 `EMOMET`

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 公式对应 |
|--------|------|-----------|----------|
| `it` | `Integer(ipr), In` | 输入 | 同位旋（1=中子, 2=质子）|
| `Eqp(mu)` | 隐式输入 | 输入 | 准粒子能量 $E_\mu$ |
| `FOPQPA(mu,nu)` | 中间量 | 中间 | $\langle 0 | \beta_\mu \beta_\nu \hat{Q}_a | 0 \rangle$ |
| `FOPQPB(mu,nu)` | 中间量 | 中间 | $\langle 0 | \beta_\mu \beta_\nu \hat{Q}_b | 0 \rangle$ |
| `EMOMET(i_coll,j_coll,K)` | 隐式输出 | 输出 | 能量矩 $M_{ab}^{(K)}$ |

##### 函数数据流

- 调用情况: `calculate_collective_mass <- energy_moments <- qp_basis_F20 / moments_computeField / moments_expectation`

---

### `qp_basis_F20`

##### 函数定义

- 函数位置: `510-544 lines`

```f90
Subroutine qp_basis_F20(nd, n_qp, Umatr, Vmatr, FOPQP, multMatElems)
```

##### 物理公式

- 物理含义: 计算单体态算符 $\hat{F} = \sum_{kl} F_{kl} a_k^\dagger a_l$ 在准粒子真空 $|0\rangle$ 中的双准粒子矩阵元：

$$
\langle 0 | \beta_\mu \beta_\nu \hat{F} | 0 \rangle
= \sum_{kl} F_{kl} \bigl( U_{k\mu}^* U_{l\nu} + V_{k\mu}^* V_{l\nu} \bigr)
$$

在实数 HFB 基下（HFBTHO 使用实系数），简化为：

$$
\mathbf{F}^{\text{qp}} = \mathbf{U}^T \mathbf{F} \, \mathbf{V} + \mathbf{V}^T \mathbf{F} \, \mathbf{U}
$$

- 实现步骤:
  1. 将线性化的 HO 矩阵元 `multMatElems` 恢复为对称矩阵 `dblmul`
  2. 利用 BLAS `dgemm` 计算 $\mathbf{V}^T \mathbf{F} \, \mathbf{U}$，累加到 `FOPQP`
  3. 利用 BLAS `dgemm` 计算 $\mathbf{U}^T \mathbf{F} \, \mathbf{V}$，累加到 `FOPQP`

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 公式对应 |
|--------|------|-----------|----------|
| `nd` | `Integer(ipr), In` | 输入 | HO 基块维度 |
| `n_qp` | `Integer(ipr), In` | 输入 | 准粒子数（$n_{\text{qp}} \le nd$）|
| `Umatr(nd,n_qp)` | `Real(pr), In` | 输入 | HFB 变换矩阵 $U_{k\mu}$ |
| `Vmatr(nd,n_qp)` | `Real(pr), In` | 输入 | HFB 变换矩阵 $V_{k\mu}$ |
| `multMatElems(:)` | `Real(pr), Inout` | 输入 | HO 基下单体态矩阵元 $F_{kl}$（下三角线性存储）|
| `FOPQP(n_qp,n_qp)` | `Real(pr), Inout` | 输出 | 准粒子基下双体态矩阵元 $\langle 0 | \beta_\mu \beta_\nu \hat{F} | 0 \rangle$ |

##### 函数数据流

- 调用情况: `energy_moments <- qp_basis_F20 <- dgemm`

---

### `INVMUL`

##### 函数定义

- 函数位置: `672-684 lines`

```f90
Subroutine INVMUL(A, B, XA, N, NX)
```

##### 物理公式

- 物理含义: 计算特殊矩阵乘积 $\mathbf{A} := \mathbf{A}^{-1} \mathbf{B} \, \mathbf{A}^{-1}$，用于将能量矩组合成惯性张量。

##### 变量对应

- 输入变量1: `A(NX,NX)`, 公式对应矩阵 $\mathbf{A}$
- 输入变量2: `B(NX,NX)`, 公式对应矩阵 $\mathbf{B}$
- 输入/输出变量1: `A(NX,NX)`, 公式对应 $\mathbf{A}^{-1} \mathbf{B} \, \mathbf{A}^{-1}$（覆盖原矩阵）
- 工作数组: `XA(NX,NX)`

##### 函数数据流

- 调用情况: `calculate_collective_mass <- INVMUL <- INVERT`

---

### `INVERT`

##### 函数定义

- 函数位置: `709-741 lines`

```f90
Subroutine INVERT(A, Y, N, NP)
```

##### 物理公式

- 物理含义: 通过 LU 分解（`LUDCMP` + `LUBKSB`）计算方阵的逆矩阵 $\mathbf{A}^{-1}$。算法基于 Crout 分解与列主元选取。

##### 变量对应

- 输入/输出变量1: `A(NP,NP)`, 公式对应矩阵 $\mathbf{A}$（输入）/ $\mathbf{A}^{-1}$（输出）
- 工作数组: `Y(NP,NP)`

##### 函数数据流

- 调用情况: `INVMUL / calculate_collective_mass <- INVERT <- LUDCMP / LUBKSB`

---

### `TRACEM`

##### 函数定义

- 函数位置: `688-702 lines`

```f90
Real(pr) Function TRACEM(A, N, NX)
```

##### 物理公式

- 物理含义: 计算矩阵的迹（trace）：

$$
\mathrm{Tr}(\mathbf{A}) = \sum_{i=1}^{N} A_{ii}
$$

##### 变量对应

- 输入变量1: `A(NX,NX)`, 公式对应矩阵 $\mathbf{A}$
- 输出变量1: `TRACEM`, 公式对应 $\mathrm{Tr}(\mathbf{A})$

##### 函数数据流

- 调用情况: `calculate_collective_mass <- TRACEM`
