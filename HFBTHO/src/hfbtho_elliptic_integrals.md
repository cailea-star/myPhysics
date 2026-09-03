# 代码文件 `hfbtho_elliptic_integrals.f90` 说明

本文档对模块 `EllipticIntegral`（460 行）中的全部子程序与函数进行说明。该模块基于 Fukushima (2009) 的快速算法，提供第二类完全椭圆积分 $E(m)$ 的双精度计算，主要用于 Vautherin 直接法求解轴对称库仑势。

---

### `CompleteEllipticFunction_2nd`

##### 函数定义

- 函数位置: `107-125 lines`

```fortran
Real(pr) Function CompleteEllipticFunction_2nd(x)
```

##### 物理公式

- 物理含义: 计算第二类完全椭圆积分 $E(m)$，定义为
$$
E(m)=\int_{0}^{\pi/2}\sqrt{1-m\sin^{2}\theta}\,\mathrm{d}\theta,
\qquad 0\le m\le 1.
$$
当 $m\to 1$ 时，$E(m)\to 1$；当 $m=0$ 时，$E(0)=\pi/2$。

- 实现步骤:
  1. 若 $x<0$ 或 $x>1$，程序终止报错。
  2. 若 $x<0.9$，直接调用 `elliptic_small_m` 用分段泰勒展开计算。
  3. 若 $x\ge 0.9$（接近 1），调用 `auxiliary` 获取互补模数 $m'=1-m$ 对应的辅助量 $E'(m')$、$K'(m')$ 和 nome $q$。
  4. 利用 Legendre 关系将结果转换回 $E(m)$：
$$
E(m)=K(m)+\frac{\pi/2-E'(m')\,K(m)}{K'(m')},
\qquad
K(m)=-\frac{\ln q}{\pi}\,K'(m').
$$

##### 变量对应

- 输入变量1: `x (Real(pr))`，模数 $m$（$0\le m\le 1$）
- 输出变量1: 函数返回值，$E(m)$

##### 函数数据流
- 调用情况: `coulom1 <- CompleteEllipticFunction_2nd <- elliptic_small_m / auxiliary / nome`
- 被 `coulom1`（`hfbtho_solver.f90`）在 Vautherin 库仑势计算中调用

---

### `auxiliary`

##### 函数定义

- 函数位置: `129-187 lines`

```fortran
Subroutine auxiliary(x, Emp, Kmp)
```

##### 物理公式

- 物理含义: 对接近 1 的模数 $m$（即 $m'=1-m$ 很小），利用互补模数 $m'$ 的泰勒展开计算辅助椭圆积分 $E'(m')$ 和 $K'(m')$：
$$
E'(m')=\sum_{i=0}^{16}c_{i}^{(E)}\,(m'-x_{0})^{i},
\qquad
K'(m')=\sum_{i=0}^{16}c_{i}^{(K)}\,(m'-x_{0})^{i},
$$
其中 $x_{0}=0.05$，系数 $c_{i}^{(E)}$、$c_{i}^{(K)}$ 为 Fukushima (2009) 预优化的 minimax 逼近系数。

- 实现步骤:
  1. 计算互补模数 $m'=1-x$。
  2. 以 $x_{0}=0.05$ 为展开中心，分别对 $E'$ 和 $K'$ 做 16 阶泰勒多项式求值。

##### 变量对应

- 输入变量1: `x (Real(pr))`，原模数 $m$
- 输出变量1: `Emp (Real(pr))`，互补模数的第二类椭圆积分 $E'(1-m)$
- 输出变量2: `Kmp (Real(pr))`，互补模数的第一类椭圆积分 $K'(1-m)$

##### 函数数据流
- 调用情况: `CompleteEllipticFunction_2nd <- auxiliary`
- 仅由 `CompleteEllipticFunction_2nd` 在 $m\ge 0.9$ 时调用

---

### `elliptic_small_m`

##### 函数定义

- 函数位置: `191-416 lines`

```fortran
Real(pr) Function elliptic_small_m(x)
```

##### 物理公式

- 物理含义: 对 $0\le m<0.9$ 的模数，采用分段中心泰勒展开计算 $E(m)$：
$$
E(m)=\sum_{i=0}^{16}c_{i}^{(E)}\,(m-x_{0})^{i},
\qquad x_{0}\in\{0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.825,0.875\}.
$$
整个区间 $[0,0.9)$ 被划分为 10 个子区间，每个子区间有独立的优化系数组，保证双精度精度。

- 实现步骤:
  1. 根据 $m$ 所在区间选择展开中心 $x_{0}$ 和对应系数数组。
  2. 执行 16 阶霍纳求值得到 $E(m)$。

##### 变量对应

- 输入变量1: `x (Real(pr))`，模数 $m$（$0\le m<0.9$）
- 输出变量1: 函数返回值，$E(m)$

##### 函数数据流
- 调用情况: `CompleteEllipticFunction_2nd <- elliptic_small_m`
- 仅由 `CompleteEllipticFunction_2nd` 在 $m<0.9$ 时直接调用

---

### `nome`

##### 函数定义

- 函数位置: `420-456 lines`

```fortran
Real(pr) Function nome(x)
```

##### 物理公式

- 物理含义: 计算椭圆积分的 nome $q(m')$，定义为
$$
q(m')=\exp\!\left[-\pi\frac{K'(m')}{K(m')}\right],
\qquad m'=1-m.
$$
对于小 $m'$，$q$ 可用级数展开：
$$
q(m')=\sum_{i=1}^{14}c_{i}^{(q)}\,(m')^{i},
\qquad c_{1}=\frac{1}{16},\;c_{2}=\frac{1}{32},\;c_{3}=\frac{21}{1024},\;\dots
$$
该级数在 $m'\ll 1$ 时快速收敛。若 $m'<10^{-14}$，直接返回 $10^{-14}$ 避免下溢。

##### 变量对应

- 输入变量1: `x (Real(pr))`，互补模数 $m'=1-m$
- 输出变量1: 函数返回值，nome $q(m')$

##### 函数数据流
- 调用情况: `CompleteEllipticFunction_2nd <- nome`
- 仅由 `CompleteEllipticFunction_2nd` 在 $m\ge 0.9$ 时调用，用于通过 $K(m)=-(\ln q/\pi)K'(m')$ 计算第一类完全椭圆积分

---

## 模块数据流总览

```
coulom1 (hfbtho_solver.f90)
  └── CompleteEllipticFunction_2nd(x)
        ├── elliptic_small_m(x)          [x < 0.9]
        └── auxiliary(x, Emp, Kmp)       [x ≥ 0.9]
              └── nome(x_eff)            [x_eff = 1-x]
```

## 参考

- Fukushima, T., "Fast Computation of Complete Elliptic Integrals and Jacobian Elliptic Functions," *Celest. Mech. Dyn. Astron.* **105**, 305–328 (2009b).
