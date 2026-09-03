# 代码文件 `hfbtho_tho.f90` 说明

本模块实现 **Transformed Harmonic Oscillator (THO)** 基函数包。THO 方法通过局部尺度变换（Local Scale Transformation, LST）将谐振子（HO）基函数映射到具有正确渐近行为的基函数，从而改善对弱束缚核（尤其是晕核）的描述。模块包含 LST 函数 $f(r)$ 及其各阶导数的计算、逆函数构造、雅可比行列式、THO 基函数的导数与拉普拉斯算子，以及球对称密度的 Legendre 展开。

---

### `f01234`

##### 函数定义

- 函数位置: `127-244 lines`

```fortran
Subroutine f01234(lpr)
```

##### 物理公式

- 物理含义: 在 Gauss 积分网格的每个点上计算局部尺度变换函数 $f(r)$ 及其 1–4 阶导数，构造雅可比行列式 $J$ 及其各阶对数导数，并导出波函数的一阶导数和 Laplacian 所需的全部辅助系数。

LST 函数满足 $f: R \mapsto f(R)$，逆关系通过 $q = \sqrt{\eta^2 + \xi^2}$ 与 $f(r) = q$ 确定。Jacobi 行列式为

$$
J = \left( \frac{f^2}{r^2} \frac{df}{dr} \right)^{1/2} = \frac{f\sqrt{f'}}{r}
$$

其对数一阶、二阶导数定义为 $fd1 = J'/J$，$fd2 = J''/J$：

$$
\begin{aligned}
fd1 &= \frac{f'}{f} - \frac{1}{r} + \frac{f''}{2f'} \\
fd2 &= 2\left(\frac{f''}{f} - \frac{fd1}{r}\right) - \frac{1}{4}\left(\frac{f''}{f'}\right)^2 + \frac{f'''}{2f'}
\end{aligned}
$$

辅助变量 $g = f/r$ 及其导数：

$$
\begin{aligned}
g_1 &= -\frac{f - f'r}{r^2} \\
g_2 &= \frac{2(f - f'r) + f''r^2}{r^3} \\
g_3 &= \frac{6(f'r - f) - 3f''r^2 + f'''r^3}{r^4}
\end{aligned}
$$

变换后的坐标与导数：

$$
\begin{aligned}
rr &= \frac{\sqrt{\eta}}{g}, \quad zz = \frac{\xi}{g} \\
z &= b_z \xi \to fh(ihli), \quad \rho = b_\perp \eta \to fl(ihli)
\end{aligned}
$$

一阶导数系数（用于 $z$ 和 $r$ 方向的偏导）：

$$
\begin{aligned}
\partial_z \phi_i &= fp1 \cdot (qh \cdot ql) + fp2 \cdot (qh_1 \cdot ql) + fp3 \cdot (qh \cdot ql_1)/(2v) \\
\partial_r \phi_i &= fp4 \cdot (qh \cdot ql) + fp5 \cdot (qh_1 \cdot ql) + fp6 \cdot (qh \cdot ql_1)/(2v)
\end{aligned}
$$

Laplacian 系数 $fs1$–$fs6$ 满足：

$$
\Delta \rho / j^2 = (c_z fs1 + c_r fs2 + fs3) \cdot qh_{ab} \cdot ql_{ab} + \cdots
$$

其中 $c_z = u^2 - (n_z^a + n_z^b + 1)$，$c_r = \frac{1}{4} - \frac{n_r^a + n_r^b + m + 1}{2v} + \left(\frac{m}{2v}\right)^2$。

权重因子：

$$
wdcor(i) = \frac{\pi \, wh(ih) \, wl(il) \, b_z \, b_\perp^2}{fd(i)}
$$

- 实现步骤:
  1. 遍历 $(ih, il)$ Gauss 网格，计算 $q = \sqrt{\xi^2 + \eta^2}$。
  2. **HO 情形** (`iLST1 == 0`): 直接取 $r = q$，$f = r$，$f' = 1$，$f'' = f''' = 0$。
  3. **THO 初始运行** (`iLST < 0`, `ih*il == 1`): 调用 `thofun(key=0, ...)` 进行全局初始化。
  4. **THO 一般运行**: 调用 `thofun(key=1, qq, f, f1, f2, f3, r, ...)` 求解逆函数 $r = f^{-1}(q)$ 及 $f$ 的各阶导数。
  5. 若渐近行为异常 (`iasswrong(3) /= 0`)，强制回退到 HO 结果。
  6. 按上述公式计算 $J$、$g$、坐标变换 $(rr, zz)$、$(r,z)$ 导数、$(u,v)$ 导数。
  7. 存储 $fp1$–$fp6$（一阶导数系数）和 $fs1$–$fs6$（Laplacian 系数）到模块数组。
  8. 计算并存储积分权重 $wdcor$ 和 $wdcori = 1/wdcor$。

##### 变量对应

- 输入变量: `lpr (Logical)`，控制是否打印 HO/THO 模式信息到输出文件
- 输出（模块公共数组）:
  - `fh(ihli) (Real(pr))`: 公式对应 $z = b_z \xi = b_z \cdot zz$，纵向坐标
  - `fl(ihli) (Real(pr))`: 公式对应 $\rho = b_\perp \eta = b_\perp \cdot rr$，径向坐标
  - `fli(ihli) (Real(pr))`: $1/fl$，径向坐标倒数
  - `fd(ihli) (Real(pr))`: 公式对应 $J^2 = f^2 f'/r^2$，Jacobi 行列式平方
  - `fp1`–`fp6(ihli) (Real(pr))`: 一阶导数展开系数
  - `fs1`–`fs6(ihli) (Real(pr))`: Laplacian 展开系数
  - `wdcor(i) (Real(pr))`: 坐标空间积分权重
  - `wdcori(i) (Real(pr))`: $1/wdcor$，权重倒数

##### 函数数据流

- 调用情况: `coordinateLST_THO <- f01234 <- thofun(key=0/1) <- densitr <- gaupolr`

---

### `thofun`

##### 函数定义

- 函数位置: `249-717 lines`

```fortran
Subroutine thofun(key,r,f,f1,f2,f3,fj,lpr,units)
```

##### 物理公式

- 物理含义: 计算 LST 函数 $f(r)$ 及其 1–3 阶导数。`key=0` 时进行全局初始化：从 HFB+HO 的 $L=0$ 球对称密度出发，构造修正密度使其具有正确的指数渐近行为，然后通过逆映射得到 $f_j(q) = f^{-1}(q)$ 并建立三次样条插值表。`key=1` 时利用插值表快速返回给定 $q$ 处的函数值。

**KEY=0 — 初始化流程**

1. **中子/质子渐近行为测试**: 对每个同位旋 `it=1,2`，计算 HFB+HO $L=0$ 密度 $\rho_{HO}^{(L=0)}(r)$ 的对数导数 $\partial_r \ln\rho$，找到第一个极小值点 $r_{mm1}$、第一个极大值点 $r_{mmmax}$、最后一个可接受点 $r_{mm2}$ 及中点 $r_{mmm}$。

2. **选择主导同位旋**: 比较两个同位旋在中子/质子密度极小值点处的密度大小，选择衰减常数 $\alpha$ 更小的（密度更大的）同位旋作为主导渐近行为。

3. **渐近密度拟合**: 在 $r > r_{mm1}$ 区域，密度对数满足

$$
\ln\rho(r) = \alpha_{ex} + z_1^{mm1} r + aaa \cdot \frac{r^{1-assm}(assm \cdot assm2 \cdot r_{mm1}^2 - 2 assm \cdot assm1 \cdot r_{mm1} r + assm2 \cdot assm1 \cdot r^2)}{assm3 \cdot assm2 \cdot assm1}
$$

其中 $assm$ 由 $r_{mmx}=r_{mmmax}$ 处的对数导数匹配条件确定。在 $r > r_{mmmax}$ 区域使用解析渐近形式：

$$
\ln\rho = bex + yyy(r), \quad yyy = -r\sqrt{\alpha^2 + \frac{4(aa + bb/r)}{r^2}} - \ln(r^2 \sqrt{\cdots}) - \frac{2bb}{\alpha}\ln\left(\frac{2bb + \alpha^2 r + \alpha r \sqrt{\cdots}}{\alpha}\right) + \cdots
$$

4. **归一化与样条**: 对修正密度 $dsy1$ 归一化至总核子数 $N_{it}$，然后对 $dsy1$ 及其积分 $dsy1i$ 分别建立三次样条 (`csplin`)。

5. **逆函数迭代求解**: 对每个 $q$ 网格点，通过割线法/牛顿法迭代求解

$$
f_j = f^{-1}(q) \quad \text{满足} \quad f(f_j) = q
$$

迭代收敛判据基于 $\ln(z_{fj1i}/z_{qqi}) \cdot z_{fj1i} / (z_{fj1} \cdot b_0 f_j^2)$ 的比值，其中 $z_{fj1i}$ 是修正密度积分在 $f_j$ 处的样条插值值，$z_{qqi}$ 是 HO 密度积分。

6. **导数关系**: 由 $f = q$，$f_j = f^{-1}(q)$，有

$$
\begin{aligned}
f' &= \frac{1}{f_j'} \\
f'' &= -\frac{f_j''}{(f_j')^3} \\
f''' &= \frac{3(f_j'')^2}{(f_j')^5} - \frac{f_j'''}{(f_j')^4}
\end{aligned}
$$

在 $q > r_{mm3}$ 的远渐近区使用二次外推：$f_j = amm3 + bmm3(q-r_{mm3}) + cmm3(q-r_{mm3})^2/2$。

**KEY=1 — 快速查表**

直接调用 `cseval` 对预计算的三次样条进行插值，得到 $f_j, f_j', f_j'', f_j'''$，再按上述导数关系反解 $f, f', f'', f'''$。

- 实现步骤:
  1. `key == 0`:
     - 分配工作数组 `dsx`, `dsy`, `dsyi`, `dsy1`, `spb0`–`spd0` 等。
     - 循环 `it=1,2` 测试渐近行为，调用 `densitr` 获取 HFB+HO $L=0$ 密度。
     - 找到对数密度的特征点（极小、极大、中点、末点）。
     - 解三次方程确定参数 $aa$，匹配中点条件。
     - 构造分段修正密度（内区=原密度，中区=抛物拟合，外区=解析渐近）。
     - 归一化后建立样条表，迭代求逆函数 $f_j(q)$。
     - 对逆函数做数值微分 (`deri`)，建立 $f_j, f_j', f_j'', f_j'''$ 的四组样条系数。
     - 若 `lpr=.True.` 输出密度到文件 `density.dat`。
  2. `key == 1`:
     - 通过 `cseval` 查表获取 $f_j$ 及其 1–3 阶导数。
     - 按反函数微分公式计算 $f, f', f'', f'''$。

##### 变量对应

- 输入变量:
  - `key (Integer(ipr))`: `0`=初始化，`1`=查表计算
  - `r (Real(pr))`: 输入/输出，初始化时输入为任意值（内部覆盖为 $qq$），查表时输入为 $q$
  - `lpr (Logical)`: 是否打印诊断信息并输出密度文件
  - `units (Logical)`: 是否使用有量纲坐标（`.True.` 时 $hhb=h$ 而非 $h \cdot b_0$）
- 输出变量:
  - `f (Real(pr))`: 公式对应 $f(r) = q$，LST 函数值
  - `f1 (Real(pr))`: 公式对应 $f'(r)$
  - `f2 (Real(pr))`: 公式对应 $f''(r)$
  - `f3 (Real(pr))`: 公式对应 $f'''(r)$
  - `fj (Real(pr))`: 公式对应 $f_j = f^{-1}(q)$，逆函数值
- 模块级输出（`key=0` 时填充）:
  - `decay (Real(pr))`: 选定的密度衰减常数 $\alpha$
  - `fdsx(iqqmax) (Real(pr))`: 逆函数样条节点 $q_i$
  - `fdsy`, `fdsy1`, `fdsy2`, `fdsy3`: 逆函数及其 1–3 阶导数的样条节点值
  - `fspb0`–`fspd3`: 四组三次样条系数数组
  - `rmm3`, `amm3`, `bmm3`, `cmm3`: 远渐近区二次外推参数

##### 函数数据流

- 调用情况: `f01234 <- thofun(key=1) <- densitr <- gaupolr`
- `thofun(key=0)` 内部调用: `densitr(itass, rmmho, sN, sP, msw)`, `deri(hhb, immho, dsy, spb0)`, `csplin(irmax, dsx, dsy1, spb0, spc0, spd0)`, `cseval(iqqmax, qq, fdsx, fdsy, fspb0, fspc0, fspd0, fj, iqqs)`

---

### `densitr`

##### 函数定义

- 函数位置: `725-820 lines`

```fortran
Subroutine densitr(it,xr,yr,yrP,msw)
```

##### 物理公式

- 物理含义: 通过 Legendre 展开计算中子或质子的球对称密度 $\rho(r)$ 和配对密度 $\tilde{\rho}(r)$ 在给定半径 $r$ 处的值。密度矩阵元来自 HFB 波函数的径向积分：

$$
\rho(r) = \int_0^{\pi/2} \rho(r, \theta) \, P_0(\cos\theta) \, \sin\theta \, d\theta
$$

其中 $\rho(r, \theta)$ 通过 HO/THO 基函数在柱坐标 $(z, \rho)$ 上的乘积展开得到。基函数在角度 $\theta$ 下的值为：

$$
\phi_{iw,ik} = qh(n_z, n_{gh}+1) \cdot ql(n_r, m, n_{gl}+1) \cdot wdcorin
$$

密度通过对密度矩阵 $R$（HFB 正则密度）和 $R^P$（配对密度）的双重求和得到：

$$
\begin{aligned}
yr &= \sum_{ib} \sum_{n_1 \geq n_2}^{nd} frit \cdot \sum_{iw=1}^{msw} yi(1,iw) \cdot \phi_{iw,ik} \cdot \phi_{iw,jk} \\
yrP &= \sum_{ib} \sum_{n_1 \geq n_2}^{nd} fritP \cdot \sum_{iw=1}^{msw} yi(1,iw) \cdot \phi_{iw,ik} \cdot \phi_{iw,jk}
\end{aligned}
$$

其中 $frit = R_{n_1n_2}^{(ib)}$，当 $n_1 \neq n_2$ 时乘以 2（对称性）。$yi(1,iw)$ 是 Legendre 多项式 $P_0(\cos\theta_{iw})$ 的离散展开系数，通过 `lingd` 求解线性方程组得到。

- 实现步骤:
  1. **初始化** (`it < 0`): 将 `it` 取反，分配角度网格 `xmw(msw)` 和 Legendre 系数矩阵 `yi(msw,msw)`。
  2. 在 $[0, \pi/2]$ 区间均匀采样 $msw$ 个角度点，计算各阶 Legendre 多项式 $yl(iw, il)$ 的值。
  3. 通过 `lingd` 求解 $yl \cdot yi = I$（单位矩阵），得到离散逆变换系数 $yi$。
  4. 对每个角度 $w = \theta_{iw}$，计算 $z = \cos(w) \cdot xr / b_z$，$t = (\sin(w) \cdot xr / b_\perp)^2$。
  5. 调用 `gaupolr(z, t)` 计算 HO 基函数值 `qh` 和 `ql`。
  6. 按块结构遍历基态，计算每个基函数的角向投影 `phy(iw, ik)`。
  7. 对密度矩阵元求和，分别累积正则密度 `yr` 和配对密度 `yrP`。

##### 变量对应

- 输入变量:
  - `it (Integer(ipr))`: 同位旋，`1`=中子，`2`=质子；首次调用时传入负值触发初始化
  - `xr (Real(pr))`: 公式对应 $r$，计算密度的空间点（单位：fm）
  - `msw (Integer(ipr))`: 角度网格点数，控制 Legendre 展开的精度（测试中 `msw=24` 已收敛）
- 输出变量:
  - `yr (Real(pr))`: 公式对应 $\rho(r)$，球对称正则密度
  - `yrP (Real(pr))`: 公式对应 $\tilde{\rho}(r)$，球对称配对密度

##### 函数数据流

- 调用情况: `thofun(key=0) <- densitr <- gaupolr`
- `densitr` 内部调用: `gaupolr(z, t)`, `lingd(msw, msw, msw, msw, yl, yi, s, il)`

---

### `gaupolr`

##### 函数定义

- 函数位置: `825-857 lines`

```fortran
Subroutine gaupolr(z,x)
```

##### 物理公式

- 物理含义: 在单点 $(z, x)$ 处计算归一化的 Gauss-Hermite 多项式 $qh(n, n_{gh}+1)$ 和 Gauss-Laguerre 多项式 $ql(n, l, n_{gl}+1)$。这些多项式是 HO 基函数的径向部分。

Gauss-Hermite 递推（归一化形式）:

$$
\begin{aligned}
q_0 &= \pi^{-1/4} e^{-z^2/2} \\
q_1 &= \sqrt{2} \, z \, q_0 \\
q_n &= \frac{1}{\sqrt{n}} \left( \sqrt{2} \, z \, q_{n-1} - \sqrt{n-1} \, q_{n-2} \right)
\end{aligned}
$$

Gauss-Laguerre 递推（对每个角动量 $l$）:

$$
\begin{aligned}
W_{00} &= \sqrt{2} \, e^{-x/2} \\
W_0^{(l)} &= W_{00} \sqrt{\frac{x^l}{2}} \quad (l > 0) \\
ql(0, l) &= wfi(l) \cdot W_0^{(l)} \\
ql(1, l) &= (l+1 - x) \cdot wfi(l+1) \cdot W_0^{(l)} \\
ql(n, l) &= \frac{(2n+l-1-x) \cdot ql(n-1, l) - \sqrt{(n-1)(n-1+l)} \cdot ql(n-2, l)}{\sqrt{n(n+l)}}
\end{aligned}
$$

其中 $wfi(l) = 1/\sqrt{l!}$ 为归一化因子。

- 实现步骤:
  1. 计算 Gauss-Hermite 基值 $q_0$ 和 $q_1$，通过递推公式计算至最大量子数 $n_{zm}$。
  2. 对每个角动量 $l = 0, \dots, n_{lm}$，计算 Laguerre 基值 $W_0^{(l)}$。
  3. 通过递推公式计算 $ql(n, l)$ 至最大径向量子数 $n_{rm}$。
  4. 结果存储到模块数组 `qh(n, ngh1)` 和 `ql(n, l, ngl1)` 的额外索引位置（`ngh1 = ngh+1`, `ngl1 = ngl+1`）。

##### 变量对应

- 输入变量:
  - `z (Real(pr))`: 公式对应 $z = \xi$（无量纲纵向坐标），用于 Hermite 多项式
  - `x (Real(pr))`: 公式对应 $x = \eta^2$（无量纲径向坐标平方），用于 Laguerre 多项式
- 输出（模块数组）:
  - `qh(n, ngh1) (Real(pr))`: 第 $n$ 阶归一化 Hermite 多项式在 $z$ 处的值
  - `ql(n, l, ngl1) (Real(pr))`: 第 $(n, l)$ 阶归一化 Laguerre 多项式在 $x$ 处的值

##### 函数数据流

- 调用情况: `densitr <- gaupolr`

---

### `coordinateLST_THO`

##### 函数定义

- 函数位置: `862-871 lines`

```fortran
Subroutine coordinateLST_THO()
```

##### 物理公式

- 物理含义: THO 坐标变换的顶层驱动子程序。依次调用 `f01234` 计算 LST 变换后的网格坐标、Jacobi 权重和辅助系数，再调用 `optHFBTHO_THO` 预计算 THO 基函数在变换后网格上的导数和 Laplacian。本程序是 HO 坐标初始化 `coordinateLST_HO` 的 THO 对应版本。

- 实现步骤:
  1. 调用 `f01234(.False.)` 完成 LST 变换和网格初始化。
  2. 检查全局错误标志 `ierror_flag`，非零则提前返回。
  3. 调用 `optHFBTHO_THO()` 构建 THO 基函数的优化数组。
  4. 再次检查错误标志。

##### 变量对应

- 输入: 无显式参数，通过模块变量 `ierror_flag` 传递状态
- 输出: 调用 `f01234` 和 `optHFBTHO_THO` 的副作用，填充模块公共数组

##### 函数数据流

- 调用情况: `coordinateLST (hfbtho_solver.f90) <- coordinateLST_THO <- f01234 / optHFBTHO_THO`

---

### `optHFBTHO_THO`

##### 函数定义

- 函数位置: `876-918 lines`

```fortran
Subroutine optHFBTHO_THO()
```

##### 物理公式

- 物理含义: 预计算 THO 基函数在变换后网格上的值及其一阶导数、Laplacian，存储到优化数组中供后续 HFB 迭代使用。这是 `optHFBTHO_HO` 的 THO 对应版本。

对于每个基态 $|n_z, n_r, m, s\rangle$，在网格点 $(ih, il)$ 上的基函数值为：

$$
qhla = qh(n_z, ih) \cdot ql(n_r, m, il)
$$

其一阶径向和纵向导数为：

$$
\begin{aligned}
fi1r &= fp4(ihli) \cdot qhla + fp5(ihli) \cdot qh1la + fp6(ihli) \cdot qhl1a \\
fi1z &= fp1(ihli) \cdot qhla + fp2(ihli) \cdot qh1la + fp3(ihli) \cdot qhl1a
\end{aligned}
$$

其中 $qh1la = qh(n_z, ih) \cdot \partial_\rho ql(n_r, m, il) / (2\eta)$，$qhl1a = \partial_z qh(n_z, ih) \cdot ql(n_r, m, il)$。

Laplacian 的离散表达式为：

$$
\begin{aligned}
fi2d &= \frac{1}{2 \cdot qhla} \Big[ fs1 \cdot (qh1la)^2 + fs2 \cdot (qhl1a)^2 + 4 fs4 \cdot qh1la \cdot qhl1a \\
&\quad + 2(fs5 \cdot qh1la + fs6 \cdot qhl1a) \cdot qhla \\
&\quad + ((u^2 - c_{nz}) \cdot fs1 + (1/4 - c_{nr} v^2 + m^2 v^4) \cdot fs2 + fs3) \cdot (qhla)^2 \\
&\quad - 2(fi1r^2 + fi1z^2) \Big]
\end{aligned}
$$

其中 $c_{nz} = 2n_z + 1$，$c_{nr} = 2n_r + m + 1$，$v^2 = 1/(2\eta)$。

- 实现步骤:
  1. 分配优化数组 `QHLA_opt(ntx, nghl)`, `FI1R_opt`, `FI1Z_opt`, `FI2D_opt(ntx, nghl)`, `y_opt(nghl)`。
  2. 遍历块结构 `ib=1, NB`，对每个块内的基态 $n_1=1, \dots, nd$。
  3. 遍历 Gauss 网格 $(ih, il)$，读取预计算的 LST 系数 `fli(ihli)` 和 `fp1`–`fp6`、`fs1`–`fs6`。
  4. 计算基函数值 `qhla`、纵向导数 `qh1la`、径向导数 `qhl1a`。
  5. 按上述公式计算 `fi1r`、`fi1z`、`fi2d`。
  6. 将结果存入优化数组对应位置。

##### 变量对应

- 输入: 无显式参数，依赖模块数组 `fli`, `fp1`–`fp6`, `fs1`–`fs6`, `qh`, `ql`, `qh1`, `ql1`, `xh`, `xl`
- 输出（模块数组）:
  - `QHLA_opt(ja, ihil) (Real(pr))`: 公式对应 $\phi_{ja}(r_{ihli})$，THO 基函数值
  - `FI1R_opt(ja, ihil) (Real(pr))`: 公式对应 $\partial_r \phi_{ja}$，径向一阶导数
  - `FI1Z_opt(ja, ihil) (Real(pr))`: 公式对应 $\partial_z \phi_{ja}$，纵向一阶导数
  - `FI2D_opt(ja, ihil) (Real(pr))`: 公式对应 $\Delta \phi_{ja} / (2\phi_{ja})$，Laplacian 相关量
  - `y_opt(ihil) (Real(pr))`: 公式对应 $1/fl(ihli)$，径向坐标倒数缓存

##### 函数数据流

- 调用情况: `coordinateLST_THO <- optHFBTHO_THO`
