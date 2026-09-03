# 代码文件 ```hfbtho_bessel.f90``` 说明

### `besei0`

##### 函数定义

- 函数位置: `121-131 lines`

```f90
Real(pr) Function besei0(x)
```

##### 物理公式

- 物理含义: 计算**指数缩放的修正贝塞尔函数** $e^{-|x|} I_0(x)$，其中 $I_0(x)$ 为第一类零阶修正贝塞尔函数。该缩放形式在 $x \gg 1$ 时可避免数值溢出。

$$
\text{besei0}(x) = e^{-|x|} I_0(x) = e^{-|x|} \sum_{m=0}^{\infty} \frac{1}{(m!)^2} \left(\frac{x}{2}\right)^{2m}
$$

##### 变量对应

- 输入变量1: `x (Real(pr))`, 公式对应 $x$（自变量）
- 输出变量1: `besei0 (Real(pr))`, 公式对应 $e^{-|x|} I_0(x)$

##### 函数数据流

- 调用情况: `coulompj / coulom <- besei0 <- calci0`

---

### `besei1`

##### 函数定义

- 函数位置: `155-165 lines`

```f90
Real(pr) Function besei1(x)
```

##### 物理公式

- 物理含义: 计算**指数缩放的修正贝塞尔函数** $e^{-|x|} I_1(x)$，其中 $I_1(x)$ 为第一类一阶修正贝塞尔函数。

$$
\text{besei1}(x) = e^{-|x|} I_1(x) = e^{-|x|} \sum_{m=0}^{\infty} \frac{1}{m! \, (m+1)!} \left(\frac{x}{2}\right)^{2m+1}
$$

##### 变量对应

- 输入变量1: `x (Real(pr))`, 公式对应 $x$（自变量）
- 输出变量1: `besei1 (Real(pr))`, 公式对应 $e^{-|x|} I_1(x)$

##### 函数数据流

- 调用情况: `上层模块 <- besei1 <- calci1`

---

### `calci0`

##### 函数定义

- 函数位置: `198-302 lines`

```f90
Subroutine calci0(arg, result, jint)
```

##### 物理公式

- 物理含义: `calci0` 为 `besei0` / `besei1` 的底层计算核心，采用分段 minimax 有理逼近与渐近展开计算 $I_0(x)$ 或其指数缩放形式。

$$
I_0(x) = \begin{cases}
1, & |x| < x_{\text{small}} \\[6pt]
\dfrac{\sum_{i=1}^{15} p_i x^{2(i-1)}}{\sum_{j=1}^{5} q_j (x^2 - 225)^{j-1}}, & x_{\text{small}} \le |x| < 15.0 \\[6pt]
\dfrac{e^{|x|}}{\sqrt{|x|}} \left[ \dfrac{\sum_{i=1}^{8} \text{pp}_i \xi^{i-1}}{\sum_{j=1}^{7} \text{qq}_j \xi^{j-1}} \right], & |x| \ge 15.0
\end{cases}
$$

其中 $\xi = |x|^{-1} - 1/15$，系数 `p`、`q`、`pp`、`qq` 为预计算的 minimax 逼近系数。

- 实现步骤:
  1. **极小参数** (`|x| < xsmall`): 直接返回 $I_0 \approx 1$
  2. **中等参数** (`xsmall <= |x| < 15`): 使用 15/5 阶有理多项式逼近计算 $I_0(x)$；若 `jint == 2`，再乘以 $e^{-x}$
  3. **大参数** (`|x| >= 15`): 使用 8/7 阶渐近展开；若 `jint == 2` 输出指数缩放形式 $\text{sump}/(\text{sumq}\sqrt{x})$，避免 $e^x$ 直接溢出
  4. **溢出保护**: 当 `jint == 1` 且 $|x| > x_{\max}$ 时，返回 `xinf`

##### 变量对应

- 输入变量1: `arg (Real(pr))`, 公式对应 $x$（自变量）
- 输入变量2: `jint (Integer(ipr))`, 公式对应功能选择标志（`1` 计算 $I_0(x)$，`2` 计算 $e^{-x} I_0(x)$）
- 输出变量1: `result (Real(pr))`, 公式对应 $I_0(x)$ 或 $e^{-x} I_0(x)$

##### 函数数据流

- 调用情况: `besei0 <- calci0`（无外部直接调用，`calci0` 为 `Private`）

---

### `calci1`

##### 函数定义

- 函数位置: `335-442 lines`

```f90
Subroutine calci1(arg, result, jint)
```

##### 物理公式

- 物理含义: 计算第一类一阶修正贝塞尔函数 $I_1(x)$ 或其指数缩放形式 $e^{-|x|} I_1(x)$，算法结构与 `calci0` 相同，但使用 $I_1$ 专用的 minimax 系数。

$$
I_1(x) = \begin{cases}
\dfrac{x}{2}, & |x| < x_{\text{small}} \\[6pt]
x \cdot \dfrac{\sum_{i=1}^{15} p_i x^{2(i-1)}}{\sum_{j=1}^{5} q_j (x^2 - 225)^{j-1}}, & x_{\text{small}} \le |x| < 15.0 \\[6pt]
\text{sgn}(x) \, \dfrac{e^{|x|}}{\sqrt{|x|}} \left[ \dfrac{\sum_{i=1}^{8} \text{pp}_i \xi^{i-1}}{\sum_{j=1}^{6} \text{qq}_j \xi^{j-1}} + \bar{p} \right], & |x| \ge 15.0
\end{cases}
$$

其中 $\bar{p} = 0.3984375$，$\xi = |x|^{-1} - 1/15$。

- 实现步骤:
  1. **极小参数**: 返回 $I_1 \approx x/2$
  2. **中等参数**: 使用 15/5 阶有理逼近，结果乘以 $x$；若 `jint == 2` 再乘以 $e^{-x}$
  3. **大参数**: 使用 8/6 阶渐近展开；若 `jint == 2` 输出指数缩放形式 $(\text{result} + \bar{p})/\sqrt{x}$
  4. **符号修正**: 若 `arg < 0`，取 `result = -result`（因为 $I_1(-x) = -I_1(x)$）

##### 变量对应

- 输入变量1: `arg (Real(pr))`, 公式对应 $x$（自变量）
- 输入变量2: `jint (Integer(ipr))`, 公式对应功能选择标志（`1` 计算 $I_1(x)$，`2` 计算 $e^{-x} I_1(x)$）
- 输出变量1: `result (Real(pr))`, 公式对应 $I_1(x)$ 或 $e^{-x} I_1(x)$

##### 函数数据流

- 调用情况: `besei1 <- calci1`（无外部直接调用，`calci1` 为 `Private`）
