# 代码文件 ```hfbtho_canonical.f90``` 说明

### `Canonical`

##### 函数定义

- 函数位置: `118-316 lines`

```f90
Subroutine Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibroib)
```

##### 物理公式

- 物理含义: 通过对角化 HFB 密度矩阵 $\rho$，构建**正则基**（canonical basis）。正则基下密度矩阵对角化，本征值给出正则占据数 $v_i^2$，本征矢给出正则单粒子波函数 $|\varphi_i\rangle$。

密度矩阵元由 HFB 系数 $V_{\alpha k}$ 构成：

$$
\rho_{ij} = \sum_{k} V_{ik}^{*} V_{jk}
$$

对角化后：

$$
\sum_j \rho_{ij} \, d_{j\alpha} = \lambda_\alpha \, d_{i\alpha}, \quad v_\alpha^2 = \frac{\lambda_\alpha}{2}
$$

其中 $d_{i\alpha}$ 为正则波函数在 HO 基下的展开系数（存储于 `ddc`），$v_\alpha^2$ 为正则占据数（存储于 `vk`）。

- 实现步骤:
  1. **无配对极限**（`Abs(ept(it)) < 0.0001` 且非有限温度）：直接取 HF 单粒子态，$v^2 = 0$ 或 $1$，无需对角化
  2. **有配对情况**：调用 LAPACK（`DSYEVD` 或 `DSYEVX`）对角化密度矩阵 `hfbcan`，得到本征值 `evvkcan` 和本征矢
  3. **占据数与弥散**：计算 $v_\alpha^2 = \text{evvkcan}(\alpha)/2$，并累加粒子数弥散 $\text{Dispersion} = \sum_\alpha 4 v_\alpha^2 (1-v_\alpha^2)$
  4. **阻塞标记**：若当前块为阻塞块且 $0.49 < v^2 \le 0.51$，记录阻塞态位置 `blocanon`
  5. **正则能量与能隙**（`icanon /= 0`）：通过二次型计算正则单粒子能量 $e_\alpha = \langle\varphi_\alpha|h|\varphi_\alpha\rangle$ 和配对能隙 $\Delta_\alpha = \langle\varphi_\alpha|\Delta|\varphi_\alpha\rangle$
  6. **阻塞候选搜索**（`icanon /= 0`）：在准粒子能量窗口 `pwiblo` 内搜索可阻塞态

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 公式对应 |
|--------|------|-----------|----------|
| `it` | `Integer(ipr), In` | 输入 | 粒子类型（1=中子, 2=质子）|
| `icanon` | `Integer(ipr), In` | 输入 | 详细输出标志（`0`=仅基础信息，`≠0`=计算能量/能隙/搜索阻塞候选）|
| `k2, k1` | `Integer(ipr), In` | 输入 | 准粒子态索引范围 |
| `nd` | `Integer(ipr), In` | 输入 | 当前 $\Omega$ 块的维度 |
| `i0` | `Integer(ipr), In` | 输入 | HO 基索引偏移 |
| `lc` | `Integer(ipr), Inout` | 输入/输出 | 全局正则态计数器 |
| `ib` | `Integer(ipr), In` | 输入 | 当前 $\Omega$ 块编号 |
| `ibiblo` | `Integer(ipr), In` | 输入 | 阻塞块编号 |
| `m` | `Integer(ipr), Inout` | 输入/输出 | LAPACK 返回的本征态数（`DSYEVX` 用）|
| `ibroib` | `Integer(ipr), In` | 输入 | HFB 矩阵的块偏移索引 |
| `ddc(i,lc,it)` | 隐式输出 | 输出 | 正则波函数展开系数 $d_{i\alpha}$ |
| `vk(lc,it)` | 隐式输出 | 输出 | 正则占据数 $v_\alpha^2$ |
| `ek(lc,it)` | 隐式输出 | 输出 | 正则单粒子能量 $e_\alpha$（`icanon≠0`）|
| `dk(lc,it)` | 隐式输出 | 输出 | 正则配对能隙 $\Delta_\alpha$（`icanon≠0`）|

##### 函数数据流

- 调用情况: `hfbdiag <- Canonical <- DSYEVD / DSYEVX`

---

### `canonical_on_grid`

##### 函数定义

- 函数位置: `321-622 lines`

```f90
Subroutine canonical_on_grid(derivatives)
```

##### 物理公式

- 物理含义: 将正则基波函数从谐振子基展开到圆柱坐标 $(r, z)$ 的高斯积分网格上，并输出到文件 `canonwf.dat`。正则波函数的坐标空间表示为：

$$
\varphi_\alpha^{(\sigma)}(r, z) = \sum_{i} d_{i\alpha} \, \phi_{n_z}^{(z)}(z) \, \phi_{n_r}^{m_l}(r) \, \chi_\sigma
$$

其中 $\phi_{n_z}^{(z)}(z)$ 和 $\phi_{n_r}^{m_l}(r)$ 分别为轴向和径向 HO 波函数（代码中对应 `qh` 和 `ql`），$\chi_\sigma$ 为自旋波函数（$\sigma = \pm 1/2$）。

- 实现步骤:
  1. 打开输出文件 `canonwf.dat`
  2. 写入全局配置（粒子数、化学势 `ala`、网格参数）
  3. 写入各正则态的量子数标签（$\tau$, $\Omega$, $k$, $e_{\text{sp}}$, $e_{\text{qp}}$, $\Delta$, $v^2$）
  4. 对每个正则态，在高斯-厄米特/拉盖尔网格上计算波函数值：
     - 自旋向上分量：`xqhla_up`，自旋向下分量：`xqhla_dw`
  5. 若 `derivatives=.True.`，额外计算径向/轴向一阶导数及拉普拉斯算符
  6. 关闭文件

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `derivatives` | `Logical, Optional, In` | 输入 | 是否输出波函数导数 |
| `ddc(i,k0,it)` | 隐式输入 | 输入 | 正则波函数 HO 展开系数 $d_{i\alpha}$ |
| `ek(k0,it)` | 隐式输入 | 输入 | 正则单粒子能量 $e_\alpha$ |
| `dk(k0,it)` | 隐式输入 | 输入 | 正则配对能隙 $\Delta_\alpha$ |
| `vk(k0,it)` | 隐式输入 | 输入 | 正则占据数 $v_\alpha^2$ |
| `ala(it)` | 隐式输入 | 输入 | 化学势 $\lambda$ |
| `bz, bp` | 隐式输入 | 输入 | 谐振子长度参数 $b_z, b_\perp$ |

##### 函数数据流

- 调用情况: `expect / RESU <- canonical_on_grid`

---

### `sample_occupations`

##### 函数定义

- 函数位置: `628-725 lines`

```f90
Subroutine sample_occupations()
```

##### 物理公式

- 物理含义: 将正则基占据数 $v_k^2$ 与准粒子占据数 $u_k^2$（或 $v_k^2$）分别写入文件 `occup_canonical.dat` 和 `occup_qp.dat`，用于后续裂变碎片粒子数统计（参见 Verrière et al. 2019 的 number-of-particles 技术）。

准粒子占据概率与正则占据数的关系：

$$
u_k^2 = \frac{1}{2}\left(1 - \frac{e_k - \lambda}{E_k}\right), \quad u_k^2 = 1 - v_k^2
$$

其中 $E_k = \sqrt{(e_k - \lambda)^2 + \Delta_k^2}$ 为准粒子能量。

- 实现步骤:
  1. 打开 `occup_canonical.dat`，写入头部信息（粒子数、颈部算符期望值）
  2. 遍历所有 $\Omega$ 块和正则态，输出 $v_k^2 > -10^{-4}$ 的态
  3. 打开 `occup_qp.dat`，写入准粒子占据信息
  4. 遍历所有准粒子态，输出 $u_k^2$（代码中变量 `uk`）

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `ek(k0,it)` | 隐式输入 | 输入 | 正则单粒子能量 $e_k$ |
| `vk(k0,it)` | 隐式输入 | 输入 | 正则占据数 $v_k^2$ |
| `uk(k,it)` | 隐式输入 | 输入 | 准粒子占据数 $u_k^2$ |
| `npr(1), npr(2)` | 隐式输入 | 输入 | 中子数 $N$、质子数 $Z$ |

##### 函数数据流

- 调用情况: `expect / RESU <- sample_occupations`
