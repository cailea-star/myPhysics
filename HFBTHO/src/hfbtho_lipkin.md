# 代码文件 ```hfbtho_lipkin.f90``` 说明

### `tracesln`

##### 函数定义

- 函数位置: `116-186 lines`

```f90
Subroutine tracesln()
```

##### 物理公式

- 物理含义: 在**正则基 (canonical basis)** 下计算 Lipkin-Nogami 参数 $\lambda_2$（代码中记为 `ala2`），并给出相应的能量修正 `etr`。该方法假设单极对力形式为

$$
\bar{v}_{ijkl} = -G \, \delta_{j\bar{i}}\delta_{l\bar{k}} \, \text{sign}(j)\text{sign}(l)
$$

正则基下的 LN 参数公式为

$$
\lambda_2 = \frac{G_{\text{eff}}}{4} \frac{\displaystyle \sum_{i>0} u_i^3 v_i \sum_{i>0} u_i v_i^3 - \sum_{i>0} (u_i v_i)^4}{\displaystyle \left( \sum_{i>0} u_i^2 v_i^2 \right)^2 - \sum_{i>0} (u_i v_i)^4}
$$

其中有效对强度 $G_{\text{eff}} = \Delta^2 / E_{\text{pair}}$（代码中 $\Delta$ 可为正，而配对能 $E_{\text{pair}} < 0$，故 `Geff` 为负值；代码通过额外负号保证 `ala2` 与 $\lambda_2$ 同号）。

Lipkin-Nogami 能量修正为

$$
E_{\text{LN}}^{(it)} = -4 \lambda_2 \sum_{i>0} u_i^2 v_i^2, \qquad
E_{\text{LN}}^{(\text{tot})} = E_{\text{LN}}^{(n)} + E_{\text{LN}}^{(p)}
$$

- 实现步骤:
  1. 遍历各块与正则态，累加五个关键组合量 $S_{u^3v}$、$S_{uv^3}$、$S_{u^2v^2}$、$S_{u^4v^4}$
  2. 计算分子 `SNtor` 与分母 `SDtor`
  3. 若 `SDtor` 非零，用 $G_{\text{eff}} = \Delta^2/E_{\text{pair}}$ 直接求 `ala2`；若 `SDtor` 趋于零（弱配对极限），改用对数/指数形式避免相消
  4. 对 `ala2` 做上限截断（`>= 10` 时置为 `4.0`），并计算能量修正 `etr`

##### 变量对应

- 输入变量1: `vk (Real(pr), from HFBTHO_canonical)`, 公式对应 $v_i^2$（正则基占据概率）
- 输入变量2: `ek (Real(pr), from HFBTHO_canonical)`, 公式对应 $e_i$（正则基单粒子能量，用于弱配对极限的 $G$ 估算）
- 输入变量3: `ala (Real(pr), from HFBTHO)`, 公式对应 $\lambda$（化学势，弱配对极限备用）
- 输入变量4: `del (Real(pr), from HFBTHO)`, 公式对应 $\Delta$（对能隙）
- 输入变量5: `ept (Real(pr), from HFBTHO)`, 公式对应 $E_{\text{pair}}$（配对能）
- 输出变量1: `ala2(2) (Real(pr), module variable)`, 公式对应 $\lambda_2$（中子/质子 LN 参数）
- 输出变量2: `Geff(2) (Real(pr), module variable)`, 公式对应 $G_{\text{eff}} = \Delta^2/E_{\text{pair}}$
- 输出变量3: `etr(3) (Real(pr), module variable)`, 公式对应 $E_{\text{LN}}$（中子、质子及总 LN 能量修正）

##### 函数数据流

- 调用情况: `hfbtho_solver.f90 / expect <- tracesln`

---

### `tracesln_qp`

##### 函数定义

- 函数位置: `204-268 lines`

```f90
Subroutine tracesln_qp()
```

##### 物理公式

- 物理含义: 在**准粒子基 (qp basis)** 下计算 Lipkin-Nogami 参数 $\lambda_2$（`ala2`）及能量修正 `etr`。与 `tracesln` 的物理结果等价，但直接用密度矩阵 $\rho$ 与配对矩阵 $\kappa$ 表示：

$$
\lambda_2 = \frac{G_{\text{eff}}}{4}
\frac{\text{Tr}^{<}\!\big[\kappa^{*}\rho\big] \; \text{Tr}^{<}\!\big[(1-\rho)\kappa\big]
- \sum_{ij}\big[\rho(1-\rho)\big]_{\bar{i}\bar{j}}\big[\rho(1-\rho)\big]_{ij}
}{\displaystyle \Big( \text{Tr}^{>}\!\big[\rho(1-\rho)\big] \Big)^2
- \text{Tr}^{>}\!\big[\rho^2(1-\rho)^2\big]}
$$

能量修正公式与 `tracesln` 相同：

$$
E_{\text{LN}}^{(it)} = -4 \lambda_2 \, \text{Tr}^{>}\!\big[\rho(1-\rho)\big], \qquad
E_{\text{LN}}^{(\text{tot})} = E_{\text{LN}}^{(n)} + E_{\text{LN}}^{(p)}
$$

- 实现步骤:
  1. 按块遍历 $\rho$（`rk`）与 $\kappa$（`ak`）矩阵元
  2. 利用对角/非对角元的对称性分别累加 $\text{Tr}[\rho(1-\rho)]$、$\text{Tr}[\rho^2(1-\rho)^2]$、$\text{Tr}[\rho\kappa]$、$\text{Tr}[(1-\rho)\kappa]$
  3. 组合得到 `SNtor`、`SDtor`，进而求得 `ala2` 与 `etr`

##### 变量对应

- 输入变量1: `rk (Real(pr), from HFBTHO)`, 公式对应 $\rho_{ij}$（密度矩阵）
- 输入变量2: `ak (Real(pr), from HFBTHO)`, 公式对应 $\kappa_{ij}$（配对矩阵）
- 输入变量3: `del (Real(pr), from HFBTHO)`, 公式对应 $\Delta$（对能隙）
- 输入变量4: `ept (Real(pr), from HFBTHO)`, 公式对应 $E_{\text{pair}}$（配对能）
- 输出变量1: `ala2(2) (Real(pr), module variable)`, 公式对应 $\lambda_2$
- 输出变量2: `Geff(2) (Real(pr), module variable)`, 公式对应 $G_{\text{eff}} = \Delta^2/E_{\text{pair}}$
- 输出变量3: `etr(3) (Real(pr), module variable)`, 公式对应 $E_{\text{LN}}$（中子、质子及总 LN 能量修正）

##### 函数数据流

- 调用情况: `hfbtho_solver.f90 / expect <- tracesln_qp`（公共接口，可供 qp 基分析调用）

---

### `densitln`

##### 函数定义

- 函数位置: `273-341 lines`

```f90
Subroutine densitln()
```

##### 物理公式

- 物理含义: 计算经 **Lipkin-Nogami 修正后的坐标空间密度**。在正则基中，原始占据数为 $v^2$，LN 修正后的有效占据数为

$$
v_{\text{LN}}^2 = v^2 + (uv)^2 \, \frac{\sigma_1}{\sigma_3}
\Big[ (2v^2 - 1)\sigma_1 - \sigma_2 \Big]
$$

其中 $u = \sqrt{1-v^2}$，$\sigma_1, \sigma_2, \sigma_3$ 为预计算的 LN 参数（存储在模块变量 `ssln(3,2)` 中）。

最终密度在 Gauss 网格点 $(i_h, i_l)$ 上累加为

$$
\rho_{it}(r_h, z_l) = \sum_{\alpha\beta} \phi_{\beta}(r_h, z_l) \phi_{\alpha}(r_h, z_l) \; v_{\text{LN},\alpha}^2 \; C_{\beta}^{(k)} C_{\alpha}^{(k)}
$$

（求和包含自旋简并因子与块结构，输出至 `ro` 后再乘以积分权重 `wdcori`。）

- 实现步骤:
  1. 遍历各块及基矢对 $(n_1, n_2)$
  2. 对每个正则态 $k$ 读取 `vk` 与 LN 参数 `ssln1, ssln2, ssln3`，按上式计算修正占据 `vvc`
  3. 用正则波函数 `ddc` 构造径向-角度乘积因子 `fr`
  4. 在 Gauss-Hermite / Gauss-Laguerre 网格上组装密度 `ro`
  5. 最后乘以体积元权重 `wdcori`

##### 变量对应

- 输入变量1: `vk (Real(pr), from HFBTHO_canonical)`, 公式对应 $v_k^2$（正则基占据数）
- 输入变量2: `ddc (Real(pr), from HFBTHO_canonical)`, 公式对应 $C_{\alpha}^{(k)}$（正则基波函数展开系数）
- 输入变量3: `ssln(3,2) (Real(pr), module variable)`, 公式对应 $\sigma_1, \sigma_2, \sigma_3$（LN 修正参数，由 `resu` 预计算）
- 输出变量1: `ro (Real(pr), from HFBTHO)`, 公式对应 $\rho_{it}(r_h, z_l)$（修正后的坐标空间密度数组）

##### 函数数据流

- 调用情况: `hfbtho_solver.f90 / expect <- densitln`
