# `src/hfbtho_projections.f90` 说明

本模块实现 HFB 态的对称性恢复，包括角动量投影（AMP）、粒子数投影（PNP）以及宇称投影。通过在高斯-勒让德网格上积分转动角度 $\beta$ 与规范角度 $\phi$，计算投影后的能量、半径、配对能隙等物理观测量。

---

### `HFBTHO_restore`

##### 函数定义

- 函数位置: `197-311 lines`

```fortran
Subroutine HFBTHO_restore()
```

##### 物理公式

- 物理含义: 对称性恢复的主驱动子程序，公开入口点。对每一个转动角度 $\beta$ 依次调用旋转矩阵构建、重叠计算、密度计算、能量计算，最后由 rank 0 完成投影与输出。

##### 实现步骤

1. 调用 `initialize_projections` 读取参数、分配全局数组、建立网格与 MPI 任务划分。
2. 调用 `initialize_wavefunctions` 构建简单旋-y 基下的 $U/V$ 矩阵。
3. 遍历本进程负责的 $\beta$ 角度（MPI 负载均衡）：
   - `initialize_angle`：分配当前角度所需的临时数组。
   - `calculate_rotation_matrix(ibet)`：计算 $R = e^{-i\beta\hat{j}_y}$。
   - `calculate_overlaps(ibet)`：计算重叠与旋转后的密度矩阵。
   - `calculate_densities(ibet)`：计算坐标空间旋转密度。
   - `calculate_energies(ibet)`：计算各能量分量。
   - `finalize_angle`：释放当前角度临时数组。
4. MPI 收集：通过 `MPI_Gatherv` 将各进程的 `all_overlaps`、`all_energies`、`beta_active` 汇总到 rank 0。
5. rank 0 调用 `project` 进行 $J$ 与 $(N,Z)$ 投影，随后调用 `print_project` 输出结果。
6. 调用 `finalize_projections` 释放全部模块级数组。

##### 变量对应

- 输入变量: 无显式参数，全部通过模块全局变量与 `HFBTHO` 模块的 namelist 参数传递。
- 输出变量: 结果写入 `projections.out` 文件，并通过 `projected_*` 数组存储投影后的物理量。

##### 函数数据流

- 调用情况: `hfbtho_library::Main_Program <- HFBTHO_restore <- initialize_projections, initialize_wavefunctions, initialize_angle, calculate_rotation_matrix, calculate_overlaps, calculate_densities, calculate_energies, finalize_angle, project, print_project, finalize_projections`

---

### `initialize_projections`

##### 函数定义

- 函数位置: `315-656 lines`

```fortran
Subroutine initialize_projections()
```

##### 物理公式

- AMP 网格：Gauss-Legendre 积分

$$
\int_0^{\pi/2} d\beta \sin\beta \; f(\beta) \quad (\text{force\_parity}=\text{.True.}) \quad \text{或} \quad \int_0^{\pi} d\beta \sin\beta \; f(\beta)
$$

通过 `gauleg` 在 $\cos\beta$ 上取点：$x_i \in [-1,1]$ 或 $[0,1]$，权重 $w_i$，然后 $\beta_i = \arccos(x_i)$。

- PNP 规范角网格：

$$
\phi_k = \frac{(k-1)\pi}{\max\phi}, \quad k=1,\dots,\max\phi
$$

- 规范相位因子：

$$
e^{\pm i\phi_k}, \qquad e^{-i\phi_k(N_0 + 2\Delta N)}
$$

其中 $N_0$ 为参考粒子数（`tz(it)`）。

##### 实现步骤

1. 读取 namelist 参数：`maxphi`, `maxbet`, `maxj`, `maxN`, `maxP`。
2. 根据 `force_parity` 与 `AMP_is_on` 调整实际计算的 `maxbet`；根据 `PNP_is_on` 调整 `maxphi`。
3. 调用 `gauleg` 建立 $\beta$ 角的 Gauss-Legendre 网格与权重。
4. 建立 $\phi$ 角等间距网格，预计算 `ephi`, `ephic`, `ephicN`, `ephicP`。
5. 分配所有全局数组（旋转矩阵、密度、能量、投影结果等）。
6. 调用 `simplex_basis(.False.)` 建立简单旋-y 基量子数。
7. 建立圆柱角 $\phi_{\text{cyl}}$ 网格与 $(r_\perp,z)\times\phi_{\text{cyl}}$ 转换索引。
8. 确定所有非退化的振子态对 (`iosc1_pair`, `iosc2_pair`)。
9. 确定 `iit-iphi` 配对列表。
10. MPI 任务划分：按 `maxbet` 均匀分配角度块给各进程，计算 `block_vec`, `beta_start`, `offset`。
11. 广播 `piu(it)`、Bogoliubov 矩阵等数据。
12. 分配本进程局部的 `all_overlaps`、`all_energies` 及 rank 0 的汇总数组。

##### 变量对应

- 输入变量: 通过 `HFBTHO` 模块全局变量传递（`number_of_gauge_points`, `number_of_rotational_angles`, `maximal_angular_momentum`, `delta_neutrons`, `delta_protons`, `force_parity`, `AMP_is_on`, `PNP_is_on` 等）。
- 输出变量: 初始化全部模块级 Save 数组与网格参数（`betabs`, `betaweight`, `phiabs`, `ephi`, `ephic`, `kdim`, `block_vec`, `offset` 等）。

##### 函数数据流

- 调用情况: `HFBTHO_restore <- initialize_projections <- simplex_basis, gauleg`

---

### `initialize_wavefunctions`

##### 函数定义

- 函数位置: `660-740 lines`

```fortran
Subroutine initialize_wavefunctions()
```

##### 物理公式

- 将实数 Bogoliubov 矩阵 `RUqpN/P`, `RVqpN/P` 转换为复数简单旋-y 基下的 $U$ 和 $V$ 矩阵：

$$
U_{\alpha k} = \begin{cases}
\text{RUqp}(k,\alpha) & n_{l,\text{sim}}(\alpha) \ge 0 \\
i\,\text{RUqp}(k,\alpha) & n_{l,\text{sim}}(\alpha) < 0
\end{cases}
$$

$$
V_{\alpha k} = -\begin{cases}
\text{RVqp}(k,\alpha) & n_{l,\text{sim}}(\alpha) \ge 0 \\
i\,\text{RVqp}(k,\alpha) & n_{l,\text{sim}}(\alpha) < 0
\end{cases}
$$

- 两套矩阵：
  - `VmatrixN1/P1`, `UmatrixN1/P1`：截断后的 qp 维度（`kdim(1:2)`）。
  - `VmatrixN2/P2`, `UmatrixN2/P2`：完整块维度（`kdim(3:4)=nt`），用于重叠计算。

##### 变量对应

- 输入变量: `RUqpN`, `RVqpN`, `RUqpP`, `RVqpP`, `KpwiN`, `KpwiP`（来自 `HFBTHO` 模块），以及 `nl_sim`（来自 `simplex_basis`）。
- 输出变量: `VmatrixN1`, `UmatrixN1`, `VmatrixN2`, `UmatrixN2`, `VmatrixP1`, `UmatrixP1`, `VmatrixP2`, `UmatrixP2`。

##### 函数数据流

- 调用情况: `HFBTHO_restore <- initialize_wavefunctions`

---

### `calculate_rotation_matrix`

##### 函数定义

- 函数位置: `794-865 lines`

```fortran
Subroutine calculate_rotation_matrix(ibet)
```

##### 物理公式

- 转动矩阵元在简单旋-y 基下的表示：

$$
R_{ab}(\beta) = \cos\frac{\beta}{2} \cdot \text{rotel}_1 + i\sin\frac{\beta}{2} \cdot \text{rotel}_3
$$

其中 `rotel_1` = $\langle a | e^{-i\beta j_y} | b \rangle_{\Lambda_b}$，`rotel_3` = $\langle a | e^{-i\beta j_y} | b \rangle_{-\Lambda_b}$，由 `calculate_ry` 计算。

- 对称关系：

$$
R_{ba}(\beta) = i^{l_a + l_b} R_{ab}(\beta)
$$

代码中通过 `iv(nla1+nla2)` 实现（`iv(n) = i^n`）。

- 行列式：

$$
\det R(\beta) = \det\bigl(R_{ab}\bigr)
$$

通过 `calculate_inverse_cmplx` 计算逆矩阵时得到行列式 `det`，存入 `detR(ibet)`。

##### 实现步骤

1. 取当前角度 $\beta = \text{betabs(ibet)}$。
2. OpenMP 并行遍历所有振子态对 `(iosc1,iosc2)`：
   - 调用 `calculate_ry` 分别计算 `nl2 = +|nl2|` 和 `nl2 = -|nl2|` 的矩阵元 `rotel1` 与 `rotel3`。
   - 按上式组装 `rotation_matrix(iosc1,iosc2)` 与 `rotation_matrix(iosc2,iosc1)`。
3. 复制 `rotation_matrix` 到 `inverse_rotation_matrix`，调用 `calculate_inverse_cmplx` 求逆并得到行列式。
4. 若行列式非零，标记 `beta_active=1`，否则标记为 0 并跳过该角度。

##### 变量对应

- 输入变量: `ibet (integer)`，对应公式中第 $i$ 个 $\beta$ 角；`betabs(ibet)` 对应 $\beta$；`bp`, `bz`（振子长度）。
- 输出变量: `rotation_matrix(nt,nt)` 对应 $R_{ab}$；`inverse_rotation_matrix(nt,nt)` 对应 $R^{-1}$；`detR(ibet)` 对应 $\det R$。

##### 函数数据流

- 调用情况: `HFBTHO_restore <- calculate_rotation_matrix <- calculate_ry, calculate_inverse_cmplx`

---

### `calculate_overlaps`

##### 函数定义

- 函数位置: `869-1127 lines`

```fortran
Subroutine calculate_overlaps(ibet)
```

##### 物理公式

**AMP only** (`PNP_is_on == 0`)：

定义旋转后的 $U$ 和 $V$：

$$
\tilde{V} = R V, \qquad \tilde{U} = R^{-1\dagger} U
$$

重叠矩阵：

$$
A = U^\dagger \tilde{U}^* + V^\dagger \tilde{V}^*
$$

代码实现：
- `aUmatrix = U^T * conjg(\tilde{U})`
- `aVmatrix = V^T * conjg(\tilde{V})`
- `OverlapMatrixTmp = aUmatrix + aVmatrix`

HFB 态的重叠：

$$
\langle \phi | \hat{R}(\beta) | \phi \rangle = \det(A) \cdot \det(R)
$$

**PNP or PNP&AMP** (`PNP_is_on > 0`)：

引入规范相位 $e^{\pm i\phi}$，按简单旋块区分：

$$
A_{++} = e^{-i\phi} U^\dagger \tilde{U}^* + e^{i\phi} V^\dagger \tilde{V}^*
$$

$$
A_{--} = e^{-i\phi} \bigl(U^\dagger \tilde{U}^*\bigr)^* + e^{i\phi} \bigl(V^\dagger \tilde{V}^*\bigr)^*
$$

代码中对应：
- `isimplex == 1`: `ephic(iphi)*aUmatrix + ephi(iphi)*aVmatrix`
- `isimplex == 2`: `ephic(iphi)*conjg(aUmatrix) + ephi(iphi)*conjg(aVmatrix)`

重叠存储为 `detA(1) * detR(ibet) * prefac`，其中 `prefac = exp(i*phi*kdim)`。

**旋转后的密度矩阵**：

$$
\rho = \tilde{V} A^{-1} V^\dagger, \qquad
\kappa = -\tilde{V} A^{-1} U^\dagger, \qquad
\kappa^C = \tilde{U} A^{-1} V^\dagger
$$

对 PNP 情况，密度额外乘以相应规范相位因子 `ephi(iphi)` 或 `ephic(iphi)`。

##### 实现步骤

1. 若 `beta_active == 0` 则直接返回。
2. 根据 AMP/PNP 分支分配局部数组（`rotated_Vmatrix`, `rotated_invUmatrix`, `InverseOverlapMatrix` 等）。
3. OpenMP 并行计算所有 `iit`（1=中子截断, 2=质子截断, 3=中子全, 4=质子全）与 `iphi` 组合：
   - `Zgemm` 计算 `\tilde{V} = R V` 与 `\tilde{U} = R^{-1\dagger} U`。
   - `Zgemm` 计算 `aUmatrix` 与 `aVmatrix`。
   - 组装 `OverlapMatrixTmp` 并调用 `calculate_inverse_cmplx` 求逆及行列式。
   - 保存 `InverseOverlapMatrix` 与重叠值。
4. OpenMP 并行计算旋转后的密度与配对张量 `\rho`, `\kappa`, `\kappa^C`（简单旋 ++ 与 --）。
5. 筛选出对密度/配对有非零贡献的振子态对，存入 `iosc1_contributing`, `iosc2_contributing`, `nt_contributing`。

##### 变量对应

- 输入变量: `ibet (integer)`，对应第 $i$ 个 $\beta$ 角；`rotation_matrix`, `inverse_rotation_matrix` 对应 $R$, $R^{-1}$；`VmatrixN1/P1/N2/P2`, `UmatrixN1/P1/N2/P2` 对应 $V$, $U$。
- 输出变量: `rotated_density`, `rotated_kappa`, `rotated_kappac` 对应 $\rho$, $\kappa$, $\kappa^C$；`all_overlaps` 存储各 `ibet`, `iphi` 下的重叠 $\det A \cdot \det R$。

##### 函数数据流

- 调用情况: `HFBTHO_restore <- calculate_overlaps <- calculate_inverse_cmplx, Zgemm`

---

### `calculate_densities`

##### 函数定义

- 函数位置: `1131-1515 lines`

```fortran
Subroutine calculate_densities(ibet)
```

##### 物理公式

- 对每个空间网格点 $(r_\perp, z, \phi_{\text{cyl}})$ 和每个规范角 `iphi`，通过振子基函数 `qhla_opt`, `fi1z_opt`, `fi1r_opt` 累加各密度分量。

- 时间偶（time-even）密度：

$$
\rho(\mathbf{r}) = \sum_{ab} \bigl(\rho^{++}_{ab} + \rho^{--}_{ab}\bigr) \, \phi_a^* \phi_b
$$

角向依赖为 $\cos\bigl((n_{l2}-n_{l1})\phi_{\text{cyl}}\bigr)$。

- 时间奇（time-odd）密度（自旋 $\mathbf{s}$、流 $\mathbf{j}$ 等）的角向依赖包含 $\sin\bigl((n_{l1}+n_{l2}+1)\phi_{\text{cyl}}\bigr)$、$\cos\bigl((n_{l1}+n_{l2}+1)\phi_{\text{cyl}}\bigr)$ 等。

- 坐标空间归一化：所有密度乘以

$$
\text{piu}(it) \times \text{wdcori}(i_{\text{hil}})
$$

其中 $\text{piu}(it) = t_z(it) / (2\sum \text{ro\_normalization})$。

- 投影密度（`PNP_is_on == 2`）：

$$
\rho_{\text{proj}}(it) = \frac{\sum_{\phi} e^{-i N_0 \phi} \langle \phi | \hat{R} | \phi \rangle \, \rho_\phi(it)}{\sum_{\phi} e^{-i N_0 \phi} \langle \phi | \hat{R} | \phi \rangle}
$$

- 旋转库仑场：

$$
V_{\text{C,rotated}}(\mathbf{r}) = \frac{1}{2\pi} \int_0^{2\pi} d\phi_{\text{cyl}} \int d\mathbf{r}' \, v_C(\mathbf{r},\mathbf{r}') \, \rho_\phi^{(p)}(\mathbf{r}')
$$

通过 `vc` 矩阵与 `zgemv` 实现。

##### 实现步骤

1. 若 `beta_active == 0` 则返回。
2. OpenMP 并行遍历所有 $(r_\perp,z)\times\phi_{\text{cyl}}$ 网格点。
3. 对每个 `iphi`，遍历 `nt_contributing` 中筛选出的振子对：
   - 提取 `densNpl`, `densNmi`, `densPpl`, `densPmi`（密度 ++/-- 组合）以及 `kappNmi`, `kappPmi`, `kappcNmi`, `kappcPmi`（配对）。
   - 由 HO 基函数乘积构建各类几何因子 `fac_ro`, `fac_tau`, `fac_dj`, `fac_dro`, `fac_aka`, `fac_sroz`, `fac_jz` 等。
   - 累加时间偶、时间奇及张量密度分量。
4. 若 `PNP_is_on == 2`，利用 `all_overlaps` 与 `ephicN/P` 计算 `ro_projected`。
5. 对所有密度乘以归一化因子 `piu(it) * wdcori(ihil)`。
6. 通过 `vc` 矩阵计算旋转后的库仑场 `cou_rotated`。

##### 变量对应

- 输入变量: `ibet (integer)`；`rotated_density`, `rotated_kappa`, `rotated_kappac` 对应旋转后的密度矩阵元；`qhla_opt`, `fi1z_opt`, `fi1r_opt` 对应 HO 基函数。
- 输出变量: `rotated_ro`, `rotated_tau`, `rotated_dj`, `rotated_dro`, `rotated_aka`, `rotated_akac` 等对应各坐标空间密度；`ro_projected` 对应投影密度；`cou_rotated` 对应旋转库仑势。

##### 函数数据流

- 调用情况: `HFBTHO_restore <- calculate_densities`

---

### `calculate_energies`

##### 函数定义

- 函数位置: `1519-1897 lines`

```fortran
Subroutine calculate_energies(ibet)
```

##### 物理公式

对每个规范角对 $(\phi_n, \phi_p)$，在各空间网格点构造同位旋标量/矢量密度后，计算 Skyrme EDF 的能量密度（与 `hfbtho_unedf.f90` 采用相同泛函形式）：

- 动能：

$$
E_{\text{kin}}^{(n)} = \frac{\hbar^2}{2m_n} \tau_n \, \text{facECM}, \qquad
E_{\text{kin}}^{(p)} = \frac{\hbar^2}{2m_p} \tau_p \, \text{facECM}
$$

- 体积能：

$$
E_{\text{vol}} = C_{\rho}^{(0)} \rho_0^2 + C_{\rho}^{(1)} \rho_1^2 + C_{\text{dr}\rho}^{(0)} \rho_0^{\sigma+2} + C_{\text{dr}\rho}^{(1)} \rho_1^{\sigma+2} + C_{\tau}^{(0)} \rho_0 \tau_0 + C_{\tau}^{(1)} \rho_1 \tau_1
$$

- 表面能：

$$
E_{\text{surf}} = C_{r\Delta r}^{(0)} \rho_0 \Delta\rho_0 + C_{r\Delta r}^{(1)} \rho_1 \Delta\rho_1
$$

- 自旋-轨道能：

$$
E_{\text{so}} = C_{r\nabla J}^{(0)} \rho_0 \nabla\cdot\mathbf{J}_0 + C_{r\nabla J}^{(1)} \rho_1 \nabla\cdot\mathbf{J}_1
$$

- 库仑能（直接 + 交换 Slater）：

$$
E_{\text{Coul,direct}} = \frac{1}{2} V_C \rho_p, \qquad
E_{\text{Coul,exch}} = -C_{\text{ExPar}} \, c_{\text{ex}} \, \rho_p^{4/3}
$$

- 配对能（中子/质子）：

$$
E_{\text{pair}}^{(it)} = C_{pV0}^{(it)} \Bigl(1 - \frac{\rho_0}{\rho_c} C_{pV1}^{(it)}\Bigr) \, \tilde{\rho}^2
$$

- 时间奇项：$\mathbf{s}\cdot\mathbf{s}$、$\mathbf{j}\cdot\mathbf{j}$、$\mathbf{s}\cdot\Delta\mathbf{s}$、$\mathbf{J}\cdot\nabla\times\mathbf{s}$ 等。
- 张量项（`use_j2terms`）：由 `SFIZ`, `SFIR`, `SZZ`, `SRZ` 等构建 `TENS` 与 `ST`。

积分方式：先对圆柱角 $\phi_{\text{cyl}}$ 调用 `integrate_complex_function` 积分，再对 $(z, r_\perp)$ 用高斯权重 `wdcor(ihil)` 积分。

结果按 `stride_size=23` 存入 `all_energies`：
1-2: 动能, 3-4: 库仑, 5-6: 体积, 7: 表面, 8: 自旋轨道, 9-10: 配对, 11-15: 时间奇, 16-17: 张量, 18-19: 粒子数, 20-21: 半径, 22-23: 配对能隙。

##### 变量对应

- 输入变量: `ibet (integer)`；各类 `rotated_*` 密度数组；`cou_rotated` 库仑势；Skyrme 耦合常数（来自 `UNEDF` 模块）。
- 输出变量: `all_energies(stride_size * betaphi_size)`，按 `(iphin, iphip)` 顺序存储 23 个能量/观测量的复数值。

##### 函数数据流

- 调用情况: `HFBTHO_restore <- calculate_energies <- integrate_complex_function`

---

### `project`

##### 函数定义

- 函数位置: `1916-2177 lines`

```fortran
Subroutine project()
```

##### 物理公式

- 角动量投影（AMP）：对给定 $J$，利用 Wigner $d$-函数 `wigner(j,0,0,beta)` 积分：

$$
\langle J | \hat{O} | J \rangle = \frac{2J+1}{2} \int_0^{\pi} d\beta \, \sin\beta \, d^J_{00}(\beta) \, \langle \phi | \hat{O} \, e^{-i\beta\hat{J}_y} | \phi \rangle
$$

代码中 `facj = half*(2*jj+1)*integration_prefactor`，`betafac = betaweight(ibet) * wigner(jj,0,0,beta)`。

- 宇称投影（`force_parity`）：

$$
\hat{P}_\pi = \frac{1}{2}\bigl(1 + \pi \hat{\Pi}\bigr)
$$

代码通过 `parity_prefac=0.5` 与 `parity_prefac_pi=0.5`（或 `force_parity` 时 1.0 / 0.0）实现；对负宇称通道用 `iv(jj) * rotated_overlap(maxbet+1-ibet, ...)`。

- 粒子数投影（PNP）：进一步对规范角积分：

$$
\langle J, N_n, N_p | \hat{O} | J, N_n, N_p \rangle = \frac{1}{(2\pi)^2} \int d\phi_n \int d\phi_p \, e^{-i(N_n\phi_n + N_p\phi_p)} \, \langle \phi | \hat{O} \, e^{-i\beta\hat{J}_y} | \phi \rangle_{\phi_n,\phi_p}
$$

代码中通过 `ephicN(iphin,nn_neut) * ephicP(iphip,nn_prot)` 实现相位，并除以 `maxphi**2`。

- 投影后的重叠（范数）：

$$
c_J = \text{facj} \times \sum_{\beta} \text{betafac} \times c_{J,\phi}
$$

其中 `projected_overlap(jj) = Dble(cj) * facj`。所有能量观测量均以 $c_J$ 为分母进行归一化。

##### 实现步骤

1. 从 `all_overlaps_gthr` 解包得到 `rotated_overlap(ibet,iphi,it)`。
2. 从 `all_energies_gthr` 解包得到各 `ibet, iphin, iphip` 下的能量分量（`ekinN_rotated` 等）。
3. 对 $J=0,\dots,J_{\max}$（步长 `jjstep`）：
   - 对每个 $\beta$ 角，计算 `betafac`。
   - 对每对 $(\phi_n, \phi_p)$，计算重叠因子 `fac_overlap` 与宇称修正 `fac_overlap_pi`。
   - 累加各能量分量乘以重叠因子，得到 `cj_phi` 与各 `*_phi`。
   - 将 `*_phi` 除以 `maxphi**2` 后，按 `betafac` 累加到 `cj` 与各总量。
   - 最后以 `facj` 归一化并存入 `projected_*` 数组。
4. 若 `PNP_is_on > 0`，进一步计算各 $(N,Z)$ 分解 `projected_NP(jj,nn_neut,nn_prot)` 与归一化 `projected_NP_norm`。

##### 变量对应

- 输入变量: `all_overlaps_gthr` 对应各角度下的重叠；`all_energies_gthr` 对应各角度下的能量；`betabs`, `betaweight`；`ephicN`, `ephicP`。
- 输出变量: `projected_overlap(0:maxj)` 对应 $c_J$；`projected_ekinN`, `projected_ekinP`, `projected_ecodi`, `projected_ecoex` 等对应各投影能量分量；`projected_NP` 对应粒子数分解；`projected_NP_norm` 对应归一化常数。

##### 函数数据流

- 调用情况: `HFBTHO_restore <- project <- wigner, iv`

---

### `print_project`

##### 函数定义

- 函数位置: `2181-2404 lines`

```fortran
Subroutine print_project()
```

##### 物理公式

- 按投影后的 $J$ 值输出各能量分量：

$$
E_{\text{even}} = E_{\text{kin}} + E_{\text{Coul}} + E_{\text{vol}} + E_{\text{surf}} + E_{\text{so}} + E_{\text{pair}}
$$

$$
E_{\text{odd}} = E_{ss} + E_{ss}^{(DD)} + E_{jj} + E_{s\Delta s} + E_{j\nabla\times s}
$$

总投影能量 $E_{\text{proj}} = E_{\text{even}} + E_{\text{odd}}$（若 `use_j2terms` 则加入张量项）。

- 电荷半径（含质子有限尺寸、中子电磁均方位移、Darwin-Foldy 修正）：

$$
r_{\text{ch}} = \sqrt{r_p^2 + 0.769 + \frac{N}{Z}(-0.1161) + 0.033}
$$

##### 实现步骤

1. 根据 `AMP_is_on` 与 `PNP_is_on` 的组合，分别输出：
   - 仅 AMP：$\beta$ 网格、重叠、$J$ 分解 $|c_J|^2$。
   - 仅 PNP：$\phi$ 网格、$(N,Z)$ 分解。
   - AMP+PNP：$\beta$ 与 $\phi$ 网格、$J$ 分解、各 $J$ 下的 $(N,Z)$ 分解。
2. 输出各 $J$（或 $J=0$）下的投影能量明细（动能、体积、表面、自旋轨道、库仑、配对、时间奇、张量）。
3. 输出均方根半径、电荷半径、配对能隙。
4. 输出各子程序的系统时钟耗时统计。

##### 变量对应

- 输入变量: 所有 `projected_*` 数组（`projected_overlap`, `projected_ekinN`, `projected_ekinP`, `projected_ecodi`, `projected_ecoex`, `projected_EVOL_rho_tau`, `projected_EVOL_rho_rho`, `projected_ESURF_rho_drho`, `projected_ESO_rho_nablaj`, `projected_eptN`, `projected_eptP`, `projected_ess0`, `projected_essD`, `projected_ecurr`, `projected_esDs`, `projected_ejNxs`, `projected_ETENS`, `projected_esT`, `projected_xn1`, `projected_xn2`, `projected_rms1`, `projected_rms2`, `projected_delN`, `projected_delP`, `projected_NP`, `projected_NP_norm`）。
- 输出变量: 写入文件 `projections.out`（单元号 `lproj`）。

##### 函数数据流

- 调用情况: `HFBTHO_restore <- print_project`

---

### `finalize_projections`

##### 函数定义

- 函数位置: `2408-2434 lines`

```fortran
Subroutine finalize_projections()
```

##### 物理公式

- 无特定物理公式，为资源清理子程序。

##### 实现步骤

1. 关闭输出文件 `lproj`。
2. 依次释放 `initialize_projections` 中分配的全部模块级数组，包括：
   - 角度网格：`betabs`, `betaweight`, `phiabs`, `ephi`, `ephic`, `ephicN`, `ephicP`
   - 重叠与旋转矩阵：`rotated_overlap`, `detR`
   - 投影结果：`projected_overlap`, `projected_ekinN`, ..., `projected_NP`, `projected_NP_norm`
   - 简单旋基：`nz_sim`, `nr_sim`, `nl_sim`
   - 空间网格转换：`phicyl`, `ihil_convert`, `iphicyl_convert`, `ihil_iphicyl_convert`
   - 振子对索引：`iosc1_pair`, `iosc2_pair`, `itiphi_pair1`, `itiphi_pair2`
   - 波函数矩阵：`VmatrixN1`, `UmatrixN1`, ..., `VmatrixP2`, `UmatrixP2`
   - MPI 相关：`all_overlaps`, `all_energies`, `beta_active`, `offset`, `block_vec`, `beta_start`, `all_overlaps_gthr`, `all_energies_gthr`, `beta_active_gthr`
   - 额外：`ro_normalization`, `xl_ihil`, `xh_ihil`

##### 变量对应

- 输入变量: 无显式参数。
- 输出变量: 无，仅释放内存。

##### 函数数据流

- 调用情况: `HFBTHO_restore <- finalize_projections`

---

### `simplex_basis`

##### 函数定义

- 函数位置: `2439-2472 lines`

```fortran
Subroutine simplex_basis(lpr)
```

##### 物理公式

- 将标准 HO 基 $|n_z, n_r, n_l, n_s\rangle$ 映射到简单旋-y 基，量子数变换为：

$$
n_{z,\text{sim}} = n_z, \qquad n_{r,\text{sim}} = n_r
$$

$$
n_{l,\text{sim}} = \begin{cases}
+n_l & n_s = +1 \\
-n_l & n_s = -1
\end{cases}
$$

代码中 `ns(iosc) == 1` 时 `nl_sim = nl`，否则 `nl_sim = -nl`。

##### 变量对应

- 输入变量: `lpr (logical)`，控制是否打印基矢对应表；`nz`, `nr`, `nl`, `ns`（来自 `HFBTHO` 模块）。
- 输出变量: `nz_sim`, `nr_sim`, `nl_sim`。

##### 函数数据流

- 调用情况: `initialize_projections <- simplex_basis`

---

### `calculate_ry`

##### 函数定义

- 函数位置: `2490-2643 lines`

```fortran
Subroutine calculate_ry(bb,bp1,bz1,bp2,bz2,nz1,nr1,nl1,nz2,nr2,nl2,rotel)
```

##### 物理公式

- 计算 HO 态之间的转动矩阵元：

$$
\langle n_z n_r \Lambda | e^{-i\beta\hat{j}_y} | n'_z n'_r \Lambda' \rangle
$$

基于 Nazmitdinov 表示（Eq. 47），要求总壳量子数守恒：

$$
N = 2n_r + |\Lambda| + n_z = 2n'_r + |\Lambda'| + n'_z
$$

- 构造辅助矩阵 $F$, $G$, $K$ 以及 $M$ 矩阵：

$$
M = \frac{1}{2}\begin{pmatrix}
\bigl(\frac{b_{p1}}{b_{p2}}\cos\beta\bigr)^2 + \bigl(\frac{b_{p1}}{b_{z2}}\sin\beta\bigr)^2 + 1 & 0 & \bigl(-\frac{b_{z1}}{b_{p2}}\frac{b_{p1}}{b_{p2}} + \frac{b_{z1}}{b_{z2}}\frac{b_{p1}}{b_{z2}}\bigr)\cos\beta\sin\beta \\
0 & \bigl(\frac{b_{p1}}{b_{p2}}\bigr)^2 + 1 & 0 \\
\bigl(-\frac{b_{z1}}{b_{p2}}\frac{b_{p1}}{b_{p2}} + \frac{b_{z1}}{b_{z2}}\frac{b_{p1}}{b_{z2}}\bigr)\cos\beta\sin\beta & 0 & \bigl(\frac{b_{z1}}{b_{z2}}\cos\beta\bigr)^2 + \bigl(\frac{b_{z1}}{b_{p2}}\sin\beta\bigr)^2 + 1
\end{pmatrix}
$$

通过 `calculate_inverse_real` 计算 $\det M$。

- 预因子：

$$
\text{prefac} = i^{n_{p1}+n_{p2}+n_{r1}+n_{r2}} \, \frac{b_{p1}}{b_{p2}}\sqrt{\frac{b_{z1}}{b_{z2}}} \, \frac{1}{\sqrt{\det M \, 2^{N+N'}}} \, \sqrt{n_{p1}!\,n_{m1}!\,n_{z1}!\,n_{p2}!\,n_{m2}!\,n_{z2}!}
$$

其中 $n_{p}, n_{m}$ 由 $\Lambda$ 的符号决定（$\Lambda\ge 0$ 时 $n_p=n_r+\Lambda$, $n_m=n_r$；$\Lambda<0$ 时 $n_p=n_r$, $n_m=n_r+|\Lambda|$）。

- 通过多重求和（`mm1-mm6`, `ii1`, `jj1`, `kk1`, `ll1`）计算 $F$、$G$、$K$ 矩阵的幂次展开贡献 `bfac`, `cfac`, `afac`，最终：

$$
\text{rotel} = \text{prefac} \times \sum_{\{m\}} b_{\text{fac}} \, a_{\text{fac}} \, c_{\text{fac}}
$$

##### 变量对应

- 输入变量:
  - `bb` 对应 $\beta$
  - `bp1`, `bz1` 对应 bra 的振子长度 $b_\perp^{(1)}$, $b_z^{(1)}$
  - `bp2`, `bz2` 对应 ket 的振子长度 $b_\perp^{(2)}$, $b_z^{(2)}$
  - `nz1`, `nr1`, `nl1` 对应 bra 的 $(n_z, n_r, \Lambda)$
  - `nz2`, `nr2`, `nl2` 对应 ket 的 $(n'_z, n'_r, \Lambda')$
- 输出变量: `rotel` 对应矩阵元 $\langle n_z n_r \Lambda | e^{-i\beta\hat{j}_y} | n'_z n'_r \Lambda' \rangle$。

##### 函数数据流

- 调用情况: `calculate_rotation_matrix <- calculate_ry <- calculate_inverse_real, iv, fak, fi`
