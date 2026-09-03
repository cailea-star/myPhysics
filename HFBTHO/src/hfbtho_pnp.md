# 代码文件 ```hfbtho_pnp.f90``` 说明

该文件实现 HFB 解在 **canonical 基** 下的 **粒子数投影 (Particle Number Projection, PNP)**。核心思想是对 HFB 准粒子真空在规范角 $\phi$ 上进行高斯积分，投影到确定的粒子数 $N$ 子空间，从而计算投影后的密度、能量及多极矩等期望值。

---

### `expectpj`

##### 函数定义

- 函数位置: `147-475 lines`

```fortran
Subroutine expectpj(lpr)
```

##### 物理公式

- **物理含义**: 粒子数投影 (PNP) 下 HFB 态的期望值计算主入口。对中子和质子分别选取规范角网格 $(\phi_n, \phi_p)$，计算投影后的能量、多极矩及均方根半径等。投影公式为

$
|N_{\text{proj}}\rangle = \int_0^{2\pi} d\phi\; e^{-iN\phi} |\text{HFB}(\phi)\rangle
$

其中 $|\text{HFB}(\phi)\rangle$ 为规范旋转后的 HFB 真空。

- **投影能量公式**: 总投影能量由规范角上的积分给出

$
E_{\text{proj}} = \frac{\displaystyle\int d\phi_n \, d\phi_p \; E(\phi_n, \phi_p) \, \langle \text{HFB}|e^{i\hat{N}_n\phi_n} e^{i\hat{N}_p\phi_p}|\text{HFB}\rangle}
                 {\displaystyle\int d\phi_n \, d\phi_p \; \langle \text{HFB}|e^{i\hat{N}_n\phi_n} e^{i\hat{N}_p\phi_p}|\text{HFB}\rangle}
$

或等价地，代码中用权重因子 `pijk = pjk(i,1)*pjk(j,2)`（即归一化后的 norm overlap）对能量密度进行加权求和。

- **投影后各能量分量**: 对每一对规范角 $(i,j)$，先计算局域能量密度贡献
  - 体积项: $E_{\text{vol}} = \int d\mathbf{r}\; \big[ (tv_1\rho_t^2 - tv_2\rho_{np}^2) + (tv_3\rho_t^2 - tv_4\rho_{np}^2)\rho_t^\sigma + tv_5\rho_t\tau + tv_6(\rho_n\tau_n + \rho_p\tau_p) \big]$
  - 表面项: $E_{\text{sur}} = \int d\mathbf{r}\; \big[ ts_1\rho_t\Delta\rho_t + ts_2(\rho_n\Delta\rho_n + \rho_p\Delta\rho_p) \big]$
  - 自旋轨道项: $E_{\text{so}} = \int d\mathbf{r}\; \big[ C_{rdJ}^{(0)}\rho_t \pm C_{rdJ}^{(1)}\rho_{np} \big] \mathbf{J}$
  - 库仑直接项: $E_{\text{Coul,d}} = \frac{1}{2}\int d\mathbf{r}\; V_C(\mathbf{r})\rho_p(\mathbf{r})$
  - 库仑交换项 (Slater): $E_{\text{Coul,e}} = -C_{ex}\int d\mathbf{r}\; \rho_p^{4/3}$
  - 对能: $E_{\text{pair}} = \int d\mathbf{r}\; \Delta^{(it)}(\mathbf{r}) \kappa^{(it)}(\mathbf{r})$
  - 动能: $E_{\text{kin}} = \int d\mathbf{r}\; \frac{\hbar^2}{2m}\tau$

  然后用 `pijk` 加权累加得到投影后的各分量。

- **实现步骤**:
  1. 调用 `allocate_arrays_PNP()` 分配 PNP 相关数组。
  2. 调用 `densitpj()` 计算复数投影密度 `ropj`, `taupj`, `dropj`, `djpj`, `akapj`。
  3. 双重循环遍历中子规范角 `ign = 1..ilpjn` 和质子规范角 `igp = 1..ilpjp`。
  4. 对每个角点对，读取投影密度并在高斯网格上积分，累加各能量项、粒子数、rms、四极矩及十六极矩。
  5. 若 `icstr /= 0`，将约束项贡献加入半投影能量 `epj`。
  6. 计算实数化的投影能量 `retotpj`、HFB 能量 `rehfbcan`、关联能 `depnp = retotpj - rehfbcan`。
  7. 若 `lpr=.true.` 则将结果输出到屏幕及 `akzout.dat`。
  8. 调用 `deallocate_arrays_PNP()` 释放数组。

##### 变量对应

- 输入变量 `lpr` (Logical): 是否打印结果到输出文件。
- 输出变量 `retotpj` (Real): 投影总能量 $E_{\text{proj}}$，对应公式中的积分加权总能量实部。
- 输出变量 `rehfbcan` (Real): HFB 在 canonical 基下的总能量（对应 $\phi_n=\phi_p=0$ 时的 `cetot` 实部）。
- 输出变量 `depnp` (Real): 投影关联能 $E_{\text{proj}} - E_{\text{HFB(can)}}$。
- 局部变量 `pijk` (Complex): 当前角点对的归一化重叠权重，公式对应
  $
  pijk = \frac{\langle\text{HFB}|e^{i\hat{N}_n\phi_i}e^{i\hat{N}_p\phi_j}|\text{HFB}\rangle}
              {\sum_{i',j'}\langle\text{HFB}|e^{i\hat{N}_n\phi_{i'}}e^{i\hat{N}_p\phi_{j'}}|\text{HFB}\rangle}
  $
- 局部变量 `ropj(ihli,i,it)` (Complex): 投影密度 $\rho_{\text{proj}}^{(it)}(\mathbf{r}_{ihli}, \phi_i)$。
- 局部变量 `taupj`, `dropj`, `djpj`, `akapj` (Complex): 投影动能密度、拉普拉斯密度、自旋密度、对密度。
- 局部变量 `coupj` (Complex): 投影库仑势 $V_C(\mathbf{r})$。
- 局部数组 `epj(ilpjmax,2)` (Complex): 半投影能量，用于 VAP (Variation After Projection) 计算。

##### 函数数据流

- 调用情况: `hfbtho_solver.f90 (line 5556) <- expectpj <- allocate_arrays_PNP, densitpj, deallocate_arrays_PNP`

---

### `densitpj`

##### 函数定义

- 函数位置: `481-773 lines`

```fortran
Subroutine densitpj()
```

##### 物理公式

- **物理含义**: 在高斯积分网格上计算复数投影密度。对每种同位旋 `it` (1=中子, 2=质子) 和每个规范角 $\phi_l$，先构造相位因子，再计算 norm overlap 和各单粒子态的投影权重，最后将波函数按权重叠加得到空间密度。

- **Norm overlap (范数重叠)**:

$
\langle\text{HFB}|e^{i\hat{N}\phi}|\text{HFB}\rangle
= \prod_k \big( u_k^2 + v_k^2 e^{2i\phi} \big) \cdot e^{-iN\phi}
$

代码中 `ppjk(ig)` 即为此重叠。归一化后

$
ppjk(ig) \leftarrow \frac{ppjk(ig)}{\sum_{ig'} ppjk(ig')}
$

- **投影单粒子权重**:

$
\rho_k(\phi) = \frac{v_k^2 e^{2i\phi}}{u_k^2 + v_k^2 e^{2i\phi}}
$

对应代码中 `prpj(kkk,ig) = v2ig/dig`。

- **投影密度**:

$
\rho_{\text{proj}}(\mathbf{r}, \phi) = \sum_k \rho_k(\phi) \, |\phi_k(\mathbf{r})|^2
$

其中 $|\phi_k(\mathbf{r})|^2$ 为 canonical 波函数的模方。类似地可定义动能密度 $\tau$、对密度 $\kappa$、自旋密度 $\mathbf{J}$ 等。

- **实现步骤**:
  1. 对每种同位旋 `it`，根据网格数 `ilpj` 生成相位因子 `exp1iphy`, `exp2iphy` 等。
  2. 初始化 `ppjk` 为 $e^{-iN\phi}$，再逐个 canonical 态乘入因子 $(u_k^2 + v_k^2 e^{2i\phi})$。
  3. 对 PWI (pairing window) 截断内的态，计算投影权重 `prpj`、`pupj`、`pcpj`、`pypj`。
  4. 在高斯网格 $(r_l, z_h)$ 上，将 canonical 波函数及其导数按权重叠加，得到投影后的 $\rho, \tau, \Delta\rho, \mathbf{J}, \kappa$ 等。
  5. 用 `dnfactor(it)` 对密度进行再归一化，使总粒子数等于 `npr(it)`。
  6. 将结果存入全局数组 `ropj`, `taupj`, `dropj`, `djpj`, `akapj` 及自旋相关数组。
  7. 调用 `coulompj()` 计算投影库仑场。

##### 变量对应

- 输入/全局变量 `vk(k,it)` (Real): canonical 基下的占据概率 $v_k^2$。
- 输出数组 `ropj(nghl, ilpjmax, 2)` (Complex): 投影密度矩阵元 $\rho_{\text{proj}}(ihli, \phi_{ig}, it)$。
- 输出数组 `taupj(nghl, ilpjmax, 2)` (Complex): 投影动能密度 $\tau_{\text{proj}}$。
- 输出数组 `dropj(nghl, ilpjmax, 2)` (Complex): 投影拉普拉斯密度 $\Delta\rho_{\text{proj}}$。
- 输出数组 `djpj(nghl, ilpjmax, 2)` (Complex): 投影自旋-轨道密度 $\mathbf{J}_{\text{proj}}$。
- 输出数组 `akapj(nghl, ilpjmax, 2)` (Complex): 投影对密度 $\kappa_{\text{proj}}$。
- 输出数组 `pjk(ilpjmax, 2)` (Complex): 归一化后的 norm overlap。
- 局部变量 `dig` (Complex): $u_k^2 + v_k^2 e^{2i\phi}$。
- 局部变量 `v2ig` (Complex): $v_k^2 e^{2i\phi}$。

##### 函数数据流

- 调用情况: `expectpj <- densitpj <- coulompj`

---

### `coulompj`

##### 函数定义

- 函数位置: `775-863 lines`

```fortran
Subroutine coulompj()
```

##### 物理公式

- **物理含义**: 利用高斯展开法计算投影密度的库仑直接势。与标准 HFB 库仑计算采用相同的数值积分技术，但作用于投影后的质子密度 `ropj(:,:,2)` 上。

- **库仑势的积分表示 (高斯展开)**:

$
\frac{1}{|\mathbf{r} - \mathbf{r}'|}= \frac{1}{\sqrt{\pi}} \int_0^D du\;\Big[ e^{-(r-r')^2 u^2} + e^{-(r+r')^2 u^2} \Big]\cdot e^{-(z-z')^2 u^2}
$

其中对 $u$ 的积分通过 Gauss-Legendre 求积实现。利用贝塞尔函数 $I_0$ 的恒等式，径向部分可写为

$
V_C(i,k)= \frac{e^2}{\sqrt{\pi}} \int_0^D du\;
I_0(2u^2 r r') \, e^{-2u^2 r r'}
\Big[ e^{-[(r-r')^2 + (z-z')^2]u^2}+ \Pi \, e^{-[(r-r')^2 + (z+z')^2]u^2} \Big]
$

其中 $I_0$ 由 `besei0`（指数缩放型）计算；$\Pi$ 为宇称因子：`Parity=.true.` 时 $\Pi=1$（对称反射贡献），`Parity=.false.` 时 $\Pi=0$（仅 $z>0$ 部分，因子 2 补偿另一半积分）。

- **投影库仑场**:

$
\text{coupj}(\mathbf{r}, k)
= \sum_i V_C(i,k) \, \rho_{\text{proj}}^{(p)}(\mathbf{r}_i, \phi_k)
$

代码中即 `coupj(:,k) = coupj(:,k) + vc(:,i) * ropj(i,k,2)`。

- **实现步骤**:
  1. 若首次调用 (`icacoupj == 0`)，预计算积分核 `vc(i,k)`：
     - 对每个高斯网格点对 $(i,k)$（即 $(r,z)$ 和 $(r',z')$），计算径向距离项 `rr2`, `zd1`, `zd2`。
     - 用 Gauss-Legendre 积分对 $u$ 求积，调用 `besei0` 计算指数缩放后的 $I_0$。
     - 将结果乘以权重 `wdcor`，对称地存入 `vc(i,k)` 和 `vc(k,i)`。
  2. 对每一投影角点 $k$，将 `vc` 与质子投影密度 `ropj(:,:,2)` 做点积，得到 `coupj`。

##### 变量对应

- 输入数组 `ropj(nghl, ilpjmax, 2)` (Complex): 投影质子密度 $\rho_{\text{proj}}^{(p)}(\mathbf{r}_i, \phi_k)$。
- 输入/预计算数组 `vc(nghl, nghl)` (Real): 库仑积分核 $V_C(i,k)$。
- 输出数组 `coupj(nghl, ilpjmax)` (Complex): 投影库仑势场。
- 局部变量 `besei0(s1)` (Real): 指数缩放的修正贝塞尔函数 $I_0(s_1) e^{-s_1}$，其中 $s_1 = 2u^2 r r'$。
- 局部变量 `fac1`, `fac2` (Real): 宇称相关的反射对称因子。

##### 函数数据流

- 调用情况: `densitpj <- coulompj <- besei0`
