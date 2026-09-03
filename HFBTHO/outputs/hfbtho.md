# `hfbtho.out` — 短输出文件说明

`hfbtho.out` 是 HFBTHO 的**精简 ASCII 输出文件**（文件句柄 `lout`），记录程序运行的关键摘要信息。它与 `thoout.dat` 共享大部分内容，但不包含详细的准粒子谱与正则单粒子能级。

> **注意**：当 `do_print = 1` 时，以下各块通过 `Do iw=lout,lfile` 同时写入 `hfbtho.out` 与 `thoout.dat`；若 `do_print = 0` 则两个文件均不生成。

---

## 1. 程序头与运行环境

由 `heading()` 子程序打印，位于文件最前端。

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| ASCII art 标题 | — | `Character` | 程序大标题（AXIALLY DEFORMED CONFIGURATIONAL...） | — | `heading()` |
| `v1.66 / v2.00d / v3.00` | — | `Character` | 版本历史与对应作者 | — | `heading()` |
| `UNEDF Module Version` | — | `Character` | UNEDF 模块版本号 | — | `heading()` |
| 日期戳 | — | `Character` | 计算时间戳 | — | `heading()` |
| `git: parent: ...` | `Version` | `Character` | 当前代码的 Git commit hash | — | `heading()` |
| `FORTRAN 95 CODE (KIND= 8)` | `pr` | `Integer` | 浮点精度 KIND 参数 | — | `heading()` |
| `Nucleus: Ca (A= 47, N= 27, Z= 20)` | `nucname, npr` | `Character/Integer` | 核素符号、质量数、中子数、质子数 | `nucleus()` | `heading()` |
| `Reflection Symmetry Imposed` | `Parity` | `Logical` | 是否强制反射对称性（$z\to -z$） | `preparer()` | `heading()` |
| `Distributed parallelism with MPI` | `HFB_cores` | `Integer` | MPI 并行任务数 | MPI 初始化 | `heading()` |
| `Multi-threading framework with OpenMP` | `numThreads` | `Integer` | OpenMP 每任务线程数 | OpenMP 查询 | `heading()` |

---

## 2. 泛函参数与核物质性质

由 `print_functional_parameters()`（`UNEDF` 模块）打印，输出当前能量密度泛函的耦合常数与饱和点性质。

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `Crho(0)`, `Crho(1)` | `Crho(0:1)` | `Real` | 密度-密度耦合常数（$t=0,1$） | `set_functional_parameters()` | `print_functional_parameters()` |
| `CDrho(0)`, `CDrho(1)` | `CDrho(0:1)` | `Real` | 密度梯度耦合常数 | `set_functional_parameters()` | `print_functional_parameters()` |
| `Ctau(0)`, `Ctau(1)` | `Ctau(0:1)` | `Real` | 密度-动能密度耦合常数 | `set_functional_parameters()` | `print_functional_parameters()` |
| `CpV0(0)`, `CpV0(1)` | `CpV0(0:1)` | `Real` | 中子/质子配对强度 $V_0$ | `initialize_HFBTHO_SOLVER()` | `print_functional_parameters()` |
| `CpV1(0)`, `CpV1(1)` | `CpV1(0:1)` | `Real` | 配对表面-体积混合参数 | `initialize_HFBTHO_SOLVER()` | `print_functional_parameters()` |
| `sigma` | `sigma` | `Real` | 密度相关项幂指数 | `set_functional_parameters()` | `print_functional_parameters()` |
| `hbzero` | `hbzero` | `Real` | $\hbar^2/2m$（MeV·fm$^2$） | `set_functional_parameters()` | `print_functional_parameters()` |
| `e^2 chrg` | `chargee2` | `Real` | 质子电荷平方 $e^2$（MeV·fm） | `preparer()` | `print_functional_parameters()` |
| `E_NM` | `E_NM` | `Real` | 核物质每核子能量（MeV） | `set_functional_parameters()` | `print_functional_parameters()` |
| `K_NM` | `K_NM` | `Real` | 核物质不可压缩模量（MeV） | `set_functional_parameters()` | `print_functional_parameters()` |
| `RHO_NM` | `RHO_NM` | `Real` | 饱和点密度（fm$^{-3}$） | `set_functional_parameters()` | `print_functional_parameters()` |
| `ASS_NM` | `ASS_NM` | `Real` | 对称能（MeV） | `set_functional_parameters()` | `print_functional_parameters()` |
| `LASS_NM` | `LASS_NM` | `Real` | 对称能斜率 $L$（MeV） | `set_functional_parameters()` | `print_functional_parameters()` |
| `SMASS_NM` | `SMASS_NM` | `Real` | 有效质量 $m^*/m$ | `set_functional_parameters()` | `print_functional_parameters()` |
| `CrDr(0)`, `CrDr(1)` | `CrDr(0:1)` | `Real` | 密度-梯度耦合常数（自然单位） | `set_functional_parameters()` | `print_functional_parameters()` |
| `CrdJ(0)`, `CrdJ(1)` | `CrdJ(0:1)` | `Real` | 自旋-轨道耦合常数 | `set_functional_parameters()` | `print_functional_parameters()` |
| `CJ(0)`, `CJ(1)` | `CJ(0:1)` | `Real` | 自旋-自旋耦合常数 | `set_functional_parameters()` | `print_functional_parameters()` |
| `CExPar` | `CExPar` | `Real` | 交换参数 | `set_functional_parameters()` | `print_functional_parameters()` |
| `functional has DME couplings` | — | `Logical` | DME 耦合标志 | `set_functional_parameters()` | `print_functional_parameters()` |
| `use 3-Nucleon DME couplings` | — | `Logical` | 三核子 DME 标志 | `set_functional_parameters()` | `print_functional_parameters()` |
| `c.m. correction` | — | `Logical` | 质心修正开关 | `set_functional_parameters()` | `print_functional_parameters()` |
| `use tensor terms` | — | `Logical` | 张量力开关 | `set_functional_parameters()` | `print_functional_parameters()` |
| `P_NM` | `P_NM` | `Real` | 饱和点压强（MeV·fm$^{-3}$） | `set_functional_parameters()` | `print_functional_parameters()` |
| `VMASS_NM` | `VMASS_NM` | `Real` | 矢量有效质量 $m_V^*/m$ | `set_functional_parameters()` | `print_functional_parameters()` |
| `t0, x0, t1, x1, ...` | `t0, x0, ...` | `Real` | 传统 Skyrme $(t,x)$ 参数 | `set_functional_parameters()` | `print_functional_parameters()` |
| `b4, b4p, te, to` | `b4, b4p, te, to` | `Real` | Skyrme 附加参数 | `set_functional_parameters()` | `print_functional_parameters()` |
| `Crho_nu` ~ `Lambda_nu` | `*_nu` 系列 | `Real` | 自然单位耦合常数（非 NEDF/非有限程力时打印） | `set_functional_parameters()` | `print_functional_parameters()` |

---

## 3. 运行特征（Characteristics of the run）

由 `preparer()` 打印，汇总输入参数与数值设置。

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `Nucleus` | `nucname, npr` | `Char/Int` | 核名称与核子数 | `nucleus()` | `preparer()` |
| `Number of HO shells` | `n00` | `Integer` | 谐振子主壳层截断 | — | `preparer()` |
| `HO length b0` | `b0` | `Real` | 谐振子长度 $b_0$（fm） | `preparer()` | `preparer()` |
| `Basis deformation` | `beta0, q` | `Real` | 基函数形变参数 | `preparer()` | `preparer()` |
| `THO basis` | `iLST` | `Integer` | THO 基开关状态 | `preparer()` | `preparer()` |
| `Maximal number of iterations` | `maxi` | `Integer` | 最大迭代次数 | — | `preparer()` |
| `Initial mixing parameter` | `xmix` | `Real` | Broyden/线性混合初值 | `initialize_HFBTHO_SOLVER()` | `preparer()` |
| `Initial wave functions` | `inin` | `Integer` | 初始波函数来源 | `start()`, `inout()` | `preparer()` |
| `Energy functional` | `skyrme` | `Character` | 当前泛函名称 | `set_functional_parameters()` | `preparer()` |
| `Pairing cutoff energy` | `pwi` | `Real` | 配对截断（MeV） | `initialize_HFBTHO_SOLVER()` | `preparer()` |
| `Lipkin-Nogami procedure` | `kindhfb` | `Integer` | LN 修正开关 | — | `preparer()` |
| `PAV procedure` | `iproj` | `Integer` | 投影后变分开关 | `preparer()` | `preparer()` |
| `Constraint calculation` | `lambda_active` | `Integer` | 约束计算开关 | `Constraint_or_not()` | `preparer()` |
| `Neutron/Proton blocking` | `blocking_mode` | `Integer` | 阻塞模式 | `handle_blocking()` | `preparer()` |
| `Temperature T` | `temper` | `Real` | 温度（MeV） | — | `preparer()` |
| `Pairing regularization` | `pairing_regularization` | `Logical` | 配对正则化开关 | — | `preparer()` |
| `Collective inertia` | `collective_inertia` | `Logical` | 集体惯性开关 | — | `preparer()` |
| `Fission fragments` | `fission_fragments` | `Logical` | 裂变碎片开关 | — | `preparer()` |
| `Restart indicator` | `inin` | `Integer` | 重启文件指示符 | `inout()` | `preparer()` |
| `Broyden mixing` | `nbroyden` | `Integer` | Broyden 历史长度 | — | `preparer()` |
| `Output file` | `lfile` | `Integer` | 长输出文件单元号 | — | `preparer()` |
| `HO: b0,1/b0,bp,bz,q` | `b0, bz, bp, q` | `Real` | 详细谐振子长度与形变参数 | `preparer()` | `preparer()` |
| `h**2/(2m_n)` | `hbzero_n` | `Real` | 中子质量参数（MeV·fm$^2$） | `preparer()` | `preparer()` |
| `h**2/(2m_p)` | `hbzero_p` | `Real` | 质子质量参数（MeV·fm$^2$） | `preparer()` | `preparer()` |
| `hom=f*41.0*A^{-1/3}` | `hom` | `Real` | HO 频率估计值及系数 | `preparer()` | `preparer()` |
| 库仑方法说明 | — | `Character` | 直接库仑（substitution）+ 交换库仑（Slater） | `preparer()` | `preparer()` |

---

## 4. 谐振子基与积分网格

由 `preparer()` / `base()` / `gausspoints()` 打印。

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `NUV, NQP` | `nuv, nqp` | `Integer` | HFB 矩阵维度与准粒子数 | `thoalloc()` | `preparer()` |
| `n00` | `n00` | `Integer` | 最大壳层数 | `base0()` | `preparer()` |
| `nbx, 2*n00+1` | `nbx` | `Integer` | 最大 K-block 数 | `base0()` | `preparer()` |
| `ntx` | `ntx` | `Integer` | 最大单粒子能级数 | `base0()` | `preparer()` |
| `ndx` | `ndx` | `Integer` | 单个 block 最大维度 | `base0()` | `preparer()` |
| `nb` | `nb` | `Integer` | 实际 block 数 | `base()` | `preparer()` |
| `nt` | `nt` | `Integer` | 实际能级数 | `base()` | `preparer()` |
| `Number of Gauss-Hermite mesh points ngh` | `ngh` | `Integer` | $z$ 方向积分点数 | `gausspoints()` | `preparer()` |
| `Number of Gauss-Laguerre mesh points ngl` | `ngl` | `Integer` | $\rho$ 方向积分点数 | `gausspoints()` | `preparer()` |
| `Number of Gauss-Legendre mesh points nleg` | `nleg` | `Integer` | 库仑/PNP 积分点数 | `gausspoints()` | `preparer()` |
| `Hermite / Laguerre / Legendre` boundaries | `xh, xl, xleg` | `Real` | 积分网格边界 | `gausspoints()` | `preparer()` |
| `Orthogonality / Normalization` | — | `Real` | 基函数正交归一性偏差 | `gausspoints()` | `preparer()` |
| `nzx` | `nzx` | `Integer` | 最大 $n_z$ 量子数 | `base0()` | `preparer()` |
| `nrx` | `nrx` | `Integer` | 最大 $n_r$ 量子数 | `base0()` | `preparer()` |
| `nlx` | `nlx` | `Integer` | 最大 $m_l$ 量子数 | `base0()` | `preparer()` |

### 4.1 实际使用的谐振子基（Actual basis used）

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `Maximal 2*omega` | `nom` | `Integer` | 实际最大 $2\omega$ | `base()` | `preparer()` |
| `Maximal nz` | `nzm` | `Integer` | 实际最大 $n_z$ | `base()` | `preparer()` |
| `Maximal nr` | `nrm` | `Integer` | 实际最大 $n_r$ | `base()` | `preparer()` |
| `Maximal ml` | `nlm` | `Integer` | 实际最大 $m_l$ | `base()` | `preparer()` |
| `Maximal N=nz+2*nr+nl` | — | `Integer` | 实际最大主量子数 | `base()` | `preparer()` |
| `2 x biggest block dim.` | — | `Integer` | 最大 block 维度 $\times$ 2 | `base()` | `preparer()` |
| `Non-zero elements of h` | — | `Integer` | HFB 矩阵非零元个数 | `base()` | `preparer()` |
| `Number of Broyden elements` | — | `Integer` | Broyden 混合向量维度 | `base()` | `preparer()` |

---

## 5. 初始 Saxon-Woods 势参数

由 `preparer()` / `start()` 打印，位于积分网格之后、迭代之前。用于生成初始单粒子波函数。

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `v0ws` | `v0ws` | `Real` | Woods-Saxon 势阱深度（MeV） | `start()` | `preparer()` |
| `kappa` | `akappa` | `Real` | 自旋-轨道耦合强度 | `start()` | `preparer()` |
| `vs0` | `vso` | `Real` | 自旋-轨道势强度（MeV） | `start()` | `preparer()` |
| `r0` | `r0` | `Real` | 势阱半径参数（fm） | `start()` | `preparer()` |
| `a` | `a` | `Real` | 势阱弥散参数（fm） | `start()` | `preparer()` |
| `r0-so` | `rso` | `Real` | 自旋-轨道半径参数（fm） | `start()` | `preparer()` |
| `a-so` | `aso` | `Real` | 自旋-轨道弥散参数（fm） | `start()` | `preparer()` |
| `b2_ws` | `b2` | `Real` | WS 势四极形变 | `start()` | `preparer()` |
| `b3_ws` | `b3` | `Real` | WS 势八极形变 | `start()` | `preparer()` |
| `b4_ws` | `b4` | `Real` | WS 势十六极形变 | `start()` | `preparer()` |

---

## 6. 迭代摘要（Iteration Summary）

由 `iter()` 打印，显示自洽迭代的收敛过程。对于阻塞计算，仅保留**最终成功阻塞**的迭代摘要；对于偶偶核，即为唯一迭代序列。

### 6.1 迭代前导言

迭代表开始前打印以下标志行：

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `Initialization for the even-even core` | `npr` | `Integer` | 偶偶核心初始化声明 | `start()` | `preparer()` |
| `Reading from wel_file: ...` | `file_wel` | `Character` | 从 wel 重启文件读取波函数 | `inout()` | `preparer()` |
| `### REGULAR STAGE (...)` | `parity` | `Logical` | 反射对称性阶段标记 | `preparer()` | `preparer()` |
| `\|HFB+HO> iterations(...)` | `b0, n00, inin, npr` | `Mixed` | 迭代阶段参数摘要 | `preparer()` | `preparer()` |

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `i` | `iiter` | `Integer` | 迭代步序号 | `iter()` | `iter()` |
| `si` | `si` | `Real` | 场矢量最大绝对变化（收敛判据） | `broyden_min()` | `iter()` |
| `mix` | `xmix` | `Real` | 当前混合参数 | `broyden_min()` | `iter()` |
| `beta` | `def(3)` | `Real` | 总四极形变 $\beta_2$ | `expect()` | `iter()` |
| `Etot` | `etot` | `Real` | 当前总能量（MeV） | `expect()` | `iter()` |
| `A` | `varmas` | `Real` | 当前粒子数（阻塞时可能为 $A\pm 1$） | `densit()` | `iter()` |
| `rn`, `rp` | `rms(1:2)` | `Real` | 中子/质子 RMS 半径（fm） | `expect()` | `iter()` |
| `En`, `Dn` | `ept(1)+frept(1)`, `delLN(1)` | `Real` | 中子配对能与平均能隙（MeV）；`delLN=del+frdel`（LN 时再加 `ala2`） | `expect()` | `iter()` |
| `Ep`, `Dp` | `ept(2)+frept(2)`, `delLN(2)` | `Real` | 质子配对能与平均能隙（MeV）；`delLN=del+frdel`（LN 时再加 `ala2`） | `expect()` | `iter()` |
| `Ln`, `Lp` | `alast(1:2)` | `Real` | 中子/质子化学势 / Fermi 能（MeV） | `hfbdiag()` | `iter()` |
| `time` | `time` | `Real` | 当前迭代步 wall time（秒） | 系统计时 | `iter()` |
| `time(Gog.)` | `wct_gogny` | `Real` | Gogny 矩阵元计算耗时（秒） | `gogny_matrix_elements()` | `iter()` |
| 迭代后缀 `B`/`L`/`N` | `bbroyden` | `Char` | 混合模式：`B`=Broyden, `L`=线性, `N`=无混合 | `broyden_min()` | `iter()` |
| `* iteration converged` | `si, iiter` | — | 收敛标志、最终步数、`ho`/`tho` 错误标志 | `iter()` | `iter()` |
| `Total CPU time` | `time5` | `Real` | 累计 wall time（分钟） | 系统计时 | `iter()` |

> 阻塞尝试的额外行：`neutron Blocking: block= X state= Y Eqp= ... Dqpe= ... Ovlp= ...` 由 `iter()` 在每次阻塞迭代时打印，标记当前阻塞的准粒子态。

---

## 7. 未投影最终结果（UNPROJECTED RESULTS）

由 `expect(.True., 1)` 打印，为 HFB 自洽收敛后的核心物理量。该块同时出现在 `hfbtho.out` 与 `thoout.dat` 末尾。

### 7.1 基本物理量

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `forces:` | `skyrme` | `Character` | 当前泛函名称（再次打印） | `set_functional_parameters()` | `expect()` |
| `NB! From expect (UNPROJECTED RESULTS)` | — | `Character` | 未投影结果块标题 | `expect()` | `expect()` |
| `pairing: CpV0,CpV1` | `CpV0, CpV1` | `Real` | 配对相互作用参数 | `initialize_HFBTHO_SOLVER()` | `expect()` |
| `Blocked neutron block` | `bloblo(1)` | `Integer` | 中子阻塞 block 编号（仅中子阻塞时出现） | `handle_blocking()` | `expect()` |
| `Blocked proton block` | `bloblo(2)` | `Integer` | 质子阻塞 block 编号（仅质子阻塞时出现） | `handle_blocking()` | `expect()` |
| `Requested part.numbs.` | `tz` | `Real` | 目标中子/质子/总核子数 | — | `expect()` |
| `UnPj(av) part.numbs` | `xn` | `Real` | HFB 平均粒子数 | `densit()` | `expect()` |
| `b0, bz, bp` | `b0, bz, bp` | `Real` | 谐振子长度参数（fm） | `preparer()` | `expect()` |
| `lambda (ala)` | `ala` | `Real` | HFB 化学势（MeV） | `hfbdiag()` | `expect()` |
| `Lambda (alast)` | `alast` | `Real` | 正则基 Fermi 能（MeV） | `hfbdiag()` | `expect()` |
| `delta(n,p), pwi` | `del, pwi` | `Real` | 平均能隙与截断 | `expect()` | `expect()` |
| `pairing energy` | `ept` / `frept` | `Real` | 中子/质子/总配对能（MeV） | `expect()` | `expect()` |
| `LN lambda_2` | `ala2` | `Real` | Lipkin-Nogami 修正参数（仅 `kindhfb<0` 时出现） | `sslipkin()` | `expect()` |
| `LN energies` | `etr` | `Real` | LN 能量修正（仅 `kindhfb<0` 时出现） | `sslipkin()` | `expect()` |
| `Geff(n,p)` | `Geff` | `Real` | LN 有效配对强度（仅 `kindhfb<0` 时出现） | `sslipkin()` | `expect()` |

### 7.2 半径与多极矩

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `rms-radius` | `rms` | `Real` | 中子/质子/总 RMS 半径（fm） | `expect()` | `expect()` |
| `charge-radius, r0` | `rc, r00` | `Real` | 电荷半径与经验半径（fm） | `expect()` | `expect()` |
| `deformation beta2` | `def` | `Real` | 中子/质子/总四极形变 | `expect()` | `expect()` |
| `dipole moment[fm]` | `qmoment(1,it)` | `Real` | 电偶极矩（fm） | `moments_computeValue()` | `expect()` |
| `quadrupole moment[b]` | `qmoment(2,it)` | `Real` | 电四极矩（barn） | `moments_computeValue()` | `expect()` |
| `octupole moment` | `qmoment(3,it)` | `Real` | 电八极矩 | `moments_computeValue()` | `expect()` |
| `hexadecapole moment` | `qmoment(4,it)` | `Real` | 电十六极矩 | `moments_computeValue()` | `expect()` |
| `q5` ~ `q8` | `qmoment(5:8,it)` | `Real` | 更高阶多极矩 | `moments_computeValue()` | `expect()` |

### 7.3 角动量

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `<Jx^2> from rho` | `jxsq_rhoN, jxsq_rhoP` | `Real` | 来自密度的 $J_x^2$ 期望值 | `expect()` | `expect()` |
| `<Jx^2> from kappa` | `jxsq_kappaN, jxsq_kappaP` | `Real` | 来自配对张量的 $J_x^2$ | `expect()` | `expect()` |
| `<J^2> = 2*<Jx^2>` | — | `Real` | 总角动量平方 $2*(<Jx^2>_{\rho} + <Jx^2>_{\kappa})$ | `expect()` | `expect()` |

### 7.4 能量分解（MeV）

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `kinetic energy` | `ekt` | `Real` | 中子/质子/总动能 | `expect()` | `expect()` |
| `volume energy` | `evol` | `Real` | 体积能（$\rho\tau + \rho\rho$ 项之和） | `expect()` | `expect()` |
| `rho_tau` | `EVOL_rho_tau` | `Real` | 密度-动能密度耦合能 | `expect()` | `expect()` |
| `rho_rho` | `EVOL_rho_rho` | `Real` | 密度-密度耦合能 | `expect()` | `expect()` |
| `surface energy` | `esurf` | `Real` | 表面能（各项之和） | `expect()` | `expect()` |
| `rho_DELTA_rho` | `ESURF_rho_DELTA_rho` | `Real` | Laplacian 密度耦合能 | `expect()` | `expect()` |
| `(NABLA_rho)^2` | `ESURF_NABLA_rho_NABLA_rho` | `Real` | 密度梯度平方耦合能 | `expect()` | `expect()` |
| `spin-orbit energy` | `eso` | `Real` | 总自旋-轨道耦合能 | `expect()` | `expect()` |
| `rho_NABLA_J` | `ESO_rho_NABLA_J` | `Real` | 密度-自旋轨道密度耦合能 | `expect()` | `expect()` |
| `NABLA_rho_J` | `ESO_NABLA_rho_J` | `Real` | 密度梯度-自旋密度耦合能 | `expect()` | `expect()` |
| `coulomb energy` | `ecodi+ecoex+coept(3)` | `Real` | 库仑总能量 | `coulom()` / `expect()` | `expect()` |
| `direct` | `ecodi` | `Real` | 直接库仑能 | `coulom()` | `expect()` |
| `exchange` | `ecoex` | `Real` | 交换库仑能（Slater） | `expect()` | `expect()` |
| `pairing` | `coept(3)` | `Real` | 库仑配对能 | `expect()` | `expect()` |
| `tensor energy` | `etens` | `Real` | 张量力能量 | `expect()` | `expect()` |
| `direct Hartree E` | `E_HARTREE_DIR` | `Real` | Hartree 直接项（有限程力） | `expect()` | `expect()` |
| `Extra E` | `EExtra` | `Real` | 额外场能量 | `expect()` | `expect()` |
| `External field E` | `E_EXT_FIELD` | `Real` | 外场能量 | `expect()` | `expect()` |
| `Entropy` | `entropy` | `Real` | 中子/质子/总熵 | `densit()` / `hfbdiag()` | `expect()` |
| `tEnergy: ehfb (qp)` | `ehfb` | `Real` | HFB 准粒子总能量（MeV） | `expect()` | `expect()` |

### 7.5 修正项（未加入总能量）

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `cmc-diagonal part` | `ECMHFB` | `Real` | 质心修正对角部分 | `expect()` | `expect()` |
| `cmc-hfb` | `ECMHFB` | `Real` | 质心修正 HFB 项 | `expect()` | `expect()` |
| `cranking rot corr` | `DEROT` | `Real` | 推转转动修正 | `expect()` | `expect()` |
| `SQUJ` | `SQUJ` | `Real` | 角动量涨落相关量 | `expect()` | `expect()` |
| `CRAN x 4` | `CRAN` | `Real` | Cranking 惯性相关量 | `expect()` | `expect()` |
| `Rigit Body` | `ERIGHFB` | `Real` | 刚体转动能量 | `expect()` | `expect()` |

---

## 8. 阻塞候选信息（Blocking Candidates）

由 `resu()` / `handle_blocking()` 在计算结束后打印，列出所有在费米面附近被考虑的阻塞态。

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `eqpmin` | `eqpmin` | `Real` | 最低准粒子能量（MeV） | `hfbdiag()` | `resu()` / `iter()` |
| `pwiblo` | `pwiblo` | `Real` | 阻塞能量窗口（MeV） | `preparer()` | `resu()` |
| `num` | `k` | `Integer` | 阻塞候选序号 | `handle_blocking()` | `resu()` |
| `block` | `ib` | `Integer` | 所属 K-block 编号 | `hfbdiag()` | `resu()` |
| `state` | `k` | `Integer` | block 内态编号 | `hfbdiag()` | `resu()` |
| `Eqp` | `eqpe` | `Real` | 准粒子能量（MeV） | `hfbdiag()` | `resu()` |
| `(1-2N)E` | `(one-two*pn)*eqpe` | `Real` | 准粒子能量加权 | `hfbdiag()` | `resu()` |
| `overl` | `ovmax` | `Real` | 与 HO 基最大重叠 | `hfbdiag()` | `resu()` |
| `labels` | `tb(imax)` | `Character` | 量子数标签 `{2ω}{π}[N,nz,nl]` | `base()` | `resu()` |
| `keyblo/blomax` | `keyblo, blomax` | `Integer` | 当前阻塞候选序号 / 总候选数 | `handle_blocking()` | `iter()` / `resu()` |
| 全局态编号 | `imax` | `Integer` | 该态在整个单粒子能级中的全局编号 | `hfbdiag()` | `resu()` |

---

## 9. 文件结束标志

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `HFBTHO_SOLVER ended without errors` | — | `Character` | 正常结束标志 | `HFBTHO_SOLVER` | `heading()` |

---

## 10. 文件生成控制

| 参数 | 类型 | 说明 |
|------|------|------|
| `do_print` | `Integer` | `1` = 生成 ASCII 输出；`0` = 不生成 |
| `lout` | `Integer` | 文件句柄，对应 `hfbtho.out` |
| `toggle_output` | `Logical` | 是否真正打开输出文件（由 `set_ASCII_filenames()` 控制） |

> 若 `USE_MPI=2` 团队并行，非 rank-0 进程设置 `do_print=0`，不生成个人输出文件。
