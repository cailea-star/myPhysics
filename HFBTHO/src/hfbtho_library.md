# 代码文件 ```hfbtho_library.f90``` 说明

### `Main_Program`

##### 函数定义

- 函数位置: `85-264 lines`

```f90
Subroutine Main_Program(filename_hfbtho, filename_unedf, &
                        my_comm_world, my_comm_team, my_n_teams, my_team_color, &
                        toggle_output, filename_output, filename_dat, filename_binary)
```

##### 物理公式

- 物理含义: HFBTHO 计算流程的顶层驱动子程序。根据编译宏配置，可执行单次 HFB 计算、质量表（Mass Table）、滴线（Drip Lines）扫描或势能面（PES）计算。

##### 实现步骤

1. **MPI 环境配置**（`130-132`）：在 `USE_MPI>0` 时调用 `configure_MPI_environment` 设置进程通信域
2. **初始化 Namelist 默认值**（`139`）：调用 `initialize_HFBTHO_NAMELIST` 预设所有输入参数的默认值
3. **输入读取与广播**（`146-222`）：
   - `USE_MPI=1`：团队 leader（`team_rank==0`）读取输入并广播至团队内所有进程
   - `USE_MPI=0` 或 `2`：进程 `0` 读取输入并广播至 `COMM_world`
   - 并行模式下调用 `Construct_Vectors` 构建可广播向量，广播后非 leader 进程调用 `Deconstruct_Vectors`
4. **基组自适应调整**（`226`）：若四极约束激活（`lambda_active(2)>0`）且 `automatic_basis=.True.`，调用 `adjust_basis` 根据期望四极矩自动优化谐振子基组
5. **分派计算模式**（`230-245`）：
   - `DO_MASSTABLE=1`：调用 `compute_mass_table`
   - `DRIP_LINES=1`：调用 `compute_driplines`
   - `DO_PES=1`：调用 `compute_PES`
   - 否则：设置 ASCII 文件名并调用 `execute_HFBTHO` 执行单次计算
6. **错误收集与收尾**（`249-262`）：在 `USE_MPI=2` 下通过 `mpi_gather` 汇总各进程错误标志，关闭输出文件

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `filename_hfbtho` | `Character(Len=256)` | 输入 | HFBTHO Namelist 输入文件名 |
| `filename_unedf` | `Character(Len=256)` | 输入 | UNEDF 泛函参数文件名 |
| `my_comm_world` | `Integer` | 输入 | MPI 全局通信域句柄 |
| `my_comm_team` | `Integer` | 输入 | MPI 团队子通信域句柄 |
| `my_n_teams` | `Integer` | 输入 | 团队总数 |
| `my_team_color` | `Integer` | 输入 | 当前团队颜色标识 |
| `toggle_output` | `Logical` | 输入 | 输出开关（`.True.` 启用文件输出）|
| `filename_output` | `Character(Len=256)` | 输入 | ASCII 短输出文件名 |
| `filename_dat` | `Character(Len=256)` | 输入 | ASCII 长输出文件名 |
| `filename_binary` | `Character(Len=256)` | 输入 | 二进制重启文件名 |

##### 函数数据流

- 调用情况: `hfbthoprog <- Main_Program <- configure_MPI_environment / initialize_HFBTHO_NAMELIST / read_HFBTHO_NAMELIST / Construct_Vectors / broadcast_vectors / Deconstruct_Vectors / adjust_basis / set_ASCII_filenames / compute_mass_table / compute_driplines / compute_PES / execute_HFBTHO`

---

### `execute_HFBTHO`

##### 函数定义

- 函数位置: `268-489 lines`

```f90
Subroutine execute_HFBTHO(filename_unedf, filename_binary)
```

##### 物理公式

- 物理含义: 将用户输入的模块级物理参数映射到求解器所需的内部 `_INI` 变量，随后调用 HFBTHO-DFT 自洽求解器并可选地执行对称性恢复与局域化分析。

- 输入参数到内部求解器变量的映射关系：

$$
\begin{aligned}
n_{00}^{\text{INI}} &= N_{\text{shells}} \\
b_0^{\text{INI}} &= b_0 \\
q^{\text{INI}} &= \beta_{\text{basis}} \\
N^{\text{INI}} &= N, \quad Z^{\text{INI}} = Z \\
\beta_2^{0} &= \beta_2^{\text{WS}}, \quad \beta_3^{0} = \beta_3^{\text{WS}}, \quad \beta_4^{0} = \beta_4^{\text{WS}} \\
N_{\text{max}}^{\text{INI}} &= N_{\text{iter}}, \quad \varepsilon^{\text{INI}} = \varepsilon_{\text{conv}} \\
V_{0n}^{\text{INI}} &= V_{\text{pair}}^{(n)}, \quad V_{0p}^{\text{INI}} = V_{\text{pair}}^{(p)} \\
T &\equiv \text{temper} \quad (\text{温度模式})
\end{aligned}
$$

其中 $N_{\text{shells}}$ 为壳层数，$b_0$ 为谐振子长度，$\beta_{\text{basis}}$ 为基组形变，$\beta_{2,3,4}^{\text{WS}}$ 为 Woods-Saxon 初始势的形变参数，$V_{\text{pair}}^{(n/p)}$ 为中子/质子对相互作用强度，$T$ 为温度。

- 温度截断：若 $T < 10^{-10}$，强制关闭温度模式

$$
\text{If } |T| \le 10^{-10} \quad \Rightarrow \quad \text{switch\_on\_temperature} = .\text{False.}
$$

- 四极矩与形变的转换关系（在 `compute_driplines` 中使用）：

$$
Q_{20} = \beta \sqrt{\frac{5}{\pi}} \frac{A^{5/3}}{100}
$$

##### 实现步骤

1. **参数映射**（`301-358`）：将模块全局变量逐一赋值给求解器内部 `_INI` 变量（壳层数、谐振子长度、基组形变、核子数、迭代参数、泛函名、库仑类型、配对参数、阻塞组态、温度、网格点数等）
2. **一致性检查**（`363`）：调用 `check_consistency` 校验 `_INI` 参数之间的兼容性；若出错则输出错误信息并返回
3. **约束条件处理**（`378-391`）：
   - 统计激活的约束个数 $\lambda_{\text{active}}$
   - 若存在奇多极矩约束（$l=1,3,5,\dots$），强制关闭宇称（`Parity_INI=.False.`）
   - 若无任何约束，关闭集体惯性计算
   - 若计算裂变碎片特性，也强制关闭宇称
4. **泛函读取与特殊处理**（`396-405`）：
   - 调用 `read_UNEDF_NAMELIST` 读取 Skyrme/EDF 泛函参数
   - 对 `SeaLL1` 或 `SLY4mod` 泛函自动开启配对正则化并将截断设为 $100.0$
   - 若无法识别泛函，自动关闭投影（`iproj_INI=0`）
5. **阻塞态设置**（`416-432`）：根据 `nkblocase` 决定是阻塞所有近费米面准粒子（能量窗 $1\,\text{MeV}$）还是阻塞单个指定态
6. **调用 DFT 求解器**（`436`）：调用 `HFBTHO_DFT_SOLVER` 执行自洽迭代
7. **对称性恢复**（`440`）：若 PNP 或 AMP 开启且非裂变碎片计算，调用 `HFBTHO_restore`
8. **局域化函数**（`444`）：若 `localization_functions=.True.`，调用 `localization`
9. **错误报告**（`449-487`）：检查 `ierror_flag`，向输出文件报告求解器结束状态或错误信息

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `filename_unedf` | `Character(Len=256)` | 输入 | UNEDF 泛函参数文件名 |
| `filename_binary` | `Character(Len=256)` | 输入 | 二进制重启文件名 |
| `number_of_shells` | `Integer` | 输入 | 谐振子壳层数 $N_{\text{shells}}$，映射到 `n00_INI` |
| `oscillator_length` | `Real(pr)` | 输入 | 谐振子长度 $b_0$，映射到 `b0_INI` |
| `basis_deformation` | `Real(pr)` | 输入 | 基组四极形变 $\beta_{\text{basis}}$，映射到 `q_INI` |
| `neutron_number` | `Integer` | 输入 | 中子数 $N$，映射到 `npr_INI(1)` |
| `proton_number` | `Integer` | 输入 | 质子数 $Z$，映射到 `npr_INI(2)` |
| `beta2_deformation` | `Real(pr)` | 输入 | WS 初始势 $\beta_2$，映射到 `b2_0` |
| `beta3_deformation` | `Real(pr)` | 输入 | WS 初始势 $\beta_3$，映射到 `b3_0` |
| `beta4_deformation` | `Real(pr)` | 输入 | WS 初始势 $\beta_4$，映射到 `b4_0` |
| `number_iterations` | `Integer` | 输入 | 最大迭代次数 $N_{\text{iter}}$，映射到 `MAX_ITER_INI` |
| `accuracy` | `Real(pr)` | 输入 | 收敛精度 $\varepsilon_{\text{conv}}$，映射到 `epsi_INI` |
| `functional` | `Character` | 输入 | EDF 泛函名称，映射到 `skyrme_INI` |
| `vpair_n` | `Real(pr)` | 输入 | 中子配对强度 $V_{\text{pair}}^{(n)}$，映射到 `V0n_INI` |
| `vpair_p` | `Real(pr)` | 输入 | 质子配对强度 $V_{\text{pair}}^{(p)}$，映射到 `V0p_INI` |
| `temperature` | `Real(pr)` | 输入 | 温度 $T$，映射到 `temper` |
| `lambda_active` | `Integer` 数组 | 输入 | 多极约束激活标志，决定集体惯性计算与宇称设置 |

##### 函数数据流

- 调用情况: `Main_Program / compute_mass_table / compute_driplines / compute_PES <- execute_HFBTHO <- check_consistency / read_UNEDF_NAMELIST / replace_functional / HFBTHO_DFT_SOLVER / HFBTHO_restore / localization`

---

### `compute_mass_table`

##### 函数定义

- 函数位置: `499-577 lines`

```f90
Subroutine compute_mass_table(filename_unedf, toggle_output, filename_output, filename_dat, filename_binary)
```

##### 物理公式

- 物理含义: 在预定义的核素网格上批量执行 HFB 计算，生成质量表。质量表中的每一行对应一组 $(Z, N)$ 及可选的基组形变 $\beta$ 或四极矩 $Q_{20}$。

- 质量数与四极矩的关系：

$$
A = Z + N
$$

$$Q_{20} = \beta \sqrt{\frac{5}{\pi}} \frac{A^{5/3}}{100}
$$

其中 $A$ 为质量数，$\beta$ 为基组形变参数，$Q_{20}$ 为四极矩期望值（单位：barn）。

##### 实现步骤

1. **输出向量分配**（`518-524`）：若 `nRows>0`，按编译条件分配 `out_vectors`（`USE_MPI=2`）或 `mass_table` 数组
2. ** bookkeeping 文件创建**（`528`）：进程 `0` 创建 `TableLog.dat` 记录所有计算核素
3. **循环遍历质量表行**（`530-563`）：
   - 读取当前行的 $Z$、$N$、$\beta$、$Q_{20}$
   - 跳过第 `0` 行（表头占位）
   - 通过 `Mod(iRow, mpi_size) /= mpi_taskid` 实现 MPI 负载均衡，每个进程仅计算属于自己的行
   - 设置 `proton_number=Z_chain`、`neutron_number=N_chain`、`expectation_values(2)=Q20`、`basis_deformation=beta_deformation`
   - 调用 `set_ASCII_filenames` 与 `execute_HFBTHO` 执行单次计算
   - 计算完成后调用 `fill_out_vectors` 或 `fill_mass_table` 收集结果
4. **结果汇总与输出**（`564-569`）：在 `USE_MPI=2` 下调用 `gather_results` 汇总，最后由 `print_mass_table` 输出完整质量表
5. **关闭 bookkeeping 文件**（`570`）：团队 leader 关闭 `TableLog.dat`

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `filename_unedf` | `Character(Len=256)` | 输入 | UNEDF 泛函参数文件名 |
| `toggle_output` | `Logical` | 输入 | 文件输出开关 |
| `filename_output` | `Character(Len=256)` | 输入 | ASCII 短输出文件名 |
| `filename_dat` | `Character(Len=256)` | 输入 | ASCII 长输出文件名 |
| `filename_binary` | `Character(Len=256)` | 输入 | 二进制重启文件名 |
| `Z_chain` | `Integer` | 内部 | 当前行质子数 $Z$ |
| `N_chain` | `Integer` | 内部 | 当前行中子数 $N$ |
| `A_chain` | `Integer` | 内部 | 当前行质量数 $A=Z+N$ |
| `beta_deformation` | `Real(pr)` | 内部 | 当前行基组形变 $\beta$ |
| `Q20` | `Real(pr)` | 内部 | 当前行四极矩约束值 $Q_{20}$ |
| `iRow` | `Integer` | 内部 | 质量表行索引 |
| `icalc` | `Integer` | 内部 | 本地成功计算计数器 |

##### 函数数据流

- 调用情况: `Main_Program <- compute_mass_table <- set_ASCII_filenames / execute_HFBTHO / allocate_out_vectors / allocate_mass_table / fill_out_vectors / fill_mass_table / gather_results / print_mass_table`

---

### `compute_driplines`

##### 函数定义

- 函数位置: `588-701 lines`

```f90
Subroutine compute_driplines(filename_unedf, toggle_output, filename_output, filename_dat, filename_binary)
```

##### 物理公式

- 物理含义: 沿给定稳定线（stable line）的同位素链或同中子素链进行扫描，通过双核子分离能 $S_{2n}$ 或 $S_{2p}$ 判定滴线位置。

- 形变步长：

$$
\Delta\beta = \frac{1}{N_{\text{defs}} - 1}
$$

其中 $N_{\text{defs}}$ 为每个核素的形变采样点数。

- 四极矩与形变转换：

$$
Q_{20}(\beta) = \beta \sqrt{\frac{5}{\pi}} \frac{A^{5/3}}{100}, \qquad A = Z + N
$$

- 双核子分离能（滴线判据）：

$$
S_{2N} = E(N-2, Z) - E(N, Z)
$$

$$
S_{2Z} = E(N, Z-2) - E(N, Z)
$$

其中 $E(N, Z)$ 为当前核素的 HFB 能量。当 $S_{2N} < 0$（或 $S_{2Z} < 0$）时，认为滴线已到达，停止该链的扫描。

##### 实现步骤

1. **形变步长计算**（`605`）：$\Delta\beta = 1/(N_{\text{defs}}-1)$
2. **团队 bookkeeping 文件**（`608-611`）：团队 leader 创建 `TeamTableXXX.dat` 并分配 `Energy_chain_gthr` 数组
3. **外层循环：稳定线上的核素**（`614-692`）：
   - 通过 `Mod(iRow, number_teams) /= team_color` 实现团队级负载均衡
   - 读取当前稳定线的 $(Z, N)$ 及扫描方向（同位素/同中子素）
4. **内层循环：沿链扫描**（`621-681`）：
   - 初始化 `Minimum_Energy_Prev = 100`
   - 对每个形变点 $\beta_i = -0.5 - \Delta\beta + i\cdot\Delta\beta$（$i=1,\dots,N_{\text{defs}}$）计算 $Q_{20}$
   - 通过 `Mod(i_deformation, team_size) /= team_rank` 实现团队内进程负载均衡
   - 设置核子数与约束值，调用 `execute_HFBTHO`
   - 根据是否使用 Lipkin-Nogami（`kindhfb_INI>0` 取 `ehfb`，否则取 `etot`）记录能量到 `Energy_chain`
5. **滴线判定**（`664-679`）：
   - 调用 `find_minimum_energy` 获取当前核素最低能量 `Minimum_Energy`
   - 计算 $S_{2N} = E_{\text{prev}}^{\min} - E_{\text{curr}}^{\min}$
   - 若 $S_{2N} < 0$，到达滴线，退出内层循环
   - 否则沿链方向增加 2 个中子或 2 个质子继续扫描（上限 $N \le 310$、$Z \le 120$）
6. **收尾**：团队 leader 输出链扫描完成信息，关闭 bookkeeping 文件

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `filename_unedf` | `Character(Len=256)` | 输入 | UNEDF 泛函参数文件名 |
| `toggle_output` | `Logical` | 输入 | 文件输出开关 |
| `filename_output` | `Character(Len=256)` | 输入 | ASCII 短输出文件名 |
| `filename_dat` | `Character(Len=256)` | 输入 | ASCII 长输出文件名 |
| `filename_binary` | `Character(Len=256)` | 输入 | 二进制重启文件名 |
| `beta_step` | `Real(pr)` | 内部 | 形变采样步长 $\Delta\beta$ |
| `Z_chain` | `Integer` | 内部 | 当前链质子数 $Z$ |
| `N_chain` | `Integer` | 内部 | 当前链中子数 $N$ |
| `A_chain` | `Integer` | 内部 | 质量数 $A=Z+N$ |
| `beta_deformation` | `Real(pr)` | 内部 | 当前形变参数 $\beta$ |
| `Q20` | `Real(pr)` | 内部 | 由 $\beta$ 转换得到的四极矩 $Q_{20}$ |
| `Energy_chain` | `Real(pr)` | 内部 | 当前核素的 HFB 能量 $E$ |
| `Minimum_Energy_Prev` | `Real(pr)` | 内部 | 前一个核素的最低能量 $E_{\text{prev}}^{\min}$ |
| `Minimum_Energy` | `Real(pr)` | 内部 | 当前核素的最低能量 $E_{\text{curr}}^{\min}$ |
| `separation_2N` | `Real(pr)` | 内部 | 双核子分离能 $S_{2N}$ 或 $S_{2Z}$ |
| `direction_sl` | `Integer` 数组 | 输入 | 扫描方向标志（1=同位素链，其他=同中子素链）|

##### 函数数据流

- 调用情况: `Main_Program <- compute_driplines <- set_ASCII_filenames / execute_HFBTHO / find_minimum_energy`

---

### `compute_PES`

##### 函数定义

- 函数位置: `711-780 lines`

```f90
Subroutine compute_PES(filename_unedf, toggle_output, filename_output, filename_dat, filename_binary)
```

##### 物理公式

- 物理含义: 在多维形变空间（由 $(Z, N)$ 与多极约束期望值 $Q_{\lambda}$ 定义）的离散网格上批量执行 HFB 计算，构建原子核的势能面（Potential Energy Surface, PES）。

- 质量数：

$$
A = Z + N
$$

- 多极约束期望值直接由输入网格给出：

$$
Q_{\lambda} = Q_{\text{PES}}^{(j)}(i_{\text{row}})
$$

其中 $\lambda$ 为多极阶次（通常 $\lambda=2,3,4$），$j$ 为约束序号，$i_{\text{row}}$ 为 PES 网格点索引。

##### 实现步骤

1. **bookkeeping 文件**（`733`）：进程 `0` 创建 `TableLog.dat` 记录所有 PES 点
2. **循环遍历 PES 网格点**（`735-772`）：
   - 读取当前点的 $(Z, N)$ 及多极约束值
   - 通过 `Mod(iRow, number_teams) /= team_color` 实现团队级负载均衡
   - 设置 `proton_number=Z_chain`、`neutron_number=N_chain`
   - 读取当前点的 Woods-Saxon 初始形变 `bet2_PES`、`bet3_PES`、`bet4_PES`（若值 $>-8.0$ 则覆盖默认值）
   - 遍历每个约束 $j$：
     - 若 $\lambda=2$ 且 `automatic_basis=.True.`，调用 `adjust_basis(Q_PES(iRow,j), .False., Z, N)` 自适应调整基组
     - 将约束值写入 `expectation_values(lambda)` 并激活 `lambda_active(lambda)=1`
   - 调用 `set_ASCII_filenames` 与 `execute_HFBTHO` 执行单次计算
   - 输出进度信息并更新计数器 `icalc`
3. **关闭文件**（`773`）：若 `do_print==1`，关闭 `TableLog.dat`

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `filename_unedf` | `Character(Len=256)` | 输入 | UNEDF 泛函参数文件名 |
| `toggle_output` | `Logical` | 输入 | 文件输出开关 |
| `filename_output` | `Character(Len=256)` | 输入 | ASCII 短输出文件名 |
| `filename_dat` | `Character(Len=256)` | 输入 | ASCII 长输出文件名 |
| `filename_binary` | `Character(Len=256)` | 输入 | 二进制重启文件名 |
| `Z_chain` | `Integer` | 内部 | 当前点质子数 $Z$ |
| `N_chain` | `Integer` | 内部 | 当前点中子数 $N$ |
| `A_chain` | `Integer` | 内部 | 质量数 $A=Z+N$ |
| `Q_PES` | `Real(pr)` 数组 | 输入 | PES 网格多极约束期望值 $Q_{\lambda}$ |
| `lambda_PES` | `Integer` 数组 | 输入 | 约束对应的多极阶次 $\lambda$ |
| `bet2_PES` | `Real(pr)` 数组 | 输入 | WS 初始势 $\beta_2$ 覆盖值 |
| `bet3_PES` | `Real(pr)` 数组 | 输入 | WS 初始势 $\beta_3$ 覆盖值 |
| `bet4_PES` | `Real(pr)` 数组 | 输入 | WS 初始势 $\beta_4$ 覆盖值 |
| `iRow` | `Integer` | 内部 | PES 网格点索引 |
| `ndefs` | `Integer` | 输入 | 每个点的约束个数 |
| `npoints` | `Integer` | 输入 | PES 网格总点数 |

##### 函数数据流

- 调用情况: `Main_Program <- compute_PES <- set_ASCII_filenames / adjust_basis / execute_HFBTHO`
