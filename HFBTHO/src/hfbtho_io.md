# 代码文件 `hfbtho_io.f90` 说明

本文档对模块 `HFBTHO_io`（1083 行）中的全部子程序与函数进行说明。该模块是 HFBTHO 的输入/输出中心，负责 ASCII 日志文件的命名与打开、二进制重启文件（`*.hel`/`*.tel`）的读写，以及新旧版本格式的兼容处理。

---

### `set_ASCII_filenames`

##### 函数定义

- 函数位置: `144-197 lines`

```fortran
Subroutine set_ASCII_filenames(toggle_output, filename_output, filename_dat)
```

##### 物理公式

- 物理含义: 无物理计算。该子程序根据 MPI 并行模式设置屏幕输出文件名并打开文件句柄 `lout`（基本输出）和 `lfile`（扩展输出）。

- 实现步骤:
  1. 根据 `USE_MPI` 编译标志和可选的用户输入确定文件名。
  2. `USE_MPI=0`：单核运行，固定名 `hfbtho.out` / `thoout.dat`。
  3. `USE_MPI=1`：多核但每核独立目录，仍用固定名。
  4. `USE_MPI=2`：团队并行，所有任务在同一目录，文件名附加 `row_string`（如 `hfbtho_000001.out`）。
  5. 仅 rank 0 打开文件并设置 `do_print=1`，其余核心 `do_print=0`。

##### 变量对应

- 输入变量1: `toggle_output (Logical)`，是否真正打开输出文件
- 输入变量2: `filename_output (Character(Len=256), Optional)`，用户指定的基本输出文件名
- 输入变量3: `filename_dat (Character(Len=256), Optional)`，用户指定的扩展输出文件名
- 输出变量1: 全局文件句柄 `lout`、`lfile` 被打开
- 输出变量2: 全局标志 `do_print` 被设置

##### 函数数据流
- 调用情况: `hfbtho_solver / hfbtho_main <- set_ASCII_filenames`
- 被 `HFBTHO_DFT_SOLVER` 在求解器初始化阶段调用

---

### `set_BINARY_filenames`

##### 函数定义

- 函数位置: `202-224 lines`

```fortran
Subroutine set_BINARY_filenames(filename_binary)
```

##### 物理公式

- 物理含义: 无物理计算。设置二进制重启/输出文件名 `welfile`。

- 实现步骤:
  1. 若用户提供文件名，直接采用。
  2. 否则根据 `iLST`（HO/THO 模式）和 `USE_MPI` 标志自动生成：
     - HO 模式: `hfbtho_output.hel`
     - THO 模式: `hfbtho_output.tel`
     - MPI=2 时附加 `row_string`。

##### 变量对应

- 输入变量1: `filename_binary (Character(Len=256), Optional)`，用户指定的二进制文件名
- 输出变量1: `welfile (Character(Len=256))`，全局二进制文件名

##### 函数数据流
- 调用情况: `inout <- set_BINARY_filenames`
- 被 `inout` 内部调用

---

### `inout`

##### 函数定义

- 函数位置: `231-317 lines`

```fortran
Subroutine inout(is, iexit, filename_binary)
```

##### 物理公式

- 物理含义: 无物理计算。作为 HFBTHO 所有二进制 IO 的统一入口，根据 `is` 的值决定读或写：
  - `is=1`：读取二进制文件作为重启初猜。
  - `is=2`：写入 HFB 对角化前的数据（密度、场、约束等）。
  - `is=3`：追加写入 HFB 对角化后的数据（阻塞、温度等）。

- 实现步骤:
  1. 调用 `set_BINARY_filenames` 确定文件名。
  2. 读模式：若 `inin>0` 则从头开始（`iexit=1`）；否则检查文件、读取版本号，调用 `read_data` 或 `read_data_old`。MPI 模式下通过 `broadcast_binary_to_team` 广播数据。
  3. 写模式：检查文件状态，`is=2` 时先写版本号，再调用 `write_data`；`is=3` 且 `VERSION_DATA==2` 时调用 `write_data_old`。

##### 变量对应

- 输入变量1: `is (Integer(ipr))`，IO 模式（1=读，2=写前，3=写后）
- 输入变量2: `filename_binary (Character(Len=256), Optional)`，可选的二进制文件名
- 输入/输出变量3: `iexit (Integer(ipr))`，退出状态（0=成功，>0=失败）
- 输入变量4: `inin`、`do_print`、`iasswrong(3)` 等全局标志

##### 函数数据流
- 调用情况: `iter / hfbtho_main <- inout <- set_BINARY_filenames / read_data / write_data / version_number / check_file`
- 被 `iter`（`hfbtho_solver.f90`，每次迭代重启）和主程序调用

---

### `write_version`

##### 函数定义

- 函数位置: `321-326 lines`

```fortran
Subroutine write_version()
```

##### 物理公式

- 物理含义: 将当前二进制文件格式版本号 `VERSION_DATA`（当前为 5）写入文件首部，供后续读取时识别格式。

##### 变量对应

- 输出变量: 向 `lwou` 写入整数 `VERSION_DATA`

##### 函数数据流
- 调用情况: `inout <- write_version`
- 仅由 `inout` 在 `is=2` 时调用

---

### `version_number`

##### 函数定义

- 函数位置: `331-341 lines`

```fortran
Integer(ipr) Function version_number()
```

##### 物理公式

- 物理含义: 从已打开的二进制文件中读取版本号。若读取失败（如旧格式无版本号），返回 $-1$。

##### 变量对应

- 输入变量: 已打开的文件句柄 `lwin`
- 输出变量: 函数返回值，版本号（或 $-1$）

##### 函数数据流
- 调用情况: `inout <- version_number`
- 仅由 `inout` 在读取分支调用

---

### `check_file`

##### 函数定义

- 函数位置: `348-391 lines`

```fortran
Integer(ipr) Function check_file(filename, action, is)
```

##### 物理公式

- 物理含义: 检查指定二进制文件的状态，并根据读/写模式打开相应句柄（`lwin` 或 `lwou`）。

- 实现步骤:
  1. 读模式：文件必须存在，以 `old`+`unformatted` 打开 `lwin`。
  2. 写模式：
     - `is=2`：覆盖已存在文件，或创建新文件。
     - `is=3`：追加到已存在文件；若不存在则报错。

##### 变量对应

- 输入变量1: `filename (Character(Len=256))`，文件名
- 输入变量2: `action (Character(Len=50))`，`'Read'` 或 `'Write'`
- 输入变量3: `is (Integer(ipr))`，IO 子模式
- 输出变量1: 函数返回值（0=成功，1=失败）
- 输出变量2: 打开全局句柄 `lwin` 或 `lwou`

##### 函数数据流
- 调用情况: `inout <- check_file`
- 被 `inout` 在每次读写前调用

---

### `read_data_old`

##### 函数定义

- 函数位置: `395-503 lines`

```fortran
Subroutine read_data_old(iexit)
```

##### 物理公式

- 物理含义: 读取 HFBTHO v2.00d 格式的旧二进制文件。该格式无关键字标记，数据按固定顺序紧密排列，包括：
  - 核子数、网格点数、壳层数 `npr, ngh, ngl, n00, nb, nt`
  - HO 参数与能量 `b0, bz, bp, beta0, si, etot, rms, ...`
  - 密度矩阵维度 `ntx, NB, nhhdim`
  - Lagrange 乘子 `multLag`
  - 分块维度 `ID_r`
  - HFB 矩阵元 `brin`（含 $h$ 与 $\Delta$）
  - 基函数量子数（旧编码格式）
  - 阻塞信息
  - THO 样条数据（若 `iLST>0`）

- 实现步骤:
  1. 顺序读取各记录，用 `counterLine` 跟踪当前记录号。
  2. 若 `n00` 或 `nb` 与当前计算不一致，报错退出。
  3. 对无 Lipkin-Nogami 情况，在配对场上加小量 $0.1$ 防止配对塌缩。
  4. 将旧编码量子数转换为 `nr, nz, nl, ns`。

##### 变量对应

- 输入/输出变量1: `iexit (Integer(ipr))`，0=成功，1=失败
- 输出变量: 填充全局变量 `brin`, `ID_r`, `nr_r`, `nz_r`, `nl_r`, `ns_r`, `bloall`, `bloblo`, `blo123`, `blok1k2`, `blomax`, `bloqpdif` 等

##### 函数数据流
- 调用情况: `inout <- read_data_old`
- 被 `inout` 当 `VERSION_READ < 3` 时调用

---

### `read_data`

##### 函数定义

- 函数位置: `512-809 lines`

```fortran
Subroutine read_data(iexit)
```

##### 物理公式

- 物理含义: 读取 HFBTHO v3/v4/v5 格式的关键字标记二进制文件。文件由若干带 `key` 的数据块组成，各块可独立读取，实现跨基组、跨约束数、跨泛函的重启。

支持的关键字块：
| 关键字 | 内容 |
|--------|------|
| `Metadata` | 质子数、中子数、各种功能开关标志 |
| `SkyFunct` | Skyrme 泛函参数、核物质性质、耦合常数、$\hbar^{2}/2m$ |
| `HO-Basis` | HO 长度 $b_0, b_z, b_\perp$、壳层数、积分网格节点权重 |
| `QuantNum` | 分块维度 `ID` 与各基函数量子数 `nr, nz, nl, ns` |
| `Various.` | 收敛残差、能量、形变、化学势、质量等 |
| `Constrai` | 约束数量、多极阶数、期望值、Lagrange 乘子、颈部约束 |
| `Densits.` | 坐标空间密度矩阵 `ro` 与配对张量 `aka` |
| `FieldsN.` | 中子平均场与配对场（10 个场数组） |
| `FieldsP.` | 质子平均场与配对场（10 个场数组） |
| `Blocking` | 阻塞候选信息 |
| `Blk-rest` | 阻塞状态标志 |
| `Temperat` | 温度、熵、温度占据数 `fn_T`, `fp_T` |
| `THObasis` | THO 样条系数与 LST 参数 |
| `Regular.` | 配对正则化有效质量 `MEFF` 与耦合 `geff_inv` |
| `GognyVNN` | Gogny 有限程矩阵 `rk`, `ak`, `vrGogny`, `vzGogny` |
| `HFBmatrX` | 完整 HFB 矩阵元 `brin` 与维度 `nhhdim` |

- 实现步骤:
  1. 循环读取 `key`，根据关键字进入对应分支。
  2. 若 `functional == 'READ'`，用文件中的泛函参数覆盖当前值；否则保留当前泛函，仅读取场和密度。
  3. 检查积分网格 `ngh, ngl` 是否一致；检查基组参数 `nb, b0, bz, bp, n00, nt` 是否一致（`different_basis` 标志）。
  4. 对约束块，按多极阶数匹配新旧 Lagrange 乘子。
  5. 读取错误时跳转至 `100` 标签，设置 `iexit=1` 并输出错误信息。

##### 变量对应

- 输入/输出变量1: `iexit (Integer(ipr))`，0=成功，1=失败
- 输入变量: `VERSION_READ`, `functional`, `ngh`, `ngl`, `nb`, `b0`, `bz`, `bp`, `n00`, `nt`
- 输出变量: 填充几乎全部 HFBTHO 全局运行变量（密度、场、阻塞、温度、THO 系数、HFB 矩阵元等）

##### 函数数据流
- 调用情况: `inout <- read_data`
- 被 `inout` 当 `VERSION_READ >= 3` 时调用

---

### `write_data`

##### 函数定义

- 函数位置: `815-974 lines`

```fortran
Subroutine write_data(is)
```

##### 物理公式

- 物理含义: 将当前 HFB 计算状态写入 v5 格式的关键字标记二进制文件，供后续重启。

- 实现步骤:
  1. `is=2`（对角化前）：写入元数据、泛函参数、基组信息、量子数、各种物理量、约束信息、密度/配对张量、中子/质子场、THO 系数、配对正则化量、Gogny 矩阵、完整 HFB 矩阵（可选）。
  2. `is=3`（对角化后）：追加写入阻塞候选（排序后）、阻塞状态、温度信息。
  3. 每个数据块前均有 8 字符关键字标识。

##### 变量对应

- 输入变量1: `is (Integer(ipr))`，2=写前，3=写后
- 输入变量: 几乎全部 HFBTHO 全局变量
- 输出变量: 向 `lwou` 写入二进制数据

##### 函数数据流
- 调用情况: `inout <- write_data`
- 被 `inout` 在写模式调用

---

### `write_data_old`

##### 函数定义

- 函数位置: `979-1049 lines`

```fortran
Subroutine write_data_old()
```

##### 物理公式

- 物理含义: 按 v2.00d 旧格式写入二进制文件，仅用于调试兼容。数据按固定顺序无关键字写入，包括核子数、HO 参数、能量、矩阵维度、Lagrange 乘子、分块信息、HFB 矩阵元、量子数（旧编码）、阻塞信息、THO 数据。

##### 变量对应

- 输入变量: 对应旧格式的全部全局变量
- 输出变量: 向 `lwou` 按旧格式写入二进制数据

##### 函数数据流
- 调用情况: `inout <- write_data_old`
- 被 `inout` 当 `is=3` 且 `VERSION_DATA==2` 时调用

---

### `blosort`

##### 函数定义

- 函数位置: `1053-1079 lines`

```fortran
Subroutine blosort(it, n)
```

##### 物理公式

- 物理含义: 无物理计算。对阻塞候选态按能量差 `bloqpdif` 进行简单选择排序，使写入文件的阻塞候选按能量从低到高排列。

##### 变量对应

- 输入变量1: `it (Integer(ipr))`，核子类型（1=中子，2=质子）
- 输入变量2: `n (Integer(ipr))`，该类型的阻塞候选数量
- 输入/输出变量: `bloqpdif(:,it)`、`bloblo(:,it)`、`blo123(:,it)`、`blok1k2(:,it)` 被原地排序

##### 函数数据流
- 调用情况: `write_data / write_data_old <- blosort`
- 被 `write_data` 和 `write_data_old` 在写入阻塞信息前调用

---

## 模块全局变量

| 变量名 | 类型 | 说明 |
|--------|------|------|
| `VERSION_DATA` | `Integer(ipr)` | 当前写入的二进制文件版本号（=5） |
| `VERSION_READ` | `Integer(ipr)` | 读取的文件版本号 |
| `welfile` | `Character(Len=256)` | 二进制文件名 |
| `Z_r, N_r` | `Integer(ipr)` | 文件中的质子数、中子数 |
| `n00_r, nb_r, nt_r` | `Integer(ipr)` | 文件中的壳层数、分块数、基函数总数 |
| `ngh_r, ngl_r, nleg_r` | `Integer(ipr)` | 文件中的积分网格点数 |
| `b0_r, bz_r, bp_r` | `Real(pr)` | 文件中的 HO 长度参数 |
| `ID_r, nr_r, nz_r, nl_r, ns_r` | 数组 | 文件中的分块维度与量子数 |
| `xh_r, xl_r, wh_r, wl_r` | 数组 | 文件中的积分节点与权重 |

---

## 模块数据流总览

```
HFBTHO_DFT_SOLVER / iter / hfbtho_main
  └── inout(is, iexit, filename_binary)
        ├── set_BINARY_filenames
        ├── check_file
        ├── write_version              [is=2]
        ├── read_data / read_data_old  [is=1]
        │     └── 填充全局密度、场、矩阵、阻塞、温度等
        └── write_data / write_data_old [is=2,3]
              └── blosort
```
