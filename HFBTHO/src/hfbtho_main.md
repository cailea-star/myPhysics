# 代码文件 ```hfbtho_main.f90``` 说明

### `Program hfbthoprog`

##### 程序定义

- 程序位置: `69-192 lines`

```f90
Program hfbthoprog
```

##### 物理含义

HFBTHO 求解器的**主程序入口**。负责初始化 MPI 并行环境、解析命令行输入参数、设置默认文件名，并调用 `Main_Program` 执行实际的 HFB 计算流程（单次计算或大规模并行质量表/势能面/滴线扫描）。

---

##### 执行流程

- 实现步骤:
  1. **MPI 环境初始化**（`106-120`）：根据编译宏 `USE_MPI` 的值，初始化 MPI 通信域 `COMM_world`、获取进程总数 `mpi_size` 与当前进程号 `mpi_taskid`
  2. **文件名预设**（`122-127`）：设置默认文件名——HFBTHO 输入文件 `hfbtho_NAMELIST.dat`、泛函参数文件 `UNEDF_NAMELIST.dat`、短输出 `hfbtho.out`、长输出 `thoout.dat`、二进制重启文件 `hfbtho_output.hel`
  3. **命令行参数解析**（`129-155`）：读取命令行参数个数 `n_args`
     - `n_args = 0`：使用全部默认文件名
     - `n_args = 1`：第 1 个参数覆盖 `filename_hfbtho`
     - `n_args = 2`：第 1 个参数覆盖 `filename_hfbtho`，第 2 个参数覆盖 `filename_unedf`
     - `n_args > 2`：输出用法说明并终止程序
  4. **MPI 团队创建**（`157-174`）：在并行模式下调用 `Create_MPI_Teams` 创建子通信域 `COMM_team`，用于大规模计算的任务分发
  5. **调用主计算程序**（`176-183`）：将文件名与 MPI 通信参数传递给 `Main_Program`，由后者完成输入读取、自洽迭代与结果输出
  6. **MPI 同步与关闭**（`185-190`）：调用 `mpi_barrier` 等待全部进程结束，随后 `mpi_finalize` 释放 MPI 环境

---

##### 变量对应

| 变量名 | 类型 | 输入/输出 | 物理含义 |
|--------|------|-----------|----------|
| `n_args` | `Integer` | 隐式输入 | 命令行参数个数 |
| `filename_hfbtho` | `Character(Len=256)` | 输入/输出 | HFBTHO Namelist 输入文件名 |
| `filename_unedf` | `Character(Len=256)` | 输入/输出 | UNEDF 泛函参数文件名 |
| `filename_output` | `Character(Len=256)` | 输出 | ASCII 短输出文件名 `hfbtho.out` |
| `filename_dat` | `Character(Len=256)` | 输出 | ASCII 长输出文件名 `thoout.dat` |
| `filename_binary` | `Character(Len=256)` | 输出 | 二进制重启文件名 `hfbtho_output.hel` |
| `toggle_output` | `Logical` | 内部 | 输出开关标志（固定为 `.True.`）|
| `my_comm_world` | `Integer` | 输出 | MPI 全局通信域句柄 |
| `my_comm_team` | `Integer` | 输出 | MPI 团队子通信域句柄 |
| `my_n_teams` | `Integer` | 输出 | 团队总数 |
| `my_team_color` | `Integer` | 输出 | 当前团队颜色标识 |

---

##### 函数数据流

- 调用情况: `hfbthoprog -> Main_Program <- hfbtho_library.f90`
- MPI 子程序调用: `hfbthoprog -> MPI_INIT / MPI_COMM_SIZE / MPI_COMM_RANK / Create_MPI_Teams / mpi_barrier / mpi_finalize`
- 系统调用: `hfbthoprog -> command_argument_count / get_command_argument`
