# `thoout.dat` — 长输出文件说明

`thoout.dat` 是 HFBTHO 的**完整 ASCII 输出文件**（文件句柄 `lfile`），包含 `hfbtho.out` 的全部内容，并额外记录了详细的准粒子能谱、正则单粒子能级以及每一次阻塞尝试的完整迭代日志。该文件是诊断计算细节与进行谱学分析的主要依据。

> **文件关系**：`hfbtho.out`（`lout`）为精简版，`thoout.dat`（`lfile`）为完整版。代码中大部分输出通过 `Do iw=lout,lfile` 同时写入两者；`thoout.dat` 独有的部分仅在 `iw=lfile` 时写入。

---

## 1. 公共内容（与 `hfbtho.out` 相同）

以下各块在 `thoout.dat` 中完全保留，其物理含义、计算来源与打印函数与 `hfbtho.out` 一致。详见 [`hfbout.md`](./hfbout.md)。

| 块编号 | 内容摘要 | 关键打印函数 |
|--------|----------|--------------|
| 1 | 程序头与运行环境（版本、Git hash、核信息、并行配置） | `heading()` |
| 2 | 泛函参数与核物质性质（Skyrme/EDF 耦合常数、饱和点） | `print_functional_parameters()` |
| 3 | 运行特征（基参数、迭代设置、泛函、库仑、配对、温度等） | `preparer()` |
| 4 | 谐振子基与积分网格（维度、量子数、积分点、正交性检验） | `preparer()`, `base()`, `gausspoints()` |
| 5 | 初始 Woods-Saxon 势参数 | `preparer()` |
| 6 | 迭代摘要（收敛步、能量、半径、形变、配对、化学势、CPU 时间） | `iter()` |
| 7 | 未投影最终结果（粒子数、半径、多极矩、角动量、能量分解、修正项） | `expect()` |
| 8 | 阻塞候选列表（费米面附近准粒子态） | `resu()` |

---

## 2. `thoout.dat` 独有内容

### 2.1 详细迭代日志（Detailed Iteration Log）

对于**奇核子阻塞计算**，HFBTHO 会在偶偶核芯基础上逐个尝试不同的阻塞候选态。`hfbtho.out` 通常仅保留最终成功阻塞的迭代摘要，而 `thoout.dat` 则保存**每一次阻塞尝试的完整迭代序列**。

#### 输出格式

每一轮阻塞尝试以如下标志开始：
```
### REGULAR STAGE (reflection symmetry imposed)
|HFB+HO> iterations(b0=  1.750, Nsh= 10, inin= -1, N= 27, Z= 20)...
```

若某候选态导致计算失败或收敛到非物理解，下一轮尝试会在 `thoout.dat` 中继续输出新的迭代表，而 `hfbtho.out` 可能已被覆盖或仅保留最后一轮。

#### 多阻塞尝试的文件结构（`blocking_mode = -1`）

当自动阻塞模式开启时，`thoout.dat` 会顺序拼接**每一次阻塞尝试的完整计算实例**。以 Ca47（6 个阻塞候选）为例，文件结构如下：

```
### REGULAR STAGE (reflection symmetry imposed)
|HFB+HO> iterations(...)...
  [完整迭代表：0N → kB → XN → 收敛脚注]
  [中子准粒子能谱表]
  [中子正则单粒子能级表]
  [质子准粒子能谱表]
  [质子正则单粒子能级表]
  NB! From expect (UNPROJECTED RESULTS)
  [未投影物理量完整输出]
  Blocking candidates are:
  [候选列表]
### REGULAR STAGE (reflection symmetry imposed)   ← 下一轮尝试
  ...（重复 6 次）
```

> **注意**：每个阻塞尝试后都会重复打印 `expect()` 的完整输出，因此 `thoout.dat` 中 `UNPROJECTED RESULTS` 会出现多次；而 `hfbtho.out` 通常只保留最后一次（即最终成功阻塞的那一轮）。

#### 逐行变量说明

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `i` | `iiter` | `Integer` | 迭代步序号 | `iter()` | `iter()` |
| `si` | `si` | `Real` | 场矢量最大绝对变化 $\max|\text{brout}-\text{brin}|$（收敛判据） | `broyden_min()` | `iter()` |
| `mix` | `xmix` | `Real` | 当前 Broyden/线性混合参数 | `broyden_min()` | `iter()` |
| `beta` | `def(3)` | `Real` | 总四极形变 $\beta_2$ | `expect()` | `iter()` |
| `Etot` | `etot` | `Real` | 当前 HFB 总能量（MeV） | `expect()` | `iter()` |
| `A` | `varmas` | `Real` | 当前粒子数（阻塞时 $A_{\text{core}}\pm 1$） | `densit()` | `iter()` |
| `rn`, `rp` | `rms(1:2)` | `Real` | 中子/质子 RMS 半径（fm） | `expect()` | `iter()` |
| `En`, `Dn` | `ept(1)+frept(1)`, `delLN(1)` | `Real` | 中子配对能与平均能隙（MeV）；`delLN=del+frdel`（LN 时再加 `ala2`） | `expect()` | `iter()` |
| `Ep`, `Dp` | `ept(2)+frept(2)`, `delLN(2)` | `Real` | 质子配对能与平均能隙（MeV）；`delLN=del+frdel`（LN 时再加 `ala2`） | `expect()` | `iter()` |
| `Ln`, `Lp` | `alast(1:2)` | `Real` | 中子/质子化学势 / Fermi 能（MeV） | `hfbdiag()` | `iter()` |
| `time` | `time` | `Real` | 当前迭代步 wall time（秒） | 系统计时 | `iter()` |
| `time(Gog.)` | — | `Real` | Gogny 矩阵元计算耗时（秒） | `gogny_matrix_elements()` | `iter()` |

迭代步标记规则（后缀字符由 `bbroyden` 控制，与核阻塞无关）：
- **`B`**：Broyden 混合生效（`broyden_min()` 中曲率条件满足）
- **`L`**：线性混合（初始迭代、Broyden 回退或 `M <= 0`）
- **`N`**：无混合（`DoMixing=.False.`，如最后一步存储迭代）

> **注意**：正常流程下第 0 步为 `0N`（`iter()` 开始前 `gamdel()` 将 `bbroyden` 置为 `'N'`）；收敛前最后一行（如 `17N`）是关闭混合的**存储迭代**，场不再更新。

阻塞信息行（紧随迭代步之后）：
```
neutron Blocking: block=  2 state=  2   Eqp=   3.39030207  Dqpe=   1.83342175  Ovlp=   0.99999976   6/  6  1-[ 3, 3, 0]
```

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `block` | `ib` | `Integer` | 阻塞准粒子所在的 K-block | `hfbdiag()` | `iter()` |
| `state` | `k0` | `Integer` | 在该 block 内的态序号（重叠最大后） | `hfbdiag()` | `iter()` |
| `Eqp` | `eqpe` | `Real` | 被阻塞准粒子的能量（MeV） | `hfbdiag()` | `iter()` |
| `Dqpe` | `eqpe - eqpmin(it)` | `Real` | 准粒子能量与最低准粒子能量之差（相对激发能） | `hfbdiag()` | `iter()` |
| `Ovlp` | `s3` | `Real` | 与 Woods-Saxon 初猜的最大重叠（阻塞跟踪重叠） | `start()` / `hfbdiag()` | `iter()` |
| `k/total` | `k, total` | `Integer` | 当前尝试的候选序号 / 总候选数 | `handle_blocking()` | `iter()` |
| `labels` | `tb(imax)` | `Character` | 量子数标签 `{2ω}{π}[N,nz,nl]` | `base()` | `iter()` |

#### 收敛/中断脚注

迭代表结束后打印以下三行之一：

```
*   iteration converged   after  17 steps   si=     0.0000091025 ho=  0 tho=  0
  Total CPU time=   0.033 minutes
```

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `converged` / `interrupted` / `iterations limit interrupt` | — | `Character` | 收敛状态：正常收敛 / 异常中断 / 达最大迭代次数未收敛 | `iter()` | `iter()` |
| `after X steps` | `iiter` | `Integer` | 实际迭代步数 | `iter()` | `iter()` |
| `si=` | `si` | `Real` | 最终收敛指标 | `broyden_min()` | `iter()` |
| `ho=` | `iError_in_HO` | `Integer` | HO 基错误标志（`0`=正常） | `base()` | `iter()` |
| `tho=` | `iError_in_THO` | `Integer` | THO 基错误标志（`0`=正常） | `base()` | `iter()` |
| `Total CPU time=` | `time5` | `Real` | 累计 wall time（分钟） | 系统计时 | `iter()` |

---

### 2.2 准粒子能谱（Quasiparticle Energies）

由 `resu()` 在 `iw=lfile` 时打印，仅出现在 `thoout.dat` 中。列出满足 `pn > -1\mathrm{d-14}` 的活跃准粒子态（几乎包含所有非零占据态，不受 `pwi` 截断限制）。

#### 输出头说明
```
 #quasiparticle energies neutrons
 --------------------------------
  eqp(k) -> q.p. energy 
  e(k)   -> referent s.p. energy 
  p(k)   -> occ.probability 
  del(k) -> referent s.p. gap 
  fermi energy alast= -10.181604177933211
```

#### 数据列说明

| 列名 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|------|----------|------|----------|--------------|--------------|
| `#k` | `k` | `Integer` | 全局准粒子编号 | `hfbdiag()` | `resu()` |
| `block#` | `ib` | `Integer` | 所属 K-block | `hfbdiag()` | `resu()` |
| `eqp(k)` | `eqpe` | `Real` | 准粒子能量 $E_k$（MeV） | `hfbdiag()` | `resu()` |
| `e(k)` | `enb` | `Real` | 参考单粒子能 $\epsilon_k = E_k(1-2v_k^2) + \lambda$（MeV） | `hfbdiag()` | `resu()` |
| `(1-2N)E` | `(one-two*pn)*eqpe` | `Real` | 准粒子能量加权 $(1-2N)E_k$ | `hfbdiag()` | `resu()` |
| `decay` | `skk` | `Real` | 准粒子衰减长度（fm$^{-1}$） | `hfbdiag()` | `resu()` |
| `p(k)` | `pn` | `Real` | 占据概率 $v_k^2$ | `hfbdiag()` | `resu()` |
| `del(k)` | `delb` | `Real` | 参考单粒子能隙 $\Delta_k = \sqrt{E_k^2 - \epsilon_k^2}$（MeV） | `hfbdiag()` | `resu()` |
| `overl` | `ovmax` | `Real` | 与主 HO 分波的最大重叠 | `hfbdiag()` | `resu()` |
| `labels` | `tb(imax)` | `Character` | 主分量量子数标签 | `base()` | `resu()` |

> **注意**：对于质子，会再打印一套完全相同的表格，标题改为 `protons`。

---

### 2.3 正则单粒子能级（Canonical Single-Particle Energies）

由 `resu()` 在 `iw=lfile` 时打印，仅出现在 `thoout.dat` 中。通过正则变换将 HFB 准粒子态转换为单粒子表象，便于与实验谱对比。

#### 输出头说明
```
 #canonical s.p. energies neutrons
 ---------------------------------
  labels -> {2*omega}{parity}[nn=nz+2*nr+nl,nz,nl]
  cqpe   -> canonical q.p. energies
  ce     -> canonical s.p. energies
  fermi energy= -10.181604177933211
  average cdelt=  1.310749530953485
```

#### 数据列说明

| 列名 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|------|----------|------|----------|--------------|--------------|
| `k0` | `k0` | `Integer` | 正则能级全局编号 | `canonical` 模块 | `resu()` |
| `ceqp` | `ekk + ala(it)` | `Real` | 正则准粒子能量 $E_{\text{c},k} = \sqrt{(\epsilon_k-\lambda)^2 + \Delta_k^2} + \lambda$（MeV） | `canonical` 模块 | `resu()` |
| `ce` | `ek0` | `Real` | 正则单粒子能量 $\epsilon_k$（MeV） | `canonical` 模块 | `resu()` |
| `v*v` | `vk0` | `Real` | 正则占据概率 $v_k^2$ | `canonical` 模块 | `resu()` |
| `u*v` | `uuvv` | `Real` | 正则配对因子 $u_k v_k$ | `canonical` 模块 | `resu()` |
| `cdel` | `dk0` | `Real` | 正则单粒子能隙 $\Delta_k$（MeV） | `canonical` 模块 | `resu()` |
| `overl` | `ovmax` | `Real` | 与 HO 基最大重叠幅度 | `canonical` 模块 | `resu()` |
| `labels` | `tb(im+imax)` | `Character` | 主分量量子数标签 | `base()` | `resu()` |

> **零占据态说明**：正则表中大量 `v*v = 0.00000000`、`u*v = 0.0000` 的能级是 Fermi 面以上的**未占据壳层**，并非数据缺失或占位符。这些态在正则变换中仍保留，但占据概率为零。

#### 尾部统计

| 输出标签 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| `#all active are` | `j, nt` | `Integer` | 活跃准粒子数 / 总基函数数 | `hfbdiag()` | `resu()` |
| `#since the cut off is pwi=` | `pwi` | `Real` | 配对窗口截断（MeV） | — | `resu()` |
| `#check: number of` | `sum` | `Real` | 积分校验粒子数 | `densit()` | `resu()` |
| `Sum canonical e_v*V^2_k=` | `two*enjacek` | `Real` | 正则占据加权和 | `canonical` 模块 | `resu()` |

> **注意**：对于质子，会再打印一套完全相同的表格，标题改为 `protons`。

---

### 2.4 正则波函数在坐标空间输出（可选）

若输入文件中设置 `canonical_wavefunctions = T`，`resu()` 会调用 `canonical_on_grid()` 将正则波函数输出到坐标空间网格。该数据量极大，仅写入 `thoout.dat`（或额外二进制文件，视编译选项而定）。

| 输出内容 | 代码变量 | 类型 | 物理含义 | 关键计算函数 | 关键打印函数 |
|----------|----------|------|----------|--------------|--------------|
| 正则波函数格点值 | `canonical_on_grid` 内部数组 | `Real` | 各正则态在 $(r_\perp, z)$ 网格上的波函数值 | `canonical_on_grid()` | `canonical_on_grid()` |

---

## 3. 文件大小差异说明

| 输出项 | `hfbtho.out` | `thoout.dat` | 来源 |
|--------|--------------|--------------|------|
| 程序头 | ✓ | ✓ | `heading()` |
| 泛函参数 | ✓ | ✓ | `print_functional_parameters()` |
| 运行特征 | ✓ | ✓ | `preparer()` |
| HO 基与网格 | ✓ | ✓ | `preparer()`, `base()`, `gausspoints()` |
| 迭代摘要（最终） | ✓ | ✓ | `iter()` |
| 迭代摘要（每次阻塞尝试） | 仅保留最后 | ✓ 全部 | `iter()` |
| 未投影物理量 | ✓（仅最后一次） | ✓（每次阻塞尝试后重复） | `expect()` |
| 阻塞候选列表 | ✓ | ✓ | `resu()` |
| **准粒子能谱** | **✗** | **✓** | **`resu()` (iw=lfile)** |
| **正则单粒子能级** | **✗** | **✓** | **`resu()` (iw=lfile)** |
| 正则波函数格点（可选） | ✗ | ✓ | `canonical_on_grid()` |

---

## 4. 文件生成控制

| 参数 | 类型 | 说明 |
|------|------|------|
| `do_print` | `Integer` | `1` = 生成 ASCII 输出；`0` = 不生成 |
| `lfile` | `Integer` | 文件句柄，对应 `thoout.dat` |
| `Print_Screen` | `Logical` | 控制是否打印准粒子谱与正则能级（`resu()` 中使用） |
| `canonical_wavefunctions` | `Logical` | 是否额外输出正则波函数格点 |
| `toggle_output` | `Logical` | 总输出开关（由 `set_ASCII_filenames()` 控制） |

> 若 `USE_MPI=2` 团队并行，非 rank-0 进程设置 `do_print=0`，不生成个人 `thoout.dat`。
