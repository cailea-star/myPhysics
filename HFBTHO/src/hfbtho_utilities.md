# 代码文件 ```hfbtho_utilities.f90``` 说明

### `ord`

##### 函数定义

- 函数位置: `116-137 lines`

```f90
Subroutine ord(n, e)
```

##### 物理公式

- 物理含义: 对一维实数向量进行选择排序（Selection Sort），按升序重新排列元素。

$$
e_{(1)} \le e_{(2)} \le \cdots \le e_{(n)}
$$

其中 $e_{(i)}$ 表示排序后第 $i$ 个位置的元素值。

- 实现步骤:
  1. 遍历位置 $i = 1, \dots, n$
  2. 在子区间 $[i, n]$ 中寻找最小值元素 $e_{(k)}$
  3. 若 $k \ne i$，交换 $e_{(i)}$ 与 $e_{(k)}$

##### 变量对应

- 输入变量1: `n (Integer(ipr))`, 公式对应 $n$（向量长度）
- 输入/输出变量1: `e(1:n) (Real(pr))`, 公式对应 $e_{(i)}$（待排序向量，InOut）

##### 函数数据流

- 调用情况: `hfbdiag <- ord`

---

### `get_CPU_time`

##### 函数定义

- 函数位置: `141-159 lines`

```f90
Subroutine get_CPU_time(subname, is)
```

##### 物理公式

- 物理含义: 利用系统时钟记录并输出指定子程序的实际运行耗时（Wall-clock time）。

$$
\Delta t = \frac{t_{2} - t_{1}}{f_{\text{countrate}}}
$$

其中 $t_{1}$ 为起始时刻的时钟计数，`is = 0` 时记录；$t_{2}$ 为结束时刻的时钟计数，`is = 1` 时记录；$f_{\text{countrate}}$ 为系统时钟频率（counts per second）。

- 实现步骤:
  1. 若 `is == 0`，调用 `system_clock` 获取起始计数 $t_{1}$ 及频率 $f_{\text{countrate}}$
  2. 若 `is == 1`，调用 `system_clock` 获取结束计数 $t_{2}$
  3. 计算耗时 $\Delta t = (t_{2} - t_{1}) / f_{\text{countrate}}$ 并输出到文件单元 `lout` 与 `lfile`

##### 变量对应

- 输入变量1: `is (Integer(ipr))`, 公式对应 $t$ 的触发标志（`0` 记录起始 $t_{1}$，`1` 记录结束 $t_{2}$）
- 输入变量2: `subname (Character*(*))`, 公式对应子程序标识符（不参与公式计算，仅用于输出标注）
- 输出（隐式）: 打印耗时 $\Delta t$（Real(pr)）到输出文件

##### 函数数据流

- 调用情况: `coulompj / hfbdiag / coulom1 / coulom / coulom_test / densit / field / gamdel / broyden / expect / getLagrange <- get_CPU_time <- system_clock`
