# 代码文件 `hfbtho_localization.f90` 说明

## 模块概述

`HFBTHO_localization` 模块提供计算核子局域化函数（Nucleon Localization Function, NLF）的接口，在 Gauss 积分网格上分别计算中子和质子的局域化函数并输出到文件。

---

### `localization`

##### 函数定义

- 函数位置: `97-129 lines`

```fortran
Subroutine localization()
```

##### 物理公式

- 物理含义: 计算核子局域化函数（NLF）。NLF 衡量给定空间点处的核子是否处于局域化状态，取值在 0 到 1 之间，越接近 1 表示局域化程度越高。计算分为两步：首先计算 Thomas-Fermi 动能密度，然后代入局域化函数公式。

Thomas-Fermi 动能密度：

$$
\tau_{\text{TF},q}(r,z) = \frac{3}{5}(6\pi^2)^{2/3} \rho_q^{5/3}(r,z)
$$

局域化函数：

$$
C_q(r,z) = \frac{1}{1 + \left( \dfrac{\rho_q \tau_q - \frac{1}{4}\left[ (\nabla_r \rho_q)^2 + (\nabla_z \rho_q)^2 \right] }{\rho_q \tau_{\text{TF},q} \left(\frac{1}{2}\right)^{2/3}} \right)^2 }
$$

其中 $q \in \{n, p\}$ 分别对应中子（代码中索引 1）和质子（代码中索引 2）。

- 实现步骤:
  1. 对每个网格点 `irz` 计算 Thomas-Fermi 动能密度 `tautf(irz,q)`
  2. 计算局域化函数 `loca(irz,q)`
  3. 将结果写入文件 `localization.dat`

##### 变量对应

- 输入变量1: `ro(nghl,2) (real(pr))`，公式对应 $\rho_q(r,z)$，粒子密度
- 输入变量2: `tau(nghl,2) (real(pr))`，公式对应 $\tau_q(r,z)$，动能密度
- 输入变量3: `NABLAR(nghl,2) (real(pr))`，公式对应 $\nabla_r \rho_q$，密度的径向梯度
- 输入变量4: `NABLAZ(nghl,2) (real(pr))`，公式对应 $\nabla_z \rho_q$，密度的轴向梯度
- 输入变量5: `fl(nghl) (real(pr))`，公式对应 $r$，径向坐标网格
- 输入变量6: `fh(nghl) (real(pr))`，公式对应 $z$，轴向坐标网格
- 输出: 写入文件 `localization.dat`，列分别为 $r$, $z$, $\rho_n/2$, $\rho_p/2$, $C_n$, $C_p$

##### 函数数据流

- 调用情况: `hfbtho_library.f90` <- `localization`
