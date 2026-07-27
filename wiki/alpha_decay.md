---
[tags]: alpha_decay
[tag-type]: phenomenon
---

# alpha_decay



### Definition & Origin

```references
  - D_S_Delion_Y.2018_J.Phys.G_Nucl.Part.Phys._Vol.45Nol.P.053001: abstract.
  - Jianmin_Dong_Y.2025_NuclearPhysicsA_Vol.1059Nol.123089P.123089: section 3.3. α-decay.
```

α 衰变态可在 α–子核势中表述为窄的出射 Gamow 共振：本征值的实部对应 \(Q_\alpha\)，虚部决定总 α 衰变宽度的一半。在超重质量区，α 衰变被视为主要衰变模式之一，\(Q_\alpha\) 与半衰期是实验可测的两个关键物理量。



### Composition & Mechanism

##### 精细结构的通道组成

```references
  - D_S_Delion_Y.2018_J.Phys.G_Nucl.Part.Phys._Vol.45Nol.P.053001: section 1. Introduction.
```

- **description**: 对形变 α 发射体，不同能量的 α 谱线组对应到子核不同能级的跃迁；这些能差与 γ 射线能量相符，由此形成 α 衰变精细结构的线状谱。

##### 势垒穿透与预形成

```references
  - D._S._Delion_Y.2006_Phys.Rev.C_Vol.73Nol.014315P.: section I. Introduction.
```

- **description**: 简单势垒图像可以描述半衰期的相对变化；要得到绝对半衰期，还需用 α 粒子谱因子或预形成概率乘以势垒穿透率。



### Quantities & Properties

```tags
[tags]: alpha_decay_energy, alpha_decay_width, alpha_decay_intensity, alpha_decay_half_life, alpha_decay_branching_ratio
```

##### α 衰变能

```references
  - D._S._Delion_Y.2006_Phys.Rev.C_Vol.73Nol.014315P.: abstract.
```

```math
E_{\mathrm{res}}=Q_\alpha^{(\mathrm{exp})}
```
- **Definition**: \(Q_\alpha\) 是实验 α 衰变能；模型通过调节排斥势深度，使共振态能量与其一致。
- **description**: 它连接实验释放能与势阱内衰变态的共振能量。

##### α 衰变宽度

```references
  - D._S._Delion_Y.2006_Phys.Rev.C_Vol.73Nol.014315P.: section II.B. Resonant states.
```

```math
\Gamma=\sum_J\Gamma_J
=\sum_J\hbar v_J\lim_{R\to\infty}|f_J(R)|^2
=\sum_J\hbar v_J|N_J|^2
```

- **Definition**: 总衰变宽度 \(\Gamma\) 是各通道分宽度 \(\Gamma_J\) 之和。
- **description**: 每个通道的贡献由通道速度与渐近振幅模平方决定。

##### α 衰变精细结构强度

```references
  - D._S._Delion_Y.2006_Phys.Rev.C_Vol.73Nol.014315P.: section II.B. Resonant states.
```

```math
I_J=\log_{10}\!\left(\frac{\Gamma_0}{\Gamma_J}\right)
```

- **Definition**: \(I_J\) 比较基态通道宽度 \(\Gamma_0\) 与第 \(J\) 个激发态通道宽度 \(\Gamma_J\)。
- **description**: 该量用于描述精细结构，并避免使用具有模型依赖性的阻碍因子。

##### α 衰变半衰期

```references
  - D._S._Delion_Y.2006_Phys.Rev.C_Vol.73Nol.014315P.: section III, Fig. 5.
```

```math
\log_{10}T
```

- **Definition**: \(T\) 表示总 α 衰变半衰期，系统学比较采用其十进对数。
- **description**: 在所考察的参数变化中，\(Q_\alpha\) 与 \(\log_{10}T\) 对淬火参数 \(v_a\) 表现出较强依赖。

##### α 衰变分支比

```references
  - D_S_Delion_Y.2018_J.Phys.G_Nucl.Part.Phys._Vol.45Nol.P.053001: section 1. Introduction.
```

```math
\mathrm{b.r.}
```

- **Definition**: 分支比表示 α 衰变进入给定转动带成员通道的相对分支。
- **description**: 每个转动带成员对应一个通道，并且分支比随成员能级升高而明显降低。



### Related Methods

##### 耦合道方法

```tags
[tags]: coupled_channels_method, rotational_model, double_folding_potential, m3y_interaction
```

```references
  - D._S._Delion_Y.2006_Phys.Rev.C_Vol.73Nol.014315P.: abstract.
  - D._S._Delion_Y.2006_Phys.Rev.C_Vol.73Nol.014315P.: section IV, Fig. 6.
  - D._S._Delion_Y.2006_Phys.Rev.C_Vol.73Nol.014315P.: section IV, Fig. 7.
```

- **description**: 稳态耦合道方法用于分析偶偶核到转动态的 α 衰变；集体激发采用刚性转子模型，α–子核相互作用由采用 M3Y 与库仑核子—核子力的双折叠势给出。
- **Outputs**: 计算的总半衰期几乎复现实验值。对 \(I_4\)，仅 \(Z=90\) 同位素链达到很好符合；\(Z=96\) 与 \(Z=98\) 链的偏差在半个数量级内，而 \(Z=94\) 附近的中心峰未能复现。



### Previous Studies

##### D. S. Delion (2006)

```tags
[tags]: alpha_decay, alpha_decay_fine_structure, alpha_decay_width, alpha_decay_half_life, coupled_channels_method, m3y_interaction, alpha_clustering, rotational_nuclei
```

```references
  - D._S._Delion_Y.2006_Phys.Rev.C_Vol.73Nol.014315P.: abstract; section I. Introduction; section V. Conclusions.
```

- **Motivation**: 检验微观计算的碎片间势能否同时描述总衰变宽度以及到转动能级的复杂衰变宽度分布。
- **Innovation or Result**: 指出 α 衰变精细结构可作为探测内部区核相互作用及其对 α 团簇依赖性的敏感工具。
- **Quantities**: 总半衰期、到 \(J=2^+\) 与 \(J=4^+\) 转动态的衰变宽度，以及 M3Y 相互作用的淬火因子 \(v_a\)。
- **Precision**: 采用 \(v_a=0.668-0.004(A-208)\) 时，总半衰期及到 \(J=2^+\) 态的衰变宽度与实验符合良好；到 \(J=4^+\) 态的结果对 \(Z=90\) 中子链符合良好，对 \(Z=92,96,98\) 链给出令人满意的描述。
- **Meaning**: 对多个同位素链的系统符合表明，到转动态的 α 衰变能够约束有效核子–核子相互作用随质量数及 α 团簇结构的变化。

##### D S Delion (2018)

```tags
[tags]: alpha_decay, alpha_decay_fine_structure, coupled_channels_method, semiclassical_method, alpha_decay_branching_ratio, alpha_decay_intensity, excitation_energy
```

```references
  - D_S_Delion_Y.2018_J.Phys.G_Nucl.Part.Phys._Vol.45Nol.P.053001: section 1. Introduction; section 7.1. Even-even emitters, Fig. 12; section 8. Conclusions.
```

- **Motivation**: 到子核激发态的 α 跃迁对核结构细节十分敏感，可用于研究低激发态结构。
- **Innovation or Result**: 耦合道计算给出精细结构强度 \(\vartheta_2\) 与子核激发能 \(E_2\) 之间的线性关联。
- **Quantities**: α 衰变分支比、精细结构强度 \(\vartheta_2\) 和激发能 \(E_2\)。
- **Precision**: 对部分偶偶 Pu、Cm 和 Cf 发射体，半经典计算将到 \(4^+\) 激发态的分支比高估约一个数量级，而耦合道结果与实验符合良好。
- **Meaning**: 激发态通道间的耦合不可忽略；耦合道方法比独立通道的半经典处理更可靠地描述 α 衰变精细结构。

##### Jianmin Dong (2025)

```tags
[tags]: alpha_decay, superheavy_nuclei, RMF, symmetry_energy, alpha_decay_energy, alpha_decay_half_life
```

```references
  - Jianmin_Dong_Y.2025_NuclearPhysicsA_Vol.1059Nol.123089P.123089: section 1. Introduction; section 3.3. α-decay, Fig. 6; section 4. Summary.
```

- **Motivation**: 超重核的 α 衰变链连接未知母核与已知子核，是鉴别新元素的重要依据。
- **Innovation or Result**: 形变 RMF–BCS 计算表明，较柔软的对称能给出更小的 α 衰变能 \(Q_\alpha\) 和更长的半衰期。
- **Quantities**: 对称能密度依赖、\(Q_\alpha\) 与 α 衰变半衰期。
- **Precision**: 对 \(^{270}\mathrm{Hs}\) 与 \(^{304}120\)，当 \(\Lambda_V\) 从 \(0\) 增至 \(0.04\) 时，半衰期增加超过两个数量级。
- **Meaning**: α 衰变能与半衰期对对称能的密度依赖具有显著敏感性，可作为超重核结构研究的观测量。



### Next Papers

##### Candidate Paper 1

```tags
[tags]: coupled_channels_method
[doi]: 10.1103/PhysRevC.68.041303
[title]: Probing mean field of neutron rich nuclei by cold fission
```

- **Citation**: Physical Review C 68, 041303(R), 2003.
- **Reason**: Delion 2006 将其作为稳态耦合道形式主要理论构件的前置来源；适合补充耦合道方法从冷裂变双精细结构研究到 α 衰变计算的沿革。

##### Candidate Paper 2

```tags
[tags]: double_folding_potential, alpha_decay
[doi]: 10.1016/0375-9474(77)90392-X
[title]: Interactions for inelastic scattering derived from realistic potentials
```

- **Citation**: Nuclear Physics A 284, 399–419, 1977.
- **Reason**: Delion 2018 将其列为构造 α–子核双折叠势所用有效相互作用的基础来源；适合补充折叠势的微观相互作用依据。

##### Candidate Paper 3

```tags
[tags]: double_folding_potential, alpha_decay
[doi]: 10.1016/0370-1573(79)90081-4
[title]: Folding model potentials from realistic interactions for heavy-ion scattering
```

- **Citation**: Physics Reports 55, 183–254, 1979.
- **Reason**: Delion 2018 将其列为 α–子核双折叠势积分的基础来源；适合补充由真实相互作用构造重离子折叠势的方法体系。

##### Candidate Paper 4

```tags
[tags]: double_folding_potential, alpha_decay
[doi]: 10.1016/0003-4916(92)90153-D
[title]: A new method of evaluating folding type integrals
```

- **Citation**: Annals of Physics 217, 279–303, 1992.
- **Reason**: Delion 2018 将其作为双折叠积分的计算方法来源；适合补充 α–子核势中折叠型积分的数值求值方法。

##### Candidate Paper 5

```tags
[tags]: excitation_energy, alpha_decay_half_life, alpha_decay_energy
[doi]: 10.1006/ndsh.1998.0009
[title]: Review of Alpha–Decay Data from Doubly–Even Nuclei
```

- **Citation**: Nuclear Data Sheets 84, 1–114, 1998.
- **Reason**: Delion 2006 使用该汇编提供实验激发能、总半衰期和 Q 值；适合补充 α 衰变系统学计算所依据的实验数据来源。

##### Candidate Paper 6

```tags
[tags]: quadrupole_deformation
[doi]: 10.1006/adnd.1995.1002
[title]: Nuclear Ground-State Masses and Deformations
```

- **Citation**: Atomic Data and Nuclear Data Tables 59, 185–381, 1995.
- **Reason**: Delion 2006 使用该文的系统学结果提供形变参数；适合补充 α 衰变耦合道计算所需核形变输入的来源。

##### Candidate Paper 7

```tags
[tags]: coupled_channels_method, double_folding_potential
[doi]: 10.1103/PhysRevC.78.034608
[title]: α-decay spectroscopy of deformed nuclei reexamined
```

- **Citation**: Physical Review C 78, 034608, 2008.
- **Reason**: Delion 2018 将其列为使用耦合道形式和双折叠势重新分析形变核 α 衰变数据的研究；适合补充该方法在转动核中的应用。

##### Candidate Paper 8

```tags
[tags]: coupled_channels_method, cluster_core_model, odd_mass_nuclei
[doi]: 10.1103/PhysRevC.86.054608
[title]: Systematic calculation of fine structure in the α decay of heavy odd-mass nuclei
```

- **Citation**: Physical Review C 86, 054608, 2012.
- **Reason**: Delion 2018 指出该研究把耦合道框架中的多通道团簇模型由偶偶转动核扩展到重奇质量核；适合补充奇质量核 α 衰变精细结构的系统计算。

##### Candidate Paper 9

```tags
[tags]: alpha_decay_half_life, alpha_decay_energy
[doi]: 10.1103/PhysRevC.77.037602
[title]: Recent α decay half-lives and analytic expression predictions including superheavy nuclei
```

- **Citation**: Physical Review C 77, 037602, 2008.
- **Reason**: Dong 2025 使用该文给出的 Royer 公式估算 α 衰变半衰期；适合补充包括超重核在内的半衰期解析表达式及其预测。

##### Candidate Paper 10

```tags
[tags]: alpha_decay_energy, alpha_decay_half_life
[doi]: 10.1006/adnd.1997.0746
[title]: Nuclear properties for astrophysical and radioactive-ion-beam applications
```

- **Citation**: Atomic Data and Nuclear Data Tables 66, 131–343, 1997.
- **Reason**: Dong 2025 引用该文说明重元素区 Qα 误差对半衰期的显著放大效应；适合补充 α 衰变能与半衰期敏感性关系的数据依据。

##### Candidate Paper 11

```tags
[tags]: double_folding_potential
[doi]: 10.1103/PhysRevC.63.034007
[title]: α-nucleus optical potential in the double-folding model
```

- **Citation**: Physical Review C 63, 034007, 2001.
- **Reason**: Delion 2006 将其作为双折叠 α–子核势中标准 α 粒子高斯密度参数的来源；适合补充 α–核光学势的微观双折叠构造及密度输入。

##### Candidate Paper 12

```tags
[tags]: pairing_correlation, BCS
[doi]: 10.1140/epja/i2006-10147-0
[title]: Improved BCS-type pairing for the relativistic mean-field theory
```

- **Citation**: European Physical Journal A 30, 519–529, 2006.
- **Reason**: Dong 2025 采用该文的密度依赖 δ 相互作用处理同位旋矢量道的 BCS 型配对关联；适合补充形变 RMF–BCS 框架中的配对项。

##### Candidate Paper 13

```tags
[tags]: FSUGarnet, symmetry_energy, density_dependent
[doi]: 10.1093/mnras/stae2170
[title]: Effects of nucleon–nucleon short-range correlation and symmetry energy on the evolution of newly born magnetars
```

- **Citation**: Monthly Notices of the Royal Astronomical Society 534, 1763–1774, 2024.
- **Reason**: Dong 2025 引用该文所用程序，通过调节 \(\Lambda_V\) 与 \(g_\rho\) 构造具有不同对称能密度依赖的 FSUGarnet 参数族；适合补充参数族的构造来源。

##### Candidate Paper 14

```tags
[tags]: two_nucleon_shell_gap, shell_closure
[doi]: 10.1103/PhysRevC.56.238
[title]: Superheavy nuclei in self-consistent nuclear calculations
```

- **Citation**: Physical Review C 56, 238–243, 1997.
- **Reason**: Dong 2025 将其作为用双核子壳隙判断壳闭合与壳稳定性的来源；适合补充超重核双中子、双质子壳隙指标及其计算式。

##### Candidate Paper 15

```tags
[tags]: RMF, symmetry_energy, superheavy_nuclei, single_particle_energy_level
[doi]: 10.1103/PhysRevC.81.044306
[title]: Effects of the density dependence of the nuclear symmetry energy on the properties of superheavy nuclei
```

- **Citation**: Physical Review C 81, 044306, 2010.
- **Reason**: Dong 2025 将其作为 RMF 框架下对称能软化影响球形超重核轨道移动的前置研究；适合补充当前形变超重核研究的比较基线。
