---
[tags]: symmetry_energy
[tag_types]: quantity
---

# symmetry_energy



### Definition & Origin

```references
  - Jianmin_Dong_Y.2025_NuclearPhysicsA_Vol.1059Nol.123089P.123089: section 2.2. Symmetry energy.
```

对于总密度 \(\rho=\rho_p+\rho_n\)、同位旋不对称度
\(\delta=(\rho_n-\rho_p)/(\rho_n+\rho_p)\) 的均匀核物质，对称能
\(E_{\mathrm{sym}}(\rho)\) 是每核子能量中相对于对称核物质的二阶不对称项系数：

```math
E(\rho,\delta)=E(\rho,0)+E_{\mathrm{sym}}(\rho)\delta^2.
```



### Composition & Mechanism

##### Nonlinear RMF expression

```references
  - Subrata_Pal_Y.2009_PhysicsLettersB_Vol.682Nol.P.23-26: main text, nonlinear RMF model; Eq. (2).
```

- **description**: 在该非线性相对论平均场拉氏量中，对称能写成两个相加项，并由费米动量 \(k_F\)、核子有效质量 \(m^*\) 和有效 \(\rho\)-介子质量 \(m_\rho^*\) 决定；非线性耦合 \(\Lambda_v\) 显式出现在 \(m_\rho^*\) 的定义中。

```math
\begin{aligned}
E_{\rm sym}(\rho)
&= \frac{k_F^2}{6E_F^*}
 + \frac{g_\rho^2}{12\pi^2m_\rho^{*2}}k_F^3,\\
E_F^* &= \sqrt{k_F^2+m^{*2}},\\
m^* &= m-g_s\phi_0,\\
m_\rho^{*2} &= m_\rho^2+2g_\rho^2(\Lambda_vg_v^2V_0^2).
\end{aligned}
```



### Quantities & Properties

```tags
[tags]: slope_parameter, density_dependent, electron_fraction, neutron_skin_thickness, ground_state_energy, single_particle_energy_level
```

##### Density-dependence parameters

```references
  - Subrata_Pal_Y.2009_PhysicsLettersB_Vol.682Nol.P.23-26: main text, Eq. (3) discussion.
```

```math
L=3\rho_0\left.\frac{\partial E_{\rm sym}(\rho)}{\partial\rho}\right|_{\rho_0},
\qquad
K_{\rm sym}=9\rho_0^2\left.\frac{\partial^2E_{\rm sym}(\rho)}{\partial\rho^2}\right|_{\rho_0}.
```

- **Definition**: 在饱和密度 \(\rho_0\)，斜率参数 \(L\) 与曲率参数 \(K_{\rm sym}\) 分别由 \(E_{\rm sym}(\rho)\) 对密度的一阶和二阶导数定义。
- **description**: 这两个参数刻画对称能在饱和密度处的密度依赖。

##### Neutron-skin constraint

```references
  - Subrata_Pal_Y.2009_PhysicsLettersB_Vol.682Nol.P.23-26: main text, Fig. 2 discussion; Table 1.
```

```math
S=R_n-R_p.
```

- **Definition**: 中子皮厚度 \(S\) 定义为中子半径 \(R_n\) 与质子半径 \(R_p\) 之差。
- **description**: 在 NL3 与 FSUGold 模型中，\(\Lambda_v\geq0.02\) 的柔软对称能参数组可以再现中子皮数据，并给出斜率参数的保守上限 \(L\sim60\,\mathrm{MeV}\)。

##### Electron fraction and neutron-skin response

```references
  - C._J._Horowitz_Y.2001_Phys.Rev.Lett._Vol.86Nol.P.5647-5650: Fig. 2 discussion.
```

- **description**: 在改变 \(\Lambda_{\rm v}\) 的模型中，电子分数 \(Y_e\) 由对称能决定，而 \(R_n-R_p\) 对对称能的密度依赖性敏感。

##### Ground-state energy response

```references
  - Jianmin_Dong_Y.2025_NuclearPhysicsA_Vol.1059Nol.123089P.123089: section 3.1, Fig. 2.
```

- **description**: 随着 \(\Lambda_V\) 增大，对称能变软且在次饱和密度下数值增大，基态能量随之增加。

##### Single-particle level response

```references
  - Jianmin_Dong_Y.2025_NuclearPhysicsA_Vol.1059Nol.123089P.123089: section 3.2, Fig. 4.
```

- **description**: 随着 \(\Lambda_V\) 增大、对称能变软，质子单粒子能级及其费米能级显著下移，而中子单粒子能级整体上移。



### Related Methods

##### Nonlinear meson coupling

```tags
[tags]: nonlinear_meson_coupling
```

```references
  - Subrata_Pal_Y.2009_PhysicsLettersB_Vol.682Nol.P.23-26: main text, model paragraph before Table 1; Fig. 1 discussion.
```

- **description**: 在 NL3 拉氏量中加入并调节同位旋矢量耦合 \(\Lambda_v\)，同时在 NL3 与 FSUGold 中改变该耦合，以生成不同的 \(E_{\rm sym}(\rho)\)。
- **Outputs**: 这些参数组显示，超饱和密度 \(\rho\geq\rho_0\) 下较刚性的对称能对应次饱和密度下较柔软的密度依赖；precision: none.

##### Nonlinear sigma–rho and omega–rho coupling

```tags
[tags]: nonlinear_meson_coupling
```

```references
  - C._J._Horowitz_Y.2001_Phys.Rev.Lett._Vol.86Nol.P.5647-5650: methods, model setup; Fig. 2 discussion.
```

- **description**: 在拉格朗日量中加入新的非线性 \(\sigma-\rho\) 与 \(\omega-\rho\) 耦合，以改变对称能的密度依赖性及其对 \(^{208}\mathrm{Pb}\) 中子皮和中子星壳层的影响。
- **Outputs**: 电子分数 \(Y_e\) 由对称能决定，而 \(R_n-R_p\) 对对称能的密度依赖性敏感；precision: none.

##### FSUGarnet interaction family

```tags
[tags]: FSUGarnet, nonlinear_meson_coupling
```

```references
  - Jianmin_Dong_Y.2025_NuclearPhysicsA_Vol.1059Nol.123089P.123089: section 2.2 Symmetry energy; Table 2.
```

- **description**: 以原始 FSUGarnet 相互作用为基础，同时调节非线性 \(\omega-\rho\) 耦合参数 \(\Lambda_V\) 与 \(\rho\)-核子耦合参数 \(g_\rho\)，构造用于研究不同密度依赖对称能的相互作用族。
- **Outputs**: 对称能在饱和密度处保持固定，其他参数不变时对称核物质状态方程也保持不变；precision: 固定点为 \(\rho=0.153\,\mathrm{fm}^{-3}\)。

##### Deformed RMF and BCS

```tags
[tags]: RMF, BCS
```

```references
  - Jianmin_Dong_Y.2025_NuclearPhysicsA_Vol.1059Nol.123089P.123089: section 1 Introduction; section 3.1, Fig. 2; section 3.2, Fig. 4.
```

- **description**: 采用形变 RMF 理论并结合 BCS 方法完成核结构计算。
- **Outputs**: 对称能变软时，基态能量增加，质子单粒子能级及其费米能级显著下移，而中子单粒子能级整体上移；precision: none.



### Previous Studies

##### Wu and Dong (2025)

```tags
[tags]: symmetry_energy, nonlinear_meson_coupling, RMF, alpha_decay_energy, alpha_decay_half_life
```

```references
  - Jianmin_Dong_Y.2025_NuclearPhysicsA_Vol.1059Nol.123089P.123089: section 1 Introduction; section 3.3, Fig. 6; section 4 Summary.
```

- **Motivation**: 不同的密度依赖对称能可能改变费米面附近的单粒子能级，进而影响超重核稳定性。
- **Innovation or Result**: 通过调节非线性 \(\omega-\rho\) 耦合 \(\Lambda_V\) 与 \(\rho\)-核子耦合 \(g_\rho\)，构造 FSUGarnet0.153 相互作用族，同时将饱和密度 \(\rho=0.153\,\mathrm{fm}^{-3}\) 处的对称能固定。
- **Quantities**: 较柔软的对称能给出较小的 \(\alpha\) 衰变能 \(Q_\alpha\) 和较长的半衰期。
- **Precision**: 对 \(^{270}\mathrm{Hs}\) 与 \(^{304}120\)，从 \(\Lambda_V=0\) 到 \(0.04\)，半衰期增加超过两个数量级。
- **Meaning**: 对称能软化时基态能量增加，而核形变几乎不变。

##### Horowitz and Piekarewicz (2001)

```tags
[tags]: symmetry_energy, nonlinear_meson_coupling, RMF
```

```references
  - C._J._Horowitz_Y.2001_Phys.Rev.Lett._Vol.86Nol.P.5647-5650: abstract; methods, model setup; conclusion, points 1–2.
```

- **Motivation**: 既有核观测量不能严格约束中子皮厚度，而中子皮厚度测量可以约束对称能的密度依赖性。
- **Innovation or Result**: 在拉格朗日量中加入新的非线性 \(\sigma-\rho\) 与 \(\omega-\rho\) 耦合，以调节对称能的密度依赖性及其对 \(^{208}\mathrm{Pb}\) 中子皮和中子星壳层的影响。
- **Quantities**: 在中子丰核物质中，预测较厚中子皮的模型同时给出较大的电子分数和较低的液—固相变密度。
- **Precision**: none.
- **Meaning**: 这些模型给出的高密度对称能随密度增长得比早期相对论平均场模型更慢。

##### Sharma and Pal (2009)

```tags
[tags]: symmetry_energy, nonlinear_meson_coupling, RMF, slope_parameter, soft
```

```references
  - Subrata_Pal_Y.2009_PhysicsLettersB_Vol.682Nol.P.23-26: Introduction; Abstract; main text, Fig. 2 discussion; Summary.
```

- **Motivation**: 对称能的密度依赖性仍缺乏约束，而中子星为约束其高密度行为提供了途径。
- **Innovation or Result**: 现有中子皮数据支持亚饱和密度下较硬、超饱和密度下较软的对称能。
- **Quantities**: 研究以中子皮厚度 \(S=R_n-R_p\) 和对称能斜率参数 \(L\) 表征密度依赖性。
- **Precision**: 在 NL3 和 FSUGold 模型中，\(\Lambda_v\geq0.02\) 的柔软对称能参数组给出保守上限 \(L\sim60\,\mathrm{MeV}\)。
- **Meaning**: 两个经过校准的 RMF 模型与稳定核中子皮数据的比较共同指向整体较柔软的对称能。



### Next Papers

##### Candidate Paper 1

```tags
[tags]: symmetry_energy, density_dependent
[doi]: 10.1093/mnras/stae2170
[title]: Effects of nucleon–nucleon short-range correlation and symmetry energy on the evolution of newly born magnetars
```

- **Citation**: Monthly Notices of the Royal Astronomical Society 534(3), 1763–1774, 2024.
- **Reason**: Dong 2025 将其列为调节 FSUGarnet 相互作用族中 \(\Lambda_V\) 与 \(g_\rho\) 的方法来源；适合补充密度依赖对称能的方法沿革。

##### Candidate Paper 2

```tags
[tags]: symmetry_energy, RMF
[doi]: 10.1103/PhysRevC.81.044306
[title]: Effects of the density dependence of the nuclear symmetry energy on the properties of superheavy nuclei
```

- **Citation**: Physical Review C 81, 044306, 2010.
- **Reason**: Dong 2025 将其列为 RMF 框架下对称能软化影响超重核球形轨道移动的原始研究；适合补充该机制与结果的一手证据。

##### Candidate Paper 3

```tags
[tags]: parity_violating_electron_scattering, Pb-208
[doi]: 10.1103/PhysRevC.63.025501
[title]: Parity violating measurements of neutron densities
```

- **Citation**: Physical Review C 63, 025501, 2001.
- **Reason**: Horowitz 2001 在 PREX 讨论中引用该文，以说明通过宇称破缺电子散射测量 \(^{208}\mathrm{Pb}\) 中子密度的方法；适合补充中子皮观测约束对称能密度依赖性的直接方法来源。

##### Candidate Paper 4

```tags
[tags]: RMF
[doi]: 10.1016/0375-9474(96)00187-X
[title]: Relativistic mean-field theory and the high-density nuclear equation of state
```

- **Citation**: Nuclear Physics A 606, 508–537, 1996.
- **Reason**: Horowitz 2001 将该文作为描述有限核并外推致密中子丰物质状态方程的相对论有效场论基础；适合补充非线性耦合扩展所依托的 RMF 模型来源。

##### Candidate Paper 5

```tags
[tags]: random_phase_approximation, transition_density
[doi]: 10.1016/0375-9474(94)00506-I
[title]: The inner boundary of a neutron-star crust
```

- **Citation**: Nuclear Physics A 584, 675–703, 1995.
- **Reason**: Horowitz 2001 在利用 RPA 估计均匀中子丰物质失稳与壳层相变密度时引用该文；适合补充中子星壳层内边界和相变密度方法的一手来源。

##### Candidate Paper 6

```tags
[tags]: random_phase_approximation, transition_density
[doi]: 10.1016/S0370-2693(00)00672-9
[title]: Inner edge of neutron-star crust with SLy effective nucleon-nucleon interactions
```

- **Citation**: Physics Letters B 485, 107–114, 2000.
- **Reason**: Horowitz 2001 在壳层相变密度计算及其下界讨论中引用该文；适合补充采用 SLy 有效核子相互作用确定中子星壳层内缘的结果。

##### Candidate Paper 7

```tags
[tags]: nuclear_matter
[doi]: 10.1016/0375-9474(81)90649-7
[title]: Hot and cold, nuclear and neutron matter
```

- **Citation**: Nuclear Physics A 361, 502–520, 1981.
- **Reason**: Horowitz 2001 引用该文的非相对论微观状态方程，并将其用于壳层相变密度比较；适合补充核物质状态方程与中子星壳层性质之间的原始模型来源。

##### Candidate Paper 8

```tags
[tags]: transition_density, neutron_star_crust
[doi]: 10.1103/PhysRevLett.70.379
[title]: Neutron star crusts
```

- **Citation**: Physical Review Letters 70, 379–382, 1993.
- **Reason**: Horowitz 2001 引用该文给出 Friedman–Pandharipande 状态方程对应的壳层相变密度 \(\rho_c=0.096\,\mathrm{fm}^{-3}\)；适合补充中子星壳层相变密度计算的一手结果。

##### Candidate Paper 9

```tags
[tags]: neutron_skin_thickness, nuclear_matter
[doi]: 10.1103/PhysRevLett.85.5296
[title]: Neutron Radii in Nuclei and the Neutron Equation of State
```

- **Citation**: Physical Review Letters 85, 5296–5299, 2000.
- **Reason**: Horowitz 2001 引用该文得到同一状态方程下 \(R_n-R_p=0.16\pm0.02\,\mathrm{fm}\)，并用于检验中子皮厚度与壳层相变密度的相关性；适合补充中子半径和中子物质状态方程之间的直接证据。

##### Candidate Paper 10

```tags
[tags]: BCS, pairing_correlation
[doi]: 10.1140/epja/i2006-10147-0
[title]: Improved BCS-type pairing for the relativistic mean-field theory
```

- **Citation**: European Physical Journal A 30, 519–529, 2006.
- **Reason**: Dong 2025 引用该文作为同位旋矢量 BCS 配对中密度依赖 \(\delta\) 相互作用的来源；适合补充形变 RMF+BCS 计算所采用的配对方法。

##### Candidate Paper 11

```tags
[tags]: two_nucleon_shell_gap, shell_closure
[doi]: 10.1103/PhysRevC.56.238
[title]: Superheavy nuclei in self-consistent nuclear calculations
```

- **Citation**: Physical Review C 56, 238–243, 1997.
- **Reason**: Dong 2025 引用该文定义双中子与双质子壳隙，并将其作为判断超重核壳层稳定性的指标；适合补充双核子壳隙方法的原始来源。

##### Candidate Paper 12

```tags
[tags]: alpha_decay_energy, alpha_decay_half_life
[doi]: 10.1006/adnd.1997.0746
[title]: Nuclear properties for astrophysical and radioactive-ion-beam applications
```

- **Citation**: Atomic Data and Nuclear Data Tables 66, 131–343, 1997.
- **Reason**: Dong 2025 引用该文说明重核区 \(Q_\alpha\) 的 \(1\,\mathrm{MeV}\) 不确定性可使半衰期产生约 \(10^3\)–\(10^5\) 倍的不确定性；适合补充 \(\alpha\) 衰变预测的误差敏感性来源。

##### Candidate Paper 13

```tags
[tags]: alpha_decay_half_life, alpha_decay_energy
[doi]: 10.1103/PhysRevC.77.037602
[title]: Recent α decay half-lives and analytic expression predictions including superheavy nuclei
```

- **Citation**: Physical Review C 77, 037602, 2008.
- **Reason**: Dong 2025 使用该文给出的 Royer 公式估算偶偶核的 \(\alpha\) 衰变半衰期；适合补充半衰期解析计算方法的原始来源。

##### Candidate Paper 14

```tags
[tags]: RMF, nonlinear_meson_coupling, direct_Urca_process
[doi]: 10.1103/PhysRevC.66.055803
[title]: Constraining URCA cooling of neutron stars from the neutron radius of 208 Pb
```

- **Citation**: Physical Review C 66, 055803, 2002.
- **Reason**: Sharma–Pal 2009 同时将其用于非线性 RMF 相互作用和直接 Urca 冷却讨论；适合补充中子半径、对称能与直接 Urca 阈值之间的联系。

##### Candidate Paper 15

```tags
[tags]: RMF, NL3, finite_nuclei
[doi]: 10.1103/PhysRevC.55.540
[title]: New parametrization for the Lagrangian density of relativistic mean field theory
```

- **Citation**: Physical Review C 55, 540–543, 1997.
- **Reason**: Sharma–Pal 2009 将其作为 NL3 参数化来源；适合补充由有限核基态性质标定 RMF 模型的原始证据。

##### Candidate Paper 16

```tags
[tags]: RMF, FSUGold, finite_nuclei
[doi]: 10.1103/PhysRevLett.95.122501
[title]: Neutron-Rich Nuclei and Neutron Stars: A New Accurately Calibrated Interaction for the Study of Neutron-Rich Matter
```

- **Citation**: Physical Review Letters 95, 122501, 2005.
- **Reason**: Sharma–Pal 2009 将其作为 FSUGold 相互作用的来源；适合补充面向有限核和中子丰物质的模型标定依据。

##### Candidate Paper 17

```tags
[tags]: direct_Urca_process
[doi]: 10.1016/j.physrep.2005.02.004
[title]: Isospin asymmetry in nuclei and neutron stars
```

- **Citation**: Physics Reports 411, 325–375, 2005.
- **Reason**: Sharma–Pal 2009 将其用于直接 Urca 快速冷却条件的背景；适合补充核与中子星同位旋不对称性的综述性来源。

##### Candidate Paper 18

```tags
[tags]: direct_Urca_process
[doi]: 10.1103/PhysRevLett.66.2701
[title]: Direct URCA process in neutron stars
```

- **Citation**: Physical Review Letters 66, 2701–2704, 1991.
- **Reason**: Sharma–Pal 2009 将其作为直接 Urca 快速中微子冷却条件的原始来源；适合补充该过程及其阈值的一手证据。

##### Candidate Paper 19

```tags
[tags]: RMF, kaon_condensation
[doi]: 10.1016/0375-9474(94)90407-3
[title]: Composition, structure and evolution of neutron stars with kaon condensates
```

- **Citation**: Nuclear Physics A 572, 693–731, 1994.
- **Reason**: Sharma–Pal 2009 在研究 \(K^-\) 凝聚对中子星半径与中子皮关联的影响时引用该文；适合补充含 kaon condensate 中子星的组成和结构。

##### Candidate Paper 20

```tags
[tags]: RMF, kaon_condensation
[doi]: 10.1103/PhysRevC.52.3470
[title]: Strangeness in hadronic stellar matter
```

- **Citation**: Physical Review C 52, 3470–3482, 1995.
- **Reason**: Sharma–Pal 2009 将其作为中子星中 \(K^-\) 凝聚研究的来源；适合补充奇异强子物质和 kaon condensation 的模型依据。

##### Candidate Paper 21

```tags
[tags]: RMF, kaon_condensation
[doi]: 10.1103/PhysRevC.53.1416
[title]: Hyperon-rich matter in neutron stars
```

- **Citation**: Physical Review C 53, 1416–1429, 1996.
- **Reason**: Sharma–Pal 2009 在构建反 kaon 最小耦合 RMF 拉氏量时引用该文；适合补充含奇异自由度中子星物质的相对论模型来源。

##### Candidate Paper 22

```tags
[tags]: RMF, kaon_condensation
[doi]: 10.1016/S0375-9474(00)00175-5
[title]: Antikaon condensation in neutron stars
```

- **Citation**: Nuclear Physics A 674, 553–577, 2000.
- **Reason**: Sharma–Pal 2009 将其作为反 kaon 最小耦合 RMF 拉氏量的来源；适合补充中子星反 kaon 凝聚模型的一手研究。
