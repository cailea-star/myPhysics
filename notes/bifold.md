---
[tags]: double_folding_potential
[tag-type]: method
[ref]: 10.1016/j.cpc.2022.108613
---

# BiFold

### Recursive Reduction of an Arbitrary-Order Folding Potential

Given $n$ density distributions $\rho_i(\mathbf r_i)$, define the $n$-th-order folding potential as

$$
V_n(\mathbf R)
=
\int
\prod_{i=1}^{n} d\mathbf r_i\,
\rho_i(\mathbf r_i)\,
v\left(
\mathbf R-\sum_{i=1}^{n}\mathbf r_i
\right).
$$

Separating the first density distribution from the integral gives

$$
V_n(\mathbf R)
=
\int d\mathbf r_1\,\rho_1(\mathbf r_1)
\left[
\int
\prod_{i=2}^{n} d\mathbf r_i\,
\rho_i(\mathbf r_i)\,
v\left(
\mathbf R-\mathbf r_1-\sum_{i=2}^{n}\mathbf r_i
\right)
\right].
$$

Define

$$
V_{n-1}(\mathbf X)
=
\int
\prod_{i=2}^{n} d\mathbf r_i\,
\rho_i(\mathbf r_i)\,
v\left(
\mathbf X-\sum_{i=2}^{n}\mathbf r_i
\right),
$$

The folding potential then satisfies the recurrence relation

$$
V_n(\mathbf R)
=
\int d\mathbf r_1\,
\rho_1(\mathbf r_1)\,
V_{n-1}(\mathbf R-\mathbf r_1)
=
(\rho_1 * V_{n-1})(\mathbf R).
$$

Thus, a higher-order folding integral can be constructed recursively, reducing one high-dimensional integral to a sequence of lower-order convolutions.

### Fourier Transform

See [Fourier Transform](mathematical_physics/fourier_transform.md).

##### Fourier Decomposition of a Density

Applying the preceding decomposition to a density $\rho(\mathbf r)$ gives

$$
\rho(\mathbf r)
=
\sum_{l,m}
\rho_{lm}(r)Y_{lm}(\hat{\mathbf r}),
\qquad
\rho_{lm}(r)
=
\int d\Omega_{\mathbf r}\,
Y_{lm}^{*}(\hat{\mathbf r})\rho(\mathbf r).
$$

Its momentum-space representation is

$$
\tilde\rho(\mathbf q)
=
\sum_{l,m}
i^{-l}\tilde\rho_{lm}(q)Y_{lm}(\hat{\mathbf q}),
$$

with

$$
\tilde\rho_{lm}(q)
=
4\pi
\int_0^\infty r^2\,dr\,
\rho_{lm}(r)j_l(qr).
$$

For an axially symmetric density, only the $m=0$ components remain:

$$
\rho(\mathbf r)
=
\sum_l \rho_{l0}(r)Y_{l0}(\hat{\mathbf r}),
\qquad
\tilde\rho(\mathbf q)
=
\sum_l i^{-l}\tilde\rho_{l0}(q)Y_{l0}(\hat{\mathbf q}).
$$

If the density is also reflection symmetric, only even values of $l$ contribute.

##### Fourier Decomposition of a Central Interaction

For a central interaction $v(\mathbf r)=v(r)$, only the $l=0$ partial wave contributes. Its Fourier transform therefore reduces to

$$
\tilde v(q)
=
4\pi
\int_0^\infty r^2\,dr\,
v(r)j_0(qr).
$$

The corresponding inverse transform is

$$
v(r)
=
\frac{1}{2\pi^2}
\int_0^\infty q^2\,dq\,
\tilde v(q)j_0(qr).
$$

A Yukawa interaction written as

$$
v(s)
=
A\frac{e^{-\mu s}}{\mu s},
$$

where $\mu$ is the inverse range, has the momentum-space representation

$$
\tilde v(q)
=
4\pi
\int_0^\infty s^2\,ds\,
A\frac{e^{-\mu s}}{\mu s}
j_0(qs)
=
\frac{4\pi A}
{\mu(q^2+\mu^2)}.
$$

More generally, for a sum of Yukawa terms,

$$
v(s)
=
\sum_i
A_i\frac{e^{-\mu_i s}}{\mu_i s},
$$

the Fourier transform is

$$
\tilde v(q)
=
4\pi
\sum_i
\frac{A_i}
{\mu_i(q^2+\mu_i^2)}.
$$

Therefore, the momentum-space forms of the M3Y-Reid and M3Y-Paris interactions can be obtained analytically by transforming each Yukawa term separately.

Thus, a central interaction remains spherically symmetric in momentum space and depends only on the magnitude $q=|\mathbf q|$.

### Local Density Approximation

The exchange potential depends on the off-diagonal one-body density matrix $\rho(\mathbf r,\mathbf r+\mathbf s)$ rather than only on the local density $\rho(\mathbf r)$. In the local-density approximation, the neighborhood around the midpoint $\mathbf r+\mathbf s/2$ is treated as locally uniform nuclear matter.

The density matrix is approximated by

$$
\rho(\mathbf r,\mathbf r+\mathbf s)
\approx
\rho\left(\mathbf r+\frac{\mathbf s}{2}\right)
\hat j_1\left[
k_F\left(\mathbf r+\frac{\mathbf s}{2}\right)s
\right],
$$

where $k_F$ is the local Fermi momentum and

$$
\hat j_1(x)
=
\frac{3j_1(x)}{x}.
$$

Because

$$
\lim_{x\to0}\hat j_1(x)=1,
$$

the approximation recovers the diagonal density in the limit $\mathbf s\to0$:

$$
\rho(\mathbf r,\mathbf r)
=
\rho(\mathbf r).
$$

The factor $\hat j_1(k_Fs)$ describes the loss of coherence between two spatial points as their separation increases. The characteristic nonlocality range is of order $1/k_F$.

This approximation is most reliable when the density varies slowly over the separation $s$. It may become less accurate in regions with large density gradients, particularly near the nuclear surface. The exchange term is not converted into a direct local interaction; its approximate nonlocal dependence is retained through $s$ and $k_F$.

### Direct Potential

The direct part of the double-folding potential is

$$
V^{\mathrm D}(\mathbf R)
=
\int d^3\mathbf r_{\mathrm t}
\int d^3\mathbf r_{\mathrm p}\,
\rho_{\mathrm t}(\mathbf r_{\mathrm t})
\rho_{\mathrm p}(\mathbf r_{\mathrm p})
v^{\mathrm D}(s),
\qquad
\mathbf s
=
\mathbf R+\mathbf r_{\mathrm p}-\mathbf r_{\mathrm t}.
$$

Here,

- $\mathbf R$ points from the center of the target to the center of the projectile;
- $\mathbf r_{\mathrm t}$ is the position of a nucleon relative to the center of the target;
- $\mathbf r_{\mathrm p}$ is the position of a nucleon relative to the center of the projectile;
- $\mathbf s$ is the relative displacement between the two interacting nucleons, with $s=|\mathbf s|$.

##### Momentum-Space Representation

Using the Fourier-transform convention defined above,

$$
v^{\mathrm D}(s)
=
\frac{1}{(2\pi)^3}
\int d^3\mathbf q\,
\tilde v^{\mathrm D}(q)
e^{i\mathbf q\cdot\mathbf s}.
$$

Recall that the relative displacement is

$$
\mathbf s
=
\mathbf R+\mathbf r_{\mathrm p}-\mathbf r_{\mathrm t}.
$$

Substitution into the coordinate-space folding integral gives

$$
V^{\mathrm D}(\mathbf R)
=
\frac{1}{(2\pi)^3}
\int d^3\mathbf q\,
e^{i\mathbf q\cdot\mathbf R}
\tilde v^{\mathrm D}(q)
\left[
\int d^3\mathbf r_{\mathrm t}\,
\rho_{\mathrm t}(\mathbf r_{\mathrm t})
e^{-i\mathbf q\cdot\mathbf r_{\mathrm t}}
\right]
\left[
\int d^3\mathbf r_{\mathrm p}\,
\rho_{\mathrm p}(\mathbf r_{\mathrm p})
e^{i\mathbf q\cdot\mathbf r_{\mathrm p}}
\right].
$$

Recognizing the two density transforms,

$$
\tilde\rho_{\mathrm t}(\mathbf q)
=
\int d^3\mathbf r_{\mathrm t}\,
\rho_{\mathrm t}(\mathbf r_{\mathrm t})
e^{-i\mathbf q\cdot\mathbf r_{\mathrm t}},
$$

$$
\tilde\rho_{\mathrm p}(-\mathbf q)
=
\int d^3\mathbf r_{\mathrm p}\,
\rho_{\mathrm p}(\mathbf r_{\mathrm p})
e^{i\mathbf q\cdot\mathbf r_{\mathrm p}},
$$

the direct potential becomes

$$
V^{\mathrm D}(\mathbf R)
=
\frac{1}{(2\pi)^3}
\int d^3\mathbf q\,
e^{i\mathbf q\cdot\mathbf R}
\tilde\rho_{\mathrm t}(\mathbf q)
\tilde\rho_{\mathrm p}(-\mathbf q)
\tilde v^{\mathrm D}(q).
$$

For a real density,

$$
\tilde\rho_{\mathrm p}(-\mathbf q)
=
\tilde\rho_{\mathrm p}(\mathbf q)^{*}.
$$

Thus, the coordinate-space double-folding integral becomes a product of the two density form factors and the interaction in momentum space.

### Exchange Potential

The finite-range knock-on exchange term is localized with a local plane-wave approximation:

$$
V^{\mathrm E}(\mathbf R)
=
\int d^3\mathbf r_{\mathrm t}
\int d^3\mathbf r_{\mathrm p}\,
\rho_{\mathrm t}
(\mathbf r_{\mathrm t},\mathbf r_{\mathrm t}+\mathbf s)
\rho_{\mathrm p}
(\mathbf r_{\mathrm p},\mathbf r_{\mathrm p}-\mathbf s)
v^{\mathrm E}(s)
\exp\left[
i\frac{\mathbf K(\mathbf R)\cdot\mathbf s}{M}
\right],
$$

where

$$
\mathbf s
=
\mathbf R+\mathbf r_{\mathrm p}-\mathbf r_{\mathrm t},
\qquad
M
=
\frac{A_{\mathrm p}A_{\mathrm t}}
{A_{\mathrm p}+A_{\mathrm t}}.
$$

Here, $\mathbf K(\mathbf R)$ is the local relative wave number and $M$ is the dimensionless reduced mass number.

Introducing the relative displacement explicitly gives

$$
V^{\mathrm E}(\mathbf R)
=
\int d^3\mathbf s\,
v^{\mathrm E}(s)
\exp\left[
i\frac{\mathbf K(\mathbf R)\cdot\mathbf s}{M}
\right]
\int d^3\mathbf r_{\mathrm t}
\int d^3\mathbf r_{\mathrm p}\,
\rho_{\mathrm t}
(\mathbf r_{\mathrm t},\mathbf r_{\mathrm t}+\mathbf s)
\rho_{\mathrm p}
(\mathbf r_{\mathrm p},\mathbf r_{\mathrm p}-\mathbf s)
\delta^{(3)}
(\mathbf R+\mathbf r_{\mathrm p}-\mathbf r_{\mathrm t}-\mathbf s).
$$

##### Applying the Local-Density Approximation

For the target and projectile density matrices, introduce the midpoint coordinates

$$
\mathbf u_{\mathrm t}
=
\mathbf r_{\mathrm t}+\frac{\mathbf s}{2},
\qquad
\mathbf u_{\mathrm p}
=
\mathbf r_{\mathrm p}-\frac{\mathbf s}{2}.
$$

The geometric constraint then becomes

$$
\mathbf R+\mathbf r_{\mathrm p}
-\mathbf r_{\mathrm t}-\mathbf s
=
\mathbf R+\mathbf u_{\mathrm p}
-\mathbf u_{\mathrm t}.
$$

The localized exchange potential is therefore

$$
V^{\mathrm E}(\mathbf R)
=
\int d^3\mathbf s\,
v^{\mathrm E}(s)
\exp\left[
i\frac{\mathbf K(\mathbf R)\cdot\mathbf s}{M}
\right]
\int d^3\mathbf u_{\mathrm t}
\int d^3\mathbf u_{\mathrm p}\,
\Big[
\rho_{\mathrm t}(\mathbf u_{\mathrm t})
\hat j_1\left[
k_{F,\mathrm t}(\mathbf u_{\mathrm t})s
\right]\Big]
\Big[
\rho_{\mathrm p}(\mathbf u_{\mathrm p})
\hat j_1\left[
k_{F,\mathrm p}(\mathbf u_{\mathrm p})s
\right]\Big]
\delta^{(3)}
(\mathbf R+\mathbf u_{\mathrm p}-\mathbf u_{\mathrm t}).
$$

##### Momentum-Space Representation

Using

$$
\delta^{(3)}
(\mathbf R+\mathbf u_{\mathrm p}-\mathbf u_{\mathrm t})
=
\frac{1}{(2\pi)^3}
\int d^3\mathbf q\,
e^{i\mathbf q\cdot
(\mathbf R+\mathbf u_{\mathrm p}-\mathbf u_{\mathrm t})},
$$

define the exchange form factors and their Fourier transforms:

$$
f_i(\mathbf u,s)
=
\rho_i(\mathbf u)
\hat j_1\left[k_{F,i}(\mathbf u)s\right],
\qquad
\tilde f_i(\mathbf q,s)
=
\int d^3\mathbf u\,
f_i(\mathbf u,s)e^{-i\mathbf q\cdot\mathbf u},
\qquad
i\in\{\mathrm t,\mathrm p\}.
$$

Then

$$
V^{\mathrm E}(\mathbf R)
=
\frac{1}{(2\pi)^3}
\int d^3\mathbf s\,
v^{\mathrm E}(s)
\exp\left[
i\frac{\mathbf K(\mathbf R)\cdot\mathbf s}{M}
\right]
\int d^3\mathbf q\,
e^{i\mathbf q\cdot\mathbf R}
\bigg[
\int d^3\mathbf u_{\mathrm t}\,
f_{\mathrm t}(\mathbf u_{\mathrm t},s)
e^{-i\mathbf q\cdot\mathbf u_{\mathrm t}}
\bigg]
\bigg[
\int d^3\mathbf u_{\mathrm p}\,
f_{\mathrm p}(\mathbf u_{\mathrm p},s)
e^{i\mathbf q\cdot\mathbf u_{\mathrm p}}
\bigg].
$$

Define

$$
G(\mathbf R,s)
=
\frac{1}{(2\pi)^3}
\int d^3\mathbf q\,
e^{i\mathbf q\cdot\mathbf R}
\tilde f_{\mathrm t}(\mathbf q,s)
\tilde f_{\mathrm p}(-\mathbf q,s).
$$

The exchange potential becomes

$$
V^{\mathrm E}(\mathbf R)
=
\int d^3\mathbf s\,
v^{\mathrm E}(s)
\exp\left[
i\frac{\mathbf K(\mathbf R)\cdot\mathbf s}{M}
\right]
G(\mathbf R,s).
$$

After performing the angular integral over $\mathbf s$,

$$
V^{\mathrm E}(\mathbf R)
=
4\pi
\int_0^\infty s^2\,ds\,
j_0\left[
\frac{K(\mathbf R)s}{M}
\right]
v^{\mathrm E}(s)
G(\mathbf R,s).
$$

##### Self-Consistency

The local relative wave number satisfies

$$
K^2(\mathbf R)
=
\frac{2\mu}{\hbar^2}
\left[
E_{\mathrm{cm}}
-
V^{\mathrm D}(\mathbf R)
-
V^{\mathrm E}(\mathbf R)
-
V^{\mathrm C}(\mathbf R)
\right],
\qquad
\mu
=
M m_N.
$$

Because $V^{\mathrm E}$ depends on $K$ and $K$ depends on $V^{\mathrm E}$, the two quantities must be determined self-consistently at every value of $\mathbf R$.

### Partial-Wave Decomposition

Assume a spherical projectile and a deformed target. The projectile density and form factors contain only the monopole component, whereas the target quantities are expanded in spherical harmonics.

##### Direct Potential

For a spherical projectile,

$$
\rho_{\mathrm p}(\mathbf r_{\mathrm p})
=
\rho_{\mathrm p}(r_{\mathrm p}),
$$

$$
\tilde\rho_{\mathrm p}(q)
=
4\pi
\int_0^\infty r_{\mathrm p}^2\,dr_{\mathrm p}\,
\rho_{\mathrm p}(r_{\mathrm p})
j_0(qr_{\mathrm p}).
$$

For a deformed target,

$$
\rho_{\mathrm t}(\mathbf r_{\mathrm t})
=
\sum_{l,m}
\rho_{\mathrm t,lm}(r_{\mathrm t})
Y_{lm}(\hat{\mathbf r}_{\mathrm t}),
\qquad
\rho_{\mathrm t,lm}(r_{\mathrm t})
=
\int d\Omega_{\mathrm t}\,
Y_{lm}^{*}(\hat{\mathbf r}_{\mathrm t})
\rho_{\mathrm t}(\mathbf r_{\mathrm t}).
$$

Its momentum-space density is

$$
\tilde\rho_{\mathrm t}(\mathbf q)
=
\sum_{l,m}
i^{-l}
\tilde\rho_{\mathrm t,lm}(q)
Y_{lm}(\hat{\mathbf q}),
\qquad
\tilde\rho_{\mathrm t,lm}(q)
=
4\pi
\int_0^\infty r_{\mathrm t}^2\,dr_{\mathrm t}\,
\rho_{\mathrm t,lm}(r_{\mathrm t})
j_l(qr_{\mathrm t}).
$$

Because the spherical projectile density is real and even,

$$
\tilde\rho_{\mathrm p}(-\mathbf q)
=
\tilde\rho_{\mathrm p}(q).
$$

The direct potential can therefore be expanded as

$$
V^{\mathrm D}(\mathbf R)
=
\sum_{l,m}
V^{\mathrm D}_{lm}(R)
Y_{lm}(\hat{\mathbf R}),
$$

with

$$
V^{\mathrm D}_{lm}(R)
=
\frac{4\pi}{(2\pi)^3}
\int_0^\infty q^2\,dq\,
j_l(qR)
\tilde\rho_{\mathrm t,lm}(q)
\tilde\rho_{\mathrm p}(q)
\tilde v^{\mathrm D}(q).
$$

##### Exchange Potential

For the spherical projectile, define

$$
f_{\mathrm p}(r_{\mathrm p},s)
=
\rho_{\mathrm p}(r_{\mathrm p})
\hat j_1\left[
k_{F,\mathrm p}(r_{\mathrm p})s
\right].
$$

Its Fourier transform is

$$
\tilde f_{\mathrm p}(q,s)
=
4\pi
\int_0^\infty r_{\mathrm p}^2\,dr_{\mathrm p}\,
f_{\mathrm p}(r_{\mathrm p},s)
j_0(qr_{\mathrm p}).
$$

For the deformed target,

$$
f_{\mathrm t}(\mathbf r_{\mathrm t},s)
=
\rho_{\mathrm t}(\mathbf r_{\mathrm t})
\hat j_1\left[
k_{F,\mathrm t}(\mathbf r_{\mathrm t})s
\right],
$$

$$
f_{\mathrm t}(\mathbf r_{\mathrm t},s)
=
\sum_{l,m}
f_{\mathrm t,lm}(r_{\mathrm t},s)
Y_{lm}(\hat{\mathbf r}_{\mathrm t}),
$$

where

$$
f_{\mathrm t,lm}(r_{\mathrm t},s)
=
\int d\Omega_{\mathrm t}\,
Y_{lm}^{*}(\hat{\mathbf r}_{\mathrm t})
f_{\mathrm t}(\mathbf r_{\mathrm t},s).
$$

The corresponding momentum-space expansion is

$$
\tilde f_{\mathrm t}(\mathbf q,s)
=
\sum_{l,m}
i^{-l}
\tilde f_{\mathrm t,lm}(q,s)
Y_{lm}(\hat{\mathbf q}),
$$

with

$$
\tilde f_{\mathrm t,lm}(q,s)
=
4\pi
\int_0^\infty r_{\mathrm t}^2\,dr_{\mathrm t}\,
f_{\mathrm t,lm}(r_{\mathrm t},s)
j_l(qr_{\mathrm t}).
$$

The exchange kernel becomes

$$
G(\mathbf R,s)
=
\sum_{l,m}
G_{lm}(R,s)
Y_{lm}(\hat{\mathbf R}),
$$

where

$$
G_{lm}(R,s)
=
\frac{4\pi}{(2\pi)^3}
\int_0^\infty q^2\,dq\,
j_l(qR)
\tilde f_{\mathrm t,lm}(q,s)
\tilde f_{\mathrm p}(q,s).
$$

The exchange potential remains

$$
V^{\mathrm E}(\mathbf R)
=
4\pi
\int_0^\infty s^2\,ds\,
j_0\left[
\frac{K(\mathbf R)s}{M}
\right]
v^{\mathrm E}(s)
G(\mathbf R,s).
$$

If the local wave number is approximated by an angle-independent monopole $K_0(R)$, the exchange multipoles decouple:

$$
V^{\mathrm E}(\mathbf R)
=
\sum_{l,m}
V^{\mathrm E}_{lm}(R)
Y_{lm}(\hat{\mathbf R}),
$$

$$
V^{\mathrm E}_{lm}(R)
=
4\pi
\int_0^\infty s^2\,ds\,
j_0\left[
\frac{K_0(R)s}{M}
\right]
v^{\mathrm E}(s)
G_{lm}(R,s).
$$

For an axially symmetric target, only $m=0$ contributes. If the target is also reflection symmetric, only even values of $l$ remain.

When the full angle-dependent $K(\mathbf R)$ is retained, the factor $j_0[K(\mathbf R)s/M]$ couples different multipoles. In that case, the $(l,m)$ components cannot be solved independently without an additional approximation or an explicit angular recoupling.


### M3Y Interactions

The M3Y effective nucleon–nucleon interaction contains separate direct and exchange components. Two commonly used parametrizations are M3Y-Reid and M3Y-Paris. The expressions below give the density-independent components $v_0^{\mathrm D}(s)$ and $v_0^{\mathrm E}(s)$ used in the density-dependent interaction.

Here, $s$ is measured in $\mathrm{fm}$, the inverse ranges are measured in $\mathrm{fm}^{-1}$, and the resulting interactions are given in $\mathrm{MeV}$.

##### M3Y-Reid

The direct component is

$$
v_{0,\mathrm{Reid}}^{\mathrm D}(s)
=
7999.0
\frac{e^{-4.0s}}{4.0s}
-
2134.25
\frac{e^{-2.5s}}{2.5s}.
$$

The finite-range exchange component is

$$
v_{0,\mathrm{Reid}}^{\mathrm E}(s)
=
4631.38
\frac{e^{-4.0s}}{4.0s}
-
1787.13
\frac{e^{-2.5s}}{2.5s}
-
7.8474
\frac{e^{-0.7072s}}{0.7072s}.
$$

##### M3Y-Paris

The direct component is

$$
v_{0,\mathrm{Paris}}^{\mathrm D}(s)
=
11061.625
\frac{e^{-4.0s}}{4.0s}
-
2537.5
\frac{e^{-2.5s}}{2.5s}.
$$

The finite-range exchange component is

$$
v_{0,\mathrm{Paris}}^{\mathrm E}(s)
=
-1524.25
\frac{e^{-4.0s}}{4.0s}
-
518.75
\frac{e^{-2.5s}}{2.5s}
-
7.8474
\frac{e^{-0.7072s}}{0.7072s}.
$$

Each component is a sum of Yukawa functions. Although the individual terms behave as $1/s$ near $s=0$, they are locally integrable in the three-dimensional folding integral because the radial measure contains $s^2\,ds$. Numerically, the grid may exclude the exact point $s=0$, evaluate the first integration cell analytically, or use an explicitly documented extrapolation. Linear extrapolation at the origin is a numerical prescription rather than part of the M3Y interaction itself.

For a density-dependent M3Y interaction, the corresponding Reid or Paris component is multiplied by the density- and energy-dependent factors:

$$
v^{X}(\rho,E,s)
=
g(E)F(\rho)v_0^{X}(s),
\qquad
X\in\{\mathrm D,\mathrm E\}.
$$

### Density-Dependent Interaction

The density-dependent effective interaction is written as

$$
v^{X}(\rho,E,s)
=
g(E)F(\rho)v_0^{X}(s),
\qquad
X\in\{\mathrm D,\mathrm E\},
$$

where $v_0^{X}(s)$ is the density-independent direct or exchange component of the underlying M3Y interaction. The density-dependent factor is parametrized as

$$
F(\rho)
=
C\left[
1+\alpha e^{-\beta\rho}-\gamma\rho^n
\right].
$$

The parameters $C$, $\alpha$, $\beta$, $\gamma$, and $n$ depend on the selected DDM3Y, BDM3Y, or CDM3Y parametrization and must be used with the corresponding M3Y-Reid or M3Y-Paris interaction.

When an explicit energy dependence is included,

$$
g(E)
=
1-\kappa\frac{E_{\mathrm{lab}}}{A_{\mathrm p}},
$$

where

$$
\kappa
=
\begin{cases}
0.002~\mathrm{MeV}^{-1},
&
\text{M3Y-Reid},
\\
0.003~\mathrm{MeV}^{-1},
&
\text{M3Y-Paris}.
\end{cases}
$$

For parametrizations without an explicit energy dependence, $g(E)=1$.

##### Frozen-Density Approximation

The density entering $F(\rho)$ is evaluated using the frozen-density approximation. The intrinsic densities of the projectile and target are assumed to remain unchanged during the collision, and the local overlap density is taken as their sum.

For the direct potential,

$$
\rho_{\mathrm{FDA}}^{\mathrm D}
(\mathbf r_{\mathrm t},\mathbf r_{\mathrm p})
=
\rho_{\mathrm t}(\mathbf r_{\mathrm t})
+
\rho_{\mathrm p}(r_{\mathrm p}).
$$

For the exchange potential, the density matrices are evaluated at the midpoint coordinates

$$
\mathbf u_{\mathrm t}
=
\mathbf r_{\mathrm t}
+
\frac{\mathbf s}{2},
\qquad
\mathbf u_{\mathrm p}
=
\mathbf r_{\mathrm p}
-
\frac{\mathbf s}{2},
$$

so the corresponding frozen density is

$$
\rho_{\mathrm{FDA}}^{\mathrm E}
(\mathbf u_{\mathrm t},\mathbf u_{\mathrm p})
=
\rho_{\mathrm t}(\mathbf u_{\mathrm t})
+
\rho_{\mathrm p}(u_{\mathrm p}).
$$

The density-dependent factor accounts phenomenologically for medium effects, particularly Pauli blocking and nuclear-matter saturation. It is distinct from the local-density approximation: the frozen-density approximation specifies the density entering the interaction, whereas the local-density approximation replaces the nonlocal one-body density matrices.

In the following algebra, $\rho_{\mathrm t}$ and $\rho_{\mathrm p}$ denote the target and projectile densities evaluated at the coordinates appropriate to the direct or exchange term defined above.

Under the frozen-density approximation,

$$
\rho_{\mathrm{FDA}}
=
\rho_{\mathrm t}
+
\rho_{\mathrm p},
$$

and the density-dependent factor becomes

$$
F(\rho_{\mathrm{FDA}})
=
C\left[
1
+
\alpha e^{-\beta\rho_{\mathrm t}}
e^{-\beta\rho_{\mathrm p}}
-
\gamma
(\rho_{\mathrm t}+\rho_{\mathrm p})^n
\right].
$$

For an integer $n$, the polynomial term can be expanded as

$$
(\rho_{\mathrm t}+\rho_{\mathrm p})^n
=
\sum_{k=0}^{n}
\binom{n}{k}
\rho_{\mathrm t}^{k}
\rho_{\mathrm p}^{n-k}.
$$

The density factor appearing in the folding integrand can therefore be written as

$$
\rho_{\mathrm t}\rho_{\mathrm p}
F(\rho_{\mathrm{FDA}})
=
C\rho_{\mathrm t}\rho_{\mathrm p}
+
C\alpha
\left[
\rho_{\mathrm t}e^{-\beta\rho_{\mathrm t}}
\right]
\left[
\rho_{\mathrm p}e^{-\beta\rho_{\mathrm p}}
\right]
-
C\gamma
\sum_{k=0}^{n}
\binom{n}{k}
\rho_{\mathrm t}^{k+1}
\rho_{\mathrm p}^{n-k+1}.
$$

Thus, the density-dependent contribution is reduced to a finite sum of separable products of target and projectile functions. Each product can be Fourier transformed independently and included as a separate momentum-space kernel. The coordinate-independent factor $g(E)$ may be applied after these spatial kernels are assembled.

For a deformed target, each nonlinear target function, such as $\rho_{\mathrm t}^{k+1}$ or $\rho_{\mathrm t}e^{-\beta\rho_{\mathrm t}}$, must be constructed in coordinate space before its spherical-harmonic decomposition. Nonlinear operations on the density generally mix different multipoles.

### Computational Workflow

The following workflow assumes a density-independent effective interaction.

In a density-dependent calculation, each separable term in $\rho_{\mathrm t}\rho_{\mathrm p}F(\rho_{\mathrm{FDA}})$ must pass through the same workflow independently, after which all contributions are summed to obtain the complete potential.

The calculation consists of two stages. The first stage prepares the momentum-space kernels. The second stage evaluates the direct and exchange potentials at a given $\mathbf R$.

##### Momentum-Space Kernels

Transform the direct interaction to momentum space:

$$
v^{\mathrm D}(s)
\longrightarrow
\tilde v^{\mathrm D}(q).
$$

For a spherical projectile and a deformed target, decompose the target density into spherical harmonics:

$$
\rho_{\mathrm p}(r),
\quad
\rho_{\mathrm t}(\mathbf r)
\longrightarrow
\rho_{\mathrm p}(r),
\quad
\rho_{\mathrm t,lm}(r).
$$

Transform the radial density components to momentum space:

$$
\rho_{\mathrm p}(r),
\quad
\rho_{\mathrm t,lm}(r)
\longrightarrow
\tilde\rho_{\mathrm p}(q),
\quad
\tilde\rho_{\mathrm t,lm}(q).
$$

Next, construct the local-density form factors:

$$
f_{\mathrm p}(r,s)
=
\rho_{\mathrm p}(r)
\hat j_1[k_{F,\mathrm p}(r)s],
\qquad
f_{\mathrm t}(\mathbf r,s)
=
\rho_{\mathrm t}(\mathbf r)
\hat j_1[k_{F,\mathrm t}(\mathbf r)s].
$$

Decompose the target form factor into spherical harmonics:

$$
f_{\mathrm p}(r,s),
\quad
f_{\mathrm t}(\mathbf r,s)
\longrightarrow
f_{\mathrm p}(r,s),
\quad
f_{\mathrm t,lm}(r,s).
$$

Transform the radial form factors to momentum space:

$$
f_{\mathrm p}(r,s),
\quad
f_{\mathrm t,lm}(r,s)
\longrightarrow
\tilde f_{\mathrm p}(q,s),
\quad
\tilde f_{\mathrm t,lm}(q,s).
$$

For each $(l,m)$ and $s$, construct the momentum-space exchange kernel:

$$
\tilde G_{lm}(q,s)
=
\tilde f_{\mathrm t,lm}(q,s)
\tilde f_{\mathrm p}(q,s).
$$

The momentum-space direct kernel is

$$
\tilde V^{\mathrm D}_{lm}(q)
=
\tilde\rho_{\mathrm t,lm}(q)
\tilde\rho_{\mathrm p}(q)
\tilde v^{\mathrm D}(q).
$$

At this point, both $\tilde V^{\mathrm D}_{lm}(q)$ and $\tilde G_{lm}(q,s)$ are ready for calculations at any $\mathbf R$.

##### Calculation at a Given Coordinate

For a selected $\mathbf R$, transform the direct kernels to the radial coordinate:

$$
\tilde V^{\mathrm D}_{lm}(q)
\longrightarrow
V^{\mathrm D}_{lm}(R).
$$

Reconstruct the direct potential at the specified orientation:

$$
V^{\mathrm D}(\mathbf R)
=
\sum_{l,m}
V^{\mathrm D}_{lm}(R)
Y_{lm}(\hat{\mathbf R}).
$$

Transform the exchange kernels to the radial coordinate:

$$
\tilde G_{lm}(q,s)
\longrightarrow
G_{lm}(R,s).
$$

Reconstruct the full exchange kernel at the specified orientation:

$$
G_{lm}(R,s)
\longrightarrow
G(\mathbf R,s)
=
\sum_{l,m}
G_{lm}(R,s)Y_{lm}(\hat{\mathbf R}).
$$

Starting from an initial $K(\mathbf R)$, calculate

$$
V^{\mathrm E}(\mathbf R)
=
4\pi
\int_0^\infty s^2\,ds\,
j_0\left[
\frac{K(\mathbf R)s}{M}
\right]
v^{\mathrm E}(s)
G(\mathbf R,s).
$$

Then update $K(\mathbf R)$ from the local-energy relation and repeat the calculation of $V^{\mathrm E}(\mathbf R)$ until convergence.
