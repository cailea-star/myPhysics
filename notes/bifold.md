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

