# Axial Gaussian Matrix Elements

The Gaussian matrix elements are evaluated in the [axial harmonic-oscillator basis](axial_basis.md), whose factors follow the canonical [coordinate-space harmonic-oscillator solutions](../../notes/quantum_mechanics/harmonic_oscillator_coordinate_space.md).

### Gaussian Kernel Separation

For $\mathbf r_a=(\boldsymbol\rho_a,z_a)$ with $\boldsymbol\rho_a=(r_a\cos\varphi_a,r_a\sin\varphi_a)$, define

$$
g_\mu(\mathbf r_1,\mathbf r_2)=\exp\left[-\frac{|\mathbf r_1-\mathbf r_2|^2}{\mu^2}\right],\qquad \mu>0.
$$

The squared distance separates as

$$
|\mathbf r_1-\mathbf r_2|^2=(z_1-z_2)^2+|\boldsymbol\rho_1-\boldsymbol\rho_2|^2.
$$

The transverse part is

$$
|\boldsymbol\rho_1-\boldsymbol\rho_2|^2=r_1^2+r_2^2-2r_1r_2\cos(\varphi_1-\varphi_2).
$$

Therefore,

$$
g_\mu(\mathbf r_1,\mathbf r_2)=g_\mu^z(z_1,z_2)g_\mu^r(\boldsymbol\rho_1,\boldsymbol\rho_2).
$$

The axial factor is

$$
g_\mu^z(z_1,z_2)=\exp\left[-\frac{(z_1-z_2)^2}{\mu^2}\right].
$$

The transverse factor is

$$
g_\mu^r(\boldsymbol\rho_1,\boldsymbol\rho_2)=\exp\left[-\frac{|\boldsymbol\rho_1-\boldsymbol\rho_2|^2}{\mu^2}\right].
$$

For the spatial label

$$
a=(n_z,n_r,\Lambda),\qquad n_z,n_r\in\mathbb N_0,\qquad \Lambda\in\mathbb Z,
$$

define

$$
\phi_a(\mathbf r)=\phi_{n_z}(z;b_z)\phi_{n_r}^{\Lambda}(r,\varphi;b_r).
$$

Using $d^3r=dz\,r\,dr\,d\varphi$, the Gaussian matrix element is

$$
G_{12,34}^{(\mu)}=\int d^3r_1d^3r_2\,\phi_1^*(\mathbf r_1)\phi_2^*(\mathbf r_2)g_\mu(\mathbf r_1,\mathbf r_2)\phi_3(\mathbf r_1)\phi_4(\mathbf r_2).
$$

Denoting its axial and transverse factors by $G^z$ and $G^r$,

$$
\boxed{G_{12,34}^{(\mu)}=G_{n_{z1}n_{z2}n_{z3}n_{z4}}^z(\mu,b_z)\,G_{n_{r1}\Lambda_1,n_{r2}\Lambda_2,n_{r3}\Lambda_3,n_{r4}\Lambda_4}^r(\mu,b_r).}
$$

### One-Dimensional Gaussian Matrix Elements

For an oscillator length $b>0$ and $n\in\mathbb N_0$, define

$$
\boxed{N_n=\frac{1}{\sqrt{2^nn!\sqrt\pi}},\qquad \phi_n(x;b)=\frac{N_n}{\sqrt b}H_n\left(\frac{x}{b}\right)e^{-x^2/(2b^2)}.}
$$

The one-dimensional Gaussian matrix element is

$$
G_{n_1n_2n_3n_4}^{1D}(\mu,b)=\int dx_1dx_2\,\phi_{n_1}(x_1;b)\phi_{n_2}(x_2;b)e^{-(x_1-x_2)^2/\mu^2}\phi_{n_3}(x_1;b)\phi_{n_4}(x_2;b).
$$

Define the Talmi-Moshinsky coefficient

$$
T(n_1,n_2;n)=\frac{\sqrt{n_1!n_2!n!}}{\left(\frac{-n_1+n_2+n}{2}\right)!\left(\frac{n_1-n_2+n}{2}\right)!\left(\frac{n_1+n_2-n}{2}\right)!}.
$$

Introduce

$$
\chi=1+\frac{\mu^2}{2b^2},\qquad \xi=\frac{p+q+n+1}{2}.
$$

The hypergeometric factor is

$$
\bar F(p,q,n;\chi)=\frac{\Gamma(\xi-p)\Gamma(\xi-q)\Gamma(\xi-n)}{\chi^\xi\sqrt{p!q!n!}}\,{}_2F_1(-p,-q;n+1-\xi;1-\chi).
$$

Since $-p$ and $-q$ are nonpositive integers, the hypergeometric series terminates:

$$
{}_2F_1(-p,-q;c;x)=\sum_{i=0}^{\min(p,q)}\frac{(-p)_i(-q)_i}{(c)_i\,i!}x^i.
$$

Here $(a)_i$ denotes the Pochhammer symbol. The matrix element is

$$
\boxed{G_{n_1n_2n_3n_4}^{1D}(\mu,b)=\frac{\mu}{\sqrt{2\pi^3}\,b}\sum_nT(n_2,n_4;n)\bar F(n_1,n_3,n;\chi).}
$$

The axial factor follows directly:

$$
\boxed{G_{n_{z1}n_{z2}n_{z3}n_{z4}}^z(\mu,b_z)=G_{n_{z1}n_{z2}n_{z3}n_{z4}}^{1D}(\mu,b_z).}
$$

### Polar-to-Cartesian Transformation

For $n_r\in\mathbb N_0$ and $\Lambda\in\mathbb Z$, define the two-dimensional polar oscillator state

$$
\phi_{n_r}^{\Lambda}(r,\varphi;b_r)=\frac{e^{i\Lambda\varphi}}{\sqrt{2\pi}}\phi_{n_r}^{|\Lambda|}(r;b_r),\qquad \nu=2n_r+|\Lambda|.
$$

It expands in the one-dimensional Cartesian basis as

$$
\boxed{\phi_{n_r}^{\Lambda}(x,y;b_r)=\sum_{n_y=0}^{\nu}i^{n_y}C_{n_y}^{n_r\Lambda}\phi_{n_x}(x;b_r)\phi_{n_y}(y;b_r),\qquad n_x=\nu-n_y.}
$$

The transformation coefficient is

$$
C_{n_y}^{n_r\Lambda}=(-1)^{n_r}2^{-n_r-|\Lambda|/2}\sqrt{\frac{(n_r+|\Lambda|)!n_r!}{n_x!n_y!}}\sum_q\binom{n_x}{k_x}\binom{n_y}{q}(-1)^{n_y-q},
$$

where

$$
k_x=n_r-q+\frac{|\Lambda|-\Lambda}{2}.
$$

The admissible values of $q$ are specified with the remaining selection rules.

### Radial Gaussian Matrix Elements

Let $\boldsymbol\rho_a=(x_a,y_a)$. The radial Gaussian matrix element is

$$
G_{12,34}^r(\mu,b_r)=\int d^2\rho_1d^2\rho_2\,\left[\phi_{n_{r1}}^{\Lambda_1}(\boldsymbol\rho_1;b_r)\right]^*\left[\phi_{n_{r2}}^{\Lambda_2}(\boldsymbol\rho_2;b_r)\right]^*e^{-|\boldsymbol\rho_1-\boldsymbol\rho_2|^2/\mu^2}\phi_{n_{r3}}^{\Lambda_3}(\boldsymbol\rho_1;b_r)\phi_{n_{r4}}^{\Lambda_4}(\boldsymbol\rho_2;b_r).
$$

The radial kernel separates in Cartesian coordinates:

$$
e^{-|\boldsymbol\rho_1-\boldsymbol\rho_2|^2/\mu^2}=e^{-(x_1-x_2)^2/\mu^2}e^{-(y_1-y_2)^2/\mu^2}.
$$

For each state $a=1,\ldots,4$, define

$$
n_{x,a}=\nu_a-n_{y,a},\qquad P_y=i^{-n_{y,1}-n_{y,2}+n_{y,3}+n_{y,4}}.
$$

Substituting the Cartesian expansions gives

$$
\boxed{G_{12,34}^r(\mu,b_r)=\sum_{n_{y,1}=0}^{\nu_1}\cdots\sum_{n_{y,4}=0}^{\nu_4}P_y\prod_{a=1}^{4}C_{n_{y,a}}^{n_{r,a}\Lambda_a}G_{n_{x,1}n_{x,2}n_{x,3}n_{x,4}}^{1D}(\mu,b_r)G_{n_{y,1}n_{y,2}n_{y,3}n_{y,4}}^{1D}(\mu,b_r).}
$$

### Symmetries and Selection Rules

- General one-dimensional parity:

  $$
  n_1+n_2+n_3+n_4\equiv0\pmod2.
  $$

- Axial parity:

  $$
  n_{z1}+n_{z2}+n_{z3}+n_{z4}\equiv0\pmod2.
  $$

- Conservation of the angular-momentum projection:

  $$
  \Lambda_1+\Lambda_2=\Lambda_3+\Lambda_4.
  $$

- Particle-exchange symmetry:

  $$
  G_{12,34}^{(\mu)}=G_{21,43}^{(\mu)}.
  $$

- Hermitian symmetry:

  $$
  G_{12,34}^{(\mu)}=\left[G_{34,12}^{(\mu)}\right]^*.
  $$

- Time-reversal symmetry, with $\bar a=(n_{z,a},n_{r,a},-\Lambda_a)$:

  $$
  G_{\bar1\bar2,\bar3\bar4}^{(\mu)}=\left[G_{12,34}^{(\mu)}\right]^*.
  $$

For the real Gaussian kernel, the retained matrix elements are real. Therefore,

$$
\boxed{G_{12,34}^{(\mu)}=G_{21,43}^{(\mu)}=G_{34,12}^{(\mu)}=G_{43,21}^{(\mu)},\qquad G_{\bar1\bar2,\bar3\bar4}^{(\mu)}=G_{12,34}^{(\mu)}.}
$$

Violation of any applicable selection rule gives a zero matrix element.
