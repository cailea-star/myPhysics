# Electrostatics of Conductors

### Electrostatic Conditions for Conductors

Consider a conductor in vacuum at electrostatic equilibrium. Let $S$ be its surface, $\hat{\mathbf n}$ point from the conductor into vacuum, $\sigma$ and $\rho$ denote the surface and volume charge densities, and $\Phi_c$ be the conductor potential.

$$
\mathbf E=-\nabla\Phi,\qquad \mathbf D_{\mathrm{out}}=\epsilon_0\mathbf E_{\mathrm{out}}.
$$

The static Maxwell equations and surface boundary conditions are

$$
\nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t}=0,\qquad \nabla\cdot\mathbf E=\frac{\rho}{\epsilon_0}.
$$

$$
\hat{\mathbf n}\times(\mathbf E_{\mathrm{out}}-\mathbf E_{\mathrm{in}})=0,\qquad \hat{\mathbf n}\cdot(\mathbf E_{\mathrm{out}}-\mathbf E_{\mathrm{in}})=\frac{\sigma}{\epsilon_0}.
$$

The field vanishes inside the conductor, so

$$
\mathbf E_{\mathrm{in}}=0\quad\Longrightarrow\quad \Phi_{\mathrm{in}}=\Phi_c,\qquad \rho_{\mathrm{in}}=0,\qquad \mathbf E_{\mathrm{out}}=\frac{\sigma}{\epsilon_0}\hat{\mathbf n}.
$$

Here $\rho_{\mathrm{in}}=0$ is the net bulk charge density, not the density of mobile carriers. At a smooth surface, the potential is continuous:

$$
\Phi_{\mathrm{out}}-\Phi_{\mathrm{in}}=-\lim_{h\to0}\int_{-h/2}^{h/2}\mathbf E\cdot\hat{\mathbf n}\,\mathrm{d}s=0.
$$

With $\partial_{\hat{\mathbf n}}\Phi=\hat{\mathbf n}\cdot\nabla\Phi$,

$$
\boxed{\Phi_{\mathrm{out}}|_S=\Phi_c,\qquad \partial_{\hat{\mathbf n}}\Phi_{\mathrm{out}}|_S=-\frac{\sigma}{\epsilon_0}}.
$$

When the conductor potential is prescribed, the solution determines the surface charge density.

### Method of Images

An image charge lies outside the solution region. If a trial potential satisfies the boundary conditions, the [Uniqueness Theorem](electrostatics_in_matter.md#uniqueness-theorem) identifies it as the solution.

Consider a spherical conducting boundary of radius $R$ in vacuum, with permittivity $\epsilon_0$. A real point charge $q$ lies at $a\hat{\mathbf z}$, and an image charge $q_{\mathrm{im}}$ lies at $a_{\mathrm{im}}\hat{\mathbf z}$. With $r=|\mathbf r|$, the trial potential is

$$
\Phi(\mathbf r)=\frac{1}{4\pi\epsilon_0}\left(\frac{q}{|\mathbf r-a\hat{\mathbf z}|}+\frac{q_{\mathrm{im}}}{|\mathbf r-a_{\mathrm{im}}\hat{\mathbf z}|}\right).
$$

Grounding the sphere requires the potential to vanish at every point on its surface:

$$
|\mathbf r-a\hat{\mathbf z}|+\frac{q}{q_{\mathrm{im}}}|\mathbf r-a_{\mathrm{im}}\hat{\mathbf z}|=0\qquad(r=R).
$$

Let $\theta$ be the angle between $\mathbf r$ and the $z$ axis. Squaring this relation and matching the constant and $\cos\theta$ terms gives

$$
a_{\mathrm{im}}=\left(\frac{q_{\mathrm{im}}}{q}\right)^2a,\qquad R^2+a_{\mathrm{im}}^2=\left(\frac{q_{\mathrm{im}}}{q}\right)^2(R^2+a^2),
$$

$$
(a-a_{\mathrm{im}})(R^2-aa_{\mathrm{im}})=0.
$$

The root $a_{\mathrm{im}}=a$ cancels the real charge and is discarded. Hence

$$
\boxed{a_{\mathrm{im}}=\frac{R^2}{a},\qquad q_{\mathrm{im}}=-\frac{R}{a}q}.
$$

Substituting these parameters into the trial potential defines the grounded solution $\Phi_g$. Below, $Q$ denotes the conductor's net charge.

##### Charge Outside a Conducting Sphere

For $a>R$, the image lies inside the sphere, and $\Phi_g$ applies outside. For a conductor held at potential $V_0$, or an isolated conductor with net charge $Q$, add a central image charge:

$$
\Phi_{V_0}=\Phi_g+\frac{V_0R}{r},\qquad \Phi_Q=\Phi_g+\frac{Q-q_{\mathrm{im}}}{4\pi\epsilon_0r}\quad(r>R).
$$

The center lies outside the solution region, so this term introduces no new source there.

##### Charge Inside a Conducting Shell

For $0<a<R$, the image lies outside the cavity. For a thin conducting spherical shell of radius $R$, $\Phi_g$ applies inside the cavity, while the exterior potential is zero when the shell is grounded. If the shell is held at $V_0$, add a constant inside:

$$
\Phi_{\mathrm{in}}=\Phi_g+V_0,\qquad \Phi_{\mathrm{out}}=\frac{V_0R}{r}.
$$

For an isolated shell with net charge $Q$, the inner and outer surfaces carry total charges $-q$ and $Q+q$, respectively:

$$
\Phi_{\mathrm{in}}=\Phi_g+\frac{Q+q}{4\pi\epsilon_0R},\qquad \Phi_{\mathrm{out}}=\frac{Q+q}{4\pi\epsilon_0r}.
$$

A central $1/r$ term would introduce an extra singularity inside the cavity.

### Green Reciprocity

Let $\Phi$ and $\Phi'$ describe two electrostatic configurations of the same $N$ conductors, with surfaces $S_i$, in a homogeneous medium of permittivity $\epsilon$. In the charge-free exterior $V$,

$$
\nabla^2\Phi=\nabla^2\Phi'=0,\qquad \Phi|_{S_i}=\Phi_i,\qquad \Phi'|_{S_i}=\Phi_i'.
$$

Both potentials vanish at infinity. Green's second identity gives

$$
0=\int_V(\Phi\nabla^2\Phi'-\Phi'\nabla^2\Phi)\,\mathrm{d}V=\oint_{\partial V}(\Phi\nabla\Phi'-\Phi'\nabla\Phi)\cdot\mathrm{d}\mathbf a.
$$

The contribution at infinity vanishes. On the boundary $\partial V_i$ of conductor $i$, the area element points outward from $V$ and into the conductor, so

$$
Q_i=\epsilon\oint_{\partial V_i}\nabla\Phi\cdot\mathrm{d}\mathbf a,\qquad Q_i'=\epsilon\oint_{\partial V_i}\nabla\Phi'\cdot\mathrm{d}\mathbf a.
$$

Since each conductor is equipotential, the remaining boundary integral becomes

$$
\oint_{\partial V}(\Phi\nabla\Phi'-\Phi'\nabla\Phi)\cdot\mathrm{d}\mathbf a=\sum_{i=1}^{N}\oint_{\partial V_i}(\Phi_i\nabla\Phi'-\Phi_i'\nabla\Phi)\cdot\mathrm{d}\mathbf a=\frac{1}{\epsilon}\sum_{i=1}^{N}(\Phi_iQ_i'-\Phi_i'Q_i).
$$

Therefore, Green reciprocity reads

$$
\boxed{\sum_{i=1}^{N}Q_i\Phi_i'=\sum_{i=1}^{N}Q_i'\Phi_i}.
$$

### Electrostatic Energy

Let $V$ be the charge-free exterior of $N$ conductors in a homogeneous linear dielectric of permittivity $\epsilon$.

$$
W=\frac{1}{2}\int_V\mathbf E\cdot\mathbf D\,\mathrm{d}V.
$$

Using $\mathbf E=-\nabla\Phi$ and $\nabla\cdot\mathbf D=0$ in $V$,

$$
\mathbf E\cdot\mathbf D=-\nabla\cdot(\Phi\mathbf D).
$$

The flux at infinity vanishes. On $\partial V_i$, $\Phi=\Phi_i$ and $\oint_{\partial V_i}\mathbf D\cdot\mathrm{d}\mathbf a=-Q_i$, hence

$$
W=-\frac{1}{2}\oint_{\partial V}\Phi\mathbf D\cdot\mathrm{d}\mathbf a=\frac{1}{2}\sum_{i=1}^{N}Q_i\Phi_i.
$$

Let $C_{ij}$ be the capacitance coefficients. [Green Reciprocity](#green-reciprocity) gives

$$
Q_i=\sum_{j=1}^{N}C_{ij}\Phi_j,\qquad C_{ij}=C_{ji}.
$$

$C_{ij}$ is the charge on conductor $i$ when conductor $j$ is held at unit potential and all other conductors at zero potential. Substituting $Q_i$ into the energy expression gives

$$
\boxed{W=\frac{1}{2}\sum_{i=1}^{N}Q_i\Phi_i=\frac{1}{2}\sum_{i,j=1}^{N}C_{ij}\Phi_i\Phi_j}.
$$

For two charge distributions with negligible rearrangement, let $\mathbf E_a,\varphi_a,\rho_a,W_a$ be their fields, potentials, charge densities, and individual energies ($a=1,2$). Then

$$
\boxed{W=W_1+W_2+W_{\mathrm{int}},\qquad W_{\mathrm{int}}=\epsilon\int_{\mathbb R^3}\mathbf E_1\cdot\mathbf E_2\,\mathrm{d}V=\int_{\mathbb R^3}\rho_1\varphi_2\,\mathrm{d}V=\int_{\mathbb R^3}\rho_2\varphi_1\,\mathrm{d}V}.
$$

### Electrostatic Stability

##### Thomson's Theorem

Keep conductor positions and total charges $Q_i$ fixed. Let $\delta\rho_i$ be a charge redistribution on conductor $i$, and let $V_i$ contain that conductor:

$$
\delta Q_i=\int_{V_i}\delta\rho_i\,\mathrm{d}V=0.
$$

The first-order energy variation is

$$
\delta W=\sum_i\int_{V_i}\Phi(\mathbf r)\,\delta\rho_i(\mathbf r)\,\mathrm{d}V.
$$

With Lagrange multipliers $\lambda_i$ for the fixed-charge constraints,

$$
0=\delta W-\sum_i\lambda_i\delta Q_i=\sum_i\int_{V_i}[\Phi(\mathbf r)-\lambda_i]\delta\rho_i(\mathbf r)\,\mathrm{d}V.
$$

Thus each conductor must be equipotential:

$$
\boxed{\Phi|_{V_i}=\lambda_i=\mathrm{constant}}.
$$

The positive quadratic field-energy term makes this stationary state an energy minimum.

##### Earnshaw's Theorem

Assume:

- Quasistatic motion: each conductor remains in electrostatic equilibrium.
- Large separation: motion causes negligible charge redistribution between conductors.

Under these assumptions, the self-energies remain effectively constant; only $W_{\mathrm{int}}$ changes with conductor positions.

Let $q_\alpha$ and $\mathbf r_\alpha$ be the charge and position of a small body, and let $\varphi_\alpha$ be the potential produced by all *other* bodies at $\mathbf r_\alpha$. Then

$$
W_{\mathrm{int}}=\frac{1}{2}\sum_\alpha q_\alpha\varphi_\alpha(\mathbf r_\alpha),\qquad \nabla_\alpha W_{\mathrm{int}}=q_\alpha\nabla_\alpha\varphi_\alpha=0
$$

at equilibrium. For a displacement $\delta\mathbf r_\alpha$, define $\mathcal H_{ij}^{(\alpha)}=q_\alpha\partial_{r_{\alpha i}}\partial_{r_{\alpha j}}\varphi_\alpha$. In a charge-free neighborhood,

$$
\delta^2W_{\mathrm{int}}=\frac{1}{2}\sum_{i,j}\mathcal H_{ij}^{(\alpha)}\delta r_{\alpha i}\delta r_{\alpha j},\qquad \operatorname{tr}\mathcal H^{(\alpha)}=q_\alpha\nabla_\alpha^2\varphi_\alpha=0.
$$

Because the symmetric Hessian has zero trace, it is either zero or has both positive and negative eigenvalues. The former gives no quadratic restoring force; the latter makes the equilibrium a saddle point. Neither permits a strict quadratic minimum.

### Force on a Conductor Surface

Let $\hat{\mathbf n}$ point outward from a conductor in vacuum, and let $\mathbf f_S$ be the force per unit area on its surface. With charge density $\rho$ and current density $\mathbf j$, [momentum conservation](electromagnetic_conservation_laws.md#momentum-conservation) gives

$$
\frac{\partial\mathbf g}{\partial t}+\nabla\cdot\hat{\mathbf T}=-(\rho\mathbf E+\mathbf j\times\mathbf B),\qquad \mathbf g=\epsilon_0\mathbf E\times\mathbf B.
$$

Here $\mathbf g$ is the field momentum density and $\hat{\mathbf T}$ is the outward momentum-flux tensor. For a stationary conductor in pure electrostatics, $\mathbf j=0$ and $\mathbf B=0$; hence $\mathbf g=0$ and $\nabla\cdot\hat{\mathbf T}=-\rho\mathbf E$.

##### Maxwell Stress Tensor

Let $\hat{\mathbf I}$ be the identity operator and $\otimes$ denote the outer product. The electrostatic tensor is

$$
\hat{\mathbf T}=\epsilon_0\left(\frac{\mathbf E\cdot\mathbf E}{2}\hat{\mathbf I}-\mathbf E\otimes\mathbf E\right).
$$

At the conductor surface, $\mathbf E_{\mathrm{in}}=0$ and $\mathbf E_{\mathrm{out}}=E_\perp\hat{\mathbf n}=(\sigma/\epsilon_0)\hat{\mathbf n}$, where $\sigma$ is the surface charge density. Thus

$$
\boxed{\mathbf f_S=-\hat{\mathbf T}\cdot\hat{\mathbf n}=\frac{\epsilon_0E_\perp^2}{2}\hat{\mathbf n}=\frac{\sigma^2}{2\epsilon_0}\hat{\mathbf n}}.
$$

##### Direct Force on the Surface Charge

Model the surface charge as a layer of thickness $l$ occupying $-l<z<0$, with uniform volume charge density $\rho_0$ and $z$ increasing along $\hat{\mathbf n}$. Within the layer,

$$
\mathbf E(z)=\frac{\rho_0(z+l)}{\epsilon_0}\hat{\mathbf n}.
$$

Integrating the $\rho\mathbf E$ force and taking $\sigma=\lim_{l\to0}\rho_0l$ gives

$$
\mathbf f_S=\int_{-l}^{0}\rho_0\mathbf E(z)\,\mathrm{d}z=\frac{(\rho_0l)^2}{2\epsilon_0}\hat{\mathbf n}\longrightarrow\frac{\sigma^2}{2\epsilon_0}\hat{\mathbf n}.
$$
