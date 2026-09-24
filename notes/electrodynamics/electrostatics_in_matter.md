# Electrostatics in Matter

### Electrostatic Potential and Interface Conditions

Let $\Phi$ be the electrostatic potential, $\rho_f$ the free volume charge density, and $\epsilon(\mathbf r)$ the permittivity of a linear, isotropic dielectric. Using $\mathbf E=-\nabla\Phi$, $\mathbf D=\epsilon\mathbf E$, and $\nabla\cdot\mathbf D=\rho_f$ gives

$$
\boxed{\nabla\cdot[\epsilon(\mathbf r)\nabla\Phi]=-\rho_f(\mathbf r)}.
$$

Within a homogeneous region, this becomes Poisson’s equation, or Laplace’s equation where $\rho_f=0$:

$$
\nabla^2\Phi=-\frac{\rho_f}{\epsilon},\qquad \rho_f=0\Longrightarrow\nabla^2\Phi=0.
$$

At an interface, let $\hat{\mathbf n}$ point from medium 1 to medium 2; $\epsilon_i$ and $\Phi_i$ denote the limiting permittivity and potential on side $i$, and $\sigma_f$ is the free surface charge density. As established in [Boundary Conditions](maxwell_equations.md#boundary-conditions), $\hat{\mathbf n}\cdot(\mathbf D_2-\mathbf D_1)=\sigma_f$.

Let $s$ measure distance along $\hat{\mathbf n}$ across a short segment of length $h$, and set $E_\perp=\mathbf E\cdot\hat{\mathbf n}$. If the one-sided electric fields remain finite, then

$$
\Phi_2-\Phi_1=-\lim_{h\to0}\int_{-h/2}^{h/2}E_\perp(s)\,\mathrm{d}s=0.
$$

Thus, defining $\partial_{\hat{\mathbf n}}\Phi=\hat{\mathbf n}\cdot\nabla\Phi$, the potential conditions are

$$
\boxed{\Phi_1=\Phi_2,\qquad \epsilon_2\partial_{\hat{\mathbf n}}\Phi_2-\epsilon_1\partial_{\hat{\mathbf n}}\Phi_1=-\sigma_f}.
$$

### Uniqueness Theorem

Let $V$ be a connected dielectric region with fixed free-charge distribution and permittivity $\epsilon(\mathbf r)>0$. For two electrostatic solutions, define

$$
\delta\Phi=\Phi'-\Phi'',\qquad \delta\mathbf E=\mathbf E'-\mathbf E''=-\nabla\delta\Phi,\qquad \delta\mathbf D=\mathbf D'-\mathbf D''=\epsilon(\mathbf r)\delta\mathbf E.
$$

Both solutions have the same $\rho_f$, so $\nabla\cdot\delta\mathbf D=0$. If they satisfy the same prescribed potential or normal displacement on each part of the boundary $\partial V$, the boundary integral below vanishes. Contributions from internal dielectric interfaces cancel because $\delta\Phi$ and the normal component of $\delta\mathbf D$ are continuous:

$$
0=\oint_{\partial V}\delta\Phi\,\delta\mathbf D\cdot\mathrm{d}\mathbf a=\int_V\nabla\cdot(\delta\Phi\,\delta\mathbf D)\,\mathrm{d}V=-\int_V\epsilon(\mathbf r)\,\delta\mathbf E\cdot\delta\mathbf E\,\mathrm{d}V.
$$

Since $\epsilon(\mathbf r)>0$,

$$
\boxed{\mathbf E'=\mathbf E'',\qquad \Phi'-\Phi''=\mathrm{constant}}.
$$

A prescribed boundary potential also fixes that constant.

### Method of Images

An image charge is placed outside the region where the potential is sought. Its potential has no source in that region; if the total potential satisfies the boundary conditions, the [Uniqueness Theorem](#uniqueness-theorem) identifies it as the solution.

##### Charge Outside a Conducting Sphere

Consider a conducting sphere of radius $R$ in vacuum. A point charge $q$ lies at $d\hat{\mathbf z}$, where $d>R$; a trial image charge $q_{\mathrm{im}}$ lies at $b\hat{\mathbf z}$ inside the sphere. For $\mathbf r$ outside the sphere, define $r=|\mathbf r|$, $r_q=|\mathbf r-d\hat{\mathbf z}|$, and $r_{\mathrm{im}}=|\mathbf r-b\hat{\mathbf z}|$. Requiring zero potential on $r=R$ for every angle gives

$$
\boxed{b=\frac{R^2}{d},\qquad q_{\mathrm{im}}=-\frac{R}{d}q,\qquad \Phi_g(\mathbf r)=\frac{1}{4\pi\epsilon_0}\left(\frac{q}{r_q}+\frac{q_{\mathrm{im}}}{r_{\mathrm{im}}}\right)\quad(r>R)}.
$$

Here $\epsilon_0$ is the vacuum permittivity and $\Phi_g$ is the grounded-sphere potential. For a sphere held at potential $V_0$, or an isolated sphere with total charge $Q$, add a charge at its center:

$$
\Phi_{V_0}=\Phi_g+\frac{V_0R}{r},\qquad \Phi_Q=\Phi_g+\frac{Q-q_{\mathrm{im}}}{4\pi\epsilon_0r}\quad(r>R).
$$

##### Charge Inside a Conducting Shell

For $0<d<R$, the same $b=R^2/d$ and $q_{\mathrm{im}}=-qR/d$ put the image *outside* the cavity. Let $\Phi_g$ denote the same charge-plus-image expression, now evaluated inside the cavity ($r<R$). For a grounded shell, the potential outside is zero. If the shell is held at $V_0$, add a **constant** inside—not a central image charge, which would lie in the region being solved:

$$
\Phi_{\mathrm{in}}=\Phi_g+V_0,\qquad \Phi_{\mathrm{out}}=\frac{V_0R}{r}.
$$

For an isolated shell with net charge $Q$, the inner and outer surfaces carry total charges $-q$ and $Q+q$, respectively. Therefore,

$$
\Phi_{\mathrm{in}}=\Phi_g+\frac{Q+q}{4\pi\epsilon_0R},\qquad \Phi_{\mathrm{out}}=\frac{Q+q}{4\pi\epsilon_0r}.
$$

### Separation of Variables

In a homogeneous, charge-free region, $\nabla^2\Phi=0$. Choose separated solutions that match the geometry, then determine their coefficients from the boundary conditions.

##### Spherical Coordinates

For an axially symmetric spherical problem, let $r$ and $\theta$ be the radial and polar coordinates. With Legendre polynomials $P_\ell$ and coefficients $a_\ell,b_\ell$,

$$
\Phi(r,\theta)=\sum_{\ell=0}^{\infty}\left(a_\ell r^\ell+b_\ell r^{-(\ell+1)}\right)P_\ell(\cos\theta).
$$

The $r^\ell$ terms are regular at the origin; the $r^{-(\ell+1)}$ terms decay at infinity.

##### Cylindrical Coordinates

For a problem independent of the cylindrical coordinate $z$, let $\varrho$ and $\phi$ be the radial and azimuthal coordinates. With coefficients $c_m,d_m,e_m,f_m$,

$$
\Phi(\varrho,\phi)=c_0+d_0\ln\varrho+\sum_{m=1}^{\infty}\left[(c_m\varrho^m+d_m\varrho^{-m})\cos(m\phi)+(e_m\varrho^m+f_m\varrho^{-m})\sin(m\phi)\right].
$$

##### Conducting Sphere in a Uniform Field

Let a grounded conducting sphere of radius $R$ be placed in vacuum under a uniform field $\mathbf E_0=E_0\hat{\mathbf z}$. Since $z=r\cos\theta$, the applied-field potential is

$$
\Phi_{\mathrm{ext}}=-\mathbf E_0\cdot\mathbf r=-E_0z=-E_0r\cos\theta.
$$

The boundary data contain only the $\ell=1$ mode, with $P_1(\cos\theta)=\cos\theta$. Write the exterior solution as

$$
\Phi_{\mathrm{out}}=\left(A_1r+\frac{B_1}{r^2}\right)\cos\theta.
$$

Matching the applied field at infinity and setting the sphere’s potential to zero determine the coefficients:

$$
\Phi_{\mathrm{out}}-\Phi_{\mathrm{ext}}\to0\ (r\to\infty)\Longrightarrow A_1=-E_0,\qquad \Phi_{\mathrm{out}}(R,\theta)=0\Longrightarrow B_1=E_0R^3.
$$

Thus the exterior potential and induced dipole moment $\mathbf p$ are

$$
\boxed{\Phi_{\mathrm{out}}=-E_0\left(r-\frac{R^3}{r^2}\right)\cos\theta,\qquad \mathbf p=4\pi\epsilon_0R^3\mathbf E_0}.
$$

##### Dielectric Sphere in a Uniform Field

Let a dielectric sphere of radius $R$ and permittivity $\epsilon_{\mathrm{in}}$ lie in a medium of permittivity $\epsilon_{\mathrm{out}}$, under $\mathbf E_0=E_0\hat{\mathbf z}$. Since $z=r\cos\theta$, the applied-field potential is

$$
\Phi_{\mathrm{ext}}=-\mathbf E_0\cdot\mathbf r=-E_0z=-E_0r\cos\theta.
$$

It contains only the $\ell=1$ mode, so regularity at the origin and decay of the induced potential at infinity give

$$
\Phi_{\mathrm{out}}=\left(-E_0r+\frac{B}{r^2}\right)\cos\theta,\qquad \Phi_{\mathrm{in}}=Ar\cos\theta,
$$

where $A$ and $B$ are coefficients. With no free surface charge, the [interface conditions](#electrostatic-potential-and-interface-conditions) at $r=R$ are

$$
\Phi_{\mathrm{out}}=\Phi_{\mathrm{in}},\qquad \epsilon_{\mathrm{out}}\partial_r\Phi_{\mathrm{out}}=\epsilon_{\mathrm{in}}\partial_r\Phi_{\mathrm{in}}.
$$

Applying these conditions gives

$$
A=-\frac{3\epsilon_{\mathrm{out}}}{\epsilon_{\mathrm{in}}+2\epsilon_{\mathrm{out}}}E_0,\qquad B=\frac{\epsilon_{\mathrm{in}}-\epsilon_{\mathrm{out}}}{\epsilon_{\mathrm{in}}+2\epsilon_{\mathrm{out}}}E_0R^3.
$$

The potentials are

$$
\boxed{\Phi_{\mathrm{out}}=-E_0r\cos\theta+\frac{\epsilon_{\mathrm{in}}-\epsilon_{\mathrm{out}}}{\epsilon_{\mathrm{in}}+2\epsilon_{\mathrm{out}}}\frac{E_0R^3}{r^2}\cos\theta,\qquad \Phi_{\mathrm{in}}=-\frac{3\epsilon_{\mathrm{out}}}{\epsilon_{\mathrm{in}}+2\epsilon_{\mathrm{out}}}E_0r\cos\theta}.
$$

The interior field is uniform:

$$
\mathbf E_{\mathrm{in}}=\frac{3\epsilon_{\mathrm{out}}}{\epsilon_{\mathrm{in}}+2\epsilon_{\mathrm{out}}}\mathbf E_0.
$$
