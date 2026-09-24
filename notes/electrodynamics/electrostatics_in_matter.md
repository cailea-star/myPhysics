# Electrostatics in Matter

### Electrostatic Potential and Interface Conditions

Let $\Phi$ be the electrostatic potential and $\epsilon(\mathbf r)$ the permittivity of a linear, isotropic dielectric. The fields are

$$
\mathbf E=-\nabla\Phi,\qquad \mathbf D=\epsilon(\mathbf r)\mathbf E=-\epsilon(\mathbf r)\nabla\Phi.
$$

In electrostatics, the relevant Maxwell equations reduce to

$$
\nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t}=0,\qquad \nabla\cdot\mathbf D=\rho_f.
$$

Substituting the field–potential relations gives

$$
\nabla\cdot\!\left[\epsilon(\mathbf r)\nabla\Phi\right]=-\rho_f.
$$

Here $\rho_f$ is the free volume charge density. At an interface, let $\hat{\mathbf n}$ point from medium 1 to medium 2, and let $\sigma_f$ be the free surface charge density. The field boundary conditions are

$$
\hat{\mathbf n}\times(\mathbf E_2-\mathbf E_1)=0,\qquad \hat{\mathbf n}\cdot(\mathbf D_2-\mathbf D_1)=\sigma_f.
$$

If $\mathbf E$ remains finite, its integral across the interface vanishes as the thickness $h\to0$:

$$
\Phi_2-\Phi_1=-\lim_{h\to0}\int_{-h/2}^{h/2}\mathbf E\cdot\hat{\mathbf n}\,\mathrm ds=0.
$$

Here $s$ is the normal coordinate. The normal-field condition then gives, with $\partial_{\hat{\mathbf n}}\Phi=\hat{\mathbf n}\cdot\nabla\Phi$ and $\epsilon_i$ the permittivity on side $i$,

$$
\boxed{\Phi_1=\Phi_2,\qquad \epsilon_2\partial_{\hat{\mathbf n}}\Phi_2-\epsilon_1\partial_{\hat{\mathbf n}}\Phi_1=-\sigma_f}.
$$

### Uniqueness Theorem

Let $V$ be a connected dielectric region with permittivity $\epsilon(\mathbf r)>0$. Two solutions with the same free-charge distribution and outer-boundary potential satisfy

$$
\nabla\cdot\mathbf D=\nabla\cdot\mathbf D'=\rho_f,\qquad \Phi|_{\partial V}=\Phi'|_{\partial V}=\Phi_{\mathrm{bound}}.
$$

Using the same constitutive and field–potential relations, define

$$
\delta\Phi=\Phi-\Phi',\qquad \delta\mathbf E=\mathbf E-\mathbf E'=-\nabla\delta\Phi,\qquad \delta\mathbf D=\mathbf D-\mathbf D'=\epsilon(\mathbf r)\delta\mathbf E.
$$

Subtracting the field equations and outer-boundary data gives

$$
\nabla\cdot\delta\mathbf D=0,\qquad \delta\Phi|_{\partial V}=0.
$$

If $V$ contains internal dielectric interfaces, the two solutions satisfy the same interface conditions, so the internal surface terms cancel. Applying the divergence theorem in each region gives

$$
0=\oint_{\partial V}\delta\Phi\,\delta\mathbf D\cdot\mathrm d\mathbf a=\int_V\nabla\cdot(\delta\Phi\,\delta\mathbf D)\,\mathrm dV=-\int_V\epsilon(\mathbf r)|\delta\mathbf E|^2\,\mathrm dV.
$$

Since $\epsilon(\mathbf r)>0$, the electric-field difference vanishes. The potential difference is constant in connected $V$ and is zero on $\partial V$:

$$
\boxed{\mathbf E=\mathbf E',\qquad \Phi=\Phi'}.
$$

For an application to a conducting sphere, see [Method of Images](electrostatics_of_conductors.md#method-of-images).

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

### Multipole Expansion

##### Far-Field Potential

Let $\rho(\mathbf r')$ be a localized charge density in vacuum, with $r=|\mathbf r|$, $r'=|\mathbf r'|$, and vacuum permittivity $\epsilon_0$. Its potential is

$$
\Phi(\mathbf r)=\frac{1}{4\pi\epsilon_0}\int\frac{\rho(\mathbf r')}{|\mathbf r-\mathbf r'|}\,\mathrm d^3r'.
$$

For $r\gg r'$, Taylor expansion in $\mathbf r'$ gives

$$
\frac{1}{|\mathbf r-\mathbf r'|}=\frac1r-\mathbf r'\cdot\nabla\left(\frac1r\right)+\frac12(\mathbf r'\otimes\mathbf r'):(\nabla\otimes\nabla)\left(\frac1r\right)+\cdots.
$$

Here $\nabla$ acts on $\mathbf r$; $\otimes$ and $:$ denote outer product and double contraction. With $\hat{\mathbf r}=\mathbf r/r$ and identity tensor $\hat{\mathbf I}$,

$$
\nabla\left(\frac1r\right)=-\frac{\mathbf r}{r^3}=-\frac{\hat{\mathbf r}}{r^2}.
$$

Differentiating this result once more and using $\nabla\otimes\mathbf r=\hat{\mathbf I}$ gives

$$
(\nabla\otimes\nabla)\left(\frac1r\right)=-\nabla\otimes\left(\frac{\mathbf r}{r^3}\right)=\frac{3\hat{\mathbf r}\otimes\hat{\mathbf r}-\hat{\mathbf I}}{r^3}.
$$

The angular tensor is traceless:

$$
\hat{\mathbf I}:(3\hat{\mathbf r}\otimes\hat{\mathbf r}-\hat{\mathbf I})=3(\hat{\mathbf r}\cdot\hat{\mathbf r})-\operatorname{tr}\hat{\mathbf I}=3-3=0.
$$

Let $\hat{\mathbf r}'=\mathbf r'/r'$ for $r'>0$. The quadratic term can therefore be written using a symmetric traceless source tensor:

$$
(\mathbf r'\otimes\mathbf r'):(3\hat{\mathbf r}\otimes\hat{\mathbf r}-\hat{\mathbf I})=\frac{r'^2}{3}(3\hat{\mathbf r}'\otimes\hat{\mathbf r}'-\hat{\mathbf I}):(3\hat{\mathbf r}\otimes\hat{\mathbf r}-\hat{\mathbf I}).
$$

Define the monopole charge, dipole moment, and quadrupole tensor by

$$
q=\int\rho(\mathbf r')\,\mathrm d^3r',\qquad \mathbf p=\int\rho(\mathbf r')\mathbf r'\,\mathrm d^3r',
$$

$$
\hat{\mathbf Q}^{(2)}=\int\rho(\mathbf r')r'^2(3\hat{\mathbf r}'\otimes\hat{\mathbf r}'-\hat{\mathbf I})\,\mathrm d^3r'.
$$

Substitution gives

$$
\Phi(\mathbf r)=\frac{1}{4\pi\epsilon_0}\left[\frac q r+\frac{\mathbf p\cdot\hat{\mathbf r}}{r^2}+\frac{\hat{\mathbf Q}^{(2)}:(3\hat{\mathbf r}\otimes\hat{\mathbf r}-\hat{\mathbf I})}{6r^3}+\cdots\right].
$$

Since $\hat{\mathbf Q}^{(2)}:\hat{\mathbf I}=0$,

$$
\boxed{\Phi(\mathbf r)=\frac{1}{4\pi\epsilon_0}\left[\frac q r+\frac{\mathbf p\cdot\hat{\mathbf r}}{r^2}+\frac{\hat{\mathbf r}\cdot\hat{\mathbf Q}^{(2)}\cdot\hat{\mathbf r}}{2r^3}+\cdots\right]}.
$$

##### Interaction with an External Field

Let $\Phi_{\mathrm{ext}}$ be a prescribed external potential that varies slowly across the localized charge distribution, with $\mathbf E_{\mathrm{ext}}=-\nabla\Phi_{\mathrm{ext}}$. The interaction energy is

$$
U_{\mathrm{int}}=\int\rho(\mathbf r')\Phi_{\mathrm{ext}}(\mathbf r')\,\mathrm d^3r'.
$$

Expanding the external potential about the origin gives

$$
\Phi_{\mathrm{ext}}(\mathbf r')=\Phi_{\mathrm{ext}}(0)+\mathbf r'\cdot\nabla\Phi_{\mathrm{ext}}(0)+\frac12(\mathbf r'\otimes\mathbf r'):(\nabla\otimes\nabla)\Phi_{\mathrm{ext}}(0)+\cdots.
$$

If the external sources lie outside the charge distribution, $\nabla^2\Phi_{\mathrm{ext}}=0$ there. Using the moments defined above,

$$
\boxed{U_{\mathrm{int}}=q\Phi_{\mathrm{ext}}(0)-\mathbf p\cdot\mathbf E_{\mathrm{ext}}(0)+\frac16\hat{\mathbf Q}^{(2)}:(\nabla\otimes\nabla)\Phi_{\mathrm{ext}}(0)+\cdots}.
$$

This is interaction energy with a prescribed external field, so it has no factor of $1/2$.

##### Dipole Force and Torque

Let a neutral, rigid dipole with moment $\mathbf p$ be centered at $\mathbf R$. Neglecting higher multipoles, its interaction energy is

$$
U_{\mathrm d}(\mathbf R)=-\mathbf p\cdot\mathbf E_{\mathrm{ext}}(\mathbf R).
$$

Hold $\mathbf p$ fixed when translating the dipole. Since the external field is electrostatic, $\nabla_{\mathbf R}\times\mathbf E_{\mathrm{ext}}=0$. With $\mathbf M$ denoting torque,

$$
\boxed{\mathbf F=-\nabla_{\mathbf R}U_{\mathrm d}=(\mathbf p\cdot\nabla_{\mathbf R})\mathbf E_{\mathrm{ext}}(\mathbf R),\qquad \mathbf M=\mathbf p\times\mathbf E_{\mathrm{ext}}(\mathbf R)}.
$$

A uniform field can exert a torque, but no net force.
