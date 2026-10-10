# Helmholtz Equation

## Definition and Boundary Conditions

Let $\Omega\subset\mathbb R^3$ be the solution domain and $\partial\Omega$ its boundary. Let $u(\mathbf r)$ be the scalar field, $f(\mathbf r)$ the prescribed source, and $k\geq0$ a real wave number. Use the sign convention

$$
\boxed{(\nabla^2+k^2)u(\mathbf r)=-f(\mathbf r),\qquad\mathbf r\in\Omega}.
$$

In a source-free region,

$$
f=0\qquad\Rightarrow\qquad(\nabla^2+k^2)u=0.
$$

For $k=0$, these reduce to the Poisson and Laplace equations, respectively.

### Boundary Conditions on a Bounded Domain

Let $\hat{\mathbf n}$ be the outward unit normal. Define the normal derivative by

$$
\frac{\partial u}{\partial n}\equiv\hat{\mathbf n}\cdot\nabla u.
$$

Let $g_{\mathrm D},g_{\mathrm N},g_{\mathrm R}$ be prescribed boundary data.

**Dirichlet boundary condition (first kind).** Specify the boundary value:

$$
u\big|_{\partial\Omega}=g_{\mathrm D}.
$$

**Neumann boundary condition (second kind).** Specify the normal derivative:

$$
\frac{\partial u}{\partial n}\bigg|_{\partial\Omega}=g_{\mathrm N}.
$$

**Robin boundary condition (third kind).** Specify a linear combination:

$$
\left(\alpha u+\beta\frac{\partial u}{\partial n}\right)\bigg|_{\partial\Omega}=g_{\mathrm R}.
$$

The prescribed coefficients $\alpha,\beta$ must not both vanish at any boundary point; all terms must have matching dimensions. Zero boundary data define homogeneous boundary conditions.

### Conditions at Infinity on an Unbounded Domain

For a localized source in three-dimensional free space, let $r=|\mathbf r|$.

For $k>0$, the outgoing solution satisfies the Sommerfeld radiation condition:

$$
\boxed{u=O(r^{-1}),\qquad\lim_{r\to\infty}r\left(\frac{\partial u}{\partial r}-iku\right)=0}.
$$

The limit holds uniformly in direction. Decay alone does not exclude incoming waves.

For $k=0$, impose the static decay condition:

$$
\lim_{r\to\infty}u(\mathbf r)=0.
$$

For an exterior domain, conditions on the finite boundary must also be specified.

## Separation of Variables

Decompose the solution into a particular solution $u_{\mathrm p}$ and a homogeneous solution $u_{\mathrm h}$:

$$
u=u_{\mathrm p}+u_{\mathrm h},\qquad(\nabla^2+k^2)u_{\mathrm p}=-f,\qquad(\nabla^2+k^2)u_{\mathrm h}=0.
$$

Choose coordinates adapted to the boundary geometry. Substitute a product ansatz, separate the coordinate-dependent terms, and superpose the separated solutions to satisfy the boundary conditions. The integration constants below are fixed by these conditions.

### Cartesian Coordinates

In Cartesian coordinates $(x,y,z)$, the Laplace operator is

$$
\nabla^2=\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}+\frac{\partial^2}{\partial z^2}.
$$

Let $k_x,k_y,k_z$ be separation parameters and take

$$
u_{\mathrm h}(x,y,z)=X(k_x;x)Y(k_y;y)Z(k_z;z).
$$

The separated equations are

$$
\frac{\mathrm d^2}{\mathrm dx^2}X(k_x;x)+k_x^2X(k_x;x)=0,\qquad\frac{\mathrm d^2}{\mathrm dy^2}Y(k_y;y)+k_y^2Y(k_y;y)=0,\qquad\frac{\mathrm d^2}{\mathrm dz^2}Z(k_z;z)+k_z^2Z(k_z;z)=0.
$$

With $\mathbf k=(k_x,k_y,k_z)$ and constant amplitude $C$, the exponential modes are

$$
\boxed{u_{\mathrm h}=Ce^{i\mathbf k\cdot\mathbf r},\qquad k_x^2+k_y^2+k_z^2=k^2}.
$$

The separation parameters may be complex. A zero parameter also admits a linear solution of the corresponding one-dimensional equation.

### Cylindrical Coordinates

In [cylindrical coordinates](coordinate_cylindrical.md) $(r_\perp,\phi,z)$, the Laplace operator is

$$
\nabla^2=\frac{1}{r_\perp}\frac{\partial}{\partial r_\perp}\left(r_\perp\frac{\partial}{\partial r_\perp}\right)+\frac{1}{r_\perp^2}\frac{\partial^2}{\partial\phi^2}+\frac{\partial^2}{\partial z^2}.
$$

For real separation parameter $k_z$ and angular index $m\in\mathbb Z$, take the product ansatz

$$
u_{\mathrm h}(r_\perp,\phi,z)=R_m(k_z;r_\perp)\Phi_m(\phi)Z(k_z;z).
$$

The separated equations are

$$
\frac{\mathrm d^2}{\mathrm dz^2}Z(k_z;z)+k_z^2Z(k_z;z)=0.
$$

$$
\frac{\mathrm d^2}{\mathrm d\phi^2}\Phi_m(\phi)+m^2\Phi_m(\phi)=0.
$$

$$
\frac{\mathrm d^2}{\mathrm dr_\perp^2}R_m(k_z;r_\perp)+\frac{1}{r_\perp}\frac{\mathrm d}{\mathrm dr_\perp}R_m(k_z;r_\perp)+\left(k^2-k_z^2-\frac{m^2}{r_\perp^2}\right)R_m(k_z;r_\perp)=0.
$$

Choose normalized angular Fourier modes:

$$
Z(k_z;z)=e^{ik_zz},\qquad\Phi_m(\phi)=\frac{1}{\sqrt{2\pi}}e^{im\phi},\qquad m\in\mathbb Z.
$$

The axial zero mode also admits

$$
Z(0;z)=A_z+B_zz.
$$

For $k_z^2<k^2$, let $k_\perp=\sqrt{k^2-k_z^2}$. Using the Bessel functions $J_{|m|},Y_{|m|}$,

$$
R_m(k_z;r_\perp)=A_mJ_{|m|}(k_\perp r_\perp)+B_mY_{|m|}(k_\perp r_\perp).
$$

For $k_z^2>k^2$, use the modified Bessel functions $I_{|m|},K_{|m|}$:

$$
R_m(k_z;r_\perp)=A_mI_{|m|}\!\left(\sqrt{k_z^2-k^2}\,r_\perp\right)+B_mK_{|m|}\!\left(\sqrt{k_z^2-k^2}\,r_\perp\right).
$$

For $k_z^2=k^2$, with reference length $r_0>0$,

$$
R_m(k_z;r_\perp)=A_mr_\perp^{|m|}+B_mr_\perp^{-|m|}\quad(m\neq0),\qquad R_0(k_z;r_\perp)=A_0+B_0\ln\frac{r_\perp}{r_0}.
$$

### Spherical Coordinates

In [spherical coordinates](coordinate_spherical.md) $(r,\theta,\phi)$, the Laplace operator is

$$
\nabla^2=\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)+\frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{r^2\sin^2\theta}\frac{\partial^2}{\partial\phi^2}.
$$

Using [scalar spherical harmonics](specialfunction_spherical_harmonics.md#scalar-spherical-harmonics), take

$$
u_{\mathrm h}(\mathbf r)=R_{lm}(r)Y_{lm}(\hat{\mathbf r}),\qquad l=0,1,\ldots,\qquad m=-l,\ldots,l.
$$

The angular equation is

$$
\left[\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{\sin^2\theta}\frac{\partial^2}{\partial\phi^2}\right]Y_{lm}(\hat{\mathbf r})=-l(l+1)Y_{lm}(\hat{\mathbf r}).
$$

The radial equation is

$$
\frac{1}{r^2}\frac{\mathrm d}{\mathrm dr}\left[r^2\frac{\mathrm d}{\mathrm dr}R_{lm}(r)\right]+\left[k^2-\frac{l(l+1)}{r^2}\right]R_{lm}(r)=0.
$$

This equation is independent of $m$; write its radial solutions as $R_l(r)$.

For $k>0$, use the spherical Bessel functions $j_l,y_l$:

$$
\boxed{R_l(r)=A_lj_l(kr)+B_ly_l(kr)}.
$$

Regularity at the origin selects $j_l$. Outgoing and incoming radial waves use the spherical Hankel functions

$$
h_l^{(1)}(kr)=j_l(kr)+iy_l(kr),\qquad h_l^{(2)}(kr)=j_l(kr)-iy_l(kr),
$$

respectively.

For $k=0$, the radial equation gives

$$
\boxed{R_l(r)=A_lr^l+B_lr^{-l-1}}.
$$

Regularity at the origin selects $r^l$; decay at infinity selects $r^{-l-1}$.

## Free-Space Green Function

Use the definition of [Green functions](green_functions.md#definition). Let $\mathbf r$ be the observation point and $\mathbf r'$ the source point. The free-space Green function satisfies

$$
(\nabla_{\mathbf r}^2+k^2)G(\mathbf r,\mathbf r';k)=-\delta^{(3)}(\mathbf r-\mathbf r'),
$$

with the outgoing radiation condition for $k>0$.

### Fundamental Solution

At fixed $k$, define the relative coordinate and its magnitude:

$$
\mathbf R\equiv\mathbf r-\mathbf r',\qquad R\equiv|\mathbf R|,\qquad G(\mathbf r,\mathbf r';k)=G(\mathbf R;k).
$$

Translation symmetry permits this representation; rotation symmetry makes the kernel depend only on $R$. Radial derivatives below are taken at fixed direction. For $R>0$,

$$
\frac{1}{R^2}\frac{\partial}{\partial R}\left(R^2\frac{\partial G(\mathbf R;k)}{\partial R}\right)+k^2G(\mathbf R;k)=0.
$$

For $k>0$, integration gives, with constants $A,B$,

$$
G(\mathbf R;k)=\frac{Ae^{ikR}+Be^{-ikR}}{R}.
$$

The outgoing condition selects $B=0$. Let $B_\epsilon$ be a ball of radius $\epsilon$ centered at $\mathbf r'$, and $\mathrm dS$ its boundary-area element. Integrating the defining equation gives

$$
-1=\lim_{\epsilon\to0}\left[\int_{\partial B_\epsilon}\frac{\partial G(\mathbf R;k)}{\partial R}\,\mathrm dS+k^2\int_{B_\epsilon}G(\mathbf R;k)\,\mathrm d^3\mathbf r\right]=-4\pi A.
$$

Therefore,

$$
\boxed{G(\mathbf R;k)=\frac{e^{ikR}}{4\pi R}}.
$$

For $k=0$, this reduces to the Poisson Green function.

### Source Integral

For a localized source, the outgoing free-space solution is

$$
\boxed{u(\mathbf r)=\int_{\mathbb R^3}G(\mathbf r-\mathbf r';k)f(\mathbf r')\,\mathrm d^3\mathbf r'=\frac{1}{4\pi}\int_{\mathbb R^3}\frac{e^{ik|\mathbf r-\mathbf r'|}}{|\mathbf r-\mathbf r'|}f(\mathbf r')\,\mathrm d^3\mathbf r'}.
$$

Indeed,

$$
(\nabla^2+k^2)u(\mathbf r)=-\int_{\mathbb R^3}\delta^{(3)}(\mathbf r-\mathbf r')f(\mathbf r')\,\mathrm d^3\mathbf r'=-f(\mathbf r).
$$

In a bounded domain $\Omega$, the free-space kernel gives a particular solution:

$$
u_{\mathrm p}(\mathbf r)=\int_\Omega G(\mathbf r-\mathbf r';k)f(\mathbf r')\,\mathrm d^3\mathbf r'.
$$

The prescribed boundary conditions are satisfied by adding a homogeneous solution:

$$
u=u_{\mathrm p}+u_{\mathrm h},\qquad(\nabla^2+k^2)u_{\mathrm h}=0.
$$

## Spherical Expansion of the Green Function

Let $r=|\mathbf r|$, $r'=|\mathbf r'|$, and define

$$
r_<\equiv\min(r,r'),\qquad r_>\equiv\max(r,r').
$$

Using [scalar spherical harmonics](specialfunction_spherical_harmonics.md#scalar-spherical-harmonics), define the source coefficients

$$
f_{lm}(r)\equiv\int\mathrm d\hat{\mathbf r}\,Y_{lm}^*(\hat{\mathbf r})f(\mathbf r),\qquad f(\mathbf r)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}f_{lm}(r)Y_{lm}(\hat{\mathbf r}).
$$

The kernel expansions below hold for $r_< < r_>$; equal radii can be treated by an appropriate limiting procedure.

### Spherical-Wave Expansion

For $k>0$, the outgoing kernel is

$$
\boxed{G(\mathbf r,\mathbf r';k)=\frac{e^{ik|\mathbf r-\mathbf r'|}}{4\pi|\mathbf r-\mathbf r'|}=ik\sum_{l=0}^{\infty}\sum_{m=-l}^{l}j_l(kr_<)h_l^{(1)}(kr_>)Y_{lm}(\hat{\mathbf r})Y_{lm}^*(\hat{\mathbf r}')}.
$$

The radial factors select the origin-regular solution $j_l$ and the outgoing solution $h_l^{(1)}$.

Substituting into the source integral and splitting at $r'=r$ gives

$$
\boxed{u(\mathbf r)=ik\sum_{l=0}^{\infty}\sum_{m=-l}^{l}Y_{lm}(\hat{\mathbf r})\left[h_l^{(1)}(kr)\int_0^r\mathrm dr'\,r'^2j_l(kr')f_{lm}(r')+j_l(kr)\int_r^\infty\mathrm dr'\,r'^2h_l^{(1)}(kr')f_{lm}(r')\right]}.
$$

Let the source be confined to a sphere of radius $a>0$. Define its wave-number-dependent multipole moments:

$$
Q_{lm}(k)\equiv\int_{\mathbb R^3}\mathrm d^3\mathbf r'\,j_l(kr')Y_{lm}^*(\hat{\mathbf r}')f(\mathbf r')=\int_0^a\mathrm dr'\,r'^2j_l(kr')f_{lm}(r').
$$

For $r>a$,

$$
\boxed{u(\mathbf r)=ik\sum_{l=0}^{\infty}\sum_{m=-l}^{l}Q_{lm}(k)h_l^{(1)}(kr)Y_{lm}(\hat{\mathbf r})}.
$$

For fixed $l$ and $kr\to\infty$,

$$
h_l^{(1)}(kr)\sim(-i)^{l+1}\frac{e^{ikr}}{kr},\qquad u(\mathbf r)\sim\frac{e^{ikr}}{r}\sum_{l=0}^{\infty}\sum_{m=-l}^{l}(-i)^lQ_{lm}(k)Y_{lm}(\hat{\mathbf r}).
$$

All outgoing multipoles have leading far-field radial decay $r^{-1}$.

### Static Limit

At fixed radii,

$$
\lim_{k\to0^+}ik\,j_l(kr_<)h_l^{(1)}(kr_>)=\frac{1}{2l+1}\frac{r_<^l}{r_>^{l+1}}.
$$

The static kernel is therefore

$$
\boxed{G(\mathbf r,\mathbf r';0)=\frac{1}{4\pi|\mathbf r-\mathbf r'|}=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\frac{1}{2l+1}\frac{r_<^l}{r_>^{l+1}}Y_{lm}(\hat{\mathbf r})Y_{lm}^*(\hat{\mathbf r}')}.
$$

The source integral becomes

$$
\boxed{u(\mathbf r)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\frac{Y_{lm}(\hat{\mathbf r})}{2l+1}\left[\frac{1}{r^{l+1}}\int_0^r\mathrm dr'\,r'^2r'^lf_{lm}(r')+r^l\int_r^\infty\mathrm dr'\,\frac{r'^2}{r'^{l+1}}f_{lm}(r')\right]}.
$$

Define the static multipole moments:

$$
Q_{lm}\equiv\int_{\mathbb R^3}\mathrm d^3\mathbf r'\,r'^lY_{lm}^*(\hat{\mathbf r}')f(\mathbf r')=\int_0^a\mathrm dr'\,r'^2r'^lf_{lm}(r').
$$

For $r>a$, the static exterior field is

$$
\boxed{u(\mathbf r)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\frac{Q_{lm}}{2l+1}\frac{Y_{lm}(\hat{\mathbf r})}{r^{l+1}}}.
$$

The orders $l=0,1,2,\ldots$ correspond to monopole, dipole, quadrupole, and higher components, with radial decay $r^{-1},r^{-2},r^{-3},\ldots$.
