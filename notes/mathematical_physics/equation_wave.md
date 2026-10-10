# Wave Equation

## Definition and Initial–Boundary Conditions

Let $\Omega\subset\mathbb R^3$ be the spatial domain, $u(\mathbf r,t)$ the scalar field, $f(\mathbf r,t)$ the prescribed source, and $v>0$ the constant wave speed. Use the sign convention

$$
\boxed{\left(\nabla^2-\frac{1}{v^2}\frac{\partial^2}{\partial t^2}\right)u(\mathbf r,t)=-f(\mathbf r,t),\qquad\mathbf r\in\Omega}.
$$

In a source-free region,

$$
f=0\qquad\Rightarrow\qquad\frac{\partial^2u}{\partial t^2}-v^2\nabla^2u=0.
$$

### Initial Conditions

At an initial time $t_0$, prescribe the field $u_0(\mathbf r)$ and its time derivative $u_1(\mathbf r)$:

$$
u(\mathbf r,t_0)=u_0(\mathbf r),\qquad\frac{\partial u}{\partial t}(\mathbf r,t_0)=u_1(\mathbf r).
$$

Two initial conditions are required because the equation is second order in time.

### Boundary Conditions on a Bounded Domain

Let $\partial\Omega$ be the boundary and $\hat{\mathbf n}$ its outward unit normal. Define

$$
\frac{\partial u}{\partial n}\equiv\hat{\mathbf n}\cdot\nabla u.
$$

Let $g_{\mathrm D},g_{\mathrm N},g_{\mathrm R}$ be prescribed boundary data, possibly depending on time.

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

The prescribed coefficients $\alpha,\beta$ must not both vanish at any boundary point; all terms must have matching dimensions. Initial and boundary data must be compatible at $t=t_0$.

### Conditions at Infinity on an Unbounded Domain

For localized initial data and sources, impose spatial decay at each finite time:

$$
\lim_{|\mathbf r|\to\infty}u(\mathbf r,t)=0.
$$

For a purely source-generated retarded solution, also exclude incoming radiation from infinity. Spatial decay alone does not select the retarded solution. For an exterior domain, conditions on the finite boundary must also be specified.

## Separation of Variables

Decompose the solution into a particular solution $u_{\mathrm p}$ and a homogeneous solution $u_{\mathrm h}$:

$$
u=u_{\mathrm p}+u_{\mathrm h},\qquad\left(\nabla^2-\frac{1}{v^2}\frac{\partial^2}{\partial t^2}\right)u_{\mathrm p}=-f,\qquad\left(\nabla^2-\frac{1}{v^2}\frac{\partial^2}{\partial t^2}\right)u_{\mathrm h}=0.
$$

### Time–Space Separation

For time-independent homogeneous boundary conditions, take a spatial factor $U(\mathbf r)$ and a temporal factor $T(t)$:

$$
u_{\mathrm h}(\mathbf r,t)=U(\mathbf r)T(t).
$$

For oscillatory modes, introduce a real separation parameter $k\geq0$:

$$
\frac{\nabla^2U(\mathbf r)}{U(\mathbf r)}=\frac{1}{v^2T(t)}\frac{\mathrm d^2T(t)}{\mathrm dt^2}=-k^2.
$$

The separated equations are

$$
(\nabla^2+k^2)U(\mathbf r)=0,\qquad\frac{\mathrm d^2T(t)}{\mathrm dt^2}+v^2k^2T(t)=0.
$$

The spatial solutions in Cartesian, cylindrical, and spherical coordinates are given in [Helmholtz Equation](equation_helmholtz.md#separation-of-variables).

For $k>0$, define the angular frequency $\omega=vk$. With complex constants $A,B$,

$$
\boxed{T(t)=Ae^{-i\omega(t-t_0)}+Be^{i\omega(t-t_0)},\qquad\omega=vk}.
$$

For a real temporal factor, $B=A^*$. For $k=0$, with constants $A_0,B_0$,

$$
T(t)=A_0+B_0(t-t_0).
$$

Initial data determine the temporal coefficients after expansion in the spatial modes.

### Frequency-Domain Equation

For a fixed angular frequency $\omega>0$, use the complex representation with field and source amplitudes $U(\mathbf r)$ and $F(\mathbf r)$:

$$
u(\mathbf r,t)=U(\mathbf r)e^{-i\omega t},\qquad f(\mathbf r,t)=F(\mathbf r)e^{-i\omega t}.
$$

Substitution gives the inhomogeneous Helmholtz equation:

$$
\boxed{(\nabla^2+k^2)U(\mathbf r)=-F(\mathbf r),\qquad k=\frac{\omega}{v}}.
$$

General time-dependent fields are constructed by superposing frequency components.

## Free-Space Green Function

Let $(\mathbf r,t)$ be the observation point and $(\mathbf r',t')$ the source point. The retarded Green function satisfies

$$
\left(\nabla_{\mathbf r}^2-\frac{1}{v^2}\frac{\partial^2}{\partial t^2}\right)G_{\mathrm{ret}}(\mathbf r,\mathbf r';t,t')=-\delta^{(3)}(\mathbf r-\mathbf r')\delta(t-t'),
$$

with no response before the source event:

$$
G_{\mathrm{ret}}(\mathbf r,\mathbf r';t,t')=0,\qquad t<t'.
$$

### Fundamental Solution

Define the relative coordinate, its magnitude, and the time difference:

$$
\mathbf R\equiv\mathbf r-\mathbf r',\qquad R\equiv|\mathbf R|,\qquad\tau\equiv t-t'.
$$

Let $\omega\in\mathbb R$ be the Fourier frequency and $\widetilde G$ the corresponding frequency component. Spatial and temporal translation symmetry give

$$
G_{\mathrm{ret}}(\mathbf r,\mathbf r';t,t')=G_{\mathrm{ret}}(\mathbf R;\tau),\qquad\widetilde G(\mathbf r,\mathbf r';\omega)=\widetilde G(\mathbf R;\omega).
$$

Expand in the time difference:

$$
G_{\mathrm{ret}}(\mathbf R;\tau)=\int_{-\infty}^{\infty}\frac{\mathrm d\omega}{2\pi}\,\widetilde G(\mathbf R;\omega)e^{-i\omega\tau}.
$$

Substitution replaces the second time derivative by $-\omega^2$, giving

$$
(\nabla_{\mathbf R}^2+k^2)\widetilde G(\mathbf R;\omega)=-\delta^{(3)}(\mathbf R),\qquad k=\frac{|\omega|}{v}.
$$

This is the [Helmholtz Green-function equation](equation_helmholtz.md#free-space-green-function). Select the outgoing kernel at positive frequencies and its complex conjugate at negative frequencies. Rotation symmetry makes the kernel depend only on $R$:

$$
\widetilde G(\mathbf R;\omega)=\frac{e^{i\omega R/v}}{4\pi R}.
$$

Substitution into the time expansion gives

$$
\boxed{G_{\mathrm{ret}}(\mathbf R;\tau)=\frac{1}{4\pi R}\int_{-\infty}^{\infty}\frac{\mathrm d\omega}{2\pi}e^{-i\omega(\tau-R/v)}=\frac{\delta(\tau-R/v)}{4\pi R}}.
$$

The response occurs only at $\tau=R/v$: the disturbance propagates from the source point to the observation point at speed $v$.

### Source Integral

The source-generated retarded particular solution is

$$
u_{\mathrm p}(\mathbf r,t)=\int_{-\infty}^{t}\mathrm dt'\int_{\mathbb R^3}\mathrm d^3\mathbf r'\,G_{\mathrm{ret}}(\mathbf r-\mathbf r';t-t')f(\mathbf r',t').
$$

The delta function fixes $t'=t-R/v$, giving

$$
\boxed{u_{\mathrm p}(\mathbf r,t)=\frac{1}{4\pi}\int_{\mathbb R^3}\mathrm d^3\mathbf r'\,\frac{f\!\left(\mathbf r',t-\frac{|\mathbf r-\mathbf r'|}{v}\right)}{|\mathbf r-\mathbf r'|}}.
$$

The field at time $t$ depends on each source point at its retarded time $t-R/v$. To satisfy prescribed initial data, add the corresponding homogeneous solution.

## Spherical Expansion of the Green Function

For a fixed frequency $\omega>0$, use the complex representation

$$
u(\mathbf r,t)=U(\mathbf r)e^{-i\omega t},\qquad k=\frac{\omega}{v}.
$$

The spherical Green-function expansion, source multipoles, exterior field, and static limit are given in [Helmholtz Equation](equation_helmholtz.md#spherical-expansion-of-the-green-function). Restore the time dependence by multiplying the spatial amplitude by $e^{-i\omega t}$.
