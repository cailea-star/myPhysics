# Diffusion Equation

## Definition, Initial and Boundary Conditions

Let $\Omega\subset\mathbb R^3$ be the solution domain and $\partial\Omega$ its boundary. Let $u(\mathbf r,t)$ be the scalar field, $f(\mathbf r,t)$ the prescribed source, and $D>0$ a constant diffusion coefficient:

$$
\boxed{\frac{\partial u(\mathbf r,t)}{\partial t}-D\nabla^2u(\mathbf r,t)=f(\mathbf r,t),\qquad\mathbf r\in\Omega,\quad t>t_0}.
$$

Here $t_0$ is the initial time, and $D$ has units of length squared per time. In a source-free region,

$$
f=0\qquad\Rightarrow\qquad\frac{\partial u}{\partial t}=D\nabla^2u.
$$

### Initial Condition

Specify the initial field $u_0(\mathbf r)$:

$$
u(\mathbf r,t_0)=u_0(\mathbf r).
$$

No independent initial value of $\partial u/\partial t$ is required.

### Boundary Conditions on a Bounded Domain

Let $\hat{\mathbf n}$ be the outward unit normal. Define the normal derivative by

$$
\frac{\partial u}{\partial n}\equiv\hat{\mathbf n}\cdot\nabla u.
$$

Let $g_{\mathrm D},g_{\mathrm N},g_{\mathrm R}$ be prescribed boundary data, which may depend on time.

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

The prescribed coefficients $\alpha,\beta$ must not both vanish at any boundary point; all terms must have matching dimensions. Zero boundary data define homogeneous boundary conditions. For a classical solution, initial and boundary data must agree at $t=t_0$.

### Conditions at Infinity on an Unbounded Domain

For localized initial data and source, impose spatial decay at each fixed time:

$$
\lim_{|\mathbf r|\to\infty}u(\mathbf r,t)=0,\qquad t\geq t_0.
$$

For an exterior domain, conditions on the finite boundary must also be specified.

## Separation of Variables

First use homogeneous, time-independent boundary conditions. Nonzero boundary data can be handled by subtracting a field satisfying those data.

### Time Separation

For the source-free equation, take

$$
u_{\mathrm h}(\mathbf r,t)=U(\mathbf r)T(t).
$$

Substitution gives, with separation constant $\lambda$,

$$
\frac{1}{DT(t)}\frac{\mathrm dT(t)}{\mathrm dt}=\frac{\nabla^2U(\mathbf r)}{U(\mathbf r)}=-\lambda.
$$

The separated equations are

$$
\frac{\mathrm dT(t)}{\mathrm dt}+D\lambda T(t)=0,\qquad(\nabla^2+\lambda)U(\mathbf r)=0.
$$

Normalize $T(t_0)=1$. Then

$$
\boxed{T(t)=e^{-D\lambda(t-t_0)}}.
$$

For $\lambda=k^2\geq0$, with wave number $k$, the spatial equation is the [homogeneous Helmholtz equation](equation_helmholtz.md#separation-of-variables). Its Cartesian, cylindrical, and spherical solutions apply directly. Positive $\lambda$ gives decay; $\lambda=0$ gives a stationary mode.

### Eigenmode Expansion

On a bounded domain, assume the homogeneous boundary conditions make $-\nabla^2$ self-adjoint and nonnegative. Let $n,m$ label the modes, and let $U_n(\mathbf r)$ be the orthonormal eigenfunctions, with eigenvalues $\lambda_n\geq0$:

$$
-\nabla^2U_n(\mathbf r)=\lambda_nU_n(\mathbf r),\qquad\int_\Omega\mathrm d^3\mathbf r\,U_n^*(\mathbf r)U_m(\mathbf r)=\delta_{nm}.
$$

Expand the field with mode amplitudes $a_n(t)$:

$$
u(\mathbf r,t)=\sum_n a_n(t)U_n(\mathbf r).
$$

Define the initial amplitudes $c_n$ and source projections $f_n(t)$:

$$
c_n\equiv\int_\Omega\mathrm d^3\mathbf r\,U_n^*(\mathbf r)u_0(\mathbf r),\qquad f_n(t)\equiv\int_\Omega\mathrm d^3\mathbf r\,U_n^*(\mathbf r)f(\mathbf r,t).
$$

Each mode satisfies

$$
\frac{\mathrm da_n(t)}{\mathrm dt}+D\lambda_n a_n(t)=f_n(t),\qquad a_n(t_0)=c_n.
$$

Therefore,

$$
\boxed{a_n(t)=c_ne^{-D\lambda_n(t-t_0)}+\int_{t_0}^{t}\mathrm dt'\,e^{-D\lambda_n(t-t')}f_n(t')}.
$$

For $f=0$, the integral term vanishes. For homogeneous Neumann conditions on a connected domain, the constant spatial mode has $\lambda_n=0$ and does not decay.

## Free-Space Green Function

Let $(\mathbf r,t)$ be the observation point and $(\mathbf r',t')$ the source point. The retarded Green function satisfies

$$
\left(\frac{\partial}{\partial t}-D\nabla_{\mathbf r}^2\right)G_{\mathrm{ret}}(\mathbf r,\mathbf r';t,t')=\delta^{(3)}(\mathbf r-\mathbf r')\delta(t-t'),
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

Let $\omega\in\mathbb R$ be the Fourier angular frequency and $\widetilde G$ the corresponding frequency component. Spatial and temporal translation symmetry give

$$
G_{\mathrm{ret}}(\mathbf r,\mathbf r';t,t')=G_{\mathrm{ret}}(\mathbf R;\tau),\qquad\widetilde G(\mathbf r,\mathbf r';\omega)=\widetilde G(\mathbf R;\omega).
$$

Expand in the time difference:

$$
G_{\mathrm{ret}}(\mathbf R;\tau)=\int_{-\infty}^{\infty}\frac{\mathrm d\omega}{2\pi}\,\widetilde G(\mathbf R;\omega)e^{-i\omega\tau}.
$$

Substitution replaces the time derivative by $-i\omega$, giving

$$
\left(\nabla_{\mathbf R}^2+\frac{i\omega}{D}\right)\widetilde G(\mathbf R;\omega)=-\frac{1}{D}\delta^{(3)}(\mathbf R).
$$

This is the complex-coefficient form of the [Helmholtz Green-function equation](equation_helmholtz.md#free-space-green-function). Rotation symmetry makes the kernel depend only on $R$. Select the spatially decaying kernel:

$$
\boxed{\widetilde G(\mathbf R;\omega)=\frac{1}{4\pi DR}\exp\!\left[-R\sqrt{-\frac{i\omega}{D}}\right]}.
$$

For $\omega\neq0$, choose the square-root branch with positive real part; at $\omega=0$, take the continuous limit.

The inverse transform gives the heat kernel:

$$
\boxed{G_{\mathrm{ret}}(\mathbf R;\tau)=\frac{1}{4\pi DR}\int_{-\infty}^{\infty}\frac{\mathrm d\omega}{2\pi}\exp\!\left[-i\omega\tau-R\sqrt{-\frac{i\omega}{D}}\right]=\frac{1}{(4\pi D\tau)^{3/2}}\exp\!\left(-\frac{R^2}{4D\tau}\right),\qquad\tau>0}.
$$

For $\tau<0$, the inverse transform vanishes. As $\tau\to0$ through positive values, the heat kernel tends to $\delta^{(3)}(\mathbf R)$ in the distributional sense.

### Initial and Source Integrals

Propagating the initial field and accumulating the subsequent source contributions gives the free-space solution:

$$
\boxed{u(\mathbf r,t)=\int_{\mathbb R^3}\mathrm d^3\mathbf r'\,G_{\mathrm{ret}}(\mathbf r-\mathbf r';t-t_0)u_0(\mathbf r')+\int_{t_0}^{t}\mathrm dt'\int_{\mathbb R^3}\mathrm d^3\mathbf r'\,G_{\mathrm{ret}}(\mathbf r-\mathbf r';t-t')f(\mathbf r',t'),\qquad t>t_0}.
$$

The first term satisfies the prescribed initial condition; the second is the source response with zero initial field. For a finite boundary, use a kernel satisfying the corresponding boundary conditions or add a boundary correction.

## Spherical Expansion of the Green Function

For a fixed frequency $\omega>0$, let $U(\mathbf r)$ and $F(\mathbf r)$ be the complex field and source amplitudes:

$$
u(\mathbf r,t)=U(\mathbf r)e^{-i\omega t},\qquad f(\mathbf r,t)=F(\mathbf r)e^{-i\omega t}.
$$

The corresponding spatial equation is

$$
\left(\nabla^2+\frac{i\omega}{D}\right)U(\mathbf r)=-\frac{F(\mathbf r)}{D}.
$$

The spherical Green-function expansion, source multipoles, and exterior field are given in [Helmholtz Equation](equation_helmholtz.md#spherical-expansion-of-the-green-function). Continue its wave-number parameter to $i\sqrt{-i\omega/D}$ and replace the source by $F/D$ to obtain the frequency-domain diffusion expansion:

$$
\boxed{\widetilde G(\mathbf r,\mathbf r';\omega)=\frac{1}{D}G\!\left(\mathbf r,\mathbf r';i\sqrt{-\frac{i\omega}{D}}\right)}.
$$

Choose the square-root branch with positive real part so that the exterior radial solution decays at infinity. Zero frequency corresponds to the static limit of the Helmholtz kernel. Restore the single-frequency time dependence by multiplying by $e^{-i\omega t}$; general time dependence follows from the inverse Fourier transform.
