# Berry Phase

### Adiabatic Evolution and Berry Phase

Let $\mathbf R(t)$ be slowly varying external parameters and $\hat H(t)=\hat H(\mathbf R(t))$ a smooth Hamiltonian with discrete, nondegenerate instantaneous eigenstates:

$$
\hat H(t)|n(t)\rangle=E_n(t)|n(t)\rangle,\qquad\langle m(t)|n(t)\rangle=\delta_{mn},
$$

where $\delta_{mn}$ is the Kronecker delta. The state satisfies the [Schrödinger equation](quantum_dynamics.md#schrödinger-equation) and may be expanded as

$$
i\hbar\frac{d}{dt}|\psi(t)\rangle=\hat H(t)|\psi(t)\rangle,\qquad|\psi(t)\rangle=\sum_n a_n(t)|n(t)\rangle.
$$

With a dot denoting $d/dt$, projection onto $\langle m(t)|$ gives

$$
i\hbar\dot a_m(t)=E_m(t)a_m(t)-i\hbar\sum_n a_n(t)\langle m(t)|\dot n(t)\rangle,
$$

Separate the dynamical phase by defining

$$
a_m(t)=c_m(t)\exp\left[-\frac{i}{\hbar}\int_{t_i}^{t}E_m(t')\,dt'\right].
$$

Then

$$
\dot c_m(t)=-\sum_n c_n(t)\langle m(t)|\dot n(t)\rangle\exp\left\{\frac{i}{\hbar}\int_{t_i}^{t}[E_m(t')-E_n(t')]\,dt'\right\}.
$$

Separating the diagonal and off-diagonal terms gives

$$
\dot c_m(t)=-c_m(t)\langle m(t)|\dot m(t)\rangle-\sum_{n\neq m}c_n(t)\langle m(t)|\dot n(t)\rangle\exp\left\{\frac{i}{\hbar}\int_{t_i}^{t}[E_m(t')-E_n(t')]\,dt'\right\}.
$$

For $m\neq n$, differentiate the instantaneous eigenvalue equation:

$$
\hat H(t)|n(t)\rangle=E_n(t)|n(t)\rangle\quad\Longrightarrow\quad\dot{\hat H}(t)|n(t)\rangle+\hat H(t)|\dot n(t)\rangle=\dot E_n(t)|n(t)\rangle+E_n(t)|\dot n(t)\rangle.
$$

Projection onto $\langle m(t)|$ gives

$$
\langle m(t)|\dot{\hat H}(t)|n(t)\rangle+E_m(t)\langle m(t)|\dot n(t)\rangle=E_n(t)\langle m(t)|\dot n(t)\rangle,
$$

and therefore

$$
\boxed{\langle m(t)|\dot n(t)\rangle=\frac{\langle m(t)|\dot{\hat H}(t)|n(t)\rangle}{E_n(t)-E_m(t)}\qquad(m\neq n)}.
$$

A standard adiabaticity criterion is

$$
\boxed{\hbar\frac{|\langle m(t)|\dot{\hat H}(t)|n(t)\rangle|}{|E_n(t)-E_m(t)|^2}\ll1\qquad(m\neq n)}.
$$

The energy gap must remain nonzero throughout the evolution. Suppose now that

$$
|\psi(t_i)\rangle=|n(t_i)\rangle.
$$

Neglecting the off-diagonal transitions gives

$$
\dot c_n(t)=-\langle n(t)|\dot n(t)\rangle c_n(t),\qquad c_n(t)=\exp\left[-\int_{t_i}^{t}\langle n(t')|\dot n(t')\rangle\,dt'\right].
$$

Therefore,

$$
\boxed{|\psi_n(t)\rangle=\exp\left[-\frac{i}{\hbar}\int_{t_i}^{t}E_n(t')\,dt'\right]\exp\left[-\int_{t_i}^{t}\langle n(t')|\dot n(t')\rangle\,dt'\right]|n(t)\rangle}.
$$

Equivalently,

$$
|\psi_n(t)\rangle=e^{i\delta_n(t)}e^{i\gamma_n(t)}|n(t)\rangle,\qquad\delta_n(t)=-\frac{1}{\hbar}\int_{t_i}^{t}E_n(t')\,dt',\qquad\gamma_n(t)=i\int_{t_i}^{t}\langle n(t')|\dot n(t')\rangle\,dt'.
$$

Normalization implies

$$
\frac{d}{dt}\langle n(t)|n(t)\rangle=0\quad\Longrightarrow\quad\langle n(t)|\dot n(t)\rangle=-\langle\dot n(t)|n(t)\rangle.
$$

Hence, $\langle n|\dot n\rangle$ is purely imaginary and $\gamma_n(t)$ is real.

### Berry Connection, Curvature, and Gauge Structure

Let $C$ be the parameter-space path from $\mathbf R_i=\mathbf R(t_i)$ to $\mathbf R_f=\mathbf R(t_f)$. Define the Berry connection $\mathbf A_n(\mathbf R)$ by

$$
\boxed{\mathbf A_n(\mathbf R)=i\langle n(\mathbf R)|\nabla_{\mathbf R}n(\mathbf R)\rangle,\qquad\gamma_n(C)=i\int_{t_i}^{t_f}dt\,\langle n(t)|\dot n(t)\rangle=\int_C\mathbf A_n(\mathbf R)\cdot d\mathbf R}.
$$

In a three-dimensional parameter space, the Berry curvature is

$$
\boxed{\mathbf B_n(\mathbf R)=\nabla_{\mathbf R}\times\mathbf A_n(\mathbf R)}.
$$

For a closed path $C$ and an oriented surface $S$ with boundary $\partial S=C$, Stokes’ theorem gives

$$
\boxed{\gamma_n(C)=\oint_C\mathbf A_n(\mathbf R)\cdot d\mathbf R=\int_S\mathbf B_n(\mathbf R)\cdot d\mathbf S\pmod{2\pi}}.
$$

Under a smooth gauge transformation,

$$
|n(\mathbf R)\rangle\longrightarrow|n'(\mathbf R)\rangle=e^{i\xi_n(\mathbf R)}|n(\mathbf R)\rangle,\qquad\mathbf A_n'=\mathbf A_n-\nabla_{\mathbf R}\xi_n,\qquad\mathbf B_n'=\mathbf B_n.
$$

For an open path,

$$
\gamma_n'(C)=\gamma_n(C)+\xi_n(\mathbf R_i)-\xi_n(\mathbf R_f).
$$

For a closed path, the Berry phase is gauge invariant modulo $2\pi$:

$$
\boxed{\gamma_n'(C)=\gamma_n(C)\pmod{2\pi}}.
$$

Away from degeneracies,

$$
\langle m(\mathbf R)|\nabla_{\mathbf R}n(\mathbf R)\rangle=\frac{\langle m(\mathbf R)|\nabla_{\mathbf R}\hat H(\mathbf R)|n(\mathbf R)\rangle}{E_n(\mathbf R)-E_m(\mathbf R)},\qquad m\neq n.
$$

Therefore,

$$
\mathbf B_n(\mathbf R)=i\sum_{m\neq n}\frac{\langle n(\mathbf R)|\nabla_{\mathbf R}\hat H(\mathbf R)|m(\mathbf R)\rangle\times\langle m(\mathbf R)|\nabla_{\mathbf R}\hat H(\mathbf R)|n(\mathbf R)\rangle}{[E_n(\mathbf R)-E_m(\mathbf R)]^2}.
$$

Degeneracies are singular points of the Berry curvature.
