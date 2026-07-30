# Quantum Dynamics

### Time-Evolution Operator

The Hilbert-space and operator conventions used below are defined in [Quantum Formalism](formalism.md).

##### Schrödinger Equation

Let $\mathcal H$ be the Hilbert space of a closed quantum system, $\hbar$ the reduced Planck constant, and $\hat H(t)$ a self-adjoint Hamiltonian on $\mathcal H$:

$$
\boxed{i\hbar\frac{\partial}{\partial t}|\psi(t)\rangle=\hat H(t)|\psi(t)\rangle}.
$$

For two solutions governed by the same Hamiltonian,

$$
\frac{d}{dt}\langle\phi(t)|\psi(t)\rangle=0\quad\Longrightarrow\quad\langle\phi(t)|\psi(t)\rangle=\langle\phi(t_0)|\psi(t_0)\rangle.
$$

Thus, inner products, norms, orthogonality, and transition probabilities are preserved.

##### Definition and Composition

Define the time-evolution operator $\hat U(t,t_0)$ and the identity operator $\hat I$ by

$$
\boxed{|\psi(t)\rangle=\hat U(t,t_0)|\psi(t_0)\rangle}.
$$

Its identity and composition laws are

$$
\boxed{\hat U(t_0,t_0)=\hat I,\qquad\hat U(t_3,t_2)\hat U(t_2,t_1)=\hat U(t_3,t_1)}.
$$

They imply

$$
\boxed{\hat U(t_1,t_2)=\hat U^{-1}(t_2,t_1)=\hat U^\dagger(t_2,t_1)}.
$$

Thus, $\hat U(t_2,t_1)$ is a [unitary operator](formalism.md#unitary-transformation). Its inverse evolves the state backward along the same dynamics and is not, in general, the time-reversal symmetry operator.

Substitution into the Schrödinger equation gives

$$
\boxed{i\hbar\frac{\partial}{\partial t}\hat U(t,t_0)=\hat H(t)\hat U(t,t_0)}.
$$

##### Dyson Series and Time Ordering

For $t\geq t_0$, integration gives

$$
\hat U(t,t_0)=\hat I-\frac{i}{\hbar}\int_{t_0}^{t}dt_1\,\hat H(t_1)\hat U(t_1,t_0).
$$

Repeated substitution yields the Dyson series:

$$
\hat U(t,t_0)=\hat I+\sum_{n=1}^{\infty}\left(-\frac{i}{\hbar}\right)^n\int_{t_0}^{t}dt_1\int_{t_0}^{t_1}dt_2\cdots\int_{t_0}^{t_{n-1}}dt_n\,\hat H(t_1)\hat H(t_2)\cdots\hat H(t_n).
$$

The nested limits impose time ordering. Let $\mathcal T$ denote the time-ordering operator and $\Theta(s)$ the Heaviside step function:

$$
\mathcal T[\hat H(t_1)\hat H(t_2)]=\Theta(t_1-t_2)\hat H(t_1)\hat H(t_2)+\Theta(t_2-t_1)\hat H(t_2)\hat H(t_1).
$$

Therefore,

$$
\boxed{\hat U(t,t_0)=\mathcal T\exp\left[-\frac{i}{\hbar}\int_{t_0}^{t}dt'\,\hat H(t')\right]}.
$$

For noncommuting Hamiltonians, this exponential is defined by the Dyson series and cannot be treated as an ordinary operator exponential.

##### Commuting and Time-Independent Hamiltonians

If the Hamiltonians commute at all times, time ordering is unnecessary:

$$
[\hat H(t_1),\hat H(t_2)]=0\quad\Longrightarrow\quad\hat U(t,t_0)=\exp\left[-\frac{i}{\hbar}\int_{t_0}^{t}dt'\,\hat H(t')\right].
$$

For a time-independent Hamiltonian $\hat H$,

$$
\boxed{\hat U(t,t_0)=\exp\left[-\frac{i}{\hbar}\hat H(t-t_0)\right]}.
$$

Let $\{|E_n\rangle\}$ be a complete orthonormal eigenbasis of $\hat H$:

$$
\hat H|E_n\rangle=E_n|E_n\rangle,\qquad\hat U(t,t_0)=\sum_n\exp\left[-\frac{iE_n(t-t_0)}{\hbar}\right]|E_n\rangle\langle E_n|.
$$

Each energy coefficient acquires only a phase:

$$
|\psi(t_0)\rangle=\sum_n c_n|E_n\rangle,\qquad|\psi(t)\rangle=\sum_n c_n\exp\left[-\frac{iE_n(t-t_0)}{\hbar}\right]|E_n\rangle,\qquad c_n=\langle E_n|\psi(t_0)\rangle.
$$

### Quantum Pictures

Let $t_0$ be the reference time, $\hat U(t,t_0)$ the [time-evolution operator](#time-evolution-operator), $|\psi_0\rangle$ the common initial state, and $\hat A_0$ the common reference observable:

$$
\boxed{|\psi_S(t_0)\rangle=|\psi_I(t_0)\rangle=|\psi_H(t_0)\rangle=|\psi_0\rangle,\qquad\hat A_S(t_0)=\hat A_I(t_0)=\hat A_H(t_0)=\hat A_0}.
$$

Unless stated otherwise, Schrödinger-picture observables have no explicit time dependence, while $\hat H_S(t)$ may depend on time.

##### Schrödinger Picture

The state evolves while the observable remains fixed:

$$
\boxed{|\psi_S(t)\rangle=\hat U(t,t_0)|\psi_0\rangle,\qquad\hat A_S(t)=\hat A_S(t_0)=\hat A_0}.
$$

$$
i\hbar\frac{\partial}{\partial t}|\psi_S(t)\rangle=\hat H_S(t)|\psi_S(t)\rangle.
$$

##### Heisenberg Picture

The state remains fixed while the observable evolves:

$$
\boxed{|\psi_H(t)\rangle=|\psi_H(t_0)\rangle=|\psi_0\rangle,\qquad\hat A_H(t)=\hat U^\dagger(t,t_0)\hat A_0\hat U(t,t_0)}.
$$

Defining $\hat H_H(t)=\hat U^\dagger(t,t_0)\hat H_S(t)\hat U(t,t_0)$ gives

$$
\boxed{\frac{d\hat A_H}{dt}=\frac{i}{\hbar}[\hat H_H(t),\hat A_H(t)]}.
$$

Thus, $\hat A_H$ is conserved if it commutes with $\hat H_H$.

##### Interaction Picture

Decompose the Hamiltonian into a time-independent solvable part $\hat H_0$ and an interaction $\hat V_S(t)$:

$$
\hat H_S(t)=\hat H_0+\hat V_S(t),\qquad\hat U_0(t,t_0)=\exp\left[-\frac{i}{\hbar}\hat H_0(t-t_0)\right].
$$

The interaction-picture quantities are

$$
|\psi_I(t)\rangle=\hat U_0^\dagger(t,t_0)|\psi_S(t)\rangle,\qquad\hat A_I(t)=\hat U_0^\dagger(t,t_0)\hat A_0\hat U_0(t,t_0),\qquad\hat V_I(t)=\hat U_0^\dagger(t,t_0)\hat V_S(t)\hat U_0(t,t_0).
$$

Therefore,

$$
\boxed{i\hbar\frac{\partial}{\partial t}|\psi_I(t)\rangle=\hat V_I(t)|\psi_I(t)\rangle}.
$$

Defining $\hat U_I(t,t_0)=\hat U_0^\dagger(t,t_0)\hat U(t,t_0)$ gives

$$
|\psi_I(t)\rangle=\hat U_I(t,t_0)|\psi_0\rangle,\qquad i\hbar\frac{\partial}{\partial t}\hat U_I(t,t_0)=\hat V_I(t)\hat U_I(t,t_0).
$$

Its Dyson series follows from the preceding section by replacing $\hat H(t)$ with $\hat V_I(t)$.

All three pictures give the same expectation value:

$$
\boxed{\langle\psi_S(t)|\hat A_0|\psi_S(t)\rangle=\langle\psi_I(t)|\hat A_I(t)|\psi_I(t)\rangle=\langle\psi_0|\hat A_H(t)|\psi_0\rangle}.
$$

##### Ehrenfest Theorem

For a normalized Schrödinger-picture state, define

$$
\langle\hat A\rangle_t=\langle\psi_S(t)|\hat A_S(t)|\psi_S(t)\rangle.
$$

Using the Schrödinger equation and its adjoint,

$$
i\hbar\frac{\partial}{\partial t}|\psi_S(t)\rangle=\hat H_S(t)|\psi_S(t)\rangle,\qquad-i\hbar\frac{\partial}{\partial t}\langle\psi_S(t)|=\langle\psi_S(t)|\hat H_S(t),
$$

gives the Ehrenfest theorem:

$$
\boxed{\frac{d}{dt}\langle\hat A\rangle_t=\left\langle\frac{\partial\hat A_S}{\partial t}\right\rangle_t+\frac{i}{\hbar}\left\langle[\hat H_S(t),\hat A_S(t)]\right\rangle_t}.
$$

For a one-dimensional particle,

$$
\hat H_S(t)=\frac{\hat p^2}{2m}+V(\hat x,t),\qquad[\hat x,\hat p]=i\hbar\hat I,\qquad[\hat H_S,\hat x]=-\frac{i\hbar}{m}\hat p,\qquad[\hat H_S,\hat p]=i\hbar\frac{\partial V}{\partial\hat x}.
$$

Hence,

$$
\boxed{\frac{d}{dt}\langle\hat x\rangle_t=\frac{\langle\hat p\rangle_t}{m},\qquad\frac{d}{dt}\langle\hat p\rangle_t=-\left\langle\frac{\partial V}{\partial\hat x}\right\rangle_t\quad\Longrightarrow\quad m\frac{d^2}{dt^2}\langle\hat x\rangle_t=-\left\langle\frac{\partial V}{\partial\hat x}\right\rangle_t}.
$$

This becomes the classical Newton equation when $V(x,t)$ is at most quadratic in $x$.

### Path Integral

##### Propagator

Let $t_i<t_f$ be the initial and final times, with corresponding coordinates $x_i$ and $x_f$. The position-space propagator is

$$
\boxed{K(x_f,t_f;x_i,t_i)=\langle x_f|\hat U(t_f,t_i)|x_i\rangle}.
$$

Writing $\psi(x,t)=\langle x|\psi(t)\rangle$ gives

$$
\boxed{\psi(x_f,t_f)=\int dx_i\,K(x_f,t_f;x_i,t_i)\psi(x_i,t_i)}.
$$

For an intermediate time $t_i<t_m<t_f$,

$$
K(x_f,t_f;x_i,t_i)=\int dx_m\,K(x_f,t_f;x_m,t_m)K(x_m,t_m;x_i,t_i),\qquad K(x_f,t_i;x_i,t_i)=\delta(x_f-x_i).
$$

##### Time Slicing

Divide $[t_i,t_f]$ into $N$ equal intervals:

$$
\Delta t=\frac{t_f-t_i}{N},\qquad t_j=t_i+j\Delta t,\qquad t_0=t_i,\qquad t_N=t_f,\qquad x_0=x_i,\qquad x_N=x_f.
$$

Inserting a position-space completeness relation at every intermediate time gives

$$
\boxed{\langle x_f|\hat U(t_f,t_i)|x_i\rangle=\int\prod_{j=1}^{N-1}dx_j\,\langle x_N|\hat U(t_N,t_{N-1})|x_{N-1}\rangle\cdots\langle x_1|\hat U(t_1,t_0)|x_0\rangle=\int\prod_{j=1}^{N-1}dx_j\prod_{j=0}^{N-1}\langle x_{j+1}|\hat U(t_{j+1},t_j)|x_j\rangle}.
$$

This decomposition is exact for finite $N$. The continuum limit is

$$
N\to\infty,\qquad\Delta t\to0,\qquad N\Delta t=t_f-t_i.
$$

##### Path-Integral Representation

Consider a one-dimensional particle of mass $m$:

$$
\hat H=\hat T+\hat V=\frac{\hat p^2}{2m}+V(\hat x),\qquad H(p,x)=\frac{p^2}{2m}+V(x).
$$

Let $|p\rangle$ be normalized by

$$
\int dp\,|p\rangle\langle p|=\hat I,\qquad \langle x|p\rangle=\frac{1}{\sqrt{2\pi\hbar}}\exp\left(\frac{ipx}{\hbar}\right).
$$

Insert the momentum completeness relation and apply the first-order Trotter decomposition:

$$
\langle x_{j+1}|\hat U(t_{j+1},t_j)|x_j\rangle\simeq\int dp_j\,\left\langle x_{j+1}\left|\exp\left(-\frac{i\Delta t}{\hbar}\hat T\right)\right|p_j\right\rangle\left\langle p_j\left|\exp\left(-\frac{i\Delta t}{\hbar}\hat V\right)\right|x_j\right\rangle.
$$

Using the momentum and position eigenvalue equations gives

$$
\langle x_{j+1}|\hat U(t_{j+1},t_j)|x_j\rangle\simeq\int\frac{dp_j}{2\pi\hbar}\exp\left[\frac{ip_j(x_{j+1}-x_j)}{\hbar}\right]\exp\left[-\frac{i\Delta t}{\hbar}H(p_j,x_j)\right].
$$

The accumulated Trotter error is $O(\Delta t)$ and vanishes in the continuum limit. Hence,

$$
\boxed{\langle x_f|\hat U(t_f,t_i)|x_i\rangle=\lim_{N\to\infty}\int\prod_{j=1}^{N-1}dx_j\prod_{j=0}^{N-1}\frac{dp_j}{2\pi\hbar}\exp\left\{\frac{i}{\hbar}\sum_{j=0}^{N-1}\left[p_j(x_{j+1}-x_j)-\Delta t\,H(p_j,x_j)\right]\right\}}.
$$

With $\dot x=dx/dt$, the continuum phase-space representation is

$$
\boxed{\langle x_f|\hat U(t_f,t_i)|x_i\rangle=\int_{x(t_i)=x_i}^{x(t_f)=x_f}\mathcal D x\,\mathcal D p\,\exp\left(\frac{i}{\hbar}S[p,x]\right)}.
$$

$$
\mathcal D x=\lim_{N\to\infty}\prod_{j=1}^{N-1}dx_j,\qquad\mathcal D p=\lim_{N\to\infty}\prod_{j=0}^{N-1}\frac{dp_j}{2\pi\hbar},\qquad S[p,x]=\int_{t_i}^{t_f}dt\,[p\dot x-H(p,x)].
$$

For the quadratic kinetic energy, each momentum integral is Gaussian:

$$
\int\frac{dp_j}{2\pi\hbar}\exp\left[\frac{ip_j(x_{j+1}-x_j)}{\hbar}\right]\exp\left[-\frac{i\Delta t}{\hbar}\left(\frac{p_j^2}{2m}+V(x_j)\right)\right]=\sqrt{\frac{m}{2\pi i\hbar\Delta t}}\exp\left\{\frac{i\Delta t}{\hbar}\left[\frac{m}{2}\left(\frac{x_{j+1}-x_j}{\Delta t}\right)^2-V(x_j)\right]\right\}.
$$

Therefore,

$$
\boxed{\langle x_f|\hat U(t_f,t_i)|x_i\rangle=\int_{x(t_i)=x_i}^{x(t_f)=x_f}\mathcal D_{\mathrm{conf}}x\,\exp\left(\frac{i}{\hbar}S[x]\right)}.
$$

$$
\mathcal D_{\mathrm{conf}}x=\lim_{N\to\infty}\left(\frac{m}{2\pi i\hbar\Delta t}\right)^{N/2}\prod_{j=1}^{N-1}dx_j,\qquad S[x]=\int_{t_i}^{t_f}dt\,L(x,\dot x),\qquad L(x,\dot x)=\frac{m}{2}\dot x^2-V(x).
$$

### Berry Phase

##### Adiabatic Evolution and Berry Phase

Let $\mathbf R(t)$ be slowly varying external parameters and $\hat H(t)=\hat H(\mathbf R(t))$ a smooth Hamiltonian with discrete, nondegenerate instantaneous eigenstates:

$$
\hat H(t)|n(t)\rangle=E_n(t)|n(t)\rangle,\qquad\langle m(t)|n(t)\rangle=\delta_{mn},
$$

where $\delta_{mn}$ is the Kronecker delta. The state satisfies the [Schrödinger equation](#schrödinger-equation) and may be expanded as

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

##### Berry Connection, Curvature, and Gauge Structure

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
