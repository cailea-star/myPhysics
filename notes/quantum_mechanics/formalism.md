# Quantum Formalism

### Axioms of Quantum Mechanics

- A pure state is a ray in a complex Hilbert space $\mathcal H$, represented by a normalized ket:

$$
\langle\psi|\psi\rangle=1, \qquad |\psi\rangle\sim e^{i\alpha}|\psi\rangle, \qquad \alpha\in\mathbb R.
$$

General states are described by [density operators](density_operators.md).

- An observable is represented by a self-adjoint operator:

$$
\hat A^\dagger=\hat A.
$$

- For the spectral projector $\hat P_a$ associated with outcome $a$, the Born probability and post-measurement state are

$$
\boxed{p(a)=\langle\psi|\hat P_a|\psi\rangle, \qquad |\psi\rangle\longrightarrow\frac{\hat P_a|\psi\rangle}{\sqrt{p(a)}}\quad\bigl(p(a)>0\bigr)}.
$$

- A closed system evolves according to the Schrödinger equation:

$$
\boxed{i\hbar\frac{\partial}{\partial t}|\psi(t)\rangle=\hat H(t)|\psi(t)\rangle}.
$$

- The Hilbert space of a composite system is the tensor product of its subsystems:

$$
\mathcal H=\bigotimes_k\mathcal H_k.
$$

For identical particles, let $\hat\Pi_{ij}$ exchange particles $i$ and $j$. Bosonic and fermionic states satisfy

$$
\hat\Pi_{ij}|\psi\rangle=|\psi\rangle\quad\text{(bosons)}, \qquad \hat\Pi_{ij}|\psi\rangle=-|\psi\rangle\quad\text{(fermions)}.
$$

### States

##### Inner Product

Let $\mathcal H$ be a complete complex inner-product space, with $|\phi\rangle,|\psi\rangle,|\chi\rangle\in\mathcal H$ and $a,b\in\mathbb C$.

- Conjugate symmetry:

$$
\langle\phi|\psi\rangle = \langle\psi|\phi\rangle^*.
$$

- Linearity in the second argument:

$$
\langle\chi|a\psi+b\phi\rangle = a\langle\chi|\psi\rangle + b\langle\chi|\phi\rangle.
$$

- Positive definiteness:

$$
\langle\phi|\phi\rangle\geq0, \qquad \langle\phi|\phi\rangle=0 \iff |\phi\rangle=0.
$$

- Norm:

$$
\lVert\phi\rVert = \sqrt{\langle\phi|\phi\rangle}.
$$

Completeness means that every Cauchy sequence under this norm converges in $\mathcal H$.

- Cauchy–Schwarz inequality:

$$
|\langle\phi|\psi\rangle|^2 \leq \langle\phi|\phi\rangle \langle\psi|\psi\rangle.
$$

##### Orthonormality

Let $\{|\phi_i\rangle\}$ be an orthonormal and complete basis, $\delta_{ij}$ the Kronecker delta, and $\hat I$ the identity operator:

$$
\langle\phi_i|\phi_j\rangle = \delta_{ij}, \qquad \sum_i |\phi_i\rangle\langle\phi_i| = \hat I.
$$

Any nonzero vector may be normalized:

$$
|\phi\rangle \longrightarrow \frac{|\phi\rangle}{\lVert\phi\rVert}.
$$

##### Dual Space

The inner product associates each ket in $\mathcal H$ with a bra in the dual space $\mathcal H^*$:

$$
\langle\phi| = \left( |\phi\rangle \right)^\dagger.
$$

For a linear operator $\hat A$, its adjoint $\hat A^\dagger$ is defined by

$$
\langle\phi|\hat A^\dagger|\psi\rangle = \langle\psi|\hat A|\phi\rangle^*.
$$

### Operators

##### Self-Adjoint Operator

A linear operator $\hat A$ on $\mathcal H$ is self-adjoint when

$$
\hat A^\dagger=\hat A.
$$

For a normalized state $|\phi\rangle$, its expectation value is real:

$$
\langle\hat A\rangle_\phi=\langle\phi|\hat A|\phi\rangle=\langle\phi|\hat A|\phi\rangle^*.
$$

Its eigenvalues are also real:

$$
\hat A|a_i\rangle=a_i|a_i\rangle, \qquad a_i=a_i^*.
$$

##### Orthonormal Eigenbasis

For a self-adjoint operator with a discrete spectrum, its eigenstates may be chosen orthonormally:

$$
\langle a_i|a_j\rangle=\delta_{ij}, \qquad \sum_i|a_i\rangle\langle a_i|=\hat I.
$$

The spectral decomposition is

$$
\boxed{\hat A=\sum_i a_i|a_i\rangle\langle a_i|}.
$$

##### Projection Operator

For the eigenstate $|a_i\rangle$, define

$$
\hat P_i=|a_i\rangle\langle a_i|.
$$

The projectors satisfy

$$
\hat P_i\hat P_j=\delta_{ij}\hat P_i, \qquad \sum_i\hat P_i=\hat I.
$$

A general orthogonal projector satisfies

$$
\hat P^\dagger=\hat P, \qquad \hat P^2=\hat P.
$$

For a normalized state $|\psi\rangle$,

$$
p_i=\langle\psi|\hat P_i|\psi\rangle, \qquad \langle\hat A\rangle=\sum_i p_i a_i.
$$

##### Operator Function

For a power series

$$
f(z)=\sum_{n=0}^{\infty}c_nz^n,
$$

the corresponding operator function is

$$
f(\hat A)=\sum_{n=0}^{\infty}c_n\hat A^n=\sum_i f(a_i)\hat P_i.
$$

### Canonical Quantization

##### Canonical Commutation Relations

For operators $\hat A$, $\hat B$, and $\hat C$, define

$$
[\hat A,\hat B]=\hat A\hat B-\hat B\hat A.
$$

A useful product rule is

$$
[\hat A\hat B,\hat C]=\hat A[\hat B,\hat C]+[\hat A,\hat C]\hat B.
$$

For Cartesian canonical coordinates,

$$
\boxed{[\hat x_i,\hat x_j]=0, \qquad [\hat p_i,\hat p_j]=0, \qquad [\hat x_i,\hat p_j]=i\hbar\delta_{ij}\hat I}.
$$

The position and momentum operators are self-adjoint:

$$
\hat x^\dagger=\hat x, \qquad \hat p^\dagger=\hat p.
$$

Their generalized eigenstates satisfy

$$
\hat x|x\rangle=x|x\rangle, \qquad \hat p|p\rangle=p|p\rangle.
$$

$$
\langle x|x'\rangle=\delta(x-x'), \qquad \int dx\,|x\rangle\langle x|=\hat I.
$$

$$
\langle p|p'\rangle=\delta(p-p'), \qquad \int dp\,|p\rangle\langle p|=\hat I.
$$

In $d$ spatial dimensions,

$$
\boxed{\langle\mathbf x|\mathbf p\rangle=\frac{1}{(2\pi\hbar)^{d/2}}\exp\left(\frac{i\mathbf p\cdot\mathbf x}{\hbar}\right)}.
$$

##### Translation Operator

The Hadamard lemma is

$$
e^{\hat A}\hat B e^{-\hat A}=\hat B+[\hat A,\hat B]+\frac{1}{2!}[\hat A,[\hat A,\hat B]]+\cdots.
$$

For a displacement $x_0$, define

$$
\hat T(x_0)=\exp\left(-\frac{ix_0\hat p}{\hbar}\right).
$$

Using the canonical commutation relation,

$$
\boxed{\hat T^\dagger(x_0)\hat x\hat T(x_0)=\hat x+x_0\hat I}.
$$

Its action on the position and momentum eigenstates is

$$
\hat T(x_0)|x\rangle=|x+x_0\rangle, \qquad \hat T(x_0)|p\rangle=e^{-ix_0p/\hbar}|p\rangle.
$$

##### Uncertainty Relation

Let $\hat A$ and $\hat B$ be self-adjoint operators and $|\psi\rangle$ a normalized state. Define

$$
\langle\hat A\rangle=\langle\psi|\hat A|\psi\rangle, \qquad \Delta A=\sqrt{\langle\psi|(\hat A-\langle\hat A\rangle)^2|\psi\rangle},
$$

and similarly for $\Delta B$. For a real parameter $\xi$, define

$$
|\Phi(\xi)\rangle=\left[\xi(\hat A-\langle\hat A\rangle)+i(\hat B-\langle\hat B\rangle)\right]|\psi\rangle.
$$

Positivity gives

$$
\langle\Phi(\xi)|\Phi(\xi)\rangle=(\Delta A)^2\xi^2+\langle\psi|i[\hat A,\hat B]|\psi\rangle\xi+(\Delta B)^2\geq0.
$$

Therefore,

$$
\boxed{\Delta A\,\Delta B\geq\frac{1}{2}\left|\langle\psi|[\hat A,\hat B]|\psi\rangle\right|}.
$$

In particular,

$$
\boxed{\Delta x\,\Delta p\geq\frac{\hbar}{2}}.
$$

### Representations and Unitary Transformations

##### Representation Transformation

A state and an operator are independent of the chosen basis. Let $\{|a\rangle\}$ and $\{|b\rangle\}$ be orthonormal and complete bases.

In the $a$ representation, the state and operator are

$$
|\psi\rangle=\sum_a|a\rangle\langle a|\psi\rangle.
$$

$$
\hat O=\sum_{a,a'}|a\rangle\langle a|\hat O|a'\rangle\langle a'|.
$$

Their components in the $b$ representation are

$$
\langle b|\psi\rangle=\sum_a\langle b|a\rangle\langle a|\psi\rangle.
$$

$$
\langle b|\hat O|b'\rangle=\sum_{a,a'}\langle b|a\rangle\langle a|\hat O|a'\rangle\langle a'|b'\rangle.
$$

For continuous bases, the sums are replaced by integrals.

##### Unitary Transformation

An operator $\hat U$ is unitary if

$$
\hat U^\dagger\hat U=\hat U\hat U^\dagger=\hat I.
$$

It preserves inner products:

$$
\langle\phi|\hat U^\dagger\hat U|\psi\rangle=\langle\phi|\psi\rangle.
$$

For two indexed orthonormal bases $\{|a_n\rangle\}$ and $\{|b_n\rangle\}$, define

$$
\hat U=\sum_n|b_n\rangle\langle a_n|.
$$

Then

$$
\hat U|a_n\rangle=|b_n\rangle, \qquad \langle b_n|\psi\rangle=\langle a_n|\hat U^\dagger|\psi\rangle.
$$

Under an active unitary transformation,

$$
\boxed{|\psi\rangle\longrightarrow\hat U|\psi\rangle, \qquad \hat O\longrightarrow\hat U\hat O\hat U^\dagger}.
$$

For a unitary operator with a discrete spectrum,

$$
\hat U|c_n\rangle=\lambda_n|c_n\rangle, \qquad |\lambda_n|=1, \qquad \lambda_n=e^{i\theta_n}, \qquad \theta_n\in\mathbb R.
$$

Define the self-adjoint generator

$$
\hat G=\sum_n\theta_n|c_n\rangle\langle c_n|, \qquad \hat G^\dagger=\hat G.
$$

Then

$$
\boxed{\hat U=\sum_n e^{i\theta_n}|c_n\rangle\langle c_n|=e^{i\hat G}}.
$$
