# Quantum Formalism

### Axioms of Quantum Mechanics

- The state of a microscopic system is described by a vector in a complex Hilbert space. The Hilbert space provides:
  - a complex linear structure for superposition;
  - a positive-definite inner product defining norms, normalization, and orthogonality;
  - completeness under the norm induced by the inner product;
  - a natural correspondence between vectors and dual vectors, together with adjoint operations on linear operators.
- An observable of a microscopic system is described by a self-adjoint operator on the Hilbert space.
- Canonical equal-time quantization is specified by the commutation relations between position and momentum.
- The system evolves according to the Schrödinger equation.
- For identical particles, the many-particle wave function must be either symmetric or antisymmetric.

### States

##### Inner Product

Let $\mathcal H$ be a complete complex inner-product space, with $|\phi\rangle,|\psi\rangle,|\chi\rangle\in\mathcal H$ and $a,b\in\mathbb C$.

- Conjugate symmetry:

$$
\langle\phi|\psi\rangle
=
\langle\psi|\phi\rangle^*.
$$

- Linearity in the second argument:

$$
\langle\chi|a\psi+b\phi\rangle
=
a\langle\chi|\psi\rangle
+
b\langle\chi|\phi\rangle.
$$

- Positive definiteness:

$$
\langle\phi|\phi\rangle\geq0,
\qquad
\langle\phi|\phi\rangle=0
\iff
|\phi\rangle=0.
$$

- Norm:

$$
\lVert\phi\rVert
=
\sqrt{\langle\phi|\phi\rangle}.
$$

Completeness means that every Cauchy sequence under this norm converges to a vector in $\mathcal H$.

- Cauchy–Schwarz inequality:

$$
|\langle\phi|\psi\rangle|^2
\leq
\langle\phi|\phi\rangle
\langle\psi|\psi\rangle.
$$

##### Orthonormality

Let $\{|\phi_i\rangle\}$ be an orthonormal and complete basis, $\delta_{ij}$ the Kronecker delta, and $\hat I$ the identity operator:

$$
\langle\phi_i|\phi_j\rangle
=
\delta_{ij},
\qquad
\sum_i
|\phi_i\rangle\langle\phi_i|
=
\hat I.
$$

Any nonzero vector may be normalized:

$$
|\phi\rangle
\longrightarrow
\frac{|\phi\rangle}{\lVert\phi\rVert}.
$$

The Gram–Schmidt procedure converts a linearly independent set into an orthonormal set.

##### Dual Space

Let $\mathcal H^*$ be the dual space of $\mathcal H$. The inner product associates each ket with a bra:

$$
\langle\phi|
=
\left(
|\phi\rangle
\right)^\dagger.
$$

Let $\hat A$ be a linear operator and $\hat A^\dagger$ its adjoint. If

$$
\hat A|\phi\rangle
=
|\eta\rangle,
$$

then

$$
\langle\phi|\hat A^\dagger
=
\langle\eta|.
$$

Equivalently,

$$
\langle\phi|\hat A^\dagger|\psi\rangle
=
\langle\psi|\hat A|\phi\rangle^*.
$$

### Operators

##### Self-Adjoint Operator

Let $\hat A$ be a linear operator on $\mathcal H$. It is self-adjoint when

$$
\hat A^\dagger=\hat A.
$$

For a normalized state $|\phi\rangle$,

$$
\langle\hat A\rangle_\phi
=
\langle\phi|\hat A|\phi\rangle
=
\langle\phi|\hat A|\phi\rangle^*.
$$

The eigenvalue equation is

$$
\hat A|a_i\rangle
=
a_i|a_i\rangle,
\qquad
a_i=a_i^*.
$$

##### Orthonormal Eigenbasis

For a discrete spectrum,

$$
\langle a_i|a_j\rangle
=
\delta_{ij},
\qquad
\sum_i|a_i\rangle\langle a_i|
=
\hat I.
$$

Within a degenerate eigenspace, the eigenstates may be chosen orthonormally and distinguished by additional commuting operators.

$$
\hat A
=
\sum_i
a_i|a_i\rangle\langle a_i|.
$$

##### Projection Operator

Define

$$
\hat P_i
=
|a_i\rangle\langle a_i|.
$$

Then

$$
\hat P_i\hat P_j
=
\delta_{ij}\hat P_i,
\qquad
\sum_i\hat P_i
=
\hat I.
$$

An orthogonal projector satisfies

$$
\hat P^\dagger=\hat P,
\qquad
\hat P^2=\hat P.
$$

For a normalized state $|\psi\rangle$,

$$
p_i
=
\langle\psi|\hat P_i|\psi\rangle,
\qquad
\langle\hat A\rangle
=
\sum_i p_i a_i.
$$

##### Operator Function

If

$$
f(z)
=
\sum_{n=0}^{\infty}c_nz^n,
$$

then

$$
f(\hat A)
=
\sum_{n=0}^{\infty}c_n\hat A^n
=
\sum_i f(a_i)\hat P_i.
$$

For an analytic $f(z)$ and a contour $\Gamma$ enclosing the spectrum of $\hat A$,

$$
f(\hat A)
=
\frac{1}{2\pi i}
\int_\Gamma
dz\,
f(z)
\left(
z\hat I-\hat A
\right)^{-1}.
$$

### Canonical Quantization

##### Canonical Commutation Relations

For operators $\hat A$, $\hat B$, and $\hat C$, define

$$
[\hat A,\hat B]
=
\hat A\hat B-\hat B\hat A,
\qquad
\{\hat A,\hat B\}
=
\hat A\hat B+\hat B\hat A.
$$

Useful identities are

$$
[\hat A\hat B,\hat C]
=
\hat A[\hat B,\hat C]
+
[\hat A,\hat C]\hat B,
$$

$$
\{\hat A\hat B,\hat C\}
=
\hat A\{\hat B,\hat C\}
-
[\hat A,\hat C]\hat B.
$$

For Cartesian canonical coordinates,

$$
[\hat x_i,\hat p_j]
=
i\hbar\delta_{ij}\hat I.
$$

The position and momentum operators are self-adjoint:

$$
\hat x|x\rangle=x|x\rangle,
\qquad
\hat p|p\rangle=p|p\rangle.
$$

Their generalized eigenstates satisfy

$$
\langle x|x'\rangle
=
\delta(x-x'),
\qquad
\int dx\,|x\rangle\langle x|
=
\hat I,
$$

$$
\langle p|p'\rangle
=
\delta(p-p'),
\qquad
\int dp\,|p\rangle\langle p|
=
\hat I.
$$

##### Commutators and Conserved Quantities

The Baker-Campbell-Hausdorff expansion gives

$$
e^{\hat A}\hat B e^{-\hat A}
=
\hat B
+
[\hat A,\hat B]
+
\frac{1}{2!}
[\hat A,[\hat A,\hat B]]
+\cdots.
$$

If $[\hat A,\hat B]=c\hat I$ for a scalar $c$, then

$$
e^{\hat A}\hat B e^{-\hat A}
=
\hat B+c\hat I.
$$

For a Hamiltonian $\hat H$ and an operator $\hat O$,

$$
\frac{d\hat O}{dt}
=
\frac{\partial\hat O}{\partial t}
+
\frac{i}{\hbar}
[\hat H,\hat O].
$$

Thus, if $\partial_t\hat O=0$ and $[\hat H,\hat O]=0$, then $\hat O$ is conserved.

##### Translation Operator

For a displacement $x_0$, define

$$
\hat T(x_0)
=
\exp\left(
-\frac{ix_0\hat p}{\hbar}
\right).
$$

Then

$$
\hat T^\dagger(x_0)
\hat x
\hat T(x_0)
=
\hat x+x_0\hat I,
$$

$$
\hat T(x_0)|x\rangle
=
|x+x_0\rangle,
\qquad
\hat T(x_0)|p\rangle
=
e^{-ix_0p/\hbar}|p\rangle.
$$

For $d$ spatial dimensions,

$$
\langle\mathbf x|\mathbf p\rangle
=
\frac{1}{(2\pi\hbar)^{d/2}}
\exp\left(
\frac{i\mathbf p\cdot\mathbf x}{\hbar}
\right).
$$

##### Uncertainty Relation

Let $\hat A$ and $\hat B$ be self-adjoint operators and $|\psi\rangle$ a normalized state. Define

$$
\langle\hat A\rangle
=
\langle\psi|\hat A|\psi\rangle,
\qquad
\Delta A
=
\sqrt{
\langle\psi|
\left(
\hat A-\langle\hat A\rangle
\right)^2
|\psi\rangle
}.
$$

The uncertainty relation is

$$
\Delta A\,\Delta B
\geq
\frac{1}{2}
\left|
\langle\psi|
i[\hat A,\hat B]
|\psi\rangle
\right|.
$$

For a real parameter $\xi$, define

$$
|\Phi(\xi)\rangle
=
\left[
\xi
\left(
\hat A-\langle\hat A\rangle
\right)
+
i
\left(
\hat B-\langle\hat B\rangle
\right)
\right]
|\psi\rangle.
$$

Then

$$
\langle\Phi(\xi)|\Phi(\xi)\rangle
=
(\Delta A)^2\xi^2
+
\langle\psi|i[\hat A,\hat B]|\psi\rangle\xi
+
(\Delta B)^2
\geq0,
$$

which gives the uncertainty relation. In particular,

$$
\Delta x\,\Delta p
\geq
\frac{\hbar}{2}.
$$

Commuting self-adjoint operators admit a common eigenbasis; see [Orthonormal Eigenbasis](#orthonormal-eigenbasis).

### Representations and Unitary Transformations

##### Representation Transformation

A state and an operator are independent of the chosen basis. Let $\{|a\rangle\}$ and $\{|b\rangle\}$ be orthonormal and complete bases.

In the $a$ representation,

$$
|\psi\rangle
=
\sum_a
|a\rangle\langle a|\psi\rangle,
$$

$$
\hat O
=
\sum_{a,a'}
|a\rangle
\langle a|\hat O|a'\rangle
\langle a'|.
$$

The state components transform as

$$
\langle b|\psi\rangle
=
\sum_a
\langle b|a\rangle
\langle a|\psi\rangle.
$$

The operator matrix elements transform as

$$
\langle b|\hat O|b'\rangle
=
\sum_{a,a'}
\langle b|a\rangle
\langle a|\hat O|a'\rangle
\langle a'|b'\rangle.
$$

For continuous bases, the sums are replaced by integrals.

##### Unitary Transformation

An operator $\hat U$ is unitary if

$$
\hat U^\dagger\hat U
=
\hat U\hat U^\dagger
=
\hat I.
$$

It preserves inner products:

$$
\langle\phi|
\hat U^\dagger\hat U
|\psi\rangle
=
\langle\phi|\psi\rangle.
$$

For two indexed orthonormal bases $\{|a_n\rangle\}$ and $\{|b_n\rangle\}$, define

$$
\hat U
=
\sum_n
|b_n\rangle\langle a_n|.
$$

Then

$$
\hat U|a_n\rangle
=
|b_n\rangle,
\qquad
\langle b_n|\psi\rangle
=
\langle a_n|
\hat U^\dagger
|\psi\rangle.
$$

Under an active unitary transformation,

$$
|\psi\rangle
\longrightarrow
\hat U|\psi\rangle,
\qquad
\hat O
\longrightarrow
\hat U\hat O\hat U^\dagger.
$$

In an infinite-dimensional space, $\hat U^\dagger\hat U=\hat I$ alone is insufficient. For

$$
\hat S
=
\sum_{n=0}^{\infty}
|n+1\rangle\langle n|,
$$

one has

$$
\hat S^\dagger\hat S
=
\hat I,
\qquad
\hat S\hat S^\dagger
=
\hat I-|0\rangle\langle0|.
$$

Thus, $\hat S$ is an isometry but not a unitary operator.

For a unitary operator with a discrete spectrum,

$$
\hat U|c_n\rangle
=
\lambda_n|c_n\rangle,
\qquad
|\lambda_n|=1,
\qquad
\lambda_n=e^{i\theta_n},
$$

where $\theta_n\in\mathbb R$. Therefore,

$$
\hat U
=
\sum_n
e^{i\theta_n}
|c_n\rangle\langle c_n|
=
e^{i\hat G},
$$

with the self-adjoint generator

$$
\hat G
=
\sum_n
\theta_n
|c_n\rangle\langle c_n|,
\qquad
\hat G^\dagger=\hat G.
$$
