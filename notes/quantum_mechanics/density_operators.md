# Density Operators

### Density Operators and Mixed States

##### Statistical Ensembles

Let $\mathcal H$ be the Hilbert space of a quantum system. Consider an ensemble in which the system occupies a normalized state $|\psi_\alpha\rangle\in\mathcal H$ with probability $p_\alpha$, where

$$
p_\alpha\geq0,
\qquad
\sum_\alpha p_\alpha=1,
\qquad
\langle\psi_\alpha|\psi_\alpha\rangle=1.
$$

The states $|\psi_\alpha\rangle$ need not be mutually orthogonal. The ensemble is represented by the density operator

$$
\hat\rho
=
\sum_\alpha
p_\alpha
|\psi_\alpha\rangle
\langle\psi_\alpha|.
$$

Let $\{|e_i\rangle\}$ be an orthonormal and complete basis of $\mathcal H$. The corresponding density matrix has elements

$$
\rho_{ij}
=
\langle e_i|\hat\rho|e_j\rangle
=
\sum_\alpha
p_\alpha
\langle e_i|\psi_\alpha\rangle
\langle\psi_\alpha|e_j\rangle.
$$

The density operator $\hat\rho$ is independent of the chosen basis, whereas the matrix elements $\rho_{ij}$ depend on that basis. Different ensembles may produce the same density operator; the decomposition into probabilities and state vectors is therefore not unique.

##### Characterization of a Density Operator

Let $|\phi\rangle\in\mathcal H$. A density operator is a positive-semidefinite trace-class operator with unit trace:

$$
\hat\rho^\dagger=\hat\rho
,\qquad
\langle\phi|\hat\rho|\phi\rangle
\geq0
,\qquad
\operatorname{Tr}\hat\rho
=
1.
$$

For the ensemble representation,

$$
\langle\phi|\hat\rho|\phi\rangle=
\sum_\alpha
p_\alpha
\left|
\langle\phi|\psi_\alpha\rangle
\right|^2
\geq0.
$$

The spectral decomposition of $\hat\rho$ is

$$
\hat\rho=\sum_n\lambda_n
|n\rangle\langle n|,
$$

where $\{|n\rangle\}$ is an orthonormal eigenbasis and

$$
\lambda_n\geq0,
\qquad
\sum_n\lambda_n=1.
$$

Thus, every density operator may be interpreted as an ensemble of mutually orthogonal states with probabilities $\lambda_n$.

##### Pure and Mixed States

A normalized state $|\psi\rangle$ defines the pure-state density operator

$$
\hat\rho_\psi
=
|\psi\rangle\langle\psi|.
$$

This is the [projection operator](formalism.md#projection-operator) onto the one-dimensional subspace spanned by $|\psi\rangle$. A density operator represents a pure state if and only if any of the following equivalent conditions holds:

$$
\hat\rho^2=\hat\rho
,\qquad
\operatorname{rank}\hat\rho=1
,\qquad
\operatorname{Tr}\hat\rho^2=1.
$$

Otherwise, it represents a mixed state and

$$
\operatorname{Tr}\hat\rho^2
<
1.
$$

The quantity

$$
\gamma
=
\operatorname{Tr}\hat\rho^2
$$

is called the purity of the state.

A coherent superposition is not the same as an incoherent mixture. Let $\{|\uparrow_z\rangle,|\downarrow_z\rangle\}$ be the spin basis along the $z$ axis, and define

$$
|{+x}\rangle
=
\frac{
|\uparrow_z\rangle
+
|\downarrow_z\rangle
}{
\sqrt{2}
}.
$$

The coherent pure state has density matrix

$$
\rho_{\mathrm{coh}}
=
|{+x}\rangle\langle{+x}|
=
\frac{1}{2}
\begin{pmatrix}
1&1\\
1&1
\end{pmatrix}.
$$

An equal incoherent mixture of the two $z$-spin states has density matrix

$$
\rho_{\mathrm{mix}}
=
\frac{1}{2}
|\uparrow_z\rangle\langle\uparrow_z|
+
\frac{1}{2}
|\downarrow_z\rangle\langle\downarrow_z|
=
\frac{1}{2}
\begin{pmatrix}
1&0\\
0&1
\end{pmatrix}.
$$

The diagonal elements give the same probabilities for a measurement along $z$, but the off-diagonal elements retain the phase coherence of the pure state.

##### Measurements and Expectation Values

Let the self-adjoint observable $\hat A$ have the spectral decomposition

$$
\hat A
=
\sum_a
a\hat P_a,
$$

where $a$ is a measurement outcome and $\hat P_a$ is the orthogonal projector onto its eigenspace. The probability of obtaining $a$ is

$$
p(a)
=
\operatorname{Tr}
\left(
\hat\rho\hat P_a
\right).
$$

The expectation value of $\hat A$ is

$$
\langle\hat A\rangle
=
\operatorname{Tr}
\left(
\hat\rho\hat A
\right)
=
\sum_a
a\,p(a).
$$

For a pure state, this reduces to

$$
\operatorname{Tr}
\left(
\hat\rho_\psi\hat A
\right)
=
\langle\psi|\hat A|\psi\rangle.
$$

##### Unitary Evolution

Let $\hat H$ be a time-independent Hamiltonian and define the time-evolution operator

$$
\hat U(t)
=
\exp\left(
-\frac{i\hat Ht}{\hbar}
\right).
$$

Under the [unitary transformation](formalism.md#unitary-transformation),

$$
\hat\rho(t)
=
\hat U(t)
\hat\rho(0)
\hat U^\dagger(t).
$$

Differentiation gives the Liouville–von Neumann equation

$$
i\hbar
\frac{d\hat\rho}{dt}
=
[\hat H,\hat\rho].
$$

Unitary evolution preserves the eigenvalues and purity of $\hat\rho$:

$$
\operatorname{Tr}\hat\rho^2(t)
=
\operatorname{Tr}\hat\rho^2(0).
$$

##### Reduced Density Operator

Let a composite system have Hilbert space

$$
\mathcal H_{AB}
=
\mathcal H_A
\otimes
\mathcal H_B,
$$

and let $\hat\rho_{AB}$ be its density operator. For an orthonormal and complete basis $\{|b_j\rangle\}$ of $\mathcal H_B$, the reduced density operator of subsystem $A$ is

$$
\hat\rho_A
=
\operatorname{Tr}_B\hat\rho_{AB}
=
\sum_j
\langle b_j|
\hat\rho_{AB}
|b_j\rangle.
$$

Let $\hat A$ act on $\mathcal H_A$, and let $\hat I_B$ be the identity operator on $\mathcal H_B$. Then

$$
\operatorname{Tr}_{AB}
\left[
\hat\rho_{AB}
\left(
\hat A\otimes\hat I_B
\right)
\right]
=
\operatorname{Tr}_A
\left(
\hat\rho_A\hat A
\right).
$$

A subsystem may therefore be mixed even when the composite system is pure. Define the entangled state

$$
|\Phi^+\rangle
=
\frac{
|0\rangle_A|0\rangle_B
+
|1\rangle_A|1\rangle_B
}{
\sqrt{2}
}.
$$

For

$$
\hat\rho_{AB}
=
|\Phi^+\rangle\langle\Phi^+|,
$$

the reduced density operators are

$$
\hat\rho_A
=
\hat\rho_B
=
\frac{1}{2}\hat I_2,
$$

where $\hat I_2$ is the identity operator on the corresponding two-dimensional subsystem.
