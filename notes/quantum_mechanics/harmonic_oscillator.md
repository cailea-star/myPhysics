# Quantum Harmonic Oscillator

### Fundamentals

##### Harmonic-Oscillator Hamiltonian

For a one-dimensional oscillator with mass $M$ and angular frequency $\omega$, let $\hat x$ and $\hat p$ be the position and momentum operators; see [Canonical Commutation Relations](formalism.md#canonical-commutation-relations).

$$
\hat H
=
\frac{\hat p^2}{2M}
+
\frac{M\omega^2}{2}\hat x^2.
$$

Define the characteristic momentum and length scales

$$
p_0
=
\sqrt{M\hbar\omega},
\qquad
x_0
=
\sqrt{\frac{\hbar}{M\omega}}.
$$

The Hamiltonian becomes

$$
\hat H
=
\frac{\hbar\omega}{2}
\left(
\frac{\hat p^2}{p_0^2}
+
\frac{\hat x^2}{x_0^2}
\right).
$$

##### Creation and Annihilation Operators

Using the characteristic scales $x_0$ and $p_0$, define the annihilation and creation operators by

$$
\hat a
=
\frac{1}{\sqrt{2}}
\left(
\frac{\hat x}{x_0}
+
i\frac{\hat p}{p_0}
\right),
$$

$$
\hat a^\dagger
=
\frac{1}{\sqrt{2}}
\left(
\frac{\hat x}{x_0}
-
i\frac{\hat p}{p_0}
\right).
$$

Let $\hat I$ be the identity operator. The canonical commutation relation gives

$$
[\hat a,\hat a^\dagger]
=
\hat I.
$$

The Hamiltonian becomes

$$
\hat H
=
\frac{\hbar\omega}{2}
\left(
\hat a\hat a^\dagger
+
\hat a^\dagger\hat a
\right)
=
\hbar\omega
\left(
\hat a^\dagger\hat a
+
\frac{1}{2}
\right).
$$

### Number Representation

##### Number Operator

Define the number operator

$$
\hat n
=
\hat a^\dagger\hat a.
$$

The Hamiltonian becomes

$$
\hat H
=
\hbar\omega
\left(
\hat n+\frac{1}{2}
\right).
$$

The number operator is self-adjoint and positive semidefinite:

$$
\hat n^\dagger
=
\hat n,
$$

$$
\langle\phi|\hat n|\phi\rangle
=
\lVert
\hat a|\phi\rangle
\rVert^2
\geq0.
$$

Its commutation relations are

$$
[\hat n,\hat a]
=
-\hat a,
\qquad
[\hat n,\hat a^\dagger]
=
\hat a^\dagger.
$$

Let $|m\rangle$ be an eigenstate of $\hat n$ with eigenvalue $m$:

$$
\hat n|m\rangle
=
m|m\rangle.
$$

Then

$$
\hat n\hat a|m\rangle
=
(m-1)\hat a|m\rangle,
$$

$$
\hat n\hat a^\dagger|m\rangle
=
(m+1)\hat a^\dagger|m\rangle.
$$

Let $n_0$ be the lowest eigenvalue. Positivity and the lowering relation require

$$
\hat a|n_0\rangle
=
0,
$$

and therefore

$$
\hat n|n_0\rangle
=\hat a^\dagger
\big[\hat a|n_0\rangle\big]
=0,\qquad
n_0=0.
$$

##### Number States

Let $|n\rangle$ be a normalized eigenstate of the number operator:

$$
\hat n|n\rangle
=
n|n\rangle,
\qquad
\langle n|n\rangle
=
1.
$$

The lowest state $|0\rangle$ is the vacuum state:

$$
\hat a|0\rangle
=
0.
$$

The commutation relations imply

$$
\hat a|n\rangle
=
c_n|n-1\rangle,
\qquad
\hat a^\dagger|n\rangle=
d_n|n+1\rangle.
$$

Their norms give

$$
|c_n|^2
=
\langle n|\hat a^\dagger\hat a|n\rangle=n,
\qquad
|d_n|^2
=
\langle n|\hat a\hat a^\dagger|n\rangle
=
n+1.
$$

Choosing the phases of the number states gives

$$
\hat a|n\rangle=
\sqrt{n}\,|n-1\rangle,
\qquad
\hat a^\dagger|n\rangle=
\sqrt{n+1}\,|n+1\rangle.
$$

Therefore,

$$
|n\rangle
=
\frac{1}{\sqrt{n}}
\hat a^\dagger|n-1\rangle,
$$

and repeated application gives

$$
|n\rangle
=
\frac{(\hat a^\dagger)^n}{\sqrt{n!}}|0\rangle,
\qquad
n\in\mathbb N_0.
$$

With $\delta_{mn}$ denoting the Kronecker delta, the number states satisfy

$$
\langle m|n\rangle
=
\delta_{mn},
\qquad
\sum_{n=0}^{\infty}
|n\rangle\langle n|
=
\hat I.
$$

The Hamiltonian is diagonal in the number representation:

$$
\hat H|n\rangle
=
E_n|n\rangle,
\qquad
E_n
=
\hbar\omega
\left(
n+\frac{1}{2}
\right).
$$

### Coherent States

##### Definition and Number-State Expansion

Let $\alpha\in\mathbb C$, with $\alpha^*$ denoting its complex conjugate. A coherent state $|\alpha\rangle$ is a normalized eigenstate of the annihilation operator:

$$
\hat a|\alpha\rangle
=
\alpha|\alpha\rangle.
$$

Expand $|\alpha\rangle$ in number states:

$$
|\alpha\rangle
=
\sum_{n=0}^{\infty}
c_n|n\rangle.
$$

Using

$$
\hat a|n\rangle
=
\sqrt n\,|n-1\rangle,
$$

the eigenvalue equation becomes

$$
\sum_{n=0}^{\infty}
\sqrt{n+1}\,c_{n+1}|n\rangle
=
\sum_{n=0}^{\infty}
\alpha c_n|n\rangle.
$$

Therefore,

$$
c_{n+1}
=
\frac{\alpha}{\sqrt{n+1}}c_n,
$$

and hence

$$
c_n
=
\frac{\alpha^n}{\sqrt{n!}}c_0.
$$

Normalization requires

$$
1
=
\langle\alpha|\alpha\rangle
=
|c_0|^2
\sum_{n=0}^{\infty}
\frac{|\alpha|^{2n}}{n!}
=
|c_0|^2e^{|\alpha|^2}.
$$

Choosing the overall phase such that $c_0$ is positive gives

$$
c_0
=
e^{-|\alpha|^2/2}.
$$

Thus,

$$
|\alpha\rangle
=
e^{-|\alpha|^2/2}
\sum_{n=0}^{\infty}
\frac{\alpha^n}{\sqrt{n!}}
|n\rangle.
$$

Using the number-state construction,

$$
|n\rangle
=
\frac{(\hat a^\dagger)^n}{\sqrt{n!}}|0\rangle,
$$

this can also be written as

$$
|\alpha\rangle=
e^{-|\alpha|^2/2}
\sum_{n=0}^{\infty}
\frac{\alpha^n}{n!}
(\hat a^\dagger)^n|0\rangle=
e^{-|\alpha|^2/2}
e^{\alpha\hat a^\dagger}
|0\rangle.
$$

Define the displacement operator

$$
\hat D(\alpha)
=
\exp
\left(
\alpha\hat a^\dagger
-
\alpha^*\hat a
\right).
$$

Since the commutator of $\alpha\hat a^\dagger$ and $-\alpha^*\hat a$ is $|\alpha|^2\hat I$, the Baker-Campbell-Hausdorff formula gives

$$
\hat D(\alpha)
=
e^{-|\alpha|^2/2}
e^{\alpha\hat a^\dagger}
e^{-\alpha^*\hat a}.
$$

Because $\hat a|0\rangle=0$,

$$
|\alpha\rangle
=
\hat D(\alpha)|0\rangle.
$$

##### Properties of Coherent States

Using the number-state expansion,

$$
|\alpha\rangle
=
e^{-|\alpha|^2/2}
\sum_{n=0}^{\infty}
\frac{\alpha^n}{\sqrt{n!}}
|n\rangle,\qquad
\langle\beta|=
e^{-|\beta|^2/2}
\sum_{m=0}^{\infty}
\frac{(\beta^*)^m}{\sqrt{m!}}
\langle m|,
$$

their overlap is

$$
\langle\beta|\alpha\rangle
=
e^{-(|\alpha|^2+|\beta|^2)/2}
\sum_{n=0}^{\infty}
\frac{(\beta^*\alpha)^n}{n!}=
e^{\beta^*\alpha-(|\alpha|^2+|\beta|^2)/2}=
e^{-|\alpha-\beta|^2}.
$$

Thus, coherent states with different eigenvalues are generally not orthogonal.

Let

$$
\alpha
=
re^{i\theta},
\qquad
d^2\alpha
=
d\operatorname{Re}\alpha\,
d\operatorname{Im}\alpha
=
r\,dr\,d\theta.
$$

Using the number-state expansion,

$$
\frac{1}{\pi}
\int_{\mathbb C}
d^2\alpha\,
|\alpha\rangle\langle\alpha|
=
\sum_{n,m=0}^{\infty}
\frac{|n\rangle\langle m|}{\sqrt{n!m!}}
\frac{1}{\pi}
\int_0^\infty
r\,dr\,e^{-r^2}r^{n+m}
\int_0^{2\pi}
d\theta\,e^{i(n-m)\theta}.
$$

The angular integral and the radial integral gives

$$
\int_0^{2\pi}
d\theta\,
e^{i(n-m)\theta}
=
2\pi\delta_{nm},
\qquad
2
\int_0^\infty
r^{2n+1}e^{-r^2}\,dr
=
n!.
$$

Hence,

$$
\frac{1}{\pi}
\int_{\mathbb C}
d^2\alpha\,
|\alpha\rangle\langle\alpha|
=
\sum_{n=0}^{\infty}
|n\rangle\langle n|
=
\hat I.
$$

Because coherent states are complete but not orthogonal, they form an overcomplete set.


##### Coherent-State Representation

Define the unnormalized coherent state

$$
|\bar\alpha\rangle
=
e^{|\alpha|^2/2}|\alpha\rangle
=
\sum_{n=0}^{\infty}
\frac{\alpha^n}{\sqrt{n!}}|n\rangle,
$$

where the bar labels the unnormalized state. Its inner product is

$$
\langle\bar\alpha|\bar\beta\rangle
=
\sum_{n=0}^{\infty}
\frac{(\alpha^*\beta)^n}{n!}
=
e^{\alpha^*\beta}.
$$

Define

$$
d\mu(\alpha)
=
\frac{d^2\alpha}{\pi}e^{-|\alpha|^2},
\qquad
\int_{\mathbb C}
d\mu(\alpha)\,
|\bar\alpha\rangle\langle\bar\alpha|
=
\hat I.
$$

For a state $|\psi\rangle$, its coherent-state representation is

$$
\psi(\alpha^*)
=
\langle\bar\alpha|\psi\rangle
=
\sum_{n=0}^{\infty}
\frac{(\alpha^*)^n}{\sqrt{n!}}
\langle n|\psi\rangle.
$$

The function $\psi(\alpha^*)$ is analytic in $\alpha^*$, and

$$
|\psi\rangle
=
\int_{\mathbb C}
d\mu(\alpha)\,
|\bar\alpha\rangle\psi(\alpha^*).
$$

Taking the adjoint of the number-state expansion gives

$$
\langle\bar\alpha|
=
\sum_{n=0}^{\infty}
\frac{(\alpha^*)^n}{\sqrt{n!}}\langle n|.
$$

The ladder relations then give

$$
\langle\bar\alpha|\hat a^\dagger
=
\alpha^*\langle\bar\alpha|,
\qquad
\langle\bar\alpha|\hat a
=
\frac{\partial}{\partial\alpha^*}
\langle\bar\alpha|.
$$

Therefore,

$$
\langle\bar\alpha|\hat a^\dagger|\psi\rangle
=
\alpha^*\psi(\alpha^*),
\qquad
\langle\bar\alpha|\hat a|\psi\rangle
=
\frac{\partial}{\partial\alpha^*}\psi(\alpha^*).
$$

Hence,

$$
\hat a^\dagger
\longrightarrow
\alpha^*,
\qquad
\hat a
\longrightarrow
\frac{\partial}{\partial\alpha^*}.
$$

For a normally ordered operator, these substitutions directly give its action on $\psi(\alpha^*)$.

##### Uncertainty Relation

The coherent-state eigenvalue equation gives

$$
\langle\alpha|\hat a|\alpha\rangle
=
\alpha,
\qquad
\langle\alpha|\hat a^\dagger|\alpha\rangle
=
\alpha^*.
$$

Therefore,

$$
\langle\alpha|\hat x|\alpha\rangle
=
\frac{x_0}{\sqrt2}
(\alpha+\alpha^*),
$$

$$
\langle\alpha|\hat p|\alpha\rangle
=
\frac{p_0}{\sqrt2 i}
(\alpha-\alpha^*).
$$

The required second moments follow from

$$
\langle\alpha|\hat a^2|\alpha\rangle
=
\alpha^2,
\qquad
\langle\alpha|\hat a^{\dagger 2}|\alpha\rangle
=
\alpha^{*2},
$$

$$
\langle\alpha|\hat a^\dagger\hat a|\alpha\rangle
=
|\alpha|^2,
\qquad
\langle\alpha|\hat a\hat a^\dagger|\alpha\rangle
=
|\alpha|^2+1.
$$

Defining the variances by

$$
(\Delta x)^2
=
\langle\hat x^2\rangle
-
\langle\hat x\rangle^2,
\qquad
(\Delta p)^2
=
\langle\hat p^2\rangle
-
\langle\hat p\rangle^2,
$$

gives

$$
(\Delta x)^2
=
\frac{x_0^2}{2},
\qquad
(\Delta p)^2
=
\frac{p_0^2}{2}.
$$

Since $x_0p_0=\hbar$,

$$
\Delta x\,\Delta p
=
\frac{\hbar}{2}.
$$

Thus, every coherent state saturates the position-momentum uncertainty relation.

##### Squeezed Coherent States

Let $\mu,\nu\in\mathbb R$ and define

$$
\hat b
=
\mu\hat a+\nu\hat a^\dagger,
\qquad
\hat b^\dagger
=
\mu\hat a^\dagger+\nu\hat a.
$$

The commutation relation is

$$
[\hat b,\hat b^\dagger]
=
(\mu^2-\nu^2)\hat I.
$$

Requiring $[\hat b,\hat b^\dagger]=\hat I$ gives

$$
\mu^2-\nu^2
=
1.
$$

The inverse transformation is

$$
\hat a
=
\mu\hat b-\nu\hat b^\dagger,
\qquad
\hat a^\dagger
=
\mu\hat b^\dagger-\nu\hat b.
$$

Therefore,

$$
\hat x
=
\frac{x_0}{\sqrt2(\mu+\nu)}
(\hat b+\hat b^\dagger),
\qquad
\hat p
=
\frac{p_0}{\sqrt2 i(\mu-\nu)}
(\hat b-\hat b^\dagger).
$$

Let $|\beta\rangle_b$ be a normalized eigenstate of $\hat b$:

$$
\hat b|\beta\rangle_b
=
\beta|\beta\rangle_b,
\qquad
\beta\in\mathbb C.
$$

Writing

$$
\langle\hat O\rangle_b
=
{}_b\langle\beta|\hat O|\beta\rangle_b,
$$

the eigenvalue equation gives

$$
\langle\hat b\rangle_b
=
\beta,
\qquad
\langle\hat b^\dagger\rangle_b
=
\beta^*,
$$

$$
\langle\hat b^\dagger\hat b\rangle_b
=
|\beta|^2,
\qquad
\langle\hat b\hat b^\dagger\rangle_b
=
|\beta|^2+1.
$$

The first moments are

$$
\langle\hat x\rangle_b
=
\frac{x_0}{\sqrt2(\mu+\nu)}
(\beta+\beta^*),
\qquad
\langle\hat p\rangle_b
=
\frac{p_0}{\sqrt2 i(\mu-\nu)}
(\beta-\beta^*).
$$

The second moments are

$$
\langle\hat x^2\rangle_b
=
\frac{x_0^2}{2(\mu+\nu)^2}
\left[
(\beta+\beta^*)^2+1
\right],
$$

$$
\langle\hat p^2\rangle_b
=
\frac{p_0^2}{2(\mu-\nu)^2}
\left[
1-(\beta-\beta^*)^2
\right].
$$

Therefore,

$$
(\Delta x)^2
=
\langle\hat x^2\rangle_b
-
\langle\hat x\rangle_b^2
=
\frac{x_0^2}{2(\mu+\nu)^2},
$$

$$
(\Delta p)^2
=
\langle\hat p^2\rangle_b
-
\langle\hat p\rangle_b^2
=
\frac{p_0^2}{2(\mu-\nu)^2}.
$$

Since $(\mu+\nu)(\mu-\nu)=1$,

$$
\Delta x\,\Delta p
=
\frac{\hbar}{2}.
$$

Introducing the real squeezing parameter $r$ by

$$
\mu
=
\cosh r,
\qquad
\nu
=
\sinh r,
$$

gives

$$
(\Delta x)^2
=
\frac{x_0^2}{2}e^{-2r},
\qquad
(\Delta p)^2
=
\frac{p_0^2}{2}e^{2r}.
$$

Thus, one quadrature is compressed while the other is enlarged, with their uncertainty product unchanged.
