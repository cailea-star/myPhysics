### Lie Groups

Let the elements of $G$ be locally parametrized by

$$
g(\boldsymbol\alpha)\in G,\qquad \boldsymbol\alpha=(\alpha_1,\ldots,\alpha_n)\in U\subseteq\mathbb R^n,\qquad g(\boldsymbol 0)=e.
$$

The number of independent parameters is

$$
n=\dim G.
$$

##### Group Parameters and Composition Functions

Define the composition function $\boldsymbol f$ by

$$
g(\boldsymbol\alpha)g(\boldsymbol\beta)=g\!\left(\boldsymbol f(\boldsymbol\alpha,\boldsymbol\beta)\right).
$$

- Closure

$$
\boldsymbol f:U\times U\to U.
$$

- Associativity

$$
\boldsymbol f\!\left(\boldsymbol\alpha,\boldsymbol f(\boldsymbol\beta,\boldsymbol\gamma)\right)=\boldsymbol f\!\left(\boldsymbol f(\boldsymbol\alpha,\boldsymbol\beta),\boldsymbol\gamma\right).
$$

- Identity

$$
\boldsymbol f(\boldsymbol\alpha,\boldsymbol 0)=\boldsymbol f(\boldsymbol 0,\boldsymbol\alpha)=\boldsymbol\alpha.
$$

- Inverse

Let $g(\overline{\boldsymbol\alpha})=g(\boldsymbol\alpha)^{-1}$. Then

$$
\boldsymbol f(\boldsymbol\alpha,\overline{\boldsymbol\alpha})=\boldsymbol f(\overline{\boldsymbol\alpha},\boldsymbol\alpha)=\boldsymbol 0.
$$

A Lie group is a smooth manifold whose multiplication and inversion maps are smooth:

$$
m(g,h)=gh,\qquad \iota(g)=g^{-1}.
$$

##### Local Structure and Infinitesimal Elements

An element $g(\boldsymbol\alpha)$ is infinitesimal when

$$
\|\boldsymbol\alpha\|\ll1.
$$

To first order,

$$
\boldsymbol f(\boldsymbol\alpha,\boldsymbol\beta)=\boldsymbol\alpha+\boldsymbol\beta+O\!\left(\|\boldsymbol\alpha\|\,\|\boldsymbol\beta\|\right),
$$

$$
\overline{\boldsymbol\alpha}=-\boldsymbol\alpha+O\!\left(\|\boldsymbol\alpha\|^2\right).
$$

Thus infinitesimal parameters add at first order; noncommutativity appears at higher order.

A smooth parameter curve through the identity,

$$
\boldsymbol\alpha(t)\in U,\qquad \boldsymbol\alpha(0)=\boldsymbol 0,
$$

defines a tangent vector

$$
X=\left.\frac{d g(\boldsymbol\alpha(t))}{dt}\right|_{t=0}\in T_eG,\qquad \dim T_eG=\dim G=n.
$$

The tangent space $T_eG$ contains all infinitesimal directions of $G$.

### Infinitesimal Operators and Generators

##### Infinitesimal Operators on Functions

Let

$$
D^{(X)}:G\to GL_d(\mathbb C)
$$

be the coordinate-space representation of $G$.

- Coordinate and scalar-function transformations

$$
\underline x'=D^{(X)}(g(\boldsymbol\alpha))\underline x,\qquad \psi'=P_{g(\boldsymbol\alpha)}\psi,\qquad \psi'(\underline x')=\psi(\underline x).
$$

Therefore,

$$
\boxed{\left(P_{g(\boldsymbol\alpha)}\psi\right)(\underline x)=\psi\!\left(D^{(X)}(g(\boldsymbol\alpha))^{-1}\underline x\right)}.
$$

- Infinitesimal differential operators

Use $a,b,c$ for group-parameter indices and $i,j,k$ for coordinate indices. Expanding at $\boldsymbol\alpha=\boldsymbol 0$,

$$
\left(P_{g(\boldsymbol\alpha)}\psi\right)(\underline x)=\psi(\underline x)-\sum_{a=1}^{n}\sum_{i,j=1}^{d}\alpha_a\left.\frac{\partial D_{ji}^{(X)}(g(\boldsymbol\alpha))}{\partial\alpha_a}\right|_{\boldsymbol\alpha=\boldsymbol 0}x_i\frac{\partial\psi}{\partial x_j}+O\!\left(\|\boldsymbol\alpha\|^2\right).
$$

Thus,

$$
\boxed{\widehat T_a=-i\sum_{i,j=1}^{d}\left.\frac{\partial D_{ji}^{(X)}(g(\boldsymbol\alpha))}{\partial\alpha_a}\right|_{\boldsymbol\alpha=\boldsymbol 0}x_i\frac{\partial}{\partial x_j},\qquad P_{g(\boldsymbol\alpha)}=I-i\sum_{a=1}^{n}\alpha_a\widehat T_a+O\!\left(\|\boldsymbol\alpha\|^2\right)}.
$$

##### Generators on Invariant Function Spaces

Let

$$
F=\operatorname{span}\{\psi_1,\ldots,\psi_m\},\qquad P_gF=F,\qquad \forall g\in G.
$$

- Matrix representation

$$
\left(P_{g(\boldsymbol\alpha)}\psi_\mu\right)(\underline x)=\sum_{\nu=1}^{m}\psi_\nu(\underline x)D_{\nu\mu}^{(F)}(g(\boldsymbol\alpha)).
$$

- Generator matrices

$$
\boxed{T_a^{(F)}=i\left.\frac{\partial D^{(F)}(g(\boldsymbol\alpha))}{\partial\alpha_a}\right|_{\boldsymbol\alpha=\boldsymbol 0},\qquad D^{(F)}(g(\boldsymbol\alpha))=I_m-i\sum_{a=1}^{n}\alpha_aT_a^{(F)}+O\!\left(\|\boldsymbol\alpha\|^2\right)}.
$$

Comparing the operator and matrix expansions gives

$$
\boxed{\widehat T_a\psi_\mu=\sum_{\nu=1}^{m}\psi_\nu[T_a^{(F)}]_{\nu\mu}}.
$$

- One-parameter subgroups

Let $\boldsymbol\alpha(t)$ describe a one-parameter subgroup:

$$
\boldsymbol\alpha(0)=\boldsymbol 0,\qquad g(\boldsymbol\alpha(t_1))g(\boldsymbol\alpha(t_2))=g(\boldsymbol\alpha(t_1+t_2)).
$$

Then

$$
P_{g(\boldsymbol\alpha(t))}=e^{-it\widehat T},\qquad D^{(F)}(g(\boldsymbol\alpha(t)))=e^{-itT^{(F)}}.
$$

### Lie Algebras

Let $G\subseteq GL(V)$ be a matrix Lie group, with elements denoted by $R,S$. Near the identity,

$$
R(\boldsymbol\alpha)=I-i\sum_{a=1}^{n}\alpha_aT_a+O\!\left(\|\boldsymbol\alpha\|^2\right).
$$

The corresponding Lie algebra is

$$
\mathfrak g=T_eG=\operatorname{span}_{\mathbb R}\{-iT_1,\ldots,-iT_n\},\qquad n=\dim G.
$$

##### Commutation Relations and Structure Constants

- Closure

$$
\boxed{[T_a,T_b]=i\sum_{c=1}^{n}f_{ab}^{\phantom{ab}c}T_c},
$$

where $f_{ab}^{\phantom{ab}c}\in\mathbb R$ are the structure constants of $\mathfrak g$.

- Antisymmetry

$$
f_{ab}^{\phantom{ab}c}=-f_{ba}^{\phantom{ba}c}.
$$

- Jacobi identity

$$
[T_a,[T_b,T_c]]+[T_b,[T_c,T_a]]+[T_c,[T_a,T_b]]=0.
$$

- Differential and matrix realizations

Every realization preserves the same commutation relations:

$$
[\widehat T_a,\widehat T_b]=i\sum_{c=1}^{n}f_{ab}^{\phantom{ab}c}\widehat T_c,\qquad [T_a^{(F)},T_b^{(F)}]=i\sum_{c=1}^{n}f_{ab}^{\phantom{ab}c}T_c^{(F)}.
$$

##### Adjoint Representation

For $R\in G$, define the conjugation map

$$
C_R(S)=RSR^{-1}.
$$

Its differential at the identity defines the adjoint action:

$$
\operatorname{Ad}_R=(dC_R)_e:\mathfrak g\to\mathfrak g.
$$

- Adjoint representation matrices

In the generator basis,

$$
\boxed{\operatorname{Ad}_R(T_b)=RT_bR^{-1}=\sum_{c=1}^{n}T_cD_{cb}^{(\mathrm{Ad})}(R)}.
$$

The generators $(T_1,\ldots,T_n)$ form a basis of the carrier space, and $D^{(\mathrm{Ad})}(R)$ is the presentation of $\operatorname{Ad}_R$ in this basis:

$$
D^{(\mathrm{Ad})}(RS)=D^{(\mathrm{Ad})}(R)D^{(\mathrm{Ad})}(S).
$$

Since the carrier space has dimension $n=\dim G$,

$$
D^{(\mathrm{Ad})}(R)\in GL_n(\mathbb R),\qquad T_a^{(\mathrm{Ad})}\in\mathbb C^{n\times n}.
$$

- Infinitesimal adjoint action

For $R(\boldsymbol\alpha)$ near the identity,

$$
R(\boldsymbol\alpha)T_bR(\boldsymbol\alpha)^{-1}=T_b-i\sum_{a=1}^{n}\alpha_a[T_a,T_b]+O\!\left(\|\boldsymbol\alpha\|^2\right)=T_b+\sum_{a,c=1}^{n}\alpha_af_{ab}^{\phantom{ab}c}T_c+O\!\left(\|\boldsymbol\alpha\|^2\right).
$$

Define the adjoint generators by

$$
D^{(\mathrm{Ad})}(R(\boldsymbol\alpha))=I_n-i\sum_{a=1}^{n}\alpha_aT_a^{(\mathrm{Ad})}+O\!\left(\|\boldsymbol\alpha\|^2\right).
$$

Therefore,

$$
\boxed{[T_a^{(\mathrm{Ad})}]_{cb}=if_{ab}^{\phantom{ab}c},\qquad -iT_a^{(\mathrm{Ad})}\in\mathbb R^{n\times n}}.
$$

### Integration on Compact Lie Groups

For a compact Lie group $G$, the finite-group average generalizes to a normalized group integral:

$$
\frac{1}{|G|}\sum_{R\in G}F(R)\quad\longrightarrow\quad\int_GdR\,F(R).
$$

If $R=R(\boldsymbol\alpha)$ is parametrized by $\boldsymbol\alpha\in U\subseteq\mathbb R^n$, then

$$
dR=W(\boldsymbol\alpha)\,d^n\boldsymbol\alpha,\qquad \int_GdR\,F(R)=\int_Ud^n\boldsymbol\alpha\,W(\boldsymbol\alpha)F(R(\boldsymbol\alpha)),
$$

where $W(\boldsymbol\alpha)$ is the weight function.

- Normalization

$$
\boxed{\int_GdR=1}.
$$

- Left and right invariance

For every $S\in G$,

$$
\boxed{\int_GdR\,F(R)=\int_GdR\,F(SR)=\int_GdR\,F(RS)}.
$$

The measure $dR$ satisfying these properties is the normalized Haar measure on $G$.

Thus, the rearrangement of a finite-group sum extends directly to compact-group integration.

### Linear Representations of Compact Lie Groups

Let

$$
D:G\to GL_m(\mathbb C)
$$

be a finite-dimensional continuous representation of a compact Lie group $G$.

- Unitary representations

Every finite-dimensional complex representation is equivalent to a unitary representation:

$$
D(R)^\dagger D(R)=I_m.
$$

Two equivalent unitary representations are related by a unitary similarity transformation. A real representation is equivalent to a real orthogonal representation.

- Orthogonality and completeness

Let $D^{(i)}$ run over all inequivalent irreducible unitary representations, with $d_i=\dim D^{(i)}$.

$$
\boxed{\int_GdR\,D_{\mu\nu}^{(i)}(R)^*D_{\rho\sigma}^{(j)}(R)=\frac{1}{d_i}\delta_{ij}\delta_{\mu\rho}\delta_{\nu\sigma}}.
$$

The completeness relation is

$$
\boxed{\sum_i d_i\sum_{\mu,\nu=1}^{d_i}D_{\mu\nu}^{(i)}(R)^*D_{\mu\nu}^{(i)}(S)=\delta_G(R^{-1}S)}.
$$

The group delta function is defined by

$$
\int_GdS\,\delta_G(R^{-1}S)F(S)=F(R).
$$

Equivalently,

$$
\sum_i d_i\chi_i(R^{-1}S)=\delta_G(R^{-1}S),\qquad \chi_i(R)=\operatorname{tr}D^{(i)}(R).
$$

Thus,

$$
\left\{\sqrt{d_i}\,D_{\mu\nu}^{(i)}\right\}_{i,\mu,\nu}
$$

is an orthonormal complete basis of $L^2(G)$.

- Characters and multiplicities

$$
\boxed{\int_GdR\,\chi_i(R)^*\chi_j(R)=\delta_{ij},\qquad \eta_i=\int_GdR\,\chi_i(R)^*\chi_D(R)}.
$$

Thus, the main results for finite-group representations extend to compact Lie groups by replacing normalized sums with Haar integrals.

### Global Properties of Lie Groups

##### Connected Components

A Lie group $G$ is connected if any two points in its group space can be joined by a continuous path lying entirely in $G$.

Let $G_0$ be the connected component containing the identity $e$. Then

$$
G\text{ is connected}\Longleftrightarrow G=G_0.
$$

- Identity component

$$
G_0\trianglelefteq G.
$$

Every other connected component is a coset of $G_0$:

$$
G=\bigsqcup_{\lambda}R_\lambda G_0,
$$

where $R_\lambda$ is a representative of each component.

- Infinitesimal elements

Every element of a connected Lie group can be written as a product of elements near the identity:

$$
R=A_1A_2\cdots A_m,\qquad A_i\simeq I.
$$

For a disconnected Lie group, one representative from each connected component is also required.

For example,

$$
SO(3),SU(2)\text{ are connected},\qquad O(3)\text{ is disconnected}.
$$

##### Path Classes and Covering Groups

Consider continuous paths from the identity $e$ to the same element $R$.

- Simply connected spaces

A connected group space is simply connected if any two such paths can be continuously deformed into each other while their endpoints remain fixed.

- Multiply connected spaces

If the paths separate into distinct deformation classes, these classes describe the connectivity of the group space.

For example,

$$
SU(2)\text{ is simply connected},\qquad SO(3)\text{ has two path classes}.
$$

- Covering groups

For a connected Lie group $G$, there exists a simply connected Lie group $\widetilde G$ and a covering homomorphism

$$
\pi:\widetilde G\to G.
$$

They have the same local structure:

$$
\operatorname{Lie}(\widetilde G)\simeq\operatorname{Lie}(G).
$$

If each element of $G$ has exactly $n$ preimages in $\widetilde G$, then $\widetilde G$ is an $n$-fold covering group of $G$.

Their global relation is

$$
\boxed{G\simeq\widetilde G/\ker\pi}.
$$

A representation $\widetilde D$ of $\widetilde G$ defines a single-valued representation of $G$ if and only if

$$
\widetilde D(k)=I,\qquad \forall k\in\ker\pi.
$$

Otherwise, it may be regarded as a multivalued representation of $G$.

- $SU(2)$ and $SO(3)$

$$
\pi:SU(2)\to SO(3),\qquad \ker\pi=\{I,-I\},
$$

$$
SO(3)\simeq SU(2)/\{I,-I\}.
$$

Thus, $SU(2)$ is the double covering group of $SO(3)$.

##### Compactness of Group Spaces

A Lie group $G$ is compact if its underlying manifold is compact.

$$
SO(3),SU(2)\text{ are compact},
$$

whereas the translation and Lorentz groups are noncompact.

For a compact Lie group, the invariant measure has finite total volume and can be normalized:

$$
\int_GdR=1.
$$
