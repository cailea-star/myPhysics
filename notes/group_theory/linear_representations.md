### Linear Spaces and Linear Operators

##### Vector Spaces

Let $V$ be an $m$-dimensional complex vector space. For $\boldsymbol u,\boldsymbol v,\boldsymbol w\in V$ and $a,b\in\mathbb C$,

- Addition

$$
\boldsymbol v+\boldsymbol w\in V,\qquad \boldsymbol v+\boldsymbol w=\boldsymbol w+\boldsymbol v,\qquad (\boldsymbol u+\boldsymbol v)+\boldsymbol w=\boldsymbol u+(\boldsymbol v+\boldsymbol w).
$$

- Zero vector and additive inverse

$$
\boldsymbol v+\boldsymbol 0=\boldsymbol v,\qquad \boldsymbol v+(-\boldsymbol v)=\boldsymbol 0.
$$

- Scalar multiplication

$$
a\boldsymbol v\in V,\qquad a(b\boldsymbol v)=(ab)\boldsymbol v,\qquad 1\boldsymbol v=\boldsymbol v.
$$

- Distributivity

$$
a(\boldsymbol v+\boldsymbol w)=a\boldsymbol v+a\boldsymbol w,\qquad (a+b)\boldsymbol v=a\boldsymbol v+b\boldsymbol v.
$$

Let $\mathcal E=(\boldsymbol e_1,\ldots,\boldsymbol e_m)$ be a basis of $V$. Unless otherwise stated, coordinates and operator matrices are written in this basis.

- Basis and coordinates

$$
\boldsymbol v=\sum_{\mu=1}^{m}\boldsymbol e_\mu v_\mu,\qquad \underline v=(v_1,\ldots,v_m)^{\mathsf T}.
$$

##### Linear Operators

A map $A:V\to V$ is linear if

$$
A(a\boldsymbol v+b\boldsymbol w)=aA(\boldsymbol v)+bA(\boldsymbol w).
$$

- Matrix representation

The matrix $[A]$ is defined by

$$
A\boldsymbol e_\mu=\sum_{\nu=1}^{m}\boldsymbol e_\nu[A]_{\nu\mu}.
$$

For $\boldsymbol b=A\boldsymbol v$,

$$
\underline b=[A]\underline v.
$$

- Operator products

$$
[AB]=[A][B].
$$

##### Similarity Transformations

Let $\mathcal E'=(\boldsymbol e'_1,\ldots,\boldsymbol e'_m)$ be another basis, related to $\mathcal E$ by an invertible matrix $S$.

- Basis transformation

$$
\boldsymbol e'_\mu=\sum_{\nu=1}^{m}\boldsymbol e_\nu S_{\nu\mu}.
$$

- Coordinate transformation

$$
\underline v'=S^{-1}\underline v.
$$

- Operator transformation

$$
\boxed{[A]'=S^{-1}[A]S}.
$$

### Linear Representations of Groups

##### General Linear Representations

Let $G$ be a group with identity $e$. A linear representation of $G$ on $V$ is a homomorphism

$$
P:G\to GL(V),\qquad g\mapsto P_g,
$$

where $GL(V)$ is the group of invertible linear operators on $V$.

- Representation law

$$
P_e=I,\qquad P_{g_1}P_{g_2}=P_{g_1g_2}.
$$

- Invariant carrier space

A representation on $V$ requires

$$
\boxed{P_gV=V,\qquad \forall g\in G}.
$$

Hence, for every basis vector $\boldsymbol e_\mu$,

$$
P_g\boldsymbol e_\mu\in V.
$$

- Basis-vector transformation

The transformed basis vector can therefore be expanded as

$$
\boxed{P_g\boldsymbol e_\mu=\sum_{\nu=1}^{m}\boldsymbol e_\nu D_{\nu\mu}(g)}.
$$

The coefficients $D_{\nu\mu}(g)$ form the matrix $D(g)$ of $P_g$ in the chosen basis.

- Coordinate transformation

For $\boldsymbol w=P_g\boldsymbol v$,

$$
\boldsymbol w=\sum_{\nu=1}^{m}\boldsymbol e_\nu\sum_{\mu=1}^{m}D_{\nu\mu}(g)v_\mu.
$$

Therefore,

$$
\boxed{\underline w=D(g)\underline v}.
$$

- Matrix representation

The operator product law gives

$$
\boxed{D(g_1g_2)=D(g_1)D(g_2)}.
$$

Thus,

$$
D:G\to GL_m(\mathbb C),\qquad g\mapsto D(g)
$$

is the matrix realization of $P$, with

$$
D(e)=I_m,\qquad D(g^{-1})=D(g)^{-1}.
$$

- Kernel and faithfulness

$$
\ker P=\{g\in G:P_g=I\}=\{g\in G:D(g)=I_m\}.
$$

$$
P\text{ is faithful}\Longleftrightarrow\ker P=\{e\}.
$$

- Trivial representation

$$
D_{\mathrm{triv}}(g)=1,\qquad \forall g\in G.
$$

##### Characters

- Character

$$
\chi_D(g)=\operatorname{tr}D(g),\qquad \chi_D(e)=m.
$$

- Conjugacy invariance

For $g'=xgx^{-1}$,

$$
D(g')=D(x)D(g)D(x)^{-1}.
$$

Therefore,

$$
\boxed{\chi_D(g')=\chi_D(g)}.
$$

##### Function-Space Realizations

Let $X=\mathbb C^n$ carry an $n$-dimensional matrix representation $D_X$ of $G$.

- Coordinate transformation

$$
\underline x'=D_X(g)\underline x.
$$

- Scalar-function transformation

For a scalar function $\psi:X\to\mathbb C$,

$$
\boxed{\psi'(\underline x')=\psi(\underline x)}.
$$

Define the induced operator by

$$
P_g\psi=\psi'.
$$

Equivalently,

$$
(P_g\psi)(\underline x)=\psi\bigl(D_X(g)^{-1}\underline x\bigr).
$$

- Invariant function space

Let $F$ be an $m$-dimensional function space with basis $(\psi_1,\ldots,\psi_m)$. If

$$
P_gF=F,\qquad \forall g\in G,
$$

then

$$
(P_g\psi_\mu)(\underline x)=\sum_{\nu=1}^{m}\psi_\nu(\underline x)[D_F(g)]_{\nu\mu}.
$$

The matrices $D_F(g)$ form the matrix representation of $P_g$ on $F$.

### Equivalent and Reducible Representations

##### Equivalent Representations

- Equivalence

Two $m$-dimensional matrix representations $D$ and $D'$ are equivalent if an invertible matrix $S$, independent of $g$, satisfies

$$
\boxed{D'(g)=S^{-1}D(g)S,\qquad \forall g\in G}.
$$

- Characters

$$
\chi_{D'}(g)=\operatorname{tr}\!\left(S^{-1}D(g)S\right)=\chi_D(g).
$$

##### Unitary Representations

- Unitary representation

A representation is unitary if

$$
P_g^\dagger=P_{g^{-1}},\qquad D(g)^\dagger D(g)=I_m,
$$

where $\dagger$ denotes the adjoint, represented by the conjugate transpose in an orthonormal basis.

- Unitarization

For a finite group, define

$$
H=\sum_{g\in G}D(g)^\dagger D(g).
$$

The matrix $H$ is positive-definite Hermitian. By the rearrangement theorem,

$$
D(s)^\dagger HD(s)=\sum_{g\in G}D(gs)^\dagger D(gs)=H.
$$

Writing $H=T^\dagger T$, define

$$
D_u(g)=TD(g)T^{-1}
$$

Then

$$
D_u(g)^\dagger D_u(g)=(T^{-1})^\dagger D(g)^\dagger HD(g)T^{-1}=I_m.
$$

Therefore,

$$
\boxed{\text{Every finite-dimensional complex representation of a finite group is equivalent to a unitary representation.}}
$$

##### Reducible Representations and Invariant Subspaces

Let $P$ act on $V$.

- Reducibility and invariant subspaces

A subspace $W\subseteq V$ is invariant if

$$
P_gW=W,\qquad \forall g\in G.
$$

The representation is reducible if it has a nonzero proper invariant subspace:

$$
\boxed{P\text{ is reducible}\Longleftrightarrow\exists W,\ \{0\}\subsetneq W\subsetneq V,\quad P_gW=W\ \forall g\in G}.
$$

If no such $W$ exists, $P$ is irreducible.

In a basis adapted to $W$, the matrices $D(g)$ have the form

$$
S^{-1}D(g)S=\begin{pmatrix}D_1(g)&M(g)\\0&D_2(g)\end{pmatrix}.
$$

- Complete reducibility

A representation is completely reducible if its matrices decompose into irreducible representations $D^{(1)},\ldots,D^{(r)}$:

$$
\boxed{S^{-1}D(g)S=D^{(1)}(g)\oplus\cdots\oplus D^{(r)}(g),\qquad r\ge1}.
$$

For a finite group, first choose an equivalent unitary realization. If $W$ is invariant, then for $\boldsymbol v\in W^\perp$ and $\boldsymbol w\in W$,

$$
\langle P_g\boldsymbol v,\boldsymbol w\rangle=\langle\boldsymbol v,P_{g^{-1}}\boldsymbol w\rangle=0.
$$

Hence,

$$
P_gW^\perp=W^\perp,\qquad V=W\oplus W^\perp.
$$

Repeating this decomposition yields the direct-sum form above.

### Regular Representations

Let $G$ be a finite group of order $n_G=|G|$.

##### Group Algebra

- Natural basis

$$
\mathbb C[G]=\left\{\sum_{h\in G}c_h\boldsymbol e_h:c_h\in\mathbb C\right\},\qquad \dim\mathbb C[G]=n_G.
$$

- Multiplication

$$
\boldsymbol e_g\boldsymbol e_h=\boldsymbol e_{gh},
$$

extended linearly to $\mathbb C[G]$.

##### Left and Right Regular Representations

For $g,h\in G$, define

$$
\boxed{L_g\boldsymbol e_h=\boldsymbol e_{gh},\qquad R_g\boldsymbol e_h=\boldsymbol e_{hg^{-1}}}.
$$

The inverse in the right action preserves the multiplication order:

$$
L_{g_1}L_{g_2}=L_{g_1g_2},\qquad R_{g_1}R_{g_2}=R_{g_1g_2}.
$$

- Matrix elements

Let $D_L(g)$ and $D_R(g)$ be the matrices of $L_g$ and $R_g$ in the natural basis. For $a,b\in G$,

$$
\boxed{[D_L(g)]_{ab}=\delta_{a,gb},\qquad [D_R(g)]_{ab}=\delta_{ag,b}},
$$

where $\delta$ is the Kronecker delta.

- Faithfulness

$$
L_g=I\Longleftrightarrow g=e,\qquad R_g=I\Longleftrightarrow g=e.
$$

- Inverse and transpose

$$
[D_L(g)^{-1}]_{ab}=[D_L(g^{-1})]_{ab}=\delta_{a,g^{-1}b}=\delta_{b,ga}=[D_L(g)]_{ba}=[D_L(g)^{\mathsf T}]_{ab}.
$$

$$
[D_R(g)^{-1}]_{ab}=[D_R(g^{-1})]_{ab}=\delta_{ag^{-1},b}=\delta_{a,bg}=\delta_{bg,a}=[D_R(g)]_{ba}=[D_R(g)^{\mathsf T}]_{ab}.
$$

Thus $D_L(g)$ and $D_R(g)$ are real orthogonal permutation matrices.

- Commuting actions

For every $\boldsymbol e_c$,

$$
L_gR_h\boldsymbol e_c=\boldsymbol e_{gch^{-1}}=R_hL_g\boldsymbol e_c.
$$

Therefore,

$$
\boxed{L_gR_h=R_hL_g}.
$$

- Character

For $g\ne e$, neither regular action fixes a basis vector. Hence

$$
\boxed{\chi_{\mathrm{reg}}(e)=n_G,\qquad \chi_{\mathrm{reg}}(g)=0\quad(g\ne e)}.
$$

##### Construction from the Cayley Table

Use the same ordering of $G$ for the Cayley table and the natural basis.

- $D_L(g)$: read row $g$. If column $b$ contains $a=gb$, place $1$ in row $a$, column $b$.

- $D_R(g)$: read column $g^{-1}$. If row $b$ contains $a=bg^{-1}$, place $1$ in row $a$, column $b$.

All remaining entries are zero.

### Orthogonality and Completeness Relations

Let $G$ be a finite group of order $n_G=|G|$. Define the group-function space

$$
\mathcal H_G=\{f:G\to\mathbb C\},\qquad \dim\mathcal H_G=n_G,\qquad \langle f_1,f_2\rangle_G=\sum_{g\in G}f_1(g)^*f_2(g).
$$

##### Schur's Lemma

Let $P^{(1)}$ and $P^{(2)}$ be irreducible complex representations on $V_1$ and $V_2$. A linear map $A:V_2\to V_1$ is an intertwining operator if

$$
P_g^{(1)}A=AP_g^{(2)},\qquad \forall g\in G.
$$

- Kernel

$$
\ker A=\{\boldsymbol v_2\in V_2:A\boldsymbol v_2=0\}.
$$

For $\boldsymbol v_2\in\ker A$,

$$
A\!\left(P_g^{(2)}\boldsymbol v_2\right)=P_g^{(1)}\!\left(A\boldsymbol v_2\right)=0,\Longrightarrow P_g^{(2)}\boldsymbol v_2\in\ker A.
$$

Hence,

$$
P_g^{(2)}\ker A=\ker A.
$$

- Image

$$
\operatorname{im}A=\{A\boldsymbol v_2:\boldsymbol v_2\in V_2\}.
$$

For $\boldsymbol v_1=A\boldsymbol v_2\in\operatorname{im}A$,

$$
P_g^{(1)}\boldsymbol v_1 = P_g^{(1)} \Big(A {\boldsymbol v_2}\Big) =A\!\Big(P_g^{(2)}\boldsymbol v_2\Big)\in\operatorname{im}A.
$$

Hence,

$$
P_g^{(1)}\operatorname{im}A=\operatorname{im}A.
$$

Irreducibility gives

$$
\ker A\in\{\{0\},V_2\},\qquad \operatorname{im}A\in\{\{0\},V_1\}.
$$

If $A\ne0$, then $A$ is an isomorphism and

$$
A^{-1}P_g^{(1)}A=P_g^{(2)}.
$$

Therefore,

$$
\boxed{P^{(1)}\text{ and }P^{(2)}\text{ inequivalent}\Longrightarrow A=0}.
$$

Let $P^{(1)}=P^{(2)}=P$ and suppose that $A$ commutes with every $P_g$. Since $V$ is finite-dimensional over $\mathbb C$, choose an eigenvalue $\lambda$ of $A$. Then $A-\lambda I$ is a singular intertwining operator and therefore must be zero. Hence,

$$
\boxed{A=\lambda I}.
$$

##### Matrix-Element Orthogonality

Let $D^{(1)},\ldots,D^{(\ell)}$ be all pairwise inequivalent irreducible unitary matrix representations of $G$, with $D^{(i)}$ of dimension $d_i$. Use Greek indices for $D^{(i)}$ and Latin indices for $D^{(j)}$:

$$
\mu,\nu,\rho=1,\ldots,d_i,\qquad a,b,c=1,\ldots,d_j.
$$

Define the $d_i\times d_j$ matrices

$$
[E_{\mu a}]_{\rho c}=\delta_{\rho\mu}\delta_{ca},\qquad Q_{\mu a}^{(i,j)}=\sum_{g\in G}D^{(i)}(g)^\dagger E_{\mu a}D^{(j)}(g).
$$

Their matrix elements are

$$
[Q_{\mu a}^{(i,j)}]_{\nu b}=\sum_{g\in G}\sum_{\rho=1}^{d_i}\sum_{c=1}^{d_j}[D^{(i)}(g)^\dagger]_{\nu\rho}\delta_{\rho\mu}\delta_{ca}D^{(j)}_{cb}(g)=\sum_{g\in G}D^{(i)}_{\mu\nu}(g)^*D^{(j)}_{ab}(g).
$$

For $s\in G$,

$$
[Q_{\mu a}^{(i,j)}D^{(j)}(s)]_{\nu b}=\sum_{g\in G}\sum_{c=1}^{d_j}D^{(i)}_{\mu\nu}(g)^*D^{(j)}_{ac}(g)D^{(j)}_{cb}(s)=\sum_{g\in G}D^{(i)}_{\mu\nu}(g)^*D^{(j)}_{ab}(gs).
$$

By the rearrangement theorem, set $h=gs$. Then

$$
=\sum_{h\in G}D^{(i)}_{\mu\nu}(hs^{-1})^*D^{(j)}_{ab}(h)=\sum_{h\in G}\sum_{\rho=1}^{d_i}D^{(i)}_{\nu\rho}(s)D^{(i)}_{\mu\rho}(h)^*D^{(j)}_{ab}(h)=[D^{(i)}(s)Q_{\mu a}^{(i,j)}]_{\nu b}.
$$

Thus, by Schur's lemma,

$$
\boxed{Q_{\mu a}^{(i,j)}D^{(j)}(s)=D^{(i)}(s)Q_{\mu a}^{(i,j)}\quad(\forall s\in G)\quad\Longrightarrow\quad[Q_{\mu a}^{(i,j)}]_{\nu b}=q_{\mu a}^{(i)}\delta_{ij}\delta_{\nu b}}.
$$

For $i=j$, cyclicity of the trace and unitarity give

$$
\operatorname{tr}Q_{\mu a}^{(i,i)}=\operatorname{tr}\left(q_{\mu a}^{(i)}I_{d_i}\right)=q_{\mu a}^{(i)}d_i=\sum_{g\in G}\operatorname{tr}\!\left(E_{\mu a}D^{(i)}(g)D^{(i)}(g)^\dagger\right)=n_G\delta_{\mu a}.
$$

Therefore,

$$
\boxed{\sum_{g\in G}D^{(i)}_{\mu\nu}(g)^*D^{(j)}_{ab}(g)=\frac{n_G}{d_i}\delta_{ij}\delta_{\mu a}\delta_{\nu b}}.
$$

Define $\chi_i(g)=\operatorname{tr}D^{(i)}(g)$. Summing over diagonal matrix elements gives

$$
\boxed{\sum_{g\in G}\chi_i(g)^*\chi_j(g)=n_G\delta_{ij}}.
$$

##### Character Decomposition

Let $D$ be a finite-dimensional complex representation, and let $\eta_i\in\mathbb N_0$ denote the multiplicity of $D^{(i)}$ in $D$. Complete reducibility gives

$$
\chi_D(g)=\sum_{i=1}^{\ell}\eta_i\chi_i(g).
$$

Character orthogonality gives

$$
\eta_i=\frac{1}{n_G}\sum_{g\in G}\chi_i(g)^*\chi_D(g).
$$

##### Completeness from the Regular Representation

For the regular representation,

$$
\eta_i=\frac{1}{n_G}\chi_i(e)^*\chi_{\mathrm{reg}}(e)=d_i.
$$

Since

$$
\dim D_{\mathrm{reg}}=\sum_{i=1}^{\ell}\eta_i d_i=n_G,
$$

we obtain

$$
\boxed{n_G=\sum_{i=1}^{\ell}d_i^2}.
$$

##### Orthogonality and Completeness of Matrix Elements

- Orthogonality over group elements

$$
\boxed{\sum_{g\in G}D^{(i)}_{\mu\nu}(g)^*D^{(j)}_{ab}(g)=\frac{n_G}{d_i}\delta_{ij}\delta_{\mu a}\delta_{\nu b}}.
$$

- Completeness over representation and matrix indices

$$
\boxed{\sum_{i=1}^{\ell}\sum_{\mu=1}^{d_i}\sum_{\nu=1}^{d_i}\frac{d_i}{n_G}D^{(i)}_{\mu\nu}(h)^*D^{(i)}_{\mu\nu}(g)=\delta_{hg}}.
$$

- Dimension

The number of irreducible matrix elements equals the dimension of the group-function space:

$$
\sum_{i=1}^{\ell}d_i^2=n_G=\dim\mathcal H_G.
$$

##### Orthogonality and Completeness of Characters

Let $D^{(i)}$ be an irreducible representation of $G$. Its character is

$$
\chi_i(g)=\operatorname{tr}D^{(i)}(g).
$$

Let $C_1,\ldots,C_k$ be the conjugacy classes of $G$, with $n_\alpha=|C_\alpha|$. Since characters are constant on conjugacy classes, write $\chi_{i\alpha}=\chi_i(g)$ for $g\in C_\alpha$.

- Orthogonality over conjugacy classes

$$
\boxed{\sum_{\alpha=1}^{k}n_\alpha\chi_{i\alpha}^*\chi_{j\alpha}=n_G\delta_{ij}}.
$$

Matrix-element completeness expands every class function in irreducible matrix elements. Conjugation invariance and Schur's lemma reduce each coefficient matrix to a scalar multiple of the identity, so the irreducible characters form a basis of the class-function space.

- Completeness over irreducible representations

$$
\boxed{\sum_{i=1}^{\ell}\frac{\sqrt{n_\alpha n_\beta}}{n_G}\chi_{i\alpha}^*\chi_{i\beta}=\delta_{\alpha\beta}}.
$$

- Dimension

The number of inequivalent irreducible representations equals the number of conjugacy classes:

$$
\ell=k.
$$

### Direct-Sum Decomposition and Tensor-Product Representations

##### Direct-Sum Decomposition of Representations

Let $D$ be an $N$-dimensional complex representation of the finite group $G$, and let $D^{(1)},\ldots,D^{(k)}$ be all pairwise inequivalent irreducible representations, with

$$
d_i=\dim D^{(i)},\qquad \chi_i(g)=\operatorname{tr}D^{(i)}(g).
$$

- Direct-sum decomposition

Let $\eta_i\in\mathbb N_0$ be the multiplicity of $D^{(i)}$ in $D$. Complete reducibility gives an invertible matrix $S$ such that

$$
\boxed{S^{-1}D(g)S=\bigoplus_{i=1}^{k}\left[I_{\eta_i}\otimes D^{(i)}(g)\right],\qquad \forall g\in G}.
$$

The dimensions satisfy

$$
N=\sum_{i=1}^{k}\eta_i d_i.
$$

- Character decomposition and multiplicity

$$
\chi_D(g)=\sum_{i=1}^{k}\eta_i\chi_i(g).
$$

If the conjugacy class $C_\alpha$ has size $n_\alpha$, character orthogonality gives

$$
\eta_i=\frac{1}{n_G}\sum_{g\in G}\chi_i(g)^*\chi_D(g)=\frac{1}{n_G}\sum_{\alpha=1}^{k}n_\alpha\chi_i(C_\alpha)^*\chi_D(C_\alpha).
$$

- Equivalence criterion

If $\eta_i'$ is the multiplicity of $D^{(i)}$ in $D'$, then

$$
D\simeq D'\Longleftrightarrow\eta_i=\eta_i'\ \forall i\Longleftrightarrow\chi_D(g)=\chi_{D'}(g)\ \forall g\in G.
$$

- Irreducibility criterion

$$
\frac{1}{n_G}\sum_{g\in G}|\chi_D(g)|^2=\sum_{i=1}^{k}\eta_i^2.
$$

Therefore,

$$
D\text{ is irreducible}\Longleftrightarrow\frac{1}{n_G}\sum_{g\in G}|\chi_D(g)|^2=1.
$$

##### Clebsch-Gordan Decomposition

Let $D^{(i)}$ and $D^{(j)}$ be irreducible representations of $G$. Their tensor-product representation satisfies

$$
\left(D^{(i)}\otimes D^{(j)}\right)(g)=D^{(i)}(g)\otimes D^{(j)}(g).
$$

$$
\dim\left(D^{(i)}\otimes D^{(j)}\right)=d_id_j,\qquad \chi_{i\otimes j}(g)=\chi_i(g)\chi_j(g).
$$

Let $\eta_r$ be the multiplicity of $D^{(r)}$. The direct-sum multiplicity formula gives

$$
\eta_r=\frac{1}{n_G}\sum_{g\in G}\chi_r(g)^*\chi_i(g)\chi_j(g),\qquad d_id_j=\sum_{r=1}^{k}\eta_rd_r.
$$

- Product and coupled bases

$$
\lvert i,\mu;j,\nu\rangle=\lvert i,\mu\rangle\otimes\lvert j,\nu\rangle,\qquad \lvert r,\tau,\kappa\rangle,
$$

where $\mu=1,\ldots,d_i$, $\nu=1,\ldots,d_j$, $\tau=1,\ldots,\eta_r$, and $\kappa=1,\ldots,d_r$.

- Clebsch-Gordan coefficients

The Clebsch-Gordan coefficients define the basis transformation:

$$
\lvert r,\tau,\kappa\rangle=\sum_{\mu=1}^{d_i}\sum_{\nu=1}^{d_j}\lvert i,\mu;j,\nu\rangle[U_{\mathrm{CG}}]_{(\mu,\nu),(r,\tau,\kappa)},\qquad [U_{\mathrm{CG}}]_{(\mu,\nu),(r,\tau,\kappa)}=\langle i,\mu;j,\nu\mid r,\tau,\kappa\rangle.
$$

- Block diagonalization

$$
\boxed{U_{\mathrm{CG}}^{-1}\left[D^{(i)}(g)\otimes D^{(j)}(g)\right]U_{\mathrm{CG}}=\bigoplus_{r=1}^{k}\left[I_{\eta_r}\otimes D^{(r)}(g)\right]}.
$$

##### Reduction by the Commutant

Let

$$
D(g)\in\mathbb C^{N\times N}
$$

be a known unitary representation of the finite group $G$.

- Step 1: Construct the commutant

Define

$$
\mathcal C(D)=\left\{A\in\mathbb C^{N\times N}:AD(g)=D(g)A,\ \forall g\in G\right\}.
$$

Complete reducibility gives a unitary matrix $S$ such that

$$
S^\dagger D(g)S=\bigoplus_{i=1}^{k}\left[I_{\eta_i}\otimes D^{(i)}(g)\right].
$$

By Schur's lemma, every $A\in\mathcal C(D)$ has the corresponding form

$$
S^\dagger AS=\bigoplus_{i=1}^{k}\left[A_i\otimes I_{d_i}\right],\qquad A_i\in\mathbb C^{\eta_i\times\eta_i}.
$$

Thus, a commuting matrix acts only on equivalent copies and as the identity within each irreducible representation space. A basis of $\mathcal C(D)$ is obtained by solving

$$
AD(g)-D(g)A=0,\qquad \forall g\in G.
$$

- Step 2: Construct a Hermitian commuting matrix

To obtain an orthogonal eigenspace decomposition, choose a generic Hermitian matrix in the commutant. Let $\{A_a\}$ be a basis of $\mathcal C(D)$, and define

$$
H=\sum_a\left(c_aA_a+c_a^*A_a^\dagger\right).
$$

Then

$$
H^\dagger=H,\qquad HD(g)=D(g)H,
$$

Since $H$ is Hermitian, it is normal and can be unitarily diagonalized. Its eigenspaces are mutually orthogonal and form a complete decomposition of the representation space.

In the irreducible basis,

$$
S^\dagger HS=\bigoplus_{i=1}^{k}\left[H_i\otimes I_{d_i}\right],\qquad H_i\in\mathbb C^{\eta_i\times\eta_i}.
$$

For generic coefficients $c_a$, each $H_i$ has distinct eigenvalues $\lambda_{i\tau}$, with no accidental coincidences between different irreducible blocks. Define

$$
E_{i\tau}=\ker\left(H-\lambda_{i\tau}I_N\right).
$$

The factor $I_{d_i}$ gives

$$
\boxed{\dim E_{i\tau}=d_i,\qquad \operatorname{rank}\left(H-\lambda_{i\tau}I_N\right)=N-d_i}.
$$

For $\boldsymbol v\in E_{i\tau}$,

$$
H\!\left[D(g)\boldsymbol v\right]=D(g)H\boldsymbol v=\lambda_{i\tau}D(g)\boldsymbol v.
$$

Hence,

$$
D(g)E_{i\tau}\subseteq E_{i\tau}.
$$

- Step 3: Construct and identify the irreducible blocks

Choose an orthonormal basis of $E_{i\tau}$ and arrange it as

$$
B_{i\tau}=\begin{pmatrix}\boldsymbol b_1&\cdots&\boldsymbol b_{d_i}\end{pmatrix}\in\mathbb C^{N\times d_i}.
$$

Then

$$
B_{i\tau}^\dagger B_{i\tau}=I_{d_i},\qquad \operatorname{rank}B_{i\tau}=d_i.
$$

The representation restricted to $E_{i\tau}$ is

$$
D^{(i,\tau)}(g)=B_{i\tau}^\dagger D(g)B_{i\tau}\in\mathbb C^{d_i\times d_i}.
$$

Define

$$
\chi_{i\tau}(g)=\operatorname{tr}D^{(i,\tau)}(g).
$$

Irreducibility is checked by

$$
\frac{1}{n_G}\sum_{g\in G}\left|\chi_{i\tau}(g)\right|^2=1.
$$

Two irreducible blocks are equivalent exactly when

$$
\chi_{i\tau}(g)=\chi_{j\sigma}(g),\qquad \forall g\in G.
$$

- Step 4: Process all eigenvalues

For every distinct eigenvalue of $H$, repeat

$$
\lambda_{i\tau}\longrightarrow E_{i\tau}\longrightarrow B_{i\tau}\longrightarrow D^{(i,\tau)}(g)\longrightarrow\chi_{i\tau}(g).
$$

Group blocks with identical characters, retain one representative $D^{(i)}$ from each equivalence class, and let $\eta_i$ be the number of equivalent copies. Arranging all $B_{i\tau}$ as the columns of $S$ gives

$$
\boxed{S^\dagger D(g)S=\bigoplus_{i=1}^{k}\left[I_{\eta_i}\otimes D^{(i)}(g)\right]}.
$$

If $D$ is the regular representation, this procedure yields all inequivalent irreducible representations of $G$.

### Projection Operators

Let $P$ be an $N$-dimensional unitary representation of $G$ on $F$, and choose an orthonormal basis $\mathcal E=(\boldsymbol e_1,\ldots,\boldsymbol e_N)$. The operator $P_g:F\to F$ and its matrix $D(g)$ are related by

$$
P_g\boldsymbol e_a=\sum_{b=1}^{N}\boldsymbol e_bD_{ba}(g),\qquad D(g)=[P_g]_{\mathcal E}\in\mathbb C^{N\times N}.
$$

Let $D^{(i)}(g)\in\mathbb C^{d_i\times d_i}$ be a known irreducible unitary matrix representation on $V_i$.

If $D^{(i)}$ occurs with multiplicity $\eta_i$, then

$$
F\simeq\bigoplus_{i=1}^{k}\left(\mathbb C^{\eta_i}\otimes V_i\right).
$$

Use the irreducible basis

$$
\lvert i,\tau,\mu\rangle=\lvert\tau\rangle\otimes\lvert i,\mu\rangle,\qquad \tau=1,\ldots,\eta_i,\qquad \mu=1,\ldots,d_i,
$$

where $i$ labels the irreducible representation, $\tau$ labels equivalent copies, and $\mu$ is the basis index within each copy.

##### Matrix-Element Projection Operators

- Definition

For $\mu,\nu=1,\ldots,d_i$, define the matrix-element projection operator

$$
\boxed{\mathcal P_{\mu\nu}^{(i)}=\frac{d_i}{n_G}\sum_{g\in G}D_{\mu\nu}^{(i)}(g)^*P_g}:F\to F.
$$

Its matrix in $\mathcal E$ is

$$
\Pi_{\mu\nu}^{(i)}=\left[\mathcal P_{\mu\nu}^{(i)}\right]_{\mathcal E}=\frac{d_i}{n_G}\sum_{g\in G}D_{\mu\nu}^{(i)}(g)^*D(g)\in\mathbb C^{N\times N}.
$$

The indices $\mu,\nu$ label matrix elements of $D^{(i)}(g)$, not rows and columns of $\Pi_{\mu\nu}^{(i)}$.

- Form in the irreducible basis

$$
\mathcal P_{\mu\nu}^{(i)}\lvert j,\tau,\rho\rangle=\delta_{ij}\delta_{\nu\rho}\lvert i,\tau,\mu\rangle.
$$

Hence, on the $D^{(i)}$-isotypic subspace $W_i$,

$$
\left.\mathcal P_{\mu\nu}^{(i)}\right|_{W_i}=I_{\eta_i}\otimes\lvert i,\mu\rangle\langle i,\nu\rvert=\sum_{\tau=1}^{\eta_i}\lvert i,\tau,\mu\rangle\langle i,\tau,\nu\rvert.
$$

The operator $\mathcal P_{\mu\nu}^{(i)}$ vanishes on $W_j$ for $j\neq i$.

- Matrix-unit relation

$$
\boxed{\mathcal P_{\mu\nu}^{(i)}\mathcal P_{\rho\sigma}^{(j)}=\delta_{ij}\delta_{\nu\rho}\mathcal P_{\mu\sigma}^{(i)}}.
$$

In particular,

$$
\left(\mathcal P_{\mu\mu}^{(i)}\right)^2=\mathcal P_{\mu\mu}^{(i)},\qquad \left(\mathcal P_{\mu\nu}^{(i)}\right)^\dagger=\mathcal P_{\nu\mu}^{(i)}.
$$

Moreover,

$$
\left.\mathcal P_{\mu\mu}^{(i)}\right|_{W_i}=I_{\eta_i}\otimes\lvert i,\mu\rangle\langle i,\mu\rvert=\sum_{\tau=1}^{\eta_i}\lvert i,\tau,\mu\rangle\langle i,\tau,\mu\rvert,
$$

and therefore

$$
\operatorname{rank}\mathcal P_{\mu\mu}^{(i)}=\operatorname{rank}\Pi_{\mu\mu}^{(i)}=\eta_i.
$$

Thus, each diagonal operator projects onto one basis component from every equivalent copy; the off-diagonal matrix units transfer that component between irreducible basis directions.

- Transformation

For any $\boldsymbol\psi\in F$,

$$
P_s\!\left(\mathcal P_{\mu\nu}^{(i)}\boldsymbol\psi\right)=\sum_{\rho=1}^{d_i}\left(\mathcal P_{\rho\nu}^{(i)}\boldsymbol\psi\right)D_{\rho\mu}^{(i)}(s).
$$

For fixed $\nu$, the vectors $\mathcal P_{\mu\nu}^{(i)}\boldsymbol\psi$ transform under $D^{(i)}$.

##### Character Projection Operators

- Definition

The character projection operator is the sum of the diagonal matrix-element projectors:

$$
\boxed{\mathcal P^{(i)}=\sum_{\mu=1}^{d_i}\mathcal P_{\mu\mu}^{(i)}=\frac{d_i}{n_G}\sum_{g\in G}\chi_i(g)^*P_g}.
$$

Its matrix in $\mathcal E$ is

$$
\Pi^{(i)}=\left[\mathcal P^{(i)}\right]_{\mathcal E}=\sum_{\mu=1}^{d_i}\Pi_{\mu\mu}^{(i)}=\frac{d_i}{n_G}\sum_{g\in G}\chi_i(g)^*D(g)\in\mathbb C^{N\times N}.
$$

- Form in the irreducible basis

On $W_i$,

$$
\left.\mathcal P^{(i)}\right|_{W_i}=I_{\eta_i}\otimes I_{d_i}=I_{\eta_i}\otimes\sum_{\mu=1}^{d_i}\lvert i,\mu\rangle\langle i,\mu\rvert=\sum_{\tau=1}^{\eta_i}\sum_{\mu=1}^{d_i}\lvert i,\tau,\mu\rangle\langle i,\tau,\mu\rvert.
$$

The operator $\mathcal P^{(i)}$ vanishes on $W_j$ for $j\neq i$.

- Orthogonality and completeness

$$
\mathcal P^{(i)}\mathcal P^{(j)}=\delta_{ij}\mathcal P^{(i)},\qquad \sum_{i=1}^{k}\mathcal P^{(i)}=I_F.
$$

Therefore,

$$
F=\bigoplus_{i=1}^{k}W_i,\qquad W_i=\operatorname{im}\mathcal P^{(i)},\qquad \dim W_i=\operatorname{rank}\mathcal P^{(i)}=\operatorname{rank}\Pi^{(i)}=\eta_i d_i.
$$

The coordinate vectors of $W_i$ form $\operatorname{col}\Pi^{(i)}\subseteq\mathbb C^N$.

Thus, $\mathcal P^{(i)}$ selects the complete $D^{(i)}$-isotypic subspace. If $\eta_i>1$, it selects all equivalent copies together and does not separate them.

##### Reduction by Projection Operators

Assume that the reducible unitary matrix representation

$$
D(g)\in\mathbb C^{N\times N}
$$

and the irreducible unitary matrix representations

$$
D^{(i)}(g)\in\mathbb C^{d_i\times d_i}
$$

are known. All operations below are performed on matrices in the chosen basis $\mathcal E$.

- Step 1: Determine the multiplicity

Compute

$$
\chi_D(g)=\operatorname{tr}D(g),\qquad \eta_i=\frac{1}{n_G}\sum_{g\in G}\chi_i(g)^*\chi_D(g).
$$

The numerical value of $\eta_i$ must agree with a nonnegative integer within tolerance. If $\eta_i=0$, skip $D^{(i)}$.

- Step 2: Construct the seed projector and extract the seeds

Fix an internal basis index

$$
\nu\in\{1,\ldots,d_i\}.
$$

Compute the seed projection matrix

$$
\Pi_{\nu\nu}^{(i)}=\frac{d_i}{n_G}\sum_{g\in G}D_{\nu\nu}^{(i)}(g)^*D(g)\in\mathbb C^{N\times N}.
$$

Its rank gives the number of equivalent copies:

$$
\operatorname{rank}\Pi_{\nu\nu}^{(i)}=\eta_i.
$$

Diagonalize

$$
\Pi_{\nu\nu}^{(i)}=V\Lambda V^\dagger.
$$

Its eigenvalues are $0$ or $1$. Select an orthonormal basis of the eigenspace with $\lambda=1$:

$$
V_{i\nu}=\begin{pmatrix}\boldsymbol v_{\nu1}&\cdots&\boldsymbol v_{\nu\eta_i}\end{pmatrix}\in\mathbb C^{N\times\eta_i},\qquad V_{i\nu}^\dagger V_{i\nu}=I_{\eta_i}.
$$

Each $\boldsymbol v_{\nu\tau}$ serves as a seed that generates one copy of $D^{(i)}$.

- Step 3: Generate the irreducible bases

For $\mu=1,\ldots,d_i$, compute

$$
\Pi_{\mu\nu}^{(i)}=\frac{d_i}{n_G}\sum_{g\in G}D_{\mu\nu}^{(i)}(g)^*D(g),
$$

and apply it to all seeds:

$$
V_{i\mu}=\Pi_{\mu\nu}^{(i)}V_{i\nu}=\begin{pmatrix}\boldsymbol v_{\mu1}&\cdots&\boldsymbol v_{\mu\eta_i}\end{pmatrix}.
$$

For fixed $\tau$, the vectors $\boldsymbol v_{\mu\tau}$ transform under $D^{(i)}$:

$$
D(g)\boldsymbol v_{\mu\tau}=\sum_{\rho=1}^{d_i}\boldsymbol v_{\rho\tau}D_{\rho\mu}^{(i)}(g).
$$

- Step 4: Assemble and verify

Arrange the basis matrices as

$$
B_{i\tau}=\begin{pmatrix}\boldsymbol v_{1\tau}&\cdots&\boldsymbol v_{d_i\tau}\end{pmatrix},\qquad B_i=\begin{pmatrix}B_{i1}&\cdots&B_{i\eta_i}\end{pmatrix},\qquad S=\begin{pmatrix}B_1&\cdots&B_k\end{pmatrix}.
$$

Each irreducible block satisfies

$$
D(g)B_{i\tau}=B_{i\tau}D^{(i)}(g).
$$

Hence,

$$
\boxed{S^\dagger D(g)S=\bigoplus_{i=1}^{k}\left[I_{\eta_i}\otimes D^{(i)}(g)\right]}.
$$

Verify

$$
\sum_{i=1}^{k}\eta_i d_i=N,\qquad \left\|S^\dagger S-I_N\right\|\approx0,
$$

$$
\max_{g\in G}\left\|D(g)S-S\left[\bigoplus_{i=1}^{k}\left(I_{\eta_i}\otimes D^{(i)}(g)\right)\right]\right\|\approx0.
$$

### Character Theory

##### Character Tables

Let $G$ be a finite group of order $n_G=|G|$. Let $C_1,\ldots,C_k$ be its conjugacy classes, with $n_\alpha=|C_\alpha|$ and $C_1=\{e\}$. Let $D^{(1)},\ldots,D^{(\ell)}$ be all inequivalent irreducible representations, with $D^{(1)}$ the trivial representation. Define

$$
d_i=\dim D^{(i)},\qquad \chi_{i\alpha}=\operatorname{tr}D^{(i)}(g)\quad(g\in C_\alpha).
$$

- Step 1: Determine the table size

The class-function space has dimension $k$, and the irreducible characters form an orthogonal complete basis. Therefore,

$$
\boxed{\ell=k}.
$$

The character table is a $k\times k$ matrix, with columns indexed by $C_\alpha$ and rows indexed by $D^{(i)}$.

- Step 2: Determine possible dimensions

$$
d_1=1,\qquad \sum_{i=1}^{k}d_i^2=n_G.
$$

Use $n_G$ and $k$ to determine the possible dimension tuples $(d_1,\ldots,d_k)$.

- Step 3: Fill the first row and column

$$
\chi_{1\alpha}=1,\qquad \chi_{i1}=d_i.
$$

- Step 4: Complete the remaining entries

$$
\sum_{\alpha=1}^{k}\frac{n_\alpha}{n_G}\chi_{i\alpha}^*\chi_{j\alpha}=\delta_{ij}.
$$

$$
\sum_{i=1}^{k}\frac{\sqrt{n_\alpha n_\beta}}{n_G}\chi_{i\alpha}^*\chi_{i\beta}=\delta_{\alpha\beta}.
$$
