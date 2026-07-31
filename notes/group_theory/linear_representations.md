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

- Character decomposition

$$
\chi_D(g)=\sum_{i=1}^{r}\chi_{D^{(i)}}(g).
$$

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

##### Completeness from the Regular Representation

Let $\eta_i$ be the multiplicity of the irreducible matrix representation $D^{(i)}$ in a completely reducible representation $D$. Then

$$
D\simeq\bigoplus_{i=1}^{\ell}\underbrace{\left(D^{(i)}\oplus\cdots\oplus D^{(i)}\right)}_{\eta_i\text{ copies}}.
$$

Therefore,

$$
\chi_D(g)=\sum_{i=1}^{\ell}\eta_i\chi_i(g).
$$

Using character orthogonality,

$$
\sum_{g\in G}\chi_i(g)^*\chi_D(g)=\sum_{g\in G}\chi_i(g)^*\left(\sum_{j=1}^{\ell}\eta_j\chi_j(g)\right)=\sum_j\eta_j\left(\sum_{g\in G}\chi_i(g)^*\chi_j(g)\right)=n_G\eta_i.
$$

Hence,

$$
\eta_i=\frac{1}{n_G}\sum_{g\in G}\chi_i(g)^*\chi_D(g).
$$

For the regular representation,

$$
\eta_i=\frac{1}{n_G}\chi_i(e)^*\chi_{\mathrm{reg}}(e)=d_i.
$$

Since

$$
\dim D=\sum_{i=1}^{\ell}\eta_i d_i,\qquad \dim D_{\mathrm{reg}}=\dim\mathbb C[G]=n_G,
$$

we obtain

$$
\boxed{n_G=\sum_{i=1}^{\ell}\eta_i d_i=\sum_{i=1}^{\ell}d_i^2}.
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
