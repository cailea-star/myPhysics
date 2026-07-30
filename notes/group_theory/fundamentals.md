### Groups and Cayley Tables

##### Definition of a Group

A group satisfies four axioms:

- Closure
- Associativity
- Identity
- Inverse

Let $G$ be a nonempty set with a binary operation $G\times G\to G$, written by juxtaposition. For all $g_1,g_2,g_3\in G$,

$$
g_1g_2\in G,\qquad (g_1g_2)g_3=g_1(g_2g_3).
$$

The identity $e\in G$ and inverse $g^{-1}\in G$ satisfy

$$
eg=ge=g,\qquad g^{-1}g=gg^{-1}=e.
$$

If $e,\tilde e$ are identities and $h,k$ are inverses of $g$, then

$$
e=e\tilde e=\tilde e,\qquad h=h(gk)=(hg)k=k.
$$

Thus the identity and the inverse of each element are unique.

##### Basic Terminology

Let $G$ be a group with identity $e$. Let $\mathbb Z_{>0}$ denote the set of positive integers.

- Commutativity

$$
gh=hg,\qquad \forall g,h\in G.
$$

The group is Abelian if this condition holds and non-Abelian otherwise.

- Group order

The order of $G$, denoted by $|G|$, is its number of elements. The group is finite if $|G|<\infty$.

- Element order

For $g\in G$, define

$$
\operatorname{ord}(g)=\min\{n\in\mathbb Z_{>0}:g^n=e\}.
$$

If no such $n$ exists, then $g$ has infinite order.

- Generating set

For $S\subseteq G$, let $\langle S\rangle$ denote the smallest subgroup of $G$ containing $S$. The set $S$ generates $G$ if

$$
\langle S\rangle=G.
$$

- Group rank

$$
\operatorname{rank}(G)=\min\{|S|:G=\langle S\rangle\}.
$$

##### Rearrangement Theorem

For $s\in G$, define

$$
sG=\{sg:g\in G\},\qquad Gs=\{gs:g\in G\}\quad\Longrightarrow\quad sG\subseteq G,\qquad Gs\subseteq G.
$$

For every $g\in G$,

$$
g=s(s^{-1}g)\in sG,\qquad g=(gs^{-1})s\in Gs\quad\Longrightarrow\quad G\subseteq sG,\qquad G\subseteq Gs.
$$

Therefore,

$$
\boxed{sG=G,\qquad Gs=G}.
$$

Left and right multiplication by $s$ permute the elements of $G$.

##### Cayley Tables

Let $G=\{g_1,\ldots,g_n\}$ be a finite group. Define the Cayley table by

$$
T_{ij}=g_ig_j,\qquad 1\leq i,j\leq n.
$$

The rearrangement theorem gives

$$
\boxed{\{T_{ij}:1\leq j\leq n\}=g_iG=G,\qquad \{T_{ij}:1\leq i\leq n\}=Gg_j=G}.
$$

Thus every element of $G$ occurs exactly once in each row and column.

The table is symmetric precisely when $G$ is Abelian:

$$
T_{ij}=T_{ji}\quad\Longleftrightarrow\quad g_ig_j=g_jg_i.
$$

### Subgroups and Cosets

##### Subgroups

Let $G$ be a group with identity $e$. A subset $H\subseteq G$ is a subgroup, denoted by $H\leq G$, if

$$
h_1h_2\in H,\qquad e\in H,\qquad h^{-1}\in H,\qquad \forall h,h_1,h_2\in H.
$$

Associativity is inherited from $G$. The trivial subgroups are

$$
\{e\}\leq G,\qquad G\leq G.
$$

##### Left and Right Cosets

Let $H\leq G$ and $g\in G$. The left and right cosets represented by $g$ are

$$
gH=\{gh:h\in H\},\qquad Hg=\{hg:h\in H\}.
$$

Multiplication by $g$ is bijective, so

$$
|gH|=|Hg|=|H|.
$$

Moreover,

$$
g\in H\quad\Longleftrightarrow\quad gH=H\quad\Longleftrightarrow\quad Hg=H.
$$

##### Coset Partitions

For $g_1,g_2\in G$,

$$
g_1H=g_2H\Longrightarrow g_2=g_1h\ (h\in H)\Longrightarrow g_1^{-1}g_2\in H.
$$

$$
g_1^{-1}g_2=h\in H\Longrightarrow g_2H=g_1hH=g_1H.
$$

Therefore,

$$
\boxed{g_1H=g_2H\quad\Longleftrightarrow\quad g_1^{-1}g_2\in H}.
$$

Similarly,

$$
Hg_1=Hg_2\quad\Longleftrightarrow\quad g_1g_2^{-1}\in H.
$$

If $x\in g_1H\cap g_2H$, then

$$
x=g_1h_1=g_2h_2\quad\Longrightarrow\quad g_1^{-1}g_2=h_1h_2^{-1}\in H\quad\Longrightarrow\quad g_1H=g_2H.
$$

Hence,

$$
\boxed{g_1H=g_2H\quad\text{or}\quad g_1H\cap g_2H=\varnothing}.
$$

The same result holds for right cosets.

For a finite group, let $r_1=e,r_2,\ldots,r_d$ represent the distinct left cosets. Then

$$
G=\bigcup_{i=1}^{d}r_iH,\qquad r_iH\cap r_jH=\varnothing\quad(i\neq j).
$$

##### Lagrange's Theorem

Let $G$ be finite and $H\leq G$. The index $[G:H]$ is the number of left cosets of $H$. Since every coset contains $|H|$ elements,

$$
\boxed{|G|=[G:H]|H|}.
$$

$$
[G:H]=\frac{|G|}{|H|}\in\mathbb Z_{>0}.
$$

For $g\in G$,

$$
\operatorname{ord}(g)=|\langle g\rangle|.
$$

$$
\frac{|G|}{\operatorname{ord}(g)}\in\mathbb Z_{>0},\qquad g^{|G|}=e.
$$

If $|G|=p$ is prime, then

$$
|H|\in\{1,p\}.
$$

Thus $G$ has no nontrivial proper subgroup, and every nonidentity element generates $G$.

### Conjugacy and Normal Subgroups

##### Conjugacy Classes

Let $G$ be a group. Elements $g,g'\in G$ are conjugate, written $g\sim g'$, if

$$
\boxed{g'=xgx^{-1},\qquad x\in G}.
$$

Conjugacy is an equivalence relation. The conjugacy class of $g$ is

$$
\mathcal C(g)=\{xgx^{-1}:x\in G\}.
$$

Hence,

$$
\mathcal C(g)=\mathcal C(g')\quad\text{or}\quad\mathcal C(g)\cap\mathcal C(g')=\varnothing.
$$

Moreover,

$$
\mathcal C(e)=\{e\}.
$$

$$
G\text{ is Abelian}\quad\Longleftrightarrow\quad\mathcal C(g)=\{g\},\qquad \forall g\in G.
$$

Conjugate elements have the same order:

$$
g'\in\mathcal C(g)\quad\Longrightarrow\quad\operatorname{ord}(g')=\operatorname{ord}(g).
$$

##### Centralizers and Class Sizes

For $g\in G$, define its centralizer by

$$
C_G(g)=\{x\in G:xg=gx\}.
$$

$$
C_G(g)\leq G.
$$

For $x,y\in G$,

$$
xgx^{-1}=ygy^{-1}\quad\Longleftrightarrow\quad x^{-1}y\in C_G(g).
$$

$$
xgx^{-1}=ygy^{-1}\quad\Longleftrightarrow\quad xC_G(g)=yC_G(g).
$$

If $G$ is finite, the conjugates of $g$ correspond to the left cosets of $C_G(g)$, so

$$
\boxed{|\mathcal C(g)|=[G:C_G(g)]=\frac{|G|}{|C_G(g)|}}.
$$

##### Normal Subgroups

Let $H\leq G$. The subgroup $H$ is normal in $G$, denoted by $H\trianglelefteq G$, if

$$
\boxed{gH=Hg,\qquad \forall g\in G}.
$$

$$
H\trianglelefteq G\quad\Longleftrightarrow\quad gHg^{-1}=H,\qquad \forall g\in G.
$$

Equivalently,

$$
H\trianglelefteq G\quad\Longleftrightarrow\quad ghg^{-1}\in H,\qquad \forall g\in G,\quad\forall h\in H.
$$

Thus,

$$
H\trianglelefteq G\quad\Longleftrightarrow\quad \mathcal C(h)\subseteq H,\qquad \forall h\in H.
$$

Hence a normal subgroup is precisely a union of complete conjugacy classes.

Two basic consequences are

$$
G\text{ is Abelian},\ H\leq G\quad\Longrightarrow\quad H\trianglelefteq G.
$$

$$
[G:H]=2\quad\Longrightarrow\quad H\trianglelefteq G.
$$

##### Quotient Groups

Let $H\trianglelefteq G$. Define

$$
G/H=\{gH:g\in G\}.
$$

For $g,k\in G$ and $h_1,h_2\in H$, normality gives $k^{-1}h_1k\in H$, and

$$
(gh_1)(kh_2)=gk(k^{-1}h_1k)h_2\in gkH.
$$

Therefore the coset product is well-defined:

$$
\boxed{(gH)(kH)=(gk)H}.
$$

The identity and inverse are

$$
H,\qquad (gH)^{-1}=g^{-1}H.
$$

For finite $G$,

$$
|G/H|=[G:H]=\frac{|G|}{|H|}.
$$

### Homomorphisms and Isomorphisms

##### Homomorphisms

Let $G$ and $G'$ be groups with identities $e$ and $e'$. A map

$$
\varphi:G\to G'
$$

is a group homomorphism if, for all $g_1,g_2\in G$,

$$
\boxed{\varphi(g_1g_2)=\varphi(g_1)\varphi(g_2)}.
$$

- Identity

$$
\varphi(e)=\varphi(ee)=\varphi(e)\varphi(e)\quad\Longrightarrow\quad\boxed{\varphi(e)=e'}.
$$

- Inverses

For $g\in G$,

$$
e'=\varphi(e)=\varphi(gg^{-1})=\varphi(g)\varphi(g^{-1}),
$$

so

$$
\boxed{\varphi(g^{-1})=\varphi(g)^{-1}}.
$$

- Injectivity

$$
\varphi(g_1)=\varphi(g_2)\quad\Longrightarrow\quad g_1=g_2.
$$

- Surjectivity

$$
\forall g'\in G',\qquad \exists g\in G:\varphi(g)=g'.
$$

A homomorphism need not be injective or surjective; distinct elements of $G$ may have the same image.

##### Images and Kernels

Let $\varphi:G\to G'$ be a homomorphism. Its image and kernel are

$$
\operatorname{Im}\varphi=\{\varphi(g):g\in G\},\qquad \ker\varphi=\{g\in G:\varphi(g)=e'\}.
$$

- Image

For $g_1,g_2\in G$,

$$
\varphi(g_1)\varphi(g_2)^{-1}=\varphi(g_1g_2^{-1})\in\operatorname{Im}\varphi.
$$

Hence,

$$
\boxed{\operatorname{Im}\varphi\leq G'}.
$$

- Kernel

For $h_1,h_2\in\ker\varphi$,

$$
\varphi(h_1h_2^{-1})=\varphi(h_1)\varphi(h_2)^{-1}=e',
$$

so $\ker\varphi\leq G$.

For $g\in G$ and $h\in\ker\varphi$,

$$
\varphi(ghg^{-1})=\varphi(g)e'\varphi(g)^{-1}=e'.
$$

Therefore,

$$
\boxed{\ker\varphi\trianglelefteq G}.
$$

- Fibers and cosets

For $g_1,g_2\in G$,

$$
\boxed{\varphi(g_1)=\varphi(g_2)\quad\Longleftrightarrow\quad g_1^{-1}g_2\in\ker\varphi\quad\Longleftrightarrow\quad g_1\ker\varphi=g_2\ker\varphi}.
$$

Thus the elements mapped to $\varphi(g)$ form one coset:

$$
\{x\in G:\varphi(x)=\varphi(g)\}=g\ker\varphi.
$$

- Injectivity and surjectivity

$$
\boxed{\varphi\text{ is injective}\quad\Longleftrightarrow\quad\ker\varphi=\{e\}}.
$$

$$
\boxed{\varphi\text{ is surjective}\quad\Longleftrightarrow\quad\operatorname{Im}\varphi=G'}.
$$

##### Isomorphisms

Let $\varphi:G\to G'$ be a homomorphism. It is an isomorphism if it is both injective and surjective.

If such a map exists, $G$ and $G'$ are isomorphic, written

$$
\boxed{G\cong G'}.
$$

Equivalently,

$$
\boxed{\varphi\text{ is an isomorphism}\quad\Longleftrightarrow\quad\ker\varphi=\{e\},\qquad\operatorname{Im}\varphi=G'}.
$$

Isomorphic groups have the same multiplication structure up to a relabeling of their elements.

- Inverse map

If $\varphi$ is an isomorphism, then $\varphi^{-1}:G'\to G$ is also an isomorphism. For $g_1'=\varphi(g_1)$ and $g_2'=\varphi(g_2)$,

$$
\varphi^{-1}(g_1'g_2')=\varphi^{-1}(\varphi(g_1g_2))=\varphi^{-1}(g_1')\varphi^{-1}(g_2').
$$

- Preserved structure

For every $g\in G$,

$$
\operatorname{ord}(\varphi(g))=\operatorname{ord}(g).
$$

Moreover,

$$
|G|=|G'|,
$$

$$
G\text{ is Abelian}\quad\Longleftrightarrow\quad G'\text{ is Abelian}.
$$

The identity map, inverse maps, and compositions of isomorphisms are isomorphisms. Hence isomorphism is an equivalence relation among groups.

##### First Isomorphism Theorem

Let $\varphi:G\to G'$ be a homomorphism. Then

$$
\boxed{G/\ker\varphi\cong\operatorname{Im}\varphi}.
$$

Let $K=\ker\varphi$. Since $K\trianglelefteq G$, define

$$
\overline{\varphi}:G/K\to\operatorname{Im}\varphi,\qquad \overline{\varphi}(gK)=\varphi(g).
$$

- Well-definedness

For $g_1,g_2\in G$,

$$
g_1K=g_2K\quad\Longleftrightarrow\quad\varphi(g_1)=\varphi(g_2).
$$

Hence $\overline{\varphi}(gK)$ does not depend on the representative $g$.

- Homomorphism

$$
\overline{\varphi}((g_1K)(g_2K))=\overline{\varphi}(g_1g_2K)=\varphi(g_1)\varphi(g_2).
$$

- Injectivity

$$
\overline{\varphi}(g_1K)=\overline{\varphi}(g_2K)\quad\Longrightarrow\quad g_1K=g_2K.
$$

- Surjectivity

For every $y\in\operatorname{Im}\varphi$, there exists $g\in G$ such that

$$
y=\varphi(g)=\overline{\varphi}(gK).
$$

Therefore, $\overline{\varphi}$ is an isomorphism.

- Finite groups

If $G$ is finite, then

$$
\boxed{|G|=|\ker\varphi|\,|\operatorname{Im}\varphi|}.
$$

- Natural projection

For $H\trianglelefteq G$, define

$$
\pi:G\to G/H,\qquad \pi(g)=gH.
$$

Then

$$
\ker\pi=H,\qquad\operatorname{Im}\pi=G/H.
$$

Thus every normal subgroup is the kernel of a surjective homomorphism.

### Direct Products

Let $G_1$ and $G_2$ be groups with identities $e_1$ and $e_2$. Their direct product is

$$
G_1\times G_2=\{(g_1,g_2):g_1\in G_1,\ g_2\in G_2\},
$$

with componentwise multiplication

$$
\boxed{(g_1,g_2)(h_1,h_2)=(g_1h_1,g_2h_2)}.
$$

The identity and inverse are

$$
(e_1,e_2),\qquad (g_1,g_2)^{-1}=(g_1^{-1},g_2^{-1}).
$$

For finite groups,

$$
\boxed{|G_1\times G_2|=|G_1||G_2|}.
$$

- Internal direct products

Let $H_1,H_2\leq G$, and define

$$
H_1H_2=\{h_1h_2:h_1\in H_1,\ h_2\in H_2\}.
$$

Suppose that

$$
H_1\cap H_2=\{e\},
$$

$$
h_1h_2=h_2h_1,\qquad \forall h_1\in H_1,\quad\forall h_2\in H_2,
$$

$$
G=H_1H_2.
$$

Every $g\in G$ then has a unique expression $g=h_1h_2$. Indeed,

$$
h_1h_2=k_1k_2\quad\Longrightarrow\quad k_1^{-1}h_1=k_2h_2^{-1}\in H_1\cap H_2=\{e\}.
$$

Define

$$
\mu:H_1\times H_2\to G,\qquad \mu(h_1,h_2)=h_1h_2.
$$

The three conditions make $\mu$ a bijective homomorphism. Therefore,

$$
\boxed{G\cong H_1\times H_2}.
$$

Moreover,

$$
H_1\trianglelefteq G,\qquad H_2\trianglelefteq G.
$$

For finite $G$,

$$
\boxed{|G|=|H_1||H_2|}.
$$
