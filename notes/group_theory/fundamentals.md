### Groups and Cayley Tables

##### Definition of a Group

A group satisfies four axioms:

- Closure
- Associativity
- Identity
- Inverse

Let $G$ be a nonempty set with multiplication $G\times G\to G$. For all $g_1,g_2,g_3\in G$,

$$
g_1g_2\in G,
\qquad
(g_1g_2)g_3=g_1(g_2g_3).
$$

The identity $e\in G$ and inverse $g^{-1}\in G$ satisfy

$$
eg=ge=g,
\qquad
g^{-1}g=gg^{-1}=e.
$$

For identities $e,e'$ and inverses $h,k$ of $g$,

$$
e=ee'=e',
\qquad
h=h(gk)=(hg)k=k.
$$

The left-sided conditions $eg=g$ and $g^{-1}g=e$ are sufficient. If $h$ is a left inverse of $g^{-1}$, then

$$
gg^{-1}
=
(hg^{-1})(gg^{-1})
=
h(g^{-1}g)g^{-1}
=
hg^{-1}
=
e,
$$

$$
ge=g(g^{-1}g)=(gg^{-1})g=g.
$$
