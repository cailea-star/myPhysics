# Wigner Symbols

> Main reference: D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, *Quantum Theory of Angular Momentum: Irreducible Tensors, Spherical Harmonics, Vector Coupling Coefficients, 3nj Symbols*, World Scientific, 1988.

### Wigner Symbols

##### Conventions

- Angular-momentum quantum numbers and dimension factor

  $$
  j=0,\frac12,1,\frac32,\ldots,\qquad m=-j,-j+1,\ldots,j,\qquad \hat j\equiv\sqrt{2j+1}.
  $$

- Phase convention

  All Clebsch–Gordan coefficients and Wigner symbols use the Condon–Shortley convention.

- Coupled basis

  $$
  |(j_1j_2)j_3m_3\rangle=\sum_{m_1,m_2}|j_1m_1\rangle\otimes|j_2m_2\rangle\langle j_1m_1;j_2m_2|j_3m_3\rangle.
  $$

- Orthogonality and completeness

  Orthogonality of the coupled basis:

  $$
  \boxed{\sum_{m_1,m_2}\langle j_3m_3|j_1m_1;j_2m_2\rangle\langle j_1m_1;j_2m_2|j_3'm_3'\rangle=\delta_{j_3j_3'}\delta_{m_3m_3'}}.
  $$

  Completeness of the coupled basis:

  $$
  \boxed{\sum_{j_3,m_3}\langle j_1m_1;j_2m_2|j_3m_3\rangle\langle j_3m_3|j_1m_1';j_2m_2'\rangle=\delta_{m_1m_1'}\delta_{m_2m_2'}}.
  $$

- General tensor-product decomposition

  See [Linear Representations](../group_theory/linear_representations.md#clebsch-gordan-coefficients).

##### Wigner 3j Symbols

- Definition

  $$
  \boxed{\begin{pmatrix}j_1&j_2&j_3\\m_1&m_2&m_3\end{pmatrix}=\frac{(-1)^{j_1-j_2-m_3}}{\hat j_3}\langle j_1m_1;j_2m_2|j_3,-m_3\rangle}.
  $$

- Selection rules

  $$
  m_1+m_2+m_3=0,\qquad |j_1-j_2|\leq j_3\leq j_1+j_2,\qquad |m_i|\leq j_i,\qquad j_1+j_2+j_3\in\mathbb Z.
  $$

- Symmetries, where $J_\Sigma=j_1+j_2+j_3$

  $$
  \begin{pmatrix}j_1&j_2&j_3\\m_1&m_2&m_3\end{pmatrix}=\begin{pmatrix}j_2&j_3&j_1\\m_2&m_3&m_1\end{pmatrix}=(-1)^{J_\Sigma}\begin{pmatrix}j_2&j_1&j_3\\m_2&m_1&m_3\end{pmatrix}=(-1)^{J_\Sigma}\begin{pmatrix}j_1&j_2&j_3\\-m_1&-m_2&-m_3\end{pmatrix}.
  $$

- Zero total angular momentum

  $$
  \langle j_1m_1;j_2m_2|00\rangle=\frac{(-1)^{j_1-m_1}}{\hat j_1}\delta_{j_1j_2}\delta_{m_1,-m_2}.
  $$

- Vanishing condition for integer $j_1,j_2,j_3$

  $$
  \begin{pmatrix}j_1&j_2&j_3\\0&0&0\end{pmatrix}=0\qquad\text{if }j_1+j_2+j_3\text{ is odd}.
  $$

##### Wigner 6j Symbols

- Definition by recoupling three angular momenta

  $$
  \boxed{|(j_1j_2)j_{12},j_3;JM\rangle=\sum_{j_{23}}(-1)^{j_1+j_2+j_3+J}\hat j_{12}\hat j_{23}\begin{Bmatrix}j_1&j_2&j_{12}\\j_3&J&j_{23}\end{Bmatrix}|j_1,(j_2j_3)j_{23};JM\rangle}.
  $$

- Selection rules

  Each triple must satisfy the triangle condition:

  $$
  (j_1,j_2,j_{12}),\qquad (j_1,J,j_{23}),\qquad (j_3,j_2,j_{23}),\qquad (j_3,J,j_{12}).
  $$

- Symmetries

  The Wigner $6j$ symbol is invariant under arbitrary permutations of its columns and simultaneous interchange of the upper and lower entries in any two columns:

  $$
  \boxed{\begin{Bmatrix}a&b&c\\d&e&f\end{Bmatrix}=\begin{Bmatrix}b&a&c\\e&d&f\end{Bmatrix}=\begin{Bmatrix}d&e&c\\a&b&f\end{Bmatrix}=\begin{Bmatrix}d&b&f\\a&e&c\end{Bmatrix}=\begin{Bmatrix}a&e&f\\d&b&c\end{Bmatrix}}.
  $$

  Combining the six column permutations with the four choices of interchanging zero or two column pairs gives $24$ equivalent forms. No phase factor is introduced.

  These transformations are the tetrahedral symmetries of the four triangle triples

  $$
  (a,b,c),\qquad(a,e,f),\qquad(d,b,f),\qquad(d,e,c).
  $$

- Zero-entry identity

  $$
  \begin{Bmatrix}a&b&0\\d&e&f\end{Bmatrix}=(-1)^{a+e+f}\frac{\delta_{ab}\delta_{de}}{\sqrt{(2a+1)(2d+1)}}.
  $$

##### Wigner 9j Symbols

- Definition by recoupling four angular momenta

  $$
  \boxed{|(j_1j_2)j_{12},(j_3j_4)j_{34};JM\rangle=\sum_{j_{13},j_{24}}\hat j_{12}\hat j_{34}\hat j_{13}\hat j_{24}\begin{Bmatrix}j_1&j_2&j_{12}\\j_3&j_4&j_{34}\\j_{13}&j_{24}&J\end{Bmatrix}|(j_1j_3)j_{13},(j_2j_4)j_{24};JM\rangle}.
  $$

- Selection rules

  Every row and column must satisfy the triangle condition.

- Symmetries, where $J_\Sigma=a+b+c+d+e+f+g+h+i$

  Transposition introduces no phase:

  $$
  \boxed{\begin{Bmatrix}a&b&c\\d&e&f\\g&h&i\end{Bmatrix}=\begin{Bmatrix}a&d&g\\b&e&h\\c&f&i\end{Bmatrix}}.
  $$

  Even permutations of rows or columns introduce no phase:

  $$
  \begin{Bmatrix}a&b&c\\d&e&f\\g&h&i\end{Bmatrix}=\begin{Bmatrix}d&e&f\\g&h&i\\a&b&c\end{Bmatrix}=\begin{Bmatrix}b&c&a\\e&f&d\\h&i&g\end{Bmatrix}.
  $$

  An odd permutation of either rows or columns contributes $(-1)^{J_\Sigma}$:

  $$
  \begin{Bmatrix}a&b&c\\d&e&f\\g&h&i\end{Bmatrix}=(-1)^{J_\Sigma}\begin{Bmatrix}d&e&f\\a&b&c\\g&h&i\end{Bmatrix}=(-1)^{J_\Sigma}\begin{Bmatrix}b&a&c\\e&d&f\\h&g&i\end{Bmatrix}.
  $$

  Simultaneous odd permutations of the rows and columns introduce no phase:

  $$
  \begin{Bmatrix}a&b&c\\d&e&f\\g&h&i\end{Bmatrix}=\begin{Bmatrix}e&d&f\\b&a&c\\h&g&i\end{Bmatrix}.
  $$

  Row permutations, column permutations, and transposition generate $3!\times3!\times2=72$ equivalent forms.

### Reduced Matrix Elements

##### Spherical Tensor Operators

Let $\hat{\mathbf J}$ be the angular-momentum operator defined in [Angular Momentum](../quantum_mechanics/angular_momentum.md), with $\hat J_\pm=\hat J_x\pm i\hat J_y$. An irreducible spherical tensor operator of rank $k$ consists of $2k+1$ components $\hat T_q^{(k)}$, where $q=-k,\ldots,k$, satisfying

$$
\boxed{[\hat J_z,\hat T_q^{(k)}]=\hbar q\hat T_q^{(k)},\qquad[\hat J_\pm,\hat T_q^{(k)}]=\hbar\sqrt{k(k+1)-q(q\pm1)}\,\hat T_{q\pm1}^{(k)}}.
$$

A rank-$0$ tensor is a rotational scalar:

$$
[\hat J_i,\hat T_0^{(0)}]=0,\qquad i\in\{x,y,z\}.
$$

For a vector operator $\hat{\mathbf A}=(\hat A_x,\hat A_y,\hat A_z)$, the rank-$1$ spherical components are

$$
\hat A_0^{(1)}=\hat A_z,\qquad\hat A_{\pm1}^{(1)}=\mp\frac{\hat A_x\pm i\hat A_y}{\sqrt2}.
$$

A Hermitian spherical tensor satisfies

$$
\left(\hat T_q^{(k)}\right)^\dagger=(-1)^q\hat T_{-q}^{(k)}.
$$

The tensor product of ranks $k_1$ and $k_2$ can be coupled to rank $k$:

$$
\boxed{\left[\hat A^{(k_1)}\otimes\hat B^{(k_2)}\right]_q^{(k)}=\sum_{q_1,q_2}\langle k_1q_1;k_2q_2|kq\rangle\hat A_{q_1}^{(k_1)}\hat B_{q_2}^{(k_2)}}.
$$

The allowed ranks are

$$
|k_1-k_2|\leq k\leq k_1+k_2,\qquad q=q_1+q_2.
$$

##### Wigner–Eckart Theorem

Let $|\alpha_i j_i m_i\rangle$, with $i\in\{1,2\}$, be angular-momentum eigenstates, where $\alpha_i$ denotes all quantum numbers other than $j_i$ and $m_i$. With the conventions above,

$$
\boxed{\langle\alpha_1j_1m_1|\hat T_q^{(k)}|\alpha_2j_2m_2\rangle=(-1)^{j_1-m_1}\begin{pmatrix}j_1&k&j_2\\-m_1&q&m_2\end{pmatrix}\langle\alpha_1j_1\|\hat T^{(k)}\|\alpha_2j_2\rangle=\frac{\langle j_2m_2;kq|j_1m_1\rangle}{\hat j_1}\langle\alpha_1j_1\|\hat T^{(k)}\|\alpha_2j_2\rangle}.
$$

The reduced matrix element $\langle\alpha_1j_1\|\hat T^{(k)}\|\alpha_2j_2\rangle$ is independent of $m_1$, $m_2$, and $q$. The angular dependence is entirely contained in the Wigner $3j$ symbol.

The reduced matrix element can be recovered from the full matrix elements:

$$
\langle\alpha_1j_1\|\hat T^{(k)}\|\alpha_2j_2\rangle=\sum_{m_1,m_2,q}(-1)^{j_1-m_1}\begin{pmatrix}j_1&k&j_2\\-m_1&q&m_2\end{pmatrix}\langle\alpha_1j_1m_1|\hat T_q^{(k)}|\alpha_2j_2m_2\rangle=\frac{1}{\hat j_1}\sum_{m_1,m_2,q}\langle j_2m_2;kq|j_1m_1\rangle\langle\alpha_1j_1m_1|\hat T_q^{(k)}|\alpha_2j_2m_2\rangle.
$$

The angular-momentum selection rules are

$$
m_1=m_2+q,\qquad |j_2-k|\leq j_1\leq j_2+k,\qquad j_1+j_2+k\in\mathbb Z.
$$

For a scalar operator,

$$
\langle\alpha_1j_1m_1|\hat T_0^{(0)}|\alpha_2j_2m_2\rangle=\frac{\delta_{j_1j_2}\delta_{m_1m_2}}{\hat j_1}\langle\alpha_1j_1\|\hat T^{(0)}\|\alpha_2j_2\rangle.
$$

### Common Coupling Formulas

The final and initial coupled states are denoted by $|(j_1j_2)JM\rangle$ and $|(j_1'j_2')J'M'\rangle$.

##### Tensor Acting on One Subsystem

Let $\hat T^{(k)}(1)$ act only on the first subsystem. Expanding the coupled states gives

$$
\langle(j_1j_2)JM|\hat T_q^{(k)}(1)|(j_1'j_2')J'M'\rangle=\sum_{m_1,m_2,m_1',m_2'}\langle JM|j_1m_1;j_2m_2\rangle\langle j_1'm_1';j_2'm_2'|J'M'\rangle\langle j_1m_1|\hat T_q^{(k)}|j_1'm_1'\rangle\langle j_2m_2|j_2'm_2'\rangle.
$$

The Wigner–Eckart theorem separates the total magnetic quantum numbers:

$$
\langle(j_1j_2)JM|\hat T_q^{(k)}(1)|(j_1'j_2')J'M'\rangle=(-1)^{J-M}\begin{pmatrix}J&k&J'\\-M&q&M'\end{pmatrix}\langle(j_1j_2)J\|\hat T^{(k)}(1)\|(j_1'j_2')J'\rangle.
$$

The reduced matrix element is

$$
\boxed{\langle(j_1j_2)J\|\hat T^{(k)}(1)\|(j_1'j_2')J'\rangle=\delta_{j_2j_2'}(-1)^{j_1+j_2+J'+k}\hat J\hat J'\begin{Bmatrix}j_1&J&j_2\\J'&j_1'&k\end{Bmatrix}\langle j_1\|\hat T^{(k)}\|j_1'\rangle}.
$$

##### Coupled Product of Two Tensor Operators

For tensors acting on different subsystems,

$$
\begin{aligned}&\langle(j_1j_2)JM|\left[\hat A^{(k_1)}(1)\otimes\hat B^{(k_2)}(2)\right]_q^{(k)}|(j_1'j_2')J'M'\rangle\\&=\sum_{m_1,m_2,m_1',m_2',q_1,q_2}\langle JM|j_1m_1;j_2m_2\rangle\langle j_1'm_1';j_2'm_2'|J'M'\rangle\langle k_1q_1;k_2q_2|kq\rangle\times\langle j_1m_1|\hat A_{q_1}^{(k_1)}|j_1'm_1'\rangle\langle j_2m_2|\hat B_{q_2}^{(k_2)}|j_2'm_2'\rangle.\end{aligned}
$$

Separating the total magnetic quantum numbers gives

$$
\langle(j_1j_2)JM|\left[\hat A^{(k_1)}(1)\otimes\hat B^{(k_2)}(2)\right]_q^{(k)}|(j_1'j_2')J'M'\rangle=(-1)^{J-M}\begin{pmatrix}J&k&J'\\-M&q&M'\end{pmatrix}\langle(j_1j_2)J\|\left[\hat A^{(k_1)}(1)\otimes\hat B^{(k_2)}(2)\right]^{(k)}\|(j_1'j_2')J'\rangle.
$$

The reduced matrix element is

$$
\boxed{\langle(j_1j_2)J\|\left[\hat A^{(k_1)}(1)\otimes\hat B^{(k_2)}(2)\right]^{(k)}\|(j_1'j_2')J'\rangle=\hat J\hat J'\hat k\begin{Bmatrix}j_1&j_2&J\\j_1'&j_2'&J'\\k_1&k_2&k\end{Bmatrix}\langle j_1\|\hat A^{(k_1)}\|j_1'\rangle\langle j_2\|\hat B^{(k_2)}\|j_2'\rangle}.
$$

##### Scalar Product of Tensor Operators

For equal integer ranks $k$, the rank-$0$ coupling has the explicit structure

$$
\begin{aligned}&\langle(j_1j_2)JM|\left[\hat A^{(k)}(1)\otimes\hat B^{(k)}(2)\right]_0^{(0)}|(j_1'j_2')J'M'\rangle\\&=\frac{1}{\hat k}\sum_{m_1,m_2,m_1',m_2',q}(-1)^{k-q}\langle JM|j_1m_1;j_2m_2\rangle\langle j_1'm_1';j_2'm_2'|J'M'\rangle\times\langle j_1m_1|\hat A_q^{(k)}|j_1'm_1'\rangle\langle j_2m_2|\hat B_{-q}^{(k)}|j_2'm_2'\rangle.\end{aligned}
$$

Reducing the subsystem matrix elements gives

$$
\boxed{\langle(j_1j_2)JM|\left[\hat A^{(k)}(1)\otimes\hat B^{(k)}(2)\right]_0^{(0)}|(j_1'j_2')J'M'\rangle=\delta_{JJ'}\delta_{MM'}\frac{(-1)^{j_1'+j_2+J+k}}{\hat k}\begin{Bmatrix}j_1&j_2&J\\j_2'&j_1'&k\end{Bmatrix}\langle j_1\|\hat A^{(k)}\|j_1'\rangle\langle j_2\|\hat B^{(k)}\|j_2'\rangle}.
$$

Using

$$
\hat A^{(k)}\cdot\hat B^{(k)}=\sum_{q=-k}^{k}(-1)^q\hat A_q^{(k)}\hat B_{-q}^{(k)}=(-1)^k\hat k\left[\hat A^{(k)}\otimes\hat B^{(k)}\right]_0^{(0)},
$$

the scalar-product matrix element becomes

$$
\boxed{\langle(j_1j_2)JM|\hat A^{(k)}(1)\cdot\hat B^{(k)}(2)|(j_1'j_2')JM\rangle=(-1)^{j_1'+j_2+J}\begin{Bmatrix}j_1&j_2&J\\j_2'&j_1'&k\end{Bmatrix}\langle j_1\|\hat A^{(k)}\|j_1'\rangle\langle j_2\|\hat B^{(k)}\|j_2'\rangle}.
$$
