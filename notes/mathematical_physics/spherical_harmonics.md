# Spherical Harmonics

The spherical-coordinate conventions, $\Omega=(\theta,\phi)$, and orbital angular-momentum operators used below are defined in [Spherical Coordinates](spherical_coordinates.md).

### Scalar Spherical Harmonics

Let $l=0,1,\ldots$, $m=-l,\ldots,l$, and let $P_l(x)$ and $P_l^m(x)$ denote the Legendre polynomial and associated Legendre function. In the Condon–Shortley convention, for $m\geq0$,

$$
\boxed{Y_{lm}(\theta,\phi)=\sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}\,P_l^m(\cos\theta)e^{im\phi}}.
$$

The spherical harmonics satisfy:

- Eigenvalue equations

  $$
  \hat{\mathbf L}^2Y_{lm}=\hbar^2l(l+1)Y_{lm},\qquad\hat L_zY_{lm}=\hbar mY_{lm}.
  $$

- Conjugation and parity

  $$
  Y_{l,-m}(\Omega)=(-1)^mY_{lm}^*(\Omega),\qquad Y_{lm}(-\hat{\mathbf r})=(-1)^lY_{lm}(\hat{\mathbf r}).
  $$

- Orthonormality and completeness

  $$
  \int d\Omega\,Y_{lm}^*(\Omega)Y_{l'm'}(\Omega)=\delta_{ll'}\delta_{mm'},\qquad\sum_{l=0}^{\infty}\sum_{m=-l}^{l}Y_{lm}(\Omega)Y_{lm}^*(\Omega')=\delta(\Omega-\Omega').
  $$

- Ladder relations, with $\hat L_\pm=\hat L_x\pm i\hat L_y$

  $$
  \hat L_\pm Y_{lm}=\hbar\sqrt{l(l+1)-m(m\pm1)}\,Y_{l,m\pm1}.
  $$

- Addition theorem, where $\cos\gamma=\hat{\mathbf r}\cdot\hat{\mathbf r}'$

  $$
  \sum_{m=-l}^{l}Y_{lm}(\Omega)Y_{lm}^*(\Omega')=\frac{2l+1}{4\pi}P_l(\cos\gamma).
  $$

- Gaunt integral

  $$
  \int d\Omega\,Y_{l_1m_1}Y_{l_2m_2}Y_{l_3m_3}=\sqrt{\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi}}\begin{pmatrix}l_1&l_2&l_3\\0&0&0\end{pmatrix}\begin{pmatrix}l_1&l_2&l_3\\m_1&m_2&m_3\end{pmatrix}.
  $$

  The matrices are Wigner $3j$ symbols. The integral vanishes unless $m_1+m_2+m_3=0$, $|l_1-l_2|\leq l_3\leq l_1+l_2$, and $l_1+l_2+l_3$ is even.

A scalar function can be expanded as

$$
f(r,\Omega)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\frac{u_{lm}(r)}{r}Y_{lm}(\Omega),\qquad u_{lm}(r)=r\int d\Omega\,Y_{lm}^*(\Omega)f(r,\Omega).
$$

### Vector Spherical Harmonics

Let $\mathbf e_0$ and $\mathbf e_\pm$ denote the spherical vector basis,

$$
\mathbf e_0=\hat{\mathbf z},\qquad \mathbf e_\pm=\mp\frac{\hat{\mathbf x}\pm i\hat{\mathbf y}}{\sqrt2}.
$$

Coupling $Y_{lm}(\Omega)$ to this spin-$1$ basis defines

$$
\boxed{\mathbf Y_{JlM}(\Omega)=\sum_{m=-l}^{l}\sum_{\mu=-1}^{1}\langle l,m;1,\mu|J,M\rangle Y_{lm}(\Omega)\mathbf e_\mu},
$$

where $\langle l,m;1,\mu|J,M\rangle$ is a Clebsch–Gordan coefficient, and $\mathbf e_\mu$ denotes $\mathbf e_0$, $\mathbf e_+$, or $\mathbf e_-$ according to $\mu=0,+1,-1$:

$$
J=|l-1|,\ldots,l+1,\qquad M=-J,\ldots,J.
$$

Let $\hat{\mathbf S}$ act on the vector index and define $\hat{\mathbf J}=\hat{\mathbf L}+\hat{\mathbf S}$. The vector spherical harmonics satisfy:

- Angular-momentum eigenvalue equations

  $$
  \hat{\mathbf J}^2\mathbf Y_{JlM}=\hbar^2J(J+1)\mathbf Y_{JlM},\qquad \hat J_z\mathbf Y_{JlM}=\hbar M\mathbf Y_{JlM}.
  $$

  $$
  \hat{\mathbf L}^2\mathbf Y_{JlM}=\hbar^2l(l+1)\mathbf Y_{JlM},\qquad \hat{\mathbf S}^2\mathbf Y_{JlM}=2\hbar^2\mathbf Y_{JlM}.
  $$

- Conjugation and parity

  $$
  \mathbf Y_{JlM}^*(\Omega)=(-1)^{J+l+M+1}\mathbf Y_{Jl,-M}(\Omega),\qquad \mathbf Y_{JlM}(-\hat{\mathbf r})=(-1)^l\mathbf Y_{JlM}(\hat{\mathbf r}).
  $$

  The second relation transforms the argument while keeping the Cartesian basis fixed. A polar vector acquires one additional minus sign under the parity operator.

- Orthonormality and completeness

  $$
  \int d\Omega\,\mathbf Y_{J'l'M'}^*(\Omega)\cdot\mathbf Y_{JlM}(\Omega)=\delta_{JJ'}\delta_{ll'}\delta_{MM'}.
  $$

  $$
  \sum_{l=0}^{\infty}\sum_{J=|l-1|}^{l+1}\sum_{M=-J}^{J}\left[\mathbf Y_{JlM}(\Omega)\right]_a\left[\mathbf Y_{JlM}(\Omega')\right]_b^*=\delta_{ab}\delta(\Omega-\Omega').
  $$

- For $J\geq1$, an alternative basis adapted to angular differential operators is

  $$
  \mathbf Y_{JM}^{(r)}=\hat{\mathbf r}Y_{JM},\qquad \mathbf\Psi_{JM}=\frac{\nabla_\Omega Y_{JM}}{\sqrt{J(J+1)}},\qquad \mathbf\Phi_{JM}=\frac{\hat{\mathbf r}\times\nabla_\Omega Y_{JM}}{\sqrt{J(J+1)}}.
  $$

For fixed $(J,M)$, the coupled and differential basis expansions are

$$
\mathbf F_{JM}(r,\Omega)=\sum_{l=|J-1|}^{J+1}F_{JlM}(r)\mathbf Y_{JlM}(\Omega).
$$

$$
\mathbf F_{JM}(r,\Omega)=F_{JM}^{(r)}(r)\mathbf Y_{JM}^{(r)}(\Omega)+F_{JM}^{(\Psi)}(r)\mathbf\Psi_{JM}(\Omega)+F_{JM}^{(\Phi)}(r)\mathbf\Phi_{JM}(\Omega).
$$

For $J\geq1$, the two bases span the same three-dimensional space and are related by a unitary change of basis. For $J=0$, only $\mathbf Y_{00}^{(r)}$ remains.
