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
  \boxed{\sum_{m=-l}^{l}Y_{lm}(\Omega)Y_{lm}^*(\Omega')=\frac{2l+1}{4\pi}P_l(\cos\gamma)}.
  $$

- Rotation law, where $D_{m'm}^{(l)}(\mathcal R)$ is the Wigner rotation matrix

  $$
  Y_{lm}(\mathcal R^{-1}\hat{\mathbf r})=\sum_{m'=-l}^{l}D_{m'm}^{(l)}(\mathcal R)Y_{lm'}(\hat{\mathbf r}).
  $$

- Gaunt integral

  $$
  \int d\Omega\,Y_{l_1m_1}Y_{l_2m_2}Y_{l_3m_3}=\sqrt{\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi}}\begin{pmatrix}l_1&l_2&l_3\\0&0&0\end{pmatrix}\begin{pmatrix}l_1&l_2&l_3\\m_1&m_2&m_3\end{pmatrix}.
  $$

  The matrices are Wigner $3j$ symbols. The integral vanishes unless $m_1+m_2+m_3=0$, $|l_1-l_2|\leq l_3\leq l_1+l_2$, and $l_1+l_2+l_3$ is even.

A scalar function can be expanded as

$$
\boxed{f(r,\Omega)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\frac{u_{lm}(r)}{r}Y_{lm}(\Omega),\qquad u_{lm}(r)=r\int d\Omega\,Y_{lm}^*(\Omega)f(r,\Omega)}.
$$
