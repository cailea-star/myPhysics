# Spherical Coordinates

### Coordinate Differentials

Let $(r,\theta,\phi)$ be spherical coordinates defined by

$$
x=r\sin\theta\cos\phi,\qquad y=r\sin\theta\sin\phi,\qquad z=r\cos\theta,
$$

where $r\geq0$, $0\leq\theta\leq\pi$, and $0\leq\phi<2\pi$. The position vector is

$$
\mathbf R=x\hat{\mathbf x}+y\hat{\mathbf y}+z\hat{\mathbf z}.
$$

The coordinate-line tangent vectors are

$$
\frac{\partial\mathbf R}{\partial r}=\sin\theta\cos\phi\,\hat{\mathbf x}+\sin\theta\sin\phi\,\hat{\mathbf y}+\cos\theta\,\hat{\mathbf z}=\hat{\mathbf r}.
$$

$$
\frac{\partial\mathbf R}{\partial\theta}=r\left(\cos\theta\cos\phi\,\hat{\mathbf x}+\cos\theta\sin\phi\,\hat{\mathbf y}-\sin\theta\,\hat{\mathbf z}\right)=r\hat{\boldsymbol\theta}.
$$

$$
\frac{\partial\mathbf R}{\partial\phi}=r\sin\theta\left(-\sin\phi\,\hat{\mathbf x}+\cos\phi\,\hat{\mathbf y}\right)=r\sin\theta\hat{\boldsymbol\phi}.
$$

The basis is right-handed:

$$
\hat{\mathbf r}\times\hat{\boldsymbol\theta}=\hat{\boldsymbol\phi}.
$$

Thus,

$$
\boxed{d\mathbf R=\hat{\mathbf r}\,dr+\hat{\boldsymbol\theta}\,r\,d\theta+\hat{\boldsymbol\phi}\,r\sin\theta\,d\phi}.
$$

With $d\Omega$ denoting the solid-angle element,

$$
ds^2=dr^2+r^2d\theta^2+r^2\sin^2\theta\,d\phi^2,\qquad d\Omega=\sin\theta\,d\theta\,d\phi,\qquad d^3\mathbf R=r^2\,dr\,d\Omega.
$$

### Differential Operators

In spherical coordinates,

$$
\nabla=\hat{\mathbf r}\frac{\partial}{\partial r}+\hat{\boldsymbol\theta}\frac{1}{r}\frac{\partial}{\partial\theta}+\hat{\boldsymbol\phi}\frac{1}{r\sin\theta}\frac{\partial}{\partial\phi}.
$$

$$
\nabla^2=\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)+\frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{r^2\sin^2\theta}\frac{\partial^2}{\partial\phi^2}.
$$

Let $\Omega=(\theta,\phi)$ denote the angular coordinates. Define

$$
\nabla_\Omega=\hat{\boldsymbol\theta}\frac{\partial}{\partial\theta}+\hat{\boldsymbol\phi}\frac{1}{\sin\theta}\frac{\partial}{\partial\phi},\qquad\boxed{\nabla=\hat{\mathbf r}\frac{\partial}{\partial r}+\frac{1}{r}\nabla_\Omega}.
$$

Let $\hbar$ be the reduced Planck constant. Then

$$
\hat{\mathbf L}=-i\hbar\,\hat{\mathbf r}\times\nabla_\Omega,\qquad\nabla_\Omega=-\frac{i}{\hbar}\hat{\mathbf r}\times\hat{\mathbf L},\qquad\hat L_z=-i\hbar\frac{\partial}{\partial\phi}.
$$

The angular Laplacian is

$$
\nabla_\Omega^2=\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{\sin^2\theta}\frac{\partial^2}{\partial\phi^2}=-\frac{\hat{\mathbf L}^2}{\hbar^2}.
$$

Hence,

$$
\boxed{\nabla^2=\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)+\frac{1}{r^2}\nabla_\Omega^2=\frac{1}{r^2}\left[\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)-\frac{\hat{\mathbf L}^2}{\hbar^2}\right]}.
$$

### Spherical Harmonics

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

### Radial Expansion

Define the reduced radial coefficients of $f(\mathbf R)$ by

$$
u_{lm}(r)=r\int d\Omega\,Y_{lm}^*(\Omega)f(\mathbf R),\qquad f(\mathbf R)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\frac{u_{lm}(r)}{r}Y_{lm}(\Omega).
$$

Then

$$
\boxed{\nabla^2\left[\frac{u_{lm}(r)}{r}Y_{lm}(\Omega)\right]=\frac{1}{r}\left[\frac{d^2u_{lm}}{dr^2}-\frac{l(l+1)}{r^2}u_{lm}(r)\right]Y_{lm}(\Omega)}.
$$
