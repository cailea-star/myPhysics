# Spherical Coordinates

Spherical coordinates specialize the general framework of [Orthogonal Curvilinear Coordinates](orthogonal_coordinates.md).

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
\nabla_\Omega=\hat{\boldsymbol\theta}\frac{\partial}{\partial\theta}+\hat{\boldsymbol\phi}\frac{1}{\sin\theta}\frac{\partial}{\partial\phi},\qquad\nabla=\hat{\mathbf r}\frac{\partial}{\partial r}+\frac{1}{r}\nabla_\Omega.
$$

Let $\hbar$ be the reduced Planck constant. Then

$$
\hat{\mathbf L}=-i\hbar\,\hat{\mathbf r}\times\nabla_\Omega,\qquad\nabla_\Omega=-\frac{i}{\hbar}\hat{\mathbf r}\times\hat{\mathbf L},\qquad\hat L_z=-i\hbar\frac{\partial}{\partial\phi}.
$$

Then the gradient operator can be expressed as:

$$
\boxed{\nabla = \hat{\mathbf r}\frac{\partial}{\partial r}+\frac{1}{r}\nabla_\Omega = \hat{\mathbf r}\frac{\partial}{\partial r}-\frac{i}{\hbar}\frac{\hat{\mathbf r}\times\hat{\mathbf L}}{r}}.
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

Let $l=0,1,\ldots$, $m=-l,\ldots,l$, and let $Y_{lm}(\Omega)$ be simultaneous eigenfunctions of $\hat{\mathbf L}^2$ and $\hat L_z$:

$$
\hat{\mathbf L}^2Y_{lm}=\hbar^2l(l+1)Y_{lm},\qquad \hat L_zY_{lm}=\hbar mY_{lm}.
$$

For a radial function $u_{lm}(r)$,

$$
\boxed{\nabla^2\left[\frac{u_{lm}(r)}{r}Y_{lm}(\Omega)\right]=\frac{1}{r}\left[\frac{d^2}{dr^2}-\frac{l(l+1)}{r^2}\right]u_{lm}(r)Y_{lm}(\Omega)}.
$$

Definitions and further properties are collected in [Spherical Harmonics](spherical_harmonics.md).

### Partial-Wave and Multipole Expansions

##### Angular Expansions

A general scalar function on the unit sphere admits the spherical-harmonic expansion

$$
\boxed{f(\Omega)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}f_{lm}Y_{lm}(\Omega),\qquad f_{lm}=\int d\Omega\,Y_{lm}^*(\Omega)f(\Omega)}.
$$

If $f$ is independent of $\phi$, only the $m=0$ components remain:

$$
\boxed{f(\theta)=\sum_{l=0}^{\infty}a_lP_l(\cos\theta),\qquad a_l=\frac{2l+1}{2}\int_0^\pi d\theta\,\sin\theta\,P_l(\cos\theta)f(\theta)}.
$$

A rotationally invariant angular kernel depends only on $\cos\gamma=\hat{\mathbf r}'\cdot\hat{\mathbf r}$. Since

$$
\int d\Omega'\,P_l(\cos\gamma)P_{l'}(\cos\gamma)=\frac{4\pi}{2l+1}\delta_{ll'}.
$$

its partial-wave expansion is

$$
\boxed{K(\Omega',\Omega)=\sum_{l=0}^{\infty}\frac{2l+1}{4\pi}K_lP_l(\cos\gamma),\qquad K_l=\int d\Omega'\,P_l(\cos\gamma)K(\cos\gamma)}.
$$

##### Plane-Wave Expansion

Let $\mathbf k=k\hat{\mathbf k}$ and $\mathbf r=r\hat{\mathbf r}$. With $j_l(x)$ denoting the spherical Bessel function and $P_l(x)$ the Legendre polynomial,

$$
\boxed{e^{i\mathbf k\cdot\mathbf r}=\sum_{l=0}^{\infty}(2l+1)i^l j_l(kr)P_l(\hat{\mathbf k}\cdot\hat{\mathbf r})=4\pi\sum_{l=0}^{\infty}\sum_{m=-l}^{l}i^l j_l(kr)Y_{lm}(\hat{\mathbf r})Y_{lm}^*(\hat{\mathbf k})}.
$$

##### Coulomb-Kernel Expansion

Let $r_i=|\mathbf r_i|$ and define

$$
r_<=\min(r_1,r_2),\qquad r_>=\max(r_1,r_2),\qquad \cos\gamma_{12}=\hat{\mathbf r}_1\cdot\hat{\mathbf r}_2.
$$

For $r_</r_><1$,

$$
\boxed{\frac{1}{|\mathbf r_1-\mathbf r_2|}=\sum_{\lambda=0}^{\infty}\frac{r_<^\lambda}{r_>^{\lambda+1}}P_\lambda(\cos\gamma_{12})=4\pi\sum_{\lambda=0}^{\infty}\sum_{\mu=-\lambda}^{\lambda}\frac{1}{2\lambda+1}\frac{r_<^\lambda}{r_>^{\lambda+1}}Y_{\lambda\mu}(\hat{\mathbf r}_1)Y_{\lambda\mu}^*(\hat{\mathbf r}_2)}.
$$
