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

Let $l=0,1,\ldots$, $m=-l,\ldots,l$, and let $Y_{lm}(\Omega)$ be simultaneous eigenfunctions of $\hat{\mathbf L}^2$ and $\hat L_z$:

$$
\hat{\mathbf L}^2Y_{lm}=\hbar^2l(l+1)Y_{lm},\qquad \hat L_zY_{lm}=\hbar mY_{lm}.
$$

For a radial function $u_{lm}(r)$,

$$
\boxed{\nabla^2\left[\frac{u_{lm}(r)}{r}Y_{lm}(\Omega)\right]=\frac{1}{r}\left[\frac{d^2}{dr^2}-\frac{l(l+1)}{r^2}\right]u_{lm}(r)Y_{lm}(\Omega)}.
$$

Definitions and further properties are collected in [Spherical Harmonics](spherical_harmonics.md).
