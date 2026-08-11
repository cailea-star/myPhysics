# Cylindrical Coordinates

Cylindrical coordinates specialize the general framework of [Orthogonal Curvilinear Coordinates](orthogonal_coordinates.md).

### Coordinate Differentials

Let $(r_\perp,\phi,z)$ be cylindrical coordinates defined by

$$
x=r_\perp\cos\phi,\qquad y=r_\perp\sin\phi,
$$

where $r_\perp\geq0$, $0\leq\phi<2\pi$, and $z\in\mathbb R$. The position vector is

$$
\mathbf R=r_\perp\cos\phi\,\hat{\mathbf x}+r_\perp\sin\phi\,\hat{\mathbf y}+z\hat{\mathbf z}=r_\perp\hat{\mathbf r}_\perp+z\hat{\mathbf z}.
$$

The local orthonormal basis is

$$
\hat{\mathbf r}_\perp=\cos\phi\,\hat{\mathbf x}+\sin\phi\,\hat{\mathbf y},\qquad \hat{\boldsymbol\phi}=-\sin\phi\,\hat{\mathbf x}+\cos\phi\,\hat{\mathbf y},\qquad \hat{\mathbf z}=\hat{\mathbf z}.
$$

The coordinate-line tangent vectors are

$$
\frac{\partial\mathbf R}{\partial r_\perp}=\hat{\mathbf r}_\perp,\qquad \frac{\partial\mathbf R}{\partial\phi}=r_\perp\hat{\boldsymbol\phi},\qquad \frac{\partial\mathbf R}{\partial z}=\hat{\mathbf z}.
$$

Only the azimuthal coordinate changes the local basis:

$$
\frac{\partial\hat{\mathbf r}_\perp}{\partial\phi}=\hat{\boldsymbol\phi},\qquad \frac{\partial\hat{\boldsymbol\phi}}{\partial\phi}=-\hat{\mathbf r}_\perp,\qquad \frac{\partial\hat{\mathbf z}}{\partial\phi}=0.
$$

The basis is right-handed:

$$
\hat{\mathbf r}_\perp\times\hat{\boldsymbol\phi}=\hat{\mathbf z}.
$$

The scale factors are

$$
h_{r_\perp}=1,\qquad h_\phi=r_\perp,\qquad h_z=1.
$$

Thus,

$$
\boxed{d\mathbf R=\hat{\mathbf r}_\perp\,dr_\perp+\hat{\boldsymbol\phi}\,r_\perp d\phi+\hat{\mathbf z}\,dz}.
$$

The metric and volume element are

$$
ds^2=dr_\perp^2+r_\perp^2d\phi^2+dz^2,\qquad d^3\mathbf R=r_\perp\,dr_\perp\,d\phi\,dz.
$$

### Differential Operators

In cylindrical coordinates,

$$
\nabla=\hat{\mathbf r}_\perp\frac{\partial}{\partial r_\perp}+\hat{\boldsymbol\phi}\frac{1}{r_\perp}\frac{\partial}{\partial\phi}+\hat{\mathbf z}\frac{\partial}{\partial z}.
$$

$$
\nabla^2=\frac{1}{r_\perp}\frac{\partial}{\partial r_\perp}\left(r_\perp\frac{\partial}{\partial r_\perp}\right)+\frac{1}{r_\perp^2}\frac{\partial^2}{\partial\phi^2}+\frac{\partial^2}{\partial z^2}.
$$

Define the azimuthal derivative operator by

$$
\nabla_\phi=\hat{\boldsymbol\phi}\frac{\partial}{\partial\phi},\qquad \nabla=\hat{\mathbf r}_\perp\frac{\partial}{\partial r_\perp}+\frac{1}{r_\perp}\nabla_\phi+\hat{\mathbf z}\frac{\partial}{\partial z}.
$$

Let $\hbar$ be the reduced Planck constant. Then

$$
\hat L_z=-i\hbar\frac{\partial}{\partial\phi},\qquad \nabla_\phi=\frac{i}{\hbar}\hat{\boldsymbol\phi}\hat L_z.
$$

Hence,

$$
\boxed{\nabla=\hat{\mathbf r}_\perp\frac{\partial}{\partial r_\perp}+\hat{\mathbf z}\frac{\partial}{\partial z}+\frac{i}{\hbar r_\perp}\hat{\boldsymbol\phi}\hat L_z}.
$$

The azimuthal Laplacian is

$$
\nabla_\phi^2=\frac{\partial^2}{\partial\phi^2}=-\frac{\hat L_z^2}{\hbar^2}.
$$

Therefore,

$$
\boxed{\nabla^2=\frac{1}{r_\perp}\frac{\partial}{\partial r_\perp}\left(r_\perp\frac{\partial}{\partial r_\perp}\right)+\frac{1}{r_\perp^2}\nabla_\phi^2+\frac{\partial^2}{\partial z^2}=\frac{1}{r_\perp}\frac{\partial}{\partial r_\perp}\left(r_\perp\frac{\partial}{\partial r_\perp}\right)+\frac{\partial^2}{\partial z^2}-\frac{\hat L_z^2}{\hbar^2r_\perp^2}}.
$$

### Azimuthal Harmonics

Let $m\in\mathbb Z$ and define the normalized azimuthal harmonic

$$
\Phi_m(\phi)=\frac{1}{\sqrt{2\pi}}e^{im\phi}.
$$

It satisfies

$$
\hat L_z\Phi_m=\hbar m\Phi_m,\qquad \nabla_\phi^2\Phi_m=-m^2\Phi_m,\qquad \int_0^{2\pi}d\phi\,\Phi_m^*(\phi)\Phi_{m'}(\phi)=\delta_{mm'}.
$$

A single-valued scalar function can be expanded as

$$
f(r_\perp,\phi,z)=\sum_{m=-\infty}^{\infty}f_m(r_\perp,z)\Phi_m(\phi),\qquad f_m(r_\perp,z)=\int_0^{2\pi}d\phi\,\Phi_m^*(\phi)f(r_\perp,\phi,z).
$$

Let $J_m(x)$ denote the Bessel function of the first kind, with $k_\perp\geq0$ and $k_z\in\mathbb R$. The cylindrical-wave modes are

$$
\Psi_{mk_\perp k_z}(r_\perp,\phi,z)=J_m(k_\perp r_\perp)\Phi_m(\phi)e^{ik_z z}.
$$

Thus,

$$
\boxed{\nabla^2\Psi_{mk_\perp k_z}=-\left(k_\perp^2+k_z^2\right)\Psi_{mk_\perp k_z}}.
$$

Let $\phi_k$ denote the azimuthal angle of the transverse wave vector. Since

$$
\mathbf k\cdot\mathbf R=k_\perp r_\perp\cos(\phi-\phi_k)+k_z z,
$$

the cylindrical partial-wave expansion is

$$
\boxed{e^{i\mathbf k\cdot\mathbf R}=e^{ik_z z}\sum_{m=-\infty}^{\infty}i^mJ_m(k_\perp r_\perp)e^{im(\phi-\phi_k)}}.
$$
