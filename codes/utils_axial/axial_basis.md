# Axial Harmonic-Oscillator Basis

The state labels and oscillator lengths follow [Axial Harmonic-Oscillator Configuration](axial_config.md). The basis combines the one-dimensional and polar [harmonic-oscillator solutions](../../notes/quantum_mechanics/harmonic_oscillator_coordinate_space.md), with [Hermite](../../notes/mathematical_physics/specialfunction_hermite_polynomials.md) and [generalized Laguerre polynomials](../../notes/mathematical_physics/specialfunction_generalized_laguerre_polynomials.md).

### Axial Harmonic-Oscillator Hamiltonian

Let $m$ be the particle mass and $\omega_z,\omega_r>0$ the axial and radial oscillator frequencies. In Cartesian coordinates,

$$
\hat H_0=-\frac{\hbar^2}{2m}\left(\partial_x^2+\partial_y^2+\partial_z^2\right)+\frac12m\omega_r^2(x^2+y^2)+\frac12m\omega_z^2z^2.
$$

Using

$$
x=r\cos\varphi,\qquad y=r\sin\varphi,\qquad \partial_x^2+\partial_y^2=\partial_r^2+\frac1r\partial_r+\frac1{r^2}\partial_\varphi^2,
$$

the Hamiltonian becomes

$$
\hat H_0=-\frac{\hbar^2}{2m}\left(\partial_z^2+\partial_r^2+\frac1r\partial_r+\frac1{r^2}\partial_\varphi^2\right)+\frac12m\omega_z^2z^2+\frac12m\omega_r^2r^2.
$$

For $n_z,n_r,\Lambda\in\mathbb N_0$, the energy spectrum is

$$
\boxed{E(n_z,n_r,\Lambda)=\hbar\omega_z\left(n_z+\frac12\right)+\hbar\omega_r(2n_r+\Lambda+1).}
$$

### Axial Basis Functions along the $z$ Direction

Following the [one-dimensional harmonic-oscillator solution](../../notes/quantum_mechanics/harmonic_oscillator_coordinate_space.md#one-dimensional-harmonic-oscillator), define

$$
b_z=\sqrt{\frac{\hbar}{m\omega_z}},\qquad \zeta=\frac{z}{b_z}.
$$

For $n_z\in\mathbb N_0$, the normalized axial basis function is

$$
\boxed{N_{n_z}=\frac{1}{\sqrt{2^{n_z}n_z!\sqrt\pi}},\qquad \phi_{n_z}(z)=\frac{N_{n_z}}{\sqrt{b_z}}H_{n_z}(\zeta)e^{-\zeta^2/2}.}
$$

It satisfies

$$
\int_{-\infty}^{\infty}\phi_{n_z}(z)\phi_{n_z'}(z)\,dz=\delta_{n_zn_z'}.
$$

The coordinate derivatives are

$$
\partial_z=\frac1{b_z}\partial_\zeta,\qquad \partial_z^2=\frac1{b_z^2}\partial_\zeta^2.
$$

The required [Hermite-polynomial](../../notes/mathematical_physics/specialfunction_hermite_polynomials.md) derivatives are

$$
\partial_\zeta H_{n_z}(\zeta)=2n_zH_{n_z-1}(\zeta),\qquad \partial_\zeta^2H_{n_z}(\zeta)=2\zeta\partial_\zeta H_{n_z}(\zeta)-2n_zH_{n_z}(\zeta).
$$

Therefore,

$$
\boxed{\partial_z\phi_{n_z}(z)=\frac{N_{n_z}}{b_z^{3/2}}\left[\partial_\zeta H_{n_z}(\zeta)-\zeta H_{n_z}(\zeta)\right]e^{-\zeta^2/2}.}
$$

$$
\boxed{\partial_z^2\phi_{n_z}(z)=\frac{N_{n_z}}{b_z^{5/2}}\left[\partial_\zeta^2H_{n_z}(\zeta)-2\zeta\partial_\zeta H_{n_z}(\zeta)+(\zeta^2-1)H_{n_z}(\zeta)\right]e^{-\zeta^2/2}.}
$$

### Radial Basis Functions along the $r$ Direction

Following the [two-dimensional polar harmonic-oscillator solution](../../notes/quantum_mechanics/harmonic_oscillator_coordinate_space.md#two-dimensional-isotropic-harmonic-oscillator), define

$$
b_r=\sqrt{\frac{\hbar}{m\omega_r}},\qquad \eta=\left(\frac{r}{b_r}\right)^2.
$$

For $n_r,\Lambda\in\mathbb N_0$, the normalized radial basis function is

$$
\boxed{N_{n_r}^{\Lambda}=\sqrt{\frac{n_r!}{(n_r+\Lambda)!}},\qquad \phi_{n_r}^{\Lambda}(r)=\frac{\sqrt2}{b_r}N_{n_r}^{\Lambda}L_{n_r}^{(\Lambda)}(\eta)\eta^{\Lambda/2}e^{-\eta/2}.}
$$

It satisfies

$$
\int_0^\infty r\,dr\,\phi_{n_r}^{\Lambda}(r)\phi_{n_r'}^{\Lambda}(r)=\delta_{n_rn_r'}.
$$

The coordinate derivatives are

$$
\partial_r=\frac{2r}{b_r^2}\partial_\eta,\qquad \partial_r^2=\frac{4r^2}{b_r^4}\partial_\eta^2+\frac{2}{b_r^2}\partial_\eta.
$$

The required [generalized Laguerre-polynomial](../../notes/mathematical_physics/specialfunction_generalized_laguerre_polynomials.md) derivatives are

$$
\partial_\eta L_{n_r}^{(\Lambda)}(\eta)=-L_{n_r-1}^{(\Lambda+1)}(\eta),\qquad \eta\partial_\eta^2L_{n_r}^{(\Lambda)}(\eta)=(\eta-\Lambda-1)\partial_\eta L_{n_r}^{(\Lambda)}(\eta)-n_rL_{n_r}^{(\Lambda)}(\eta).
$$

Define

$$
u_{n_r}^{\Lambda}(\eta)=L_{n_r}^{(\Lambda)}(\eta)\eta^{\Lambda/2}e^{-\eta/2}.
$$

Its derivatives are

$$
\partial_\eta u_{n_r}^{\Lambda}=\left[\partial_\eta L_{n_r}^{(\Lambda)}+\frac{\Lambda-\eta}{2\eta}L_{n_r}^{(\Lambda)}\right]\eta^{\Lambda/2}e^{-\eta/2}.
$$

$$
\partial_\eta^2u_{n_r}^{\Lambda}=\left[\partial_\eta^2L_{n_r}^{(\Lambda)}+\frac{\Lambda-\eta}{\eta}\partial_\eta L_{n_r}^{(\Lambda)}+\frac{(\Lambda-\eta)^2-2\Lambda}{4\eta^2}L_{n_r}^{(\Lambda)}\right]\eta^{\Lambda/2}e^{-\eta/2}.
$$

Therefore,

$$
\boxed{\partial_r\phi_{n_r}^{\Lambda}(r)=\frac{\sqrt2N_{n_r}^{\Lambda}}{b_r}\frac{2r}{b_r^2}\partial_\eta u_{n_r}^{\Lambda}(\eta).}
$$

$$
\boxed{\partial_r^2\phi_{n_r}^{\Lambda}(r)=\frac{\sqrt2N_{n_r}^{\Lambda}}{b_r}\left[\frac{4r^2}{b_r^4}\partial_\eta^2u_{n_r}^{\Lambda}(\eta)+\frac{2}{b_r^2}\partial_\eta u_{n_r}^{\Lambda}(\eta)\right].}
$$

### Quadrature Weights

The Gauss-Hermite and Gauss-Laguerre quadratures are

$$
\int_{-\infty}^{\infty}f(\zeta)e^{-\zeta^2}\,d\zeta\approx\sum_i w_i^{\rm GH}f(\zeta_i),\qquad \int_0^\infty f(\eta)e^{-\eta}\,d\eta\approx\sum_j w_j^{\rm GL}f(\eta_j).
$$

Thus, $w_i^{\rm GH}$ and $w_j^{\rm GL}$ already represent the exponential measures $e^{-\zeta^2}d\zeta$ and $e^{-\eta}d\eta$. However, the basis-function products already contain

$$
\phi_{n_z}(z)\phi_{n_z'}(z)\propto e^{-\zeta^2},\qquad \phi_{n_r}^{\Lambda}(r)\phi_{n_r'}^{\Lambda}(r)\propto e^{-\eta}.
$$

The exponential factors must therefore be removed from the quadrature weights:

$$
w_i^{\rm GH}\rightarrow w_i^{\rm GH}e^{\zeta_i^2},\qquad w_j^{\rm GL}\rightarrow w_j^{\rm GL}e^{\eta_j}.
$$

The coordinate Jacobians are

$$
dz=b_z\,d\zeta,\qquad r\,dr=\frac{b_r^2}{2}\,d\eta.
$$

Therefore, the physical integration weights are

$$
\boxed{w_{z,i}=b_zw_i^{\rm GH}e^{\zeta_i^2},\qquad w_{r,j}=\frac{b_r^2}{2}w_j^{\rm GL}e^{\eta_j}.}
$$

### Reflection Symmetry

Since the Hermite polynomials satisfy

$$
H_{n_z}(-\zeta)=(-1)^{n_z}H_{n_z}(\zeta),
$$

the axial basis functions obey

$$
\phi_{n_z}(-z)=(-1)^{n_z}\phi_{n_z}(z),\qquad \Phi_\alpha(-z,r,\varphi,s)=(-1)^{n_z}\Phi_\alpha(z,r,\varphi,s).
$$

Their derivatives have the parities

$$
\partial_z\phi_{n_z}(-z)=(-1)^{n_z+1}\partial_z\phi_{n_z}(z),\qquad \partial_z^2\phi_{n_z}(-z)=(-1)^{n_z}\partial_z^2\phi_{n_z}(z).
$$

For an even integrand $F(-z)=F(z)$,

$$
\int_{-\infty}^{\infty}F(z)\,dz=2\int_0^\infty F(z)\,dz.
$$

For an odd integrand $F(-z)=-F(z)$,

$$
\int_{-\infty}^{\infty}F(z)\,dz=0.
$$

The half-space quadrature therefore retains the positive Gauss-Hermite nodes and doubles their physical weights:

$$
\boxed{z_i>0,\qquad w_{z,i}^{(+)}=2w_{z,i}.}
$$
