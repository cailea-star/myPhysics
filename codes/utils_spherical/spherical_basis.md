# Spherical Harmonic-Oscillator Basis

### Spherical Harmonic-Oscillator Hamiltonian

Let $M$ be the particle mass and $\omega>0$ the oscillator frequency. In Cartesian coordinates,

$$
\hat H_0=-\frac{\hbar^2}{2M}\left(\partial_x^2+\partial_y^2+\partial_z^2\right)+\frac12M\omega^2(x^2+y^2+z^2).
$$

Using spherical coordinates,

$$
x=r\sin\theta\cos\varphi,\qquad y=r\sin\theta\sin\varphi,\qquad z=r\cos\theta,
$$

and the orbital angular-momentum operator $\hat{\mathbf L}$,

$$
\partial_x^2+\partial_y^2+\partial_z^2=\partial_r^2+\frac2r\partial_r-\frac{\hat{\mathbf L}^2}{\hbar^2r^2},
$$

the Hamiltonian becomes

$$
\hat H_0=-\frac{\hbar^2}{2M}\left(\partial_r^2+\frac2r\partial_r-\frac{\hat{\mathbf L}^2}{\hbar^2r^2}\right)+\frac12M\omega^2r^2.
$$

The single-particle label and its opposite-projection partner follow [Spherical Harmonic-Oscillator Setting](spherical_setting.md):

$$
\alpha=(n,l,j,m),\qquad\bar\alpha=(n,l,j,-m).
$$

Let $\phi_{nl}(r)$ be the real radial amplitude and $\mathcal Y_{jm}^{l}(\theta,\varphi)$ the [spinor spherical harmonic](../../notes/mathematical_physics/specialfunction_spherical_harmonics.md#spinor-spherical-harmonics). The basis spinor and its opposite-projection partner are

$$
\boxed{\phi_\alpha(r,\theta,\varphi)=\phi_{nl}(r)\mathcal Y_{jm}^{l}(\theta,\varphi),\qquad \phi_{\bar\alpha}(r,\theta,\varphi)=\phi_{nl}(r)\mathcal Y_{j,-m}^{l}(\theta,\varphi).}
$$

The spin-independent Hamiltonian satisfies

$$
\hat H_0\phi_\alpha=E_\alpha\phi_\alpha,\qquad \hat H_0\phi_{\bar\alpha}=E_\alpha\phi_{\bar\alpha},
$$

with

$$
\boxed{E_\alpha=\hbar\omega\left(2n+l+\frac32\right).}
$$

### Radial Basis Functions

Following the [three-dimensional harmonic-oscillator solution](../../notes/quantum_mechanics/harmonic_oscillator_coordinate_space.md#three-dimensional-isotropic-harmonic-oscillator), define

$$
b=\sqrt{\frac{\hbar}{M\omega}},\qquad \eta=\left(\frac{r}{b}\right)^2.
$$

For $n,l\in\mathbb N_0$, the normalized radial basis function is

$$
\boxed{\mathcal N_{nl}=\sqrt{\frac{2n!}{\Gamma(n+l+3/2)}},\qquad \phi_{nl}(r)=\frac{\mathcal N_{nl}}{b^{3/2}}L_n^{(l+1/2)}(\eta)\eta^{l/2}e^{-\eta/2}.}
$$

It satisfies

$$
\int_0^\infty r^2\,dr\,\phi_{nl}(r)\phi_{n'l}(r)=\delta_{nn'}.
$$

The coordinate derivatives are

$$
\partial_r=\frac{2r}{b^2}\partial_\eta,\qquad \partial_r^2=\frac{4r^2}{b^4}\partial_\eta^2+\frac2{b^2}\partial_\eta.
$$

The required [generalized Laguerre-polynomial](../../notes/mathematical_physics/specialfunction_generalized_laguerre_polynomials.md) derivatives are

$$
\partial_\eta L_n^{(l+1/2)}(\eta)=-L_{n-1}^{(l+3/2)}(\eta),\qquad \eta\partial_\eta^2L_n^{(l+1/2)}(\eta)=\left(\eta-l-\frac32\right)\partial_\eta L_n^{(l+1/2)}(\eta)-nL_n^{(l+1/2)}(\eta).
$$

For $n=0$, the first derivative vanishes; set $L_{-1}^{(l+3/2)}=0$.

Define the dimensionless radial amplitude

$$
u_{nl}(\eta)=L_n^{(l+1/2)}(\eta)\eta^{l/2}e^{-\eta/2}.
$$

Its derivatives are

$$
\partial_\eta u_{nl}=\left[\partial_\eta L_n^{(l+1/2)}+\frac{l-\eta}{2\eta}L_n^{(l+1/2)}\right]\eta^{l/2}e^{-\eta/2}.
$$

$$
\partial_\eta^2u_{nl}=\left[\partial_\eta^2L_n^{(l+1/2)}+\frac{l-\eta}{\eta}\partial_\eta L_n^{(l+1/2)}+\frac{(l-\eta)^2-2l}{4\eta^2}L_n^{(l+1/2)}\right]\eta^{l/2}e^{-\eta/2}.
$$

Therefore,

$$
\boxed{\partial_r\phi_{nl}(r)=\frac{\mathcal N_{nl}}{b^{3/2}}\frac{2r}{b^2}\partial_\eta u_{nl}(\eta).}
$$

$$
\boxed{\partial_r^2\phi_{nl}(r)=\frac{\mathcal N_{nl}}{b^{3/2}}\left[\frac{4r^2}{b^4}\partial_\eta^2u_{nl}(\eta)+\frac2{b^2}\partial_\eta u_{nl}(\eta)\right].}
$$

The expressions containing $1/\eta$ apply for $r>0$; values at the origin follow from the $r\to0$ limits of $\phi_{nl}(r)$.

### Spinor Spherical Harmonics

Following the [spinor spherical-harmonic convention](../../notes/mathematical_physics/specialfunction_spherical_harmonics.md#spinor-spherical-harmonics), let $Y_{lm_l}$ be the normalized scalar spherical harmonic and $\chi_\mu$ the normalized spin state with projection $\mu=\pm\frac12$. Coupling orbital angular momentum $l$ to spin $\frac12$ gives

$$
\boxed{\mathcal Y_{jm}^{l}(\theta,\varphi)=\sum_{\mu=\pm1/2}\left\langle l,m-\mu;\frac12,\mu\middle|j,m\right\rangle Y_{l,m-\mu}(\theta,\varphi)\chi_\mu.}
$$

Only terms with $|m-\mu|\leq l$ contribute. The Clebsch–Gordan coefficients use the Condon–Shortley convention.

The angular-momentum operators satisfy

$$
\hat{\mathbf L}^2\mathcal Y_{jm}^{l}=\hbar^2l(l+1)\mathcal Y_{jm}^{l},\qquad \hat{\mathbf J}^2\mathcal Y_{jm}^{l}=\hbar^2j(j+1)\mathcal Y_{jm}^{l},\qquad \hat J_z\mathcal Y_{jm}^{l}=\hbar m\mathcal Y_{jm}^{l},
$$

where $\hat{\mathbf J}=\hat{\mathbf L}+\hat{\mathbf S}$ and $\hat{\mathbf S}$ is the spin operator.

With $d\Omega=\sin\theta\,d\theta\,d\varphi$, angular orthonormality is

$$
\int d\Omega\,\left[\mathcal Y_{jm}^{l}(\theta,\varphi)\right]^\dagger\mathcal Y_{j'm'}^{l'}(\theta,\varphi)=\delta_{ll'}\delta_{jj'}\delta_{mm'}.
$$

Together with radial orthonormality, the full basis satisfies

$$
\boxed{\int_0^\infty r^2\,dr\int d\Omega\,\phi_\alpha^\dagger(r,\theta,\varphi)\phi_{\alpha'}(r,\theta,\varphi)=\delta_{nn'}\delta_{ll'}\delta_{jj'}\delta_{mm'}.}
$$

### Quadrature Weights

Using $\eta=(r/b)^2$, generalized Gauss–Laguerre quadrature with parameter $\frac12$ gives

$$
\int_0^\infty \eta^{1/2}e^{-\eta}f(\eta)\,d\eta\approx\sum_{i=1}^{N_r}w_i^{\rm GL}f(\eta_i).
$$

The weights $w_i^{\rm GL}$ already represent the measure $\eta^{1/2}e^{-\eta}d\eta$. However, the basis-function products already contain

$$
\phi_{nl}(r)\phi_{n'l'}(r)\propto e^{-\eta}.
$$

The exponential factor must therefore be removed from the quadrature weights:

$$
w_i^{\rm GL}\rightarrow w_i^{\rm GL}e^{\eta_i}.
$$

The coordinate Jacobian is

$$
r^2\,dr=\frac{b^3}{2}\eta^{1/2}\,d\eta.
$$

Therefore, the physical radial nodes and integration weights are

$$
\boxed{r_i=b\sqrt{\eta_i},\qquad w_{r,i}=\frac{b^3}{2}w_i^{\rm GL}e^{\eta_i}.}
$$

For a radial integrand $F(r)$,

$$
\int_0^\infty r^2\,dr\,F(r)\approx\sum_{i=1}^{N_r}w_{r,i}F(r_i).
$$
