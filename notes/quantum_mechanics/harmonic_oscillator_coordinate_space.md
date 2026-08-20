# Harmonic-Oscillator Coordinate-Space Solutions

### One-Dimensional Harmonic Oscillator

The operator construction and number representation are given in [Quantum Harmonic Oscillator](harmonic_oscillator.md). Here the eigenvalue problem is solved directly in coordinate space.

##### Dimensionless Schrödinger Equation

For a one-dimensional oscillator with mass $M>0$ and angular frequency $\omega>0$, the coordinate-space Hamiltonian is

$$
\hat H=-\frac{\hbar^2}{2M}\frac{d^2}{dx^2}+\frac12M\omega^2x^2.
$$

Define the characteristic length $x_0$, dimensionless coordinate $\zeta$, and dimensionless energy $\epsilon$ by

$$
\boxed{x_0=\sqrt{\frac{\hbar}{M\omega}},\qquad \zeta=\frac{x}{x_0},\qquad \epsilon=\frac{E}{\hbar\omega}}.
$$

Write the normalized coordinate wave function as

$$
\phi(x)=\frac{1}{\sqrt{x_0}}\widetilde{\phi}(\zeta),\qquad \int_{-\infty}^{\infty}|\widetilde{\phi}(\zeta)|^2\,d\zeta=1.
$$

The stationary Schrödinger equation becomes

$$
\left[-\frac12\frac{d^2}{d\zeta^2}+\frac12\zeta^2\right]\widetilde{\phi}(\zeta)=\epsilon\widetilde{\phi}(\zeta).
$$

##### Hermite Solution and Energy Spectrum

Normalizability requires Gaussian decay. Write

$$
\widetilde{\phi}(\zeta)=e^{-\zeta^2/2}h(\zeta).
$$

The reduced function satisfies

$$
h''(\zeta)-2\zeta h'(\zeta)+(2\epsilon-1)h(\zeta)=0.
$$

The polynomial solution is the physicists' [Hermite polynomial](../mathematical_physics/specialfunction_hermite_polynomials.md):

$$
H_n''(\zeta)-2\zeta H_n'(\zeta)+2nH_n(\zeta)=0,\qquad n=0,1,2,\ldots.
$$

Polynomial termination therefore gives

$$
\epsilon(n)=n+\frac12,\qquad n=0,1,2,\ldots.
$$

The Hermite normalization coefficient and dimensionless eigenfunctions are

$$
\boxed{N_n=\frac{1}{\sqrt{2^n n!\sqrt\pi}},\qquad \widetilde{\phi}_n(\zeta)=N_nH_n(\zeta)e^{-\zeta^2/2}}.
$$

The physical eigenfunctions and energy spectrum are

$$
\boxed{\phi_n(x)=\frac{N_n}{\sqrt{x_0}}H_n\left(\frac{x}{x_0}\right)e^{-x^2/(2x_0^2)},\qquad E(n)=\hbar\omega\left(n+\frac12\right)}.
$$

They satisfy

$$
\int_{-\infty}^{\infty}\phi_{n'}^*(x)\phi_n(x)\,dx=\delta_{n'n},\qquad \sum_{n=0}^{\infty}\phi_n(x)\phi_n^*(x')=\delta(x-x').
$$

A normalized coordinate wave function can therefore be expanded as

$$
\psi(x)=\sum_{n=0}^{\infty}c_n\phi_n(x),\qquad c_n=\int_{-\infty}^{\infty}\phi_n^*(x)\psi(x)\,dx,\qquad \sum_{n=0}^{\infty}|c_n|^2=1.
$$

Since $H_n(-\zeta)=(-1)^nH_n(\zeta)$,

$$
\boxed{\phi_n(-x)=(-1)^n\phi_n(x)}.
$$

Each one-dimensional energy level is nondegenerate.

### Two-Dimensional Isotropic Harmonic Oscillator

For a two-dimensional oscillator with mass $M>0$ and angular frequency $\omega>0$, the Cartesian Hamiltonian is

$$
\hat H=-\frac{\hbar^2}{2M}\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right)+\frac12M\omega^2(x^2+y^2).
$$

Using the one-dimensional eigenfunctions above,

$$
\phi_{n_x,n_y}(x,y)=\phi_{n_x}(x)\phi_{n_y}(y),\qquad E(n_x,n_y)=\hbar\omega\left(n_x+n_y+1\right).
$$

##### Polar Dimensionless Schrödinger Equation

The same system can be expressed in polar coordinates, $x=r\cos\phi$ and $y=r\sin\phi$. The Hamiltonian then becomes

$$
\hat H=-\frac{\hbar^2}{2M}\left(\frac{\partial^2}{\partial r^2}+\frac{1}{r}\frac{\partial}{\partial r}+\frac{1}{r^2}\frac{\partial^2}{\partial\phi^2}\right)+\frac12M\omega^2r^2.
$$

A normalized coordinate wave function admits the angular expansion

$$
\psi(r,\phi)=\sum_{m=-\infty}^{\infty}\frac{e^{im\phi}}{\sqrt{2\pi}}\psi^m(r),\qquad \sum_{m=-\infty}^{\infty}\int_0^\infty r\,dr\,|\psi^m(r)|^2=1.
$$

Define the characteristic length $r_0$, dimensionless radius $\rho$, and dimensionless energy $\epsilon$ by

$$
\boxed{r_0=\sqrt{\frac{\hbar}{M\omega}},\qquad \rho=\frac{r}{r_0},\qquad \epsilon=\frac{E}{\hbar\omega}}.
$$

For a fixed $m$, the physical and dimensionless radial basis functions satisfy

$$
\phi^m(r)=\frac1{r_0}\widetilde{\phi}^m(\rho),\qquad \int_0^\infty\rho\,d\rho\,|\widetilde{\phi}^m(\rho)|^2=1.
$$

The dimensionless radial basis function satisfies

$$
\left[-\frac12\left(\frac{d^2}{d\rho^2}+\frac1\rho\frac{d}{d\rho}-\frac{m^2}{\rho^2}\right)+\frac12\rho^2\right]\widetilde{\phi}^m(\rho)=\epsilon\widetilde{\phi}^m(\rho).
$$

##### Laguerre Solution and Energy Spectrum

Since the radial equation depends only on $|m|$, regularity at $\rho=0$ and normalizability require

$$
\widetilde{\phi}^{|m|}(\rho)=\rho^{|m|}e^{-\eta/2}l(\eta),\qquad \eta=\rho^2.
$$

The reduced function satisfies

$$
\eta\frac{d^2}{d\eta^2}l(\eta)+(|m|+1-\eta)\frac{d}{d\eta}l(\eta)+\frac{\epsilon-|m|-1}{2}l(\eta)=0.
$$

The polynomial solution is the [generalized Laguerre polynomial](../mathematical_physics/specialfunction_generalized_laguerre_polynomials.md):

$$
\eta\frac{d^2}{d\eta^2}L_{n_r}^{(|m|)}(\eta)+(|m|+1-\eta)\frac{d}{d\eta}L_{n_r}^{(|m|)}(\eta)+n_rL_{n_r}^{(|m|)}(\eta)=0.
$$

Polynomial termination therefore gives

$$
\epsilon(n_r,m)=2n_r+|m|+1,\qquad n_r\in\mathbb N_0.
$$

The Laguerre normalization coefficient and dimensionless radial eigenfunctions are

$$
\boxed{N_{n_r}^{(|m|)}=\sqrt{\frac{n_r!}{(n_r+|m|)!}},\qquad \widetilde{\phi}_{n_r}^{|m|}(\rho)=\sqrt2N_{n_r}^{(|m|)}\rho^{|m|}e^{-\rho^2/2}L_{n_r}^{(|m|)}(\rho^2)}.
$$

The physical eigenfunctions and energy spectrum are

$$
\boxed{\phi_{n_r}^{|m|}(r)=\frac{\sqrt2N_{n_r}^{(|m|)}}{r_0}\left(\frac{r}{r_0}\right)^{|m|}e^{-r^2/(2r_0^2)}L_{n_r}^{(|m|)}\left(\frac{r^2}{r_0^2}\right),\qquad E(n_r,m)=\hbar\omega(2n_r+|m|+1)}.
$$

The full eigenfunctions are

$$
\phi_{n_r}^{m}(r,\phi)=\frac{e^{im\phi}}{\sqrt{2\pi}}\phi_{n_r}^{|m|}(r).
$$

The orthogonality relation is

$$
\int_0^{2\pi}d\phi\int_0^\infty r\,dr\,\left[\phi_{n_r'}^{m'}(r,\phi)\right]^*\phi_{n_r}^{m}(r,\phi)=\delta_{n_r'n_r}\delta_{m'm}.
$$

The completeness relation is

$$
\sum_{m=-\infty}^{\infty}\sum_{n_r=0}^{\infty}\phi_{n_r}^{m}(r,\phi)\left[\phi_{n_r}^{m}(r',\phi')\right]^*=\frac{\delta(r-r')\delta(\phi-\phi')}{r}.
$$

A normalized coordinate wave function can therefore be expanded as

$$
\psi(r,\phi)=\sum_{m=-\infty}^{\infty}\sum_{n_r=0}^{\infty}c_{n_r}^m\phi_{n_r}^{m}(r,\phi),\qquad c_{n_r}^m=\int_0^{2\pi}d\phi\int_0^\infty r\,dr\,\left[\phi_{n_r}^{m}(r,\phi)\right]^*\psi(r,\phi),\qquad \sum_{m=-\infty}^{\infty}\sum_{n_r=0}^{\infty}|c_{n_r}^m|^2=1.
$$

Under spatial inversion, $(r,\phi)\rightarrow(r,\phi+\pi)$, the parity is

$$
\boxed{\phi_{n_r}^{m}(r,\phi+\pi)=(-1)^m\phi_{n_r}^{m}(r,\phi)}.
$$

### Three-Dimensional Isotropic Harmonic Oscillator

For a three-dimensional oscillator with mass $M>0$ and angular frequency $\omega>0$, the Cartesian Hamiltonian is

$$
\hat H=-\frac{\hbar^2}{2M}\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}+\frac{\partial^2}{\partial z^2}\right)+\frac12M\omega^2(x^2+y^2+z^2).
$$

Using the one-dimensional eigenfunctions above,

$$
\phi_{n_x,n_y,n_z}(x,y,z)=\phi_{n_x}(x)\phi_{n_y}(y)\phi_{n_z}(z),\qquad E(n_x,n_y,n_z)=\hbar\omega\left(n_x+n_y+n_z+\frac32\right).
$$

##### Spherical Dimensionless Schrödinger Equation

Using the conventions and differential operators defined in [Spherical Coordinates](../mathematical_physics/coordinate_spherical.md#differential-operators), the Hamiltonian becomes

$$
\hat H=-\frac{\hbar^2}{2M}\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)+\frac{\hat L^2}{2Mr^2}+\frac12M\omega^2r^2.
$$

Using the normalized [spherical harmonics](../mathematical_physics/specialfunction_spherical_harmonics.md), a normalized coordinate wave function admits the angular expansion

$$
\psi(r,\theta,\phi)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}Y_{lm}(\theta,\phi)\psi_{lm}(r),\qquad \sum_{l=0}^{\infty}\sum_{m=-l}^{l}\int_0^\infty r^2\,dr\,|\psi_{lm}(r)|^2=1.
$$

Define the characteristic length $r_0$, dimensionless radius $\rho$, and dimensionless energy $\epsilon$ by

$$
\boxed{r_0=\sqrt{\frac{\hbar}{M\omega}},\qquad \rho=\frac{r}{r_0},\qquad \epsilon=\frac{E}{\hbar\omega}}.
$$

For fixed $l$ and $m$, the physical and dimensionless radial basis functions satisfy

$$
\phi_{lm}(r)=\frac1{r_0^{3/2}}\widetilde{\phi}_{lm}(\rho),\qquad \int_0^\infty\rho^2\,d\rho\,|\widetilde{\phi}_{lm}(\rho)|^2=1.
$$

The dimensionless radial basis function satisfies

$$
\left[-\frac12\left(\frac{d^2}{d\rho^2}+\frac2\rho\frac{d}{d\rho}-\frac{l(l+1)}{\rho^2}\right)+\frac12\rho^2\right]\widetilde{\phi}_{lm}(\rho)=\epsilon\widetilde{\phi}_{lm}(\rho).
$$

##### Laguerre Solution and Energy Spectrum

Since the radial equation depends only on $l$, regularity at $\rho=0$ and normalizability require

$$
\widetilde{\phi}_l(\rho)=\rho^le^{-z/2}g(z),\qquad z=\rho^2.
$$

The reduced function satisfies

$$
z\frac{d^2}{dz^2}g(z)+\left(l+\frac32-z\right)\frac{d}{dz}g(z)+\frac{\epsilon-l-\frac32}{2}g(z)=0.
$$

The polynomial solution is the [generalized Laguerre polynomial](../mathematical_physics/specialfunction_generalized_laguerre_polynomials.md):

$$
z\frac{d^2}{dz^2}L_n^{(l+1/2)}(z)+\left(l+\frac32-z\right)\frac{d}{dz}L_n^{(l+1/2)}(z)+nL_n^{(l+1/2)}(z)=0.
$$

Polynomial termination therefore gives

$$
\epsilon(n,l)=2n+l+\frac32,\qquad n\in\mathbb N_0.
$$

The Laguerre normalization coefficient and dimensionless radial eigenfunctions are

$$
\boxed{N_n^{(l+1/2)}=\sqrt{\frac{2n!}{\Gamma(n+l+3/2)}},\qquad \widetilde{\phi}_{nl}(\rho)=N_n^{(l+1/2)}\rho^le^{-\rho^2/2}L_n^{(l+1/2)}(\rho^2)}.
$$

The physical eigenfunctions and energy spectrum are

$$
\boxed{\phi_{nl}(r)=\frac{N_n^{(l+1/2)}}{r_0^{3/2}}\left(\frac{r}{r_0}\right)^le^{-r^2/(2r_0^2)}L_n^{(l+1/2)}\left(\frac{r^2}{r_0^2}\right),\qquad E(n,l)=\hbar\omega\left(2n+l+\frac32\right)}.
$$

The full eigenfunctions are

$$
\phi_{nlm}(r,\theta,\phi)=\phi_{nl}(r)Y_{lm}(\theta,\phi).
$$

The orthogonality relation is

$$
\int_0^\infty r^2\,dr\int_0^\pi\sin\theta\,d\theta\int_0^{2\pi}d\phi\,\phi_{n'l'm'}^*(r,\theta,\phi)\phi_{nlm}(r,\theta,\phi)=\delta_{n'n}\delta_{l'l}\delta_{m'm}.
$$

The completeness relation is

$$
\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\sum_{n=0}^{\infty}\phi_{nlm}(r,\theta,\phi)\phi_{nlm}^*(r',\theta',\phi')=\frac{\delta(r-r')\delta(\theta-\theta')\delta(\phi-\phi')}{r^2\sin\theta}.
$$

A normalized coordinate wave function can therefore be expanded as

$$
\psi(r,\theta,\phi)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\sum_{n=0}^{\infty}c_{nlm}\phi_{nlm}(r,\theta,\phi),\qquad c_{nlm}=\int_0^\infty r^2\,dr\int_0^\pi\sin\theta\,d\theta\int_0^{2\pi}d\phi\,\phi_{nlm}^*(r,\theta,\phi)\psi(r,\theta,\phi),\qquad \sum_{l=0}^{\infty}\sum_{m=-l}^{l}\sum_{n=0}^{\infty}|c_{nlm}|^2=1.
$$

Under spatial inversion, $(r,\theta,\phi)\rightarrow(r,\pi-\theta,\phi+\pi)$, the parity is

$$
\boxed{\phi_{nlm}(r,\pi-\theta,\phi+\pi)=(-1)^l\phi_{nlm}(r,\theta,\phi)}.
$$
