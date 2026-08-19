# Harmonic-Oscillator Coordinate-Space Solutions

### One-Dimensional Harmonic Oscillator

The operator construction and number representation are given in [Quantum Harmonic Oscillator](harmonic_oscillator.md). Here the eigenvalue problem is solved directly in coordinate space.

##### Dimensionless Schrödinger Equation

For a one-dimensional oscillator with mass $M>0$ and angular frequency $\omega>0$, the coordinate-space Hamiltonian is

$$
\hat H=-\frac{\hbar^2}{2M}\frac{d^2}{dx^2}+\frac12M\omega^2x^2.
$$

Define the characteristic length $x_0$, dimensionless coordinate $\xi$, and dimensionless energy $\epsilon$ by

$$
\boxed{x_0=\sqrt{\frac{\hbar}{M\omega}},\qquad \xi=\frac{x}{x_0},\qquad \epsilon=\frac{E}{\hbar\omega}}.
$$

Write the normalized coordinate wave function as

$$
\psi(x)=\frac{1}{\sqrt{x_0}}\widetilde{\psi}(\xi),\qquad \int_{-\infty}^{\infty}|\widetilde{\psi}(\xi)|^2\,d\xi=1.
$$

The stationary Schrödinger equation becomes

$$
\left[-\frac12\frac{d^2}{d\xi^2}+\frac12\xi^2\right]\widetilde{\psi}(\xi)=\epsilon\widetilde{\psi}(\xi).
$$

##### Hermite Solution and Energy Spectrum

Normalizability requires Gaussian decay. Write

$$
\widetilde{\psi}(\xi)=e^{-\xi^2/2}h(\xi).
$$

The reduced function satisfies

$$
h''(\xi)-2\xi h'(\xi)+(2\epsilon-1)h(\xi)=0.
$$

The polynomial solution is the physicists' [Hermite polynomial](../mathematical_physics/hermite_polynomials.md):

$$
H_n''(\xi)-2\xi H_n'(\xi)+2nH_n(\xi)=0,\qquad n=0,1,2,\ldots.
$$

Polynomial termination therefore gives

$$
\epsilon_n=n+\frac12,\qquad n=0,1,2,\ldots.
$$

The Hermite normalization coefficient and dimensionless eigenfunctions are

$$
\boxed{N_n=\frac{1}{\sqrt{2^n n!\sqrt\pi}},\qquad \widetilde{\psi}_n(\xi)=N_nH_n(\xi)e^{-\xi^2/2}}.
$$

The physical eigenfunctions and energy spectrum are

$$
\boxed{\psi_n(x)=\frac{N_n}{\sqrt{x_0}}H_n\left(\frac{x}{x_0}\right)e^{-x^2/(2x_0^2)},\qquad E_n=\hbar\omega\left(n+\frac12\right)}.
$$

They satisfy

$$
\int_{-\infty}^{\infty}\psi_m^*(x)\psi_n(x)\,dx=\delta_{mn},\qquad \sum_{n=0}^{\infty}\psi_n(x)\psi_n^*(x')=\delta(x-x').
$$

Since $H_n(-\xi)=(-1)^nH_n(\xi)$,

$$
\boxed{\psi_n(-x)=(-1)^n\psi_n(x)}.
$$

Each one-dimensional energy level is nondegenerate.

### Two-Dimensional Isotropic Harmonic Oscillator

For a two-dimensional oscillator with mass $M>0$ and angular frequency $\omega>0$, the Cartesian Hamiltonian is

$$
\hat H=-\frac{\hbar^2}{2M}\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right)+\frac12M\omega^2(x^2+y^2).
$$

Using the one-dimensional eigenfunctions above,

$$
\psi_{n_x,n_y}(x,y)=\psi_{n_x}(x)\psi_{n_y}(y),\qquad E_{n_x,n_y}=\hbar\omega\left(n_x+n_y+1\right).
$$

##### Polar Dimensionless Schrödinger Equation

The same system can be expressed in polar coordinates, $x=r\cos\phi$ and $y=r\sin\phi$. The Hamiltonian then becomes

$$
\hat H=-\frac{\hbar^2}{2M}\left(\frac{\partial^2}{\partial r^2}+\frac{1}{r}\frac{\partial}{\partial r}+\frac{1}{r^2}\frac{\partial^2}{\partial\phi^2}\right)+\frac12M\omega^2r^2.
$$

Define the characteristic length $r_0$, dimensionless radius $\rho$, and dimensionless energy $\epsilon$ by

$$
\boxed{r_0=\sqrt{\frac{\hbar}{M\omega}},\qquad \rho=\frac{r}{r_0},\qquad \epsilon=\frac{E}{\hbar\omega}}.
$$

Write the normalized coordinate wave function as

$$
\psi(r,\phi)=\frac{1}{r_0}\widetilde{\psi}(\rho,\phi),\qquad \int_0^{2\pi}d\phi\int_0^\infty \rho\,d\rho\,|\widetilde{\psi}(\rho,\phi)|^2=1.
$$

The stationary Schrödinger equation becomes

$$
\left[-\frac12\left(\frac{\partial^2}{\partial\rho^2}+\frac{1}{\rho}\frac{\partial}{\partial\rho}+\frac{1}{\rho^2}\frac{\partial^2}{\partial\phi^2}\right)+\frac12\rho^2\right]\widetilde{\psi}(\rho,\phi)=\epsilon\widetilde{\psi}(\rho,\phi).
$$

##### Laguerre Solution and Energy Spectrum

For $m\in\mathbb Z$, separate the angular dependence as

$$
\widetilde{\psi}(\rho,\phi)=\frac{1}{\sqrt{2\pi}}R_m(\rho)e^{im\phi},\qquad \int_0^\infty \rho|R_m(\rho)|^2\,d\rho=1.
$$

The radial function satisfies

$$
\left[-\frac12\left(\frac{d^2}{d\rho^2}+\frac{1}{\rho}\frac{d}{d\rho}-\frac{m^2}{\rho^2}\right)+\frac12\rho^2\right]R_m(\rho)=\epsilon R_m(\rho).
$$

Regularity at $\rho=0$ and normalizability require

$$
R_m(\rho)=\rho^{|m|}e^{-\rho^2/2}l(z),\qquad z=\rho^2.
$$

The reduced function satisfies

$$
z\frac{d^2}{dz^2}l(z)+(|m|+1-z)\frac{d}{dz}l(z)+\frac{\epsilon-|m|-1}{2}l(z)=0.
$$

The polynomial solution is the [generalized Laguerre polynomial](../mathematical_physics/generalized_laguerre_polynomials.md):

$$
z\frac{d^2}{dz^2}L_{n_r}^{(|m|)}(z)+(|m|+1-z)\frac{d}{dz}L_{n_r}^{(|m|)}(z)+n_rL_{n_r}^{(|m|)}(z)=0.
$$

Polynomial termination therefore gives

$$
\epsilon_{n_rm}=2n_r+|m|+1,\qquad n_r\in\mathbb N_0.
$$

The Laguerre normalization coefficient and dimensionless eigenfunctions are

$$
\boxed{N_{n_r}^{(|m|)}=\sqrt{\frac{n_r!}{(n_r+|m|)!}},\qquad \widetilde{\psi}_{n_rm}(\rho,\phi)=\frac{N_{n_r}^{(|m|)}}{\sqrt\pi}\rho^{|m|}e^{-\rho^2/2}L_{n_r}^{(|m|)}(\rho^2)e^{im\phi}}.
$$

The physical eigenfunctions and energy spectrum are

$$
\boxed{\psi_{n_rm}(r,\phi)=\frac{N_{n_r}^{(|m|)}}{r_0\sqrt\pi}\left(\frac{r}{r_0}\right)^{|m|}e^{-r^2/(2r_0^2)}L_{n_r}^{(|m|)}\left(\frac{r^2}{r_0^2}\right)e^{im\phi},\qquad E_{n_rm}=\hbar\omega(2n_r+|m|+1)}.
$$

The orthogonality relation is

$$
\int_0^{2\pi}d\phi\int_0^\infty r\,dr\,\psi_{n_rm}^*(r,\phi)\psi_{n_r'm'}(r,\phi)=\delta_{n_rn_r'}\delta_{mm'}.
$$

The completeness relation is

$$
\sum_{m=-\infty}^{\infty}\sum_{n_r=0}^{\infty}\psi_{n_rm}(r,\phi)\psi_{n_rm}^*(r',\phi')=\frac{\delta(r-r')\delta(\phi-\phi')}{r}.
$$

Under spatial inversion, $(r,\phi)\rightarrow(r,\phi+\pi)$, the parity is

$$
\boxed{\psi_{n_rm}(r,\phi+\pi)=(-1)^m\psi_{n_rm}(r,\phi)}.
$$

### Three-Dimensional Isotropic Harmonic Oscillator

For a three-dimensional oscillator with mass $M>0$ and angular frequency $\omega>0$, the Cartesian Hamiltonian is

$$
\hat H=-\frac{\hbar^2}{2M}\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}+\frac{\partial^2}{\partial z^2}\right)+\frac12M\omega^2(x^2+y^2+z^2).
$$

Using the one-dimensional eigenfunctions above,

$$
\psi_{n_xn_yn_z}(x,y,z)=\psi_{n_x}(x)\psi_{n_y}(y)\psi_{n_z}(z),\qquad E_{n_xn_yn_z}=\hbar\omega\left(n_x+n_y+n_z+\frac32\right).
$$

##### Spherical Dimensionless Schrödinger Equation

Using the conventions and differential operators defined in [Spherical Coordinates](../mathematical_physics/spherical_coordinates.md#differential-operators), the Hamiltonian becomes

$$
\hat H=-\frac{\hbar^2}{2M}\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)+\frac{\hat L^2}{2Mr^2}+\frac12M\omega^2r^2.
$$

Define the characteristic length $r_0$, dimensionless radius $\rho$, and dimensionless energy $\epsilon$ by

$$
\boxed{r_0=\sqrt{\frac{\hbar}{M\omega}},\qquad \rho=\frac{r}{r_0},\qquad \epsilon=\frac{E}{\hbar\omega}}.
$$

Write the normalized coordinate wave function as

$$
\psi(r,\theta,\phi)=\frac{1}{r_0^{3/2}}\widetilde{\psi}(\rho,\theta,\phi),\qquad \int_0^\infty\rho^2\,d\rho\int_0^\pi\sin\theta\,d\theta\int_0^{2\pi}d\phi\,|\widetilde{\psi}(\rho,\theta,\phi)|^2=1.
$$

For $l=0,1,2,\ldots$ and $m=-l,\ldots,l$, separate the angular dependence as

$$
\widetilde{\psi}(\rho,\theta,\phi)=R_l(\rho)Y_l^m(\theta,\phi),\qquad \hat L^2Y_l^m=\hbar^2l(l+1)Y_l^m.
$$

The radial function satisfies

$$
\left[-\frac12\left(\frac{d^2}{d\rho^2}+\frac{2}{\rho}\frac{d}{d\rho}-\frac{l(l+1)}{\rho^2}\right)+\frac12\rho^2\right]R_l(\rho)=\epsilon R_l(\rho).
$$

Define the reduced radial function by

$$
R_l(\rho)=\frac{u_l(\rho)}{\rho},\qquad \int_0^\infty|u_l(\rho)|^2\,d\rho=1.
$$

It satisfies

$$
\left[-\frac12\frac{d^2}{d\rho^2}+\frac{l(l+1)}{2\rho^2}+\frac12\rho^2\right]u_l(\rho)=\epsilon u_l(\rho).
$$

##### Laguerre Solution and Energy Spectrum

Regularity at $\rho=0$ and normalizability require

$$
u_l(\rho)=\rho^{l+1}e^{-\rho^2/2}g(z),\qquad z=\rho^2.
$$

The reduced function satisfies

$$
z\frac{d^2}{dz^2}g(z)+\left(l+\frac32-z\right)\frac{d}{dz}g(z)+\frac{\epsilon-l-\frac32}{2}g(z)=0.
$$

The polynomial solution is the [generalized Laguerre polynomial](../mathematical_physics/generalized_laguerre_polynomials.md):

$$
z\frac{d^2}{dz^2}L_n^{(l+1/2)}(z)+\left(l+\frac32-z\right)\frac{d}{dz}L_n^{(l+1/2)}(z)+nL_n^{(l+1/2)}(z)=0.
$$

Polynomial termination therefore gives

$$
\epsilon_{nl}=2n+l+\frac32,\qquad n\in\mathbb N_0.
$$

The Laguerre normalization coefficient and dimensionless eigenfunctions are

$$
\boxed{N_n^{(l+1/2)}=\sqrt{\frac{2n!}{\Gamma(n+l+3/2)}},\qquad \widetilde{\psi}_{nlm}(\rho,\theta,\phi)=N_n^{(l+1/2)}\rho^le^{-\rho^2/2}L_n^{(l+1/2)}(\rho^2)Y_l^m(\theta,\phi)}.
$$

The physical eigenfunctions and energy spectrum are

$$
\boxed{\psi_{nlm}(r,\theta,\phi)=\frac{N_n^{(l+1/2)}}{r_0^{3/2}}\left(\frac{r}{r_0}\right)^le^{-r^2/(2r_0^2)}L_n^{(l+1/2)}\left(\frac{r^2}{r_0^2}\right)Y_l^m(\theta,\phi),\qquad E_{nl}=\hbar\omega\left(2n+l+\frac32\right)}.
$$

The orthogonality relation is

$$
\int_0^\infty r^2\,dr\int_0^\pi\sin\theta\,d\theta\int_0^{2\pi}d\phi\,\psi_{nlm}^*(r,\theta,\phi)\psi_{n'l'm'}(r,\theta,\phi)=\delta_{nn'}\delta_{ll'}\delta_{mm'}.
$$

The completeness relation is

$$
\sum_{n=0}^{\infty}\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\psi_{nlm}(r,\theta,\phi)\psi_{nlm}^*(r',\theta',\phi')=\frac{\delta(r-r')\delta(\theta-\theta')\delta(\phi-\phi')}{r^2\sin\theta}.
$$

Under spatial inversion, $(r,\theta,\phi)\rightarrow(r,\pi-\theta,\phi+\pi)$, the parity is

$$
\boxed{\psi_{nlm}(r,\pi-\theta,\phi+\pi)=(-1)^l\psi_{nlm}(r,\theta,\phi)}.
$$
