# Spin–Orbit Coupling

### Spherical Preliminaries

##### Spinor Spherical Harmonics

Let $\Omega=(\theta,\phi)$ and let $Y_{l m_l}(\Omega)$ be a scalar spherical harmonic. The spin-$\frac12$ basis is

$$
\chi_{1/2}=\begin{pmatrix}1\\0\end{pmatrix},\qquad \chi_{-1/2}=\begin{pmatrix}0\\1\end{pmatrix},\qquad \hat S_z\chi_\mu=\hbar\mu\chi_\mu.
$$

With $\hat{\mathbf J}=\hat{\mathbf L}+\hat{\mathbf S}$,

$$
\boxed{\mathcal Y_{jm}^{l}(\Omega)=\sum_{m_l=-l}^{l}\sum_{\mu=\pm1/2}\left\langle l,m_l;\frac12,\mu\middle|j,m\right\rangle Y_{l m_l}(\Omega)\chi_\mu}.
$$

The coupled quantum numbers satisfy

$$
j=l\pm\frac12,\qquad m=m_l+\mu,\qquad m=-j,-j+1,\ldots,j,
$$

with only $j=\frac12$ allowed for $l=0$. The spinor spherical harmonics obey

$$
\hat J^2\mathcal Y_{jm}^{l}=\hbar^2j(j+1)\mathcal Y_{jm}^{l},\qquad \hat J_z\mathcal Y_{jm}^{l}=\hbar m\mathcal Y_{jm}^{l}.
$$

$$
\hat L^2\mathcal Y_{jm}^{l}=\hbar^2l(l+1)\mathcal Y_{jm}^{l},\qquad \hat S^2\mathcal Y_{jm}^{l}=\frac34\hbar^2\mathcal Y_{jm}^{l}.
$$

Definitions and further properties are collected in [Spherical Harmonics](../mathematical_physics/spherical_harmonics.md#spinor-spherical-harmonics).

##### Spherical Gradient

Let $\hat{\mathbf r}$ be the radial unit vector and let $\nabla_\Omega$ act only on the angular coordinates. The gradient separates as

$$
\nabla=\hat{\mathbf r}\frac{\partial}{\partial r}+\frac{1}{r}\nabla_\Omega.
$$

The orbital angular momentum is

$$
\hat{\mathbf L}=-i\hbar\,\hat{\mathbf r}\times\nabla_\Omega.
$$

Since $\hat{\mathbf r}\cdot\nabla_\Omega=0$,

$$
\hat{\mathbf r}\times\hat{\mathbf L}=i\hbar\nabla_\Omega,\qquad \nabla_\Omega=-\frac{i}{\hbar}\hat{\mathbf r}\times\hat{\mathbf L}.
$$

Therefore,

$$
\boxed{\nabla=\hat{\mathbf r}\frac{\partial}{\partial r}-\frac{i}{\hbar r}\hat{\mathbf r}\times\hat{\mathbf L}}.
$$

The spherical-coordinate conventions and differential operators are collected in [Spherical Coordinates](../mathematical_physics/spherical_coordinates.md#differential-operators).

### Spherical Spin–Momentum Operator

##### Factorization of $\boldsymbol\sigma\cdot\nabla$

The Pauli vector and its product identity are defined in [Spin-$1/2$ in Matrix Form](angular_momentum.md#spin-12-in-matrix-form). Since $\hat{\mathbf r}\cdot\hat{\mathbf L}=0$,

$$
(\boldsymbol\sigma\cdot\hat{\mathbf r})(\boldsymbol\sigma\cdot\hat{\mathbf L})=i\boldsymbol\sigma\cdot(\hat{\mathbf r}\times\hat{\mathbf L}).
$$

Using the spherical gradient,

$$
\boxed{\boldsymbol\sigma\cdot\nabla=(\boldsymbol\sigma\cdot\hat{\mathbf r})\frac{\partial}{\partial r}-\frac{i}{\hbar r}\boldsymbol\sigma\cdot(\hat{\mathbf r}\times\hat{\mathbf L})=(\boldsymbol\sigma\cdot\hat{\mathbf r})\left[\frac{\partial}{\partial r}-\frac{1}{\hbar r}\boldsymbol\sigma\cdot\hat{\mathbf L}\right]}.
$$

With $\hat{\mathbf p}=-i\hbar\nabla$,

$$
\boxed{\boldsymbol\sigma\cdot\hat{\mathbf p}=-i\hbar(\boldsymbol\sigma\cdot\hat{\mathbf r})\left[\frac{\partial}{\partial r}-\frac{1}{\hbar r}\boldsymbol\sigma\cdot\hat{\mathbf L}\right]}.
$$

##### Operator $\boldsymbol\sigma\cdot\hat{\mathbf L}$

For spin $\frac12$,

$$
\hat{\mathbf S}=\frac{\hbar}{2}\boldsymbol\sigma,\qquad \hat{\mathbf J}=\hat{\mathbf L}+\hat{\mathbf S}.
$$

The scalar coupling identity gives

$$
\boldsymbol\sigma\cdot\hat{\mathbf L}=\frac{2}{\hbar}\hat{\mathbf S}\cdot\hat{\mathbf L}=\frac{1}{\hbar}\left(\hat J^2-\hat L^2-\hat S^2\right).
$$

Hence,

$$
\boldsymbol\sigma\cdot\hat{\mathbf L}\,\mathcal Y_{jm}^{l}=\hbar\left[j(j+1)-l(l+1)-\frac34\right]\mathcal Y_{jm}^{l}.
$$

For the two coupled branches,

$$
\boldsymbol\sigma\cdot\hat{\mathbf L}\,\mathcal Y_{l+1/2,m}^{l}=\hbar l\,\mathcal Y_{l+1/2,m}^{l},\qquad \boldsymbol\sigma\cdot\hat{\mathbf L}\,\mathcal Y_{l-1/2,m}^{l}=-\hbar(l+1)\mathcal Y_{l-1/2,m}^{l}.
$$

Define the Dirac quantum number $\kappa$ by

$$
\kappa=-\left(j+\frac12\right)=-(l+1)\quad\text{for }j=l+\frac12,\qquad \kappa=j+\frac12=l\quad\text{for }j=l-\frac12.
$$

Writing $\mathcal Y_{\kappa m}\equiv\mathcal Y_{jm}^{l}$,

$$
\boxed{\boldsymbol\sigma\cdot\hat{\mathbf L}\,\mathcal Y_{\kappa m}=-\hbar(\kappa+1)\mathcal Y_{\kappa m}}.
$$

The general scalar-coupling identity is derived in [Angular Momentum](angular_momentum.md#angular-momentum-coupling).

##### Operator $\boldsymbol\sigma\cdot\hat{\mathbf r}$

The operator $\boldsymbol\sigma\cdot\hat{\mathbf r}$ preserves $j$ and $m$ while exchanging the two orbital branches:

$$
[\hat J^2,\boldsymbol\sigma\cdot\hat{\mathbf r}]=0,\qquad [\hat J_z,\boldsymbol\sigma\cdot\hat{\mathbf r}]=0,\qquad l'=2j-l.
$$

With the Condon–Shortley phases used above,

$$
\boxed{(\boldsymbol\sigma\cdot\hat{\mathbf r})\mathcal Y_{jm}^{l}=-\mathcal Y_{jm}^{l'},\qquad (\boldsymbol\sigma\cdot\hat{\mathbf r})\mathcal Y_{\kappa m}=-\mathcal Y_{-\kappa,m}}.
$$

##### Radial Action

Let $G(r)$ be a radial function. Since

$$
\left[\frac{\partial}{\partial r}-\frac{1}{\hbar r}\boldsymbol\sigma\cdot\hat{\mathbf L}\right]\frac{G(r)}{r}\mathcal Y_{\kappa m}=\frac{1}{r}\left[\frac{dG}{dr}+\frac{\kappa}{r}G(r)\right]\mathcal Y_{\kappa m},
$$

the angular factor is exchanged by $\boldsymbol\sigma\cdot\hat{\mathbf r}$:

$$
\boxed{\boldsymbol\sigma\cdot\nabla\left[\frac{G(r)}{r}\mathcal Y_{\kappa m}\right]=-\frac{1}{r}\left[\frac{dG}{dr}+\frac{\kappa}{r}G(r)\right]\mathcal Y_{-\kappa,m}}.
$$

Therefore,

$$
\boxed{\boldsymbol\sigma\cdot\hat{\mathbf p}\left[\frac{G(r)}{r}\mathcal Y_{\kappa m}\right]=\frac{i\hbar}{r}\left[\frac{dG}{dr}+\frac{\kappa}{r}G(r)\right]\mathcal Y_{-\kappa,m}}.
$$

### Cylindrical Preliminaries

##### Spinor Azimuthal Harmonics

Using the [azimuthal harmonics](../mathematical_physics/cylindrical_coordinates.md#azimuthal-harmonics),

$$
\Phi_m(\phi)=\frac{1}{\sqrt{2\pi}}e^{im\phi},\qquad m\in\mathbb Z.
$$

The spin basis and Pauli matrices are defined in [Spin-$1/2$ in Matrix Form](angular_momentum.md#spin-12-in-matrix-form). The total angular-momentum projection is

$$
\hat J_z=\hat L_z+\hat S_z=-i\hbar\frac{\partial}{\partial\phi}+\frac{\hbar}{2}\sigma_z.
$$

For $\Omega\in\mathbb Z+\frac12$, define

$$
\Xi_{\Omega,+}(\phi)=\Phi_{\Omega-\frac12}(\phi)\chi_+=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}e^{i(\Omega-\frac12)\phi}\\0\end{pmatrix},\qquad \Xi_{\Omega,-}(\phi)=\Phi_{\Omega+\frac12}(\phi)\chi_-=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}0\\e^{i(\Omega+\frac12)\phi}\end{pmatrix}.
$$

They satisfy

$$
\hat J_z\Xi_{\Omega,\pm}=\hbar\Omega\Xi_{\Omega,\pm},\qquad \int_0^{2\pi}d\phi\,\Xi_{\Omega,s}^\dagger\Xi_{\Omega',s'}=\delta_{\Omega\Omega'}\delta_{ss'},\qquad s,s'\in\{+,-\}.
$$

A spinor with fixed $\Omega$ is

$$
\Psi_\Omega(r_\perp,\phi,z)=f_+(r_\perp,z)\Xi_{\Omega,+}(\phi)+f_-(r_\perp,z)\Xi_{\Omega,-}(\phi)=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}f_+(r_\perp,z)e^{i(\Omega-\frac12)\phi}\\f_-(r_\perp,z)e^{i(\Omega+\frac12)\phi}\end{pmatrix}.
$$

##### Cylindrical Gradient

The [cylindrical gradient](../mathematical_physics/cylindrical_coordinates.md#differential-operators) can be written as

$$
\nabla=\hat{\mathbf r}_\perp\frac{\partial}{\partial r_\perp}+\hat{\mathbf z}\frac{\partial}{\partial z}+\frac{i}{\hbar r_\perp}\hat{\boldsymbol\phi}\hat L_z,\qquad \hat L_z=-i\hbar\frac{\partial}{\partial\phi}.
$$

The operators $\sigma_{r_\perp}$, $\sigma_\phi$, and $\sigma_z$ are defined in [Pauli Operators in Local Frames](angular_momentum.md#pauli-operators-in-local-frames). Hence

$$
\boxed{\boldsymbol{\sigma}\cdot\nabla=\sigma_{r_\perp}\frac{\partial}{\partial r_\perp}+\sigma_z\frac{\partial}{\partial z}+\frac{i}{\hbar r_\perp}\sigma_\phi\hat L_z}.
$$

### Cylindrical Spin–Momentum Operator

With $\hat{\mathbf p}=-i\hbar\nabla$,

$$
\boxed{\boldsymbol\sigma\cdot\hat{\mathbf p}=-i\hbar\left(\sigma_{r_\perp}\frac{\partial}{\partial r_\perp}+\sigma_z\frac{\partial}{\partial z}+\frac{i}{\hbar r_\perp}\sigma_\phi\hat L_z\right)}.
$$

##### Radial Operator $\sigma_{r_\perp}$

The radial Pauli operator exchanges the two angular spinors:

$$
\sigma_{r_\perp}\Xi_{\Omega,+}=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}0&e^{-i\phi}\\e^{i\phi}&0\end{pmatrix}\begin{pmatrix}e^{i(\Omega-\frac12)\phi}\\0\end{pmatrix}=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}0\\e^{i(\Omega+\frac12)\phi}\end{pmatrix}=\Xi_{\Omega,-}.
$$

$$
\sigma_{r_\perp}\Xi_{\Omega,-}=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}0&e^{-i\phi}\\e^{i\phi}&0\end{pmatrix}\begin{pmatrix}0\\e^{i(\Omega+\frac12)\phi}\end{pmatrix}=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}e^{i(\Omega-\frac12)\phi}\\0\end{pmatrix}=\Xi_{\Omega,+}.
$$

Therefore,

$$
\sigma_{r_\perp}\frac{\partial\Psi_\Omega}{\partial r_\perp}=\sigma_{r_\perp}\frac{\partial}{\partial r_\perp}\left(f_+\Xi_{\Omega,+}+f_-\Xi_{\Omega,-}\right)=\frac{\partial f_-}{\partial r_\perp}\Xi_{\Omega,+}+\frac{\partial f_+}{\partial r_\perp}\Xi_{\Omega,-}.
$$

##### Axial Operator $\sigma_z$

The axial Pauli operator acts diagonally on the two angular spinors:

$$
\sigma_z\Xi_{\Omega,+}=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}1&0\\0&-1\end{pmatrix}\begin{pmatrix}e^{i(\Omega-\frac12)\phi}\\0\end{pmatrix}=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}e^{i(\Omega-\frac12)\phi}\\0\end{pmatrix}=\Xi_{\Omega,+}.
$$

$$
\sigma_z\Xi_{\Omega,-}=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}1&0\\0&-1\end{pmatrix}\begin{pmatrix}0\\e^{i(\Omega+\frac12)\phi}\end{pmatrix}=-\frac{1}{\sqrt{2\pi}}\begin{pmatrix}0\\e^{i(\Omega+\frac12)\phi}\end{pmatrix}=-\Xi_{\Omega,-}.
$$

Therefore,

$$
\sigma_z\frac{\partial\Psi_\Omega}{\partial z}=\sigma_z\frac{\partial}{\partial z}\left(f_+\Xi_{\Omega,+}+f_-\Xi_{\Omega,-}\right)=\frac{\partial f_+}{\partial z}\Xi_{\Omega,+}-\frac{\partial f_-}{\partial z}\Xi_{\Omega,-}.
$$

##### Azimuthal Operator $\sigma_\phi\hat L_z$

The orbital angular momentum acts diagonally on the angular spinors:

$$
\hat L_z\Xi_{\Omega,+}=\hbar\left(\Omega-\frac12\right)\Xi_{\Omega,+},\qquad \hat L_z\Xi_{\Omega,-}=\hbar\left(\Omega+\frac12\right)\Xi_{\Omega,-}.
$$

The combined action on $\Xi_{\Omega,+}$ is

$$
\frac{i}{\hbar r_\perp}\sigma_\phi\hat L_z\Xi_{\Omega,+}=\frac{i(\Omega-\frac12)}{r_\perp\sqrt{2\pi}}\begin{pmatrix}0&-ie^{-i\phi}\\ie^{i\phi}&0\end{pmatrix}\begin{pmatrix}e^{i(\Omega-\frac12)\phi}\\0\end{pmatrix}=-\frac{\Omega-\frac12}{r_\perp\sqrt{2\pi}}\begin{pmatrix}0\\e^{i(\Omega+\frac12)\phi}\end{pmatrix}=-\frac{\Omega-\frac12}{r_\perp}\Xi_{\Omega,-}.
$$

The combined action on $\Xi_{\Omega,-}$ is

$$
\frac{i}{\hbar r_\perp}\sigma_\phi\hat L_z\Xi_{\Omega,-}=\frac{i(\Omega+\frac12)}{r_\perp\sqrt{2\pi}}\begin{pmatrix}0&-ie^{-i\phi}\\ie^{i\phi}&0\end{pmatrix}\begin{pmatrix}0\\e^{i(\Omega+\frac12)\phi}\end{pmatrix}=\frac{\Omega+\frac12}{r_\perp\sqrt{2\pi}}\begin{pmatrix}e^{i(\Omega-\frac12)\phi}\\0\end{pmatrix}=\frac{\Omega+\frac12}{r_\perp}\Xi_{\Omega,+}.
$$

Therefore,

$$
\frac{i}{\hbar r_\perp}\sigma_\phi\hat L_z\Psi_\Omega=\frac{i}{\hbar r_\perp}\sigma_\phi\hat L_z\left(f_+\Xi_{\Omega,+}+f_-\Xi_{\Omega,-}\right)=\frac{\Omega+\frac12}{r_\perp}f_-\Xi_{\Omega,+}-\frac{\Omega-\frac12}{r_\perp}f_+\Xi_{\Omega,-}.
$$

##### Combined Action

Combining the radial, axial, and azimuthal contributions,

$$
\boldsymbol\sigma\cdot\nabla\,\Psi_\Omega=\left(\frac{\partial f_-}{\partial r_\perp}\Xi_{\Omega,+}+\frac{\partial f_+}{\partial r_\perp}\Xi_{\Omega,-}\right)+\left(\frac{\partial f_+}{\partial z}\Xi_{\Omega,+}-\frac{\partial f_-}{\partial z}\Xi_{\Omega,-}\right)+\left(\frac{\Omega+\frac12}{r_\perp}f_-\Xi_{\Omega,+}-\frac{\Omega-\frac12}{r_\perp}f_+\Xi_{\Omega,-}\right).
$$

Collecting the coefficients of $\Xi_{\Omega,+}$ and $\Xi_{\Omega,-}$ gives

$$
\boxed{\boldsymbol\sigma\cdot\nabla\,\Psi_\Omega=\left[\frac{\partial f_+}{\partial z}+\left(\frac{\partial}{\partial r_\perp}+\frac{\Omega+\frac12}{r_\perp}\right)f_-\right]\Xi_{\Omega,+}+\left[\left(\frac{\partial}{\partial r_\perp}-\frac{\Omega-\frac12}{r_\perp}\right)f_+-\frac{\partial f_-}{\partial z}\right]\Xi_{\Omega,-}}.
$$

Therefore,

$$
\boxed{\boldsymbol\sigma\cdot\hat{\mathbf p}\,\Psi_\Omega=-i\hbar\left\{\left[\frac{\partial f_+}{\partial z}+\left(\frac{\partial}{\partial r_\perp}+\frac{\Omega+\frac12}{r_\perp}\right)f_-\right]\Xi_{\Omega,+}+\left[\left(\frac{\partial}{\partial r_\perp}-\frac{\Omega-\frac12}{r_\perp}\right)f_+-\frac{\partial f_-}{\partial z}\right]\Xi_{\Omega,-}\right\}}.
$$
