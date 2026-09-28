# Electromagnetic Wave Propagation

### Wave Equations from Maxwell's Equations

Start from the [Maxwell equations in matter](electromagnetic_maxwell_equations.md#maxwell-equations-in-matter), where $\rho_{\mathrm f}$ and $\mathbf j_{\mathrm f}$ are the free charge and current densities:

$$
\nabla\cdot\mathbf D=\rho_{\mathrm f},\qquad \nabla\cdot\mathbf B=0.
$$

$$
\nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t},\qquad \nabla\times\mathbf H=\mathbf j_{\mathrm f}+\frac{\partial\mathbf D}{\partial t}.
$$

For a homogeneous, isotropic, linear, nondispersive medium, let $\epsilon$ and $\mu$ be constant scalars:

$$
\mathbf D=\epsilon\mathbf E,\qquad \mathbf B=\mu\mathbf H.
$$

Taking another curl of Faraday's and the Ampère–Maxwell equations gives

$$
\nabla\times(\nabla\times\mathbf E)=-\mu\frac{\partial\mathbf j_{\mathrm f}}{\partial t}-\mu\epsilon\frac{\partial^2\mathbf E}{\partial t^2},\qquad \nabla\times(\nabla\times\mathbf H)=\nabla\times\mathbf j_{\mathrm f}-\mu\epsilon\frac{\partial^2\mathbf H}{\partial t^2}.
$$

The Gauss laws become

$$
\nabla\cdot\mathbf E=\frac{\rho_{\mathrm f}}{\epsilon},\qquad \nabla\cdot\mathbf H=0.
$$

Using the vector identity

$$
\nabla\times(\nabla\times\mathbf F)=\nabla(\nabla\cdot\mathbf F)-\nabla^2\mathbf F,
$$

the wave equations with sources follow:

$$
\boxed{\left(\nabla^2-\mu\epsilon\frac{\partial^2}{\partial t^2}\right)\mathbf E=\frac{1}{\epsilon}\nabla\rho_{\mathrm f}+\mu\frac{\partial\mathbf j_{\mathrm f}}{\partial t},\qquad \left(\nabla^2-\mu\epsilon\frac{\partial^2}{\partial t^2}\right)\mathbf H=-\nabla\times\mathbf j_{\mathrm f}}.
$$

For a monochromatic field with time dependence $e^{-i\omega t}$, the electric-field equation also takes the frequency-domain curl–curl form

$$
\nabla\times(\nabla\times\mathbf E(\omega))-\mu\epsilon\omega^2\mathbf E(\omega)=i\omega\mu\mathbf j_{\mathrm f}(\omega).
$$

### Plane Waves in Source-Free Media

Consider a homogeneous, isotropic, linear medium with real constant $\epsilon,\mu>0$ and no free charge or current:

$$
\rho_{\mathrm f}=0,\qquad \mathbf j_{\mathrm f}=0,\qquad \mathbf D=\epsilon\mathbf E,\qquad \mathbf B=\mu\mathbf H.
$$

The wave equations derived above reduce to

$$
\boxed{\left(\nabla^2-\mu\epsilon\frac{\partial^2}{\partial t^2}\right)\mathbf E=0,\qquad \left(\nabla^2-\mu\epsilon\frac{\partial^2}{\partial t^2}\right)\mathbf H=0}.
$$

Let $\mathbf E_0,\mathbf H_0$ be complex amplitudes, $\mathbf k=k\hat{\mathbf k}$ the wave vector, and $\omega$ the angular frequency. A plane-wave solution is

$$
\mathbf E=\Re\!\left[\mathbf E_0e^{i(\mathbf k\cdot\mathbf r-\omega t)}\right],\qquad \mathbf H=\Re\!\left[\mathbf H_0e^{i(\mathbf k\cdot\mathbf r-\omega t)}\right].
$$

The Maxwell equations require

$$
\mathbf k\cdot\mathbf E_0=\mathbf k\cdot\mathbf H_0=0,\qquad \mathbf k\times\mathbf E_0=\omega\mu\mathbf H_0,\qquad \mathbf k\times\mathbf H_0=-\omega\epsilon\mathbf E_0.
$$

Define the medium impedance $Z=\sqrt{\mu/\epsilon}$. Then

$$
\boxed{k^2=\mu\epsilon\omega^2,\qquad v_{\mathrm p}=\frac{\omega}{k}=\frac{1}{\sqrt{\mu\epsilon}},\qquad \mathbf H_0=\frac{1}{Z}\hat{\mathbf k}\times\mathbf E_0}.
$$

With $c=(\mu_0\epsilon_0)^{-1/2}$, $\epsilon_r=\epsilon/\epsilon_0$, and $\mu_r=\mu/\mu_0$, the refractive index is $n=\sqrt{\epsilon_r\mu_r}=c/v_{\mathrm p}$.

Using the [energy density and Poynting vector in matter](electromagnetic_conservation_laws.md#energy-and-momentum-conservation-in-matter), their cycle averages are

$$
\langle u\rangle=\frac{\epsilon|\mathbf E_0|^2}{2},\qquad \boxed{\langle\mathbf S_P\rangle=\frac{|\mathbf E_0|^2}{2Z}\hat{\mathbf k}=v_{\mathrm p}\langle u\rangle\hat{\mathbf k}}.
$$

### Polarization

Take the propagation direction along $+z$. Let $E_{x0},E_{y0}\geq0$ be the amplitudes of the two transverse components and $\phi_x,\phi_y$ their constant phases. The complex electric-field amplitude is

$$
\mathbf k=k\hat{\mathbf z},\qquad \mathbf E_0=E_{x0}e^{i\phi_x}\hat{\mathbf x}+E_{y0}e^{i\phi_y}\hat{\mathbf y}.
$$

Define the phase of the $y$ component relative to the $x$ component and factor out the common phase:

$$
\Delta\phi=\phi_y-\phi_x,\qquad \mathbf E_0=e^{i\phi_x}\left(E_{x0}\hat{\mathbf x}+E_{y0}e^{i\Delta\phi}\hat{\mathbf y}\right).
$$

The physical electric field is

$$
\mathbf E(x,y,z,t)=\Re\!\left[\mathbf E_0e^{i(kz-\omega t)}\right].
$$

Only the two amplitudes and their phase difference determine the polarization:

- **Linear polarization:** One amplitude vanishes, or $\Delta\phi=0\pmod\pi$.
- **Elliptical polarization:** Both amplitudes are nonzero and $\Delta\phi\ne0\pmod\pi$.
- **Circular polarization:** The special elliptical case $E_{x0}=E_{y0}$ and $\Delta\phi=\pm\pi/2$.

With the $e^{i(kz-\omega t)}$ convention from the preceding section, define the circular polarization basis vectors

$$
\boxed{\mathbf e_{\mathrm R}=\frac{\hat{\mathbf x}-i\hat{\mathbf y}}{\sqrt2},\qquad \mathbf e_{\mathrm L}=\frac{\hat{\mathbf x}+i\hat{\mathbf y}}{\sqrt2}}.
$$

Here $\Delta\phi=-\pi/2$ gives right-handed polarization $\mathbf e_{\mathrm R}$, while $\Delta\phi=+\pi/2$ gives left-handed polarization $\mathbf e_{\mathrm L}$. Viewed from the $+z$ side toward the origin, the physical electric field rotates clockwise and counterclockwise, respectively.

### Plane Waves in a Free-Electron Medium

Consider a homogeneous, nonmagnetic free-electron medium with vacuum background. For a transverse wave in its bulk,

$$
\mathbf D=\epsilon_0\mathbf E,\qquad \mathbf B=\mu_0\mathbf H,\qquad \rho_{\mathrm f}=0.
$$

For a field of angular frequency $\omega$, take the time dependence $e^{-i\omega t}$. With electron number density $n_e$, charge magnitude $e>0$, mean drift velocity $\mathbf v$, and complex conductivity $\sigma(\omega)$,

$$
\mathbf j_{\mathrm f}(\omega)=\sigma(\omega)\mathbf E(\omega)=-n_e e\mathbf v(\omega).
$$

Let $m$ be the electron mass and $\tau$ the momentum relaxation time. With no applied static magnetic field, the wave magnetic force is second order in the field amplitude:

$$
m\frac{\mathrm d\mathbf v}{\mathrm dt}=-e(\mathbf E+\mathbf v\times\mathbf B)-\frac{m\mathbf v}{\tau}\simeq-e\mathbf E-\frac{m\mathbf v}{\tau}.
$$

Hence

$$
\boxed{\sigma(\omega)=\frac{n_e e^2\tau}{m(1-i\omega\tau)}}.
$$

The current can be absorbed into an effective permittivity:

$$
\nabla\times\mathbf B(\omega)=\mu_0\mathbf j_{\mathrm f}(\omega)-i\omega\mu_0\epsilon_0\mathbf E(\omega)=-i\omega\mu_0\epsilon_{\mathrm{eff}}(\omega)\mathbf E(\omega).
$$

Defining the plasma frequency $\omega_p$ gives

$$
\boxed{\epsilon_{\mathrm{eff}}(\omega)=\epsilon_0+\frac{i\sigma(\omega)}{\omega}=\epsilon_0\left[1-\frac{\omega_p^2}{\omega(\omega+i/\tau)}\right],\qquad \omega_p^2=\frac{n_e e^2}{m\epsilon_0}}.
$$

For a wave entering $z>0$, let $\mathbf E_0$ be its amplitude at $z=0$ and $k=\beta+i\alpha$, where $\beta$ is the phase constant and $\alpha\geq0$ the attenuation rate:

$$
\mathbf k=k\hat{\mathbf z},\qquad \mathbf E(\mathbf r,t)=\Re\!\left[\mathbf E_0e^{i(\mathbf k\cdot\mathbf r-\omega t)}\right]=\Re\!\left[\mathbf E_0e^{-\alpha z}e^{i(\beta z-\omega t)}\right].
$$

The frequency-domain wave equation from the first section gives

$$
\boxed{k^2=\mu_0\epsilon_{\mathrm{eff}}(\omega)\omega^2}.
$$

- **Low-frequency good conductor:** For $\omega\tau\ll1$,

  $$
  \sigma(\omega)\simeq\sigma_{\mathrm{dc}}=\frac{n_e e^2\tau}{m}.
  $$

  For a good conductor with $\sigma_{\mathrm{dc}}\gg\omega\epsilon_0$,

  $$
  \epsilon_{\mathrm{eff}}(\omega)\simeq\frac{i\sigma_{\mathrm{dc}}}{\omega},\qquad k\simeq\frac{1+i}{\delta},\qquad \delta=\sqrt{\frac{2}{\mu_0\sigma_{\mathrm{dc}}\omega}}.
  $$

- **High frequency with negligible scattering loss:** If $\omega\tau\gg1$,

  $$
  \frac{\epsilon_{\mathrm{eff}}}{\epsilon_0}\simeq1-\frac{\omega_p^2}{\omega^2},\qquad k^2\simeq\frac{\omega^2-\omega_p^2}{c^2}.
  $$

  For $\omega<\omega_p$, the wave is evanescent, with $k=i/\delta$ and $\delta=c/\sqrt{\omega_p^2-\omega^2}$; for $\omega>\omega_p$, it can propagate.

### Plane Waves in a Magnetized Free-Electron Medium

Consider a homogeneous free-electron medium with vacuum background. For a transverse wave in its bulk,

$$
\mathbf D=\epsilon_0\mathbf E,\qquad \mathbf B=\mu_0\mathbf H,\qquad \rho_{\mathrm f}=0.
$$

For a field of angular frequency $\omega$ with time dependence $e^{-i\omega t}$, the current response is tensorial:

$$
\mathbf j_{\mathrm f}(\omega)=\hat{\boldsymbol{\sigma}}(\omega)\cdot\mathbf E(\omega).
$$

Apply a static magnetic field along $+z$, with $B_{\mathrm{ext}}>0$:

$$
\mathbf B_{\mathrm{ext}}=B_{\mathrm{ext}}\hat{\mathbf z}.
$$

In collisionless linear response, neglect the force from the wave's magnetic field. With electron number density $n_e$, charge magnitude $e>0$, mass $m$, and drift velocity $\mathbf v$,

$$
m\frac{\mathrm d\mathbf v}{\mathrm dt}=-e(\mathbf E+\mathbf v\times\mathbf B_{\mathrm{ext}}),\qquad \mathbf j_{\mathrm f}=-n_e e\mathbf v.
$$

The conductivity tensor and cyclotron frequency are

$$
\boxed{\hat{\boldsymbol{\sigma}}(\omega)=\frac{\epsilon_0\omega_p^2}{\omega^2-\omega_B^2}\begin{pmatrix}i\omega&\omega_B&0\\-\omega_B&i\omega&0\\0&0&i(\omega^2-\omega_B^2)/\omega\end{pmatrix},\qquad \omega_B=\frac{eB_{\mathrm{ext}}}{m}}.
$$

With $\hat{\mathbf I}$ the identity tensor, the effective relative permittivity is

$$
\boxed{\hat{\boldsymbol{\epsilon}}_r(\omega)=\hat{\mathbf I}+\frac{i\hat{\boldsymbol{\sigma}}(\omega)}{\epsilon_0\omega}=\begin{pmatrix}\epsilon_1&i\epsilon_2&0\\-i\epsilon_2&\epsilon_1&0\\0&0&\epsilon_3\end{pmatrix}}.
$$

$$
\epsilon_1=1-\frac{\omega_p^2}{\omega^2-\omega_B^2},\qquad \epsilon_2=\frac{\omega_p^2\omega_B}{\omega(\omega^2-\omega_B^2)},\qquad \epsilon_3=1-\frac{\omega_p^2}{\omega^2}.
$$

Consider a plane wave propagating parallel to the static field:

$$
\mathbf k=k\hat{\mathbf z},\qquad \mathbf E(\mathbf r,t)=\Re\!\left[\mathbf E_0e^{i(kz-\omega t)}\right].
$$

Since $\rho_{\mathrm f}=0$, Gauss's law gives

$$
i\mathbf k\cdot\mathbf E_0=\frac{\rho_{\mathrm f}(\omega)}{\epsilon_0}=0\qquad\Longrightarrow\qquad E_{0z}=0.
$$

Substitution into the frequency-domain curl–curl equation gives

$$
\mathbf k\times(\mathbf k\times\mathbf E_0)+\frac{\omega^2}{c^2}\hat{\boldsymbol{\epsilon}}_r\cdot\mathbf E_0=0.
$$

Because $\mathbf k\cdot\mathbf E_0=0$,

$$
\mathbf k\times(\mathbf k\times\mathbf E_0)=\mathbf k(\mathbf k\cdot\mathbf E_0)-k^2\mathbf E_0=-k^2\mathbf E_0,
$$

and therefore

$$
k^2\mathbf E_0=\frac{\omega^2}{c^2}\hat{\boldsymbol{\epsilon}}_r\cdot\mathbf E_0.
$$

In the circular basis defined earlier, the eigenmodes are

$$
\boxed{\mathbf E_{0,\mathrm R}\parallel\mathbf e_{\mathrm R},\quad k_{\mathrm R}=\frac{\omega}{c}\sqrt{\epsilon_1+\epsilon_2};\qquad \mathbf E_{0,\mathrm L}\parallel\mathbf e_{\mathrm L},\quad k_{\mathrm L}=\frac{\omega}{c}\sqrt{\epsilon_1-\epsilon_2}.}
$$

For an initially $x$-polarized wave of amplitude $E_0$, assume equal excitation of the two modes. Define $\bar k=(k_{\mathrm R}+k_{\mathrm L})/2$ and $\Delta k=k_{\mathrm R}-k_{\mathrm L}$. Its complex amplitude is

$$
\mathbf E(\omega,z)=\frac{E_0}{\sqrt2}\left(\mathbf e_{\mathrm R}e^{ik_{\mathrm R}z}+\mathbf e_{\mathrm L}e^{ik_{\mathrm L}z}\right)=E_0e^{i\bar kz}\left[\cos\!\left(\frac{\Delta kz}{2}\right)\hat{\mathbf x}+\sin\!\left(\frac{\Delta kz}{2}\right)\hat{\mathbf y}\right].
$$

For weak gyrotropy, $|\epsilon_2|\ll\epsilon_1$, the polarization rotates after traversing a thickness $d$ by

$$
\boxed{\theta_{\mathrm F}=\frac{\Delta k\,d}{2}\simeq\frac{\omega d\,\epsilon_2}{2c\sqrt{\epsilon_1}}\propto B_{\mathrm{ext}}d.}
$$

### Plane Waves at a Material Interface

##### Snell's Law

Let the interface be $z=0$, with $+\hat{\mathbf z}$ pointing upward. Lossless, isotropic, nondispersive media 1 and 2 occupy $z>0$ (above) and $z<0$ (below), respectively.

$$
n_1=c\sqrt{\epsilon_1\mu_1},\qquad Z_1=\sqrt{\frac{\mu_1}{\epsilon_1}};\qquad n_2=c\sqrt{\epsilon_2\mu_2},\qquad Z_2=\sqrt{\frac{\mu_2}{\epsilon_2}}.
$$

At a specified angular frequency $\omega$, the wave numbers are

$$
k_1=\frac{n_1\omega}{c},\qquad k_2=\frac{n_2\omega}{c}.
$$

A wave is incident from medium 1 in the $xz$ plane. The angles $\theta_{\mathrm i},\theta_{\mathrm t}$ are measured from $-\hat{\mathbf z}$, and $\theta_{\mathrm r}$ from $+\hat{\mathbf z}$:

$$
\mathbf k_{\mathrm i}=k_1(\sin\theta_{\mathrm i}\hat{\mathbf x}-\cos\theta_{\mathrm i}\hat{\mathbf z}),\qquad \mathbf k_{\mathrm r}=k_1(\sin\theta_{\mathrm r}\hat{\mathbf x}+\cos\theta_{\mathrm r}\hat{\mathbf z}),\qquad \mathbf k_{\mathrm t}=k_2(\sin\theta_{\mathrm t}\hat{\mathbf x}-\cos\theta_{\mathrm t}\hat{\mathbf z}).
$$

For a stationary planar interface, time-translation symmetry fixes the frequency, while translation symmetry parallel to the interface fixes the tangential wave vector:

$$
\omega_{\mathrm i}=\omega_{\mathrm r}=\omega_{\mathrm t},\qquad \mathbf k_{\mathrm i,\parallel}=\mathbf k_{\mathrm r,\parallel}=\mathbf k_{\mathrm t,\parallel}.
$$

Therefore,

$$
\boxed{\theta_{\mathrm r}=\theta_{\mathrm i},\qquad n_1\sin\theta_{\mathrm i}=n_2\sin\theta_{\mathrm t}}.
$$

The following assumes a propagating transmitted wave, so $\theta_{\mathrm t}$ is real.

##### Fresnel Equations

With no free surface current, the [tangential boundary conditions](electromagnetic_maxwell_equations.md#boundary-conditions) are

$$
\hat{\mathbf z}\times(\mathbf E_{\mathrm i}+\mathbf E_{\mathrm r}-\mathbf E_{\mathrm t})=0,\qquad \hat{\mathbf z}\times(\mathbf H_{\mathrm i}+\mathbf H_{\mathrm r}-\mathbf H_{\mathrm t})=0.
$$

- **TE polarization:** At the interface, let the electric-field amplitudes be

  $$
  \mathbf E_{\mathrm i}=E_{\mathrm i}\hat{\mathbf y},\qquad \mathbf E_{\mathrm r}=E_{\mathrm r}\hat{\mathbf y},\qquad \mathbf E_{\mathrm t}=E_{\mathrm t}\hat{\mathbf y}.
  $$

  The magnetic-field amplitude of each wave follows from

  $$
  \mathbf H=\frac{1}{Z}\hat{\mathbf k}\times\mathbf E.
  $$

  Thus, for the wave-vector directions defined above,

  $$
  \mathbf H_{\mathrm i}=\frac{E_{\mathrm i}}{Z_1}(\cos\theta_{\mathrm i}\hat{\mathbf x}+\sin\theta_{\mathrm i}\hat{\mathbf z}),\qquad \mathbf H_{\mathrm r}=\frac{E_{\mathrm r}}{Z_1}(-\cos\theta_{\mathrm r}\hat{\mathbf x}+\sin\theta_{\mathrm r}\hat{\mathbf z}),\qquad \mathbf H_{\mathrm t}=\frac{E_{\mathrm t}}{Z_2}(\cos\theta_{\mathrm t}\hat{\mathbf x}+\sin\theta_{\mathrm t}\hat{\mathbf z}).
  $$

  Using $\theta_{\mathrm r}=\theta_{\mathrm i}$, the tangential boundary conditions give

  $$
  \boxed{E_{\mathrm i}+E_{\mathrm r}=E_{\mathrm t},\qquad \frac{\cos\theta_{\mathrm i}}{Z_1}(E_{\mathrm i}-E_{\mathrm r})=\frac{\cos\theta_{\mathrm t}}{Z_2}E_{\mathrm t}}.
  $$

  Hence

  $$
  \boxed{r_{\mathrm{TE}}=\frac{E_{\mathrm r}}{E_{\mathrm i}}=\frac{Z_2\cos\theta_{\mathrm i}-Z_1\cos\theta_{\mathrm t}}{Z_2\cos\theta_{\mathrm i}+Z_1\cos\theta_{\mathrm t}},\qquad t_{\mathrm{TE}}=\frac{E_{\mathrm t}}{E_{\mathrm i}}=\frac{2Z_2\cos\theta_{\mathrm i}}{Z_2\cos\theta_{\mathrm i}+Z_1\cos\theta_{\mathrm t}}}.
  $$

- **TM polarization:** At the interface, let

  $$
  \mathbf H_{\mathrm i}=H_{\mathrm i}\hat{\mathbf y},\qquad \mathbf H_{\mathrm r}=H_{\mathrm r}\hat{\mathbf y},\qquad \mathbf H_{\mathrm t}=H_{\mathrm t}\hat{\mathbf y}.
  $$

  The electric-field amplitudes follow from

  $$
  \mathbf E=-Z\hat{\mathbf k}\times\mathbf H.
  $$

  Thus

  $$
  \mathbf E_{\mathrm i}=-Z_1H_{\mathrm i}(\cos\theta_{\mathrm i}\hat{\mathbf x}+\sin\theta_{\mathrm i}\hat{\mathbf z}),\qquad \mathbf E_{\mathrm r}=Z_1H_{\mathrm r}(\cos\theta_{\mathrm r}\hat{\mathbf x}-\sin\theta_{\mathrm r}\hat{\mathbf z}),\qquad \mathbf E_{\mathrm t}=-Z_2H_{\mathrm t}(\cos\theta_{\mathrm t}\hat{\mathbf x}+\sin\theta_{\mathrm t}\hat{\mathbf z}).
  $$

  Using $\theta_{\mathrm r}=\theta_{\mathrm i}$, the tangential boundary conditions give

  $$
  \boxed{H_{\mathrm i}+H_{\mathrm r}=H_{\mathrm t},\qquad Z_1\cos\theta_{\mathrm i}(H_{\mathrm i}-H_{\mathrm r})=Z_2\cos\theta_{\mathrm t}H_{\mathrm t}}.
  $$

  Defining the coefficients using magnetic-field amplitudes,

  $$
  \boxed{r_{\mathrm{TM}}=\frac{H_{\mathrm r}}{H_{\mathrm i}}=\frac{Z_1\cos\theta_{\mathrm i}-Z_2\cos\theta_{\mathrm t}}{Z_1\cos\theta_{\mathrm i}+Z_2\cos\theta_{\mathrm t}},\qquad t_{\mathrm{TM}}=\frac{H_{\mathrm t}}{H_{\mathrm i}}=\frac{2Z_1\cos\theta_{\mathrm i}}{Z_1\cos\theta_{\mathrm i}+Z_2\cos\theta_{\mathrm t}}}.
  $$

Reflectance and transmittance are ratios of cycle-averaged energy flux *normal* to the interface:

$$
R_{\mathrm{TE}}=|r_{\mathrm{TE}}|^2,\qquad T_{\mathrm{TE}}=\frac{Z_1\cos\theta_{\mathrm t}}{Z_2\cos\theta_{\mathrm i}}|t_{\mathrm{TE}}|^2.
$$

$$
R_{\mathrm{TM}}=|r_{\mathrm{TM}}|^2,\qquad T_{\mathrm{TM}}=\frac{Z_2\cos\theta_{\mathrm t}}{Z_1\cos\theta_{\mathrm i}}|t_{\mathrm{TM}}|^2.
$$

For lossless media, $R_{\mathrm{TE}}+T_{\mathrm{TE}}=R_{\mathrm{TM}}+T_{\mathrm{TM}}=1$. At normal incidence,

$$
\boxed{R=\left|\frac{Z_2-Z_1}{Z_2+Z_1}\right|^2}.
$$

If $\mu_1=\mu_2$, TM reflection vanishes at the Brewster angle:

$$
\boxed{\tan\theta_{\mathrm B}=\frac{n_2}{n_1},\qquad R_{\mathrm{TM}}(\theta_{\mathrm B})=0}.
$$
