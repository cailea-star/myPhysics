# Quasistatic Fields

### Quasistatic Limits of Maxwell's Equations

Let $\rho_f$ and $\mathbf j_f$ be the free charge and current densities. The [Maxwell equations in matter](maxwell_equations.md#maxwell-equations-in-matter) are

$$
\nabla\cdot\mathbf D=\rho_f,\qquad \nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t},
$$

$$
\nabla\cdot\mathbf B=0,\qquad \nabla\times\mathbf H=\mathbf j_f+\frac{\partial\mathbf D}{\partial t}.
$$

For a homogeneous linear medium, let $\epsilon$ and $\mu$ denote permittivity and permeability:

$$
\mathbf D=\epsilon\mathbf E,\qquad \mathbf B=\mu\mathbf H.
$$

An Ohmic conductor additionally obeys the local Ohm law, with electrical conductivity $\sigma$:

$$
\mathbf j_f=\sigma\mathbf E.
$$

The two quasistatic limits discard different terms:

Electric quasistatics neglect Faraday induction at leading order:

$$
\boxed{\nabla\times\mathbf E\simeq0}.
$$

Magnetic quasistatics neglect the displacement current but retain Faraday induction:

$$
\boxed{\nabla\times\mathbf H\simeq\mathbf j_f,\qquad \nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t}}.
$$

### Field Diffusion in Conductors

Inside a homogeneous Ohmic conductor, use the magnetic quasistatic approximation. With constant $\mu$ and $\sigma$, the constitutive relations give

$$
\nabla\times\mathbf E=-\mu\frac{\partial\mathbf H}{\partial t},\qquad \nabla\times\mathbf H=\sigma\mathbf E.
$$

Within this approximation, Gauss's law for magnetism and the divergence of Ampère's law imply

$$
\nabla\cdot\mathbf E=0,\qquad \nabla\cdot\mathbf H=0.
$$

For any vector field $\mathbf F$, use $\nabla\times(\nabla\times\mathbf F)=\nabla(\nabla\cdot\mathbf F)-\nabla^2\mathbf F$. Taking the curl of each field equation then gives

$$
\nabla\times(\nabla\times\mathbf E)=-\nabla^2\mathbf E=-\mu\frac{\partial}{\partial t}(\nabla\times\mathbf H)=-\mu\sigma\frac{\partial\mathbf E}{\partial t},
$$

$$
\nabla\times(\nabla\times\mathbf H)=-\nabla^2\mathbf H=\sigma\nabla\times\mathbf E=-\mu\sigma\frac{\partial\mathbf H}{\partial t}.
$$

Thus the electric field and magnetic field intensity obey diffusion equations:

$$
\boxed{\frac{\partial\mathbf E}{\partial t}=\frac{1}{\mu\sigma}\nabla^2\mathbf E,\qquad \frac{\partial\mathbf H}{\partial t}=\frac{1}{\mu\sigma}\nabla^2\mathbf H}.
$$

Their diffusion coefficient is $1/(\mu\sigma)$.

### Skin Effect

Let a homogeneous conductor occupy $z>0$, with its surface at $z=0$. Suppose the electric field is polarized along $x$ and oscillates with angular frequency $\omega$. Use a complex field for the calculation; the physical field is its real part:

$$
\widetilde{\mathbf E}(z,t)=\mathcal E(z)e^{-i\omega t}\hat{\mathbf x}.
$$

Substitution into the diffusion equation gives

$$
\frac{\mathrm d^2\mathcal E}{\mathrm dz^2}+i\mu\sigma\omega\mathcal E=0.
$$

Let $E_0$ be the surface amplitude and try $\mathcal E(z)=E_0e^{-pz}$, where $p$ is a complex inverse length. Then

$$
p^2=-i\mu\sigma\omega,\qquad p=\alpha(1-i),\qquad \alpha=\sqrt{\frac{\mu\sigma\omega}{2}}.
$$

Choose the root with positive real part so that the field decays as $z\to\infty$. Thus

$$
\widetilde{\mathbf E}(z,t)=E_0e^{-\alpha z}e^{-i(\omega t-\alpha z)}\hat{\mathbf x}.
$$

If $E_0$ is real, the physical electric field is

$$
\mathbf E(z,t)=E_0e^{-\alpha z}\cos(\omega t-\alpha z)\hat{\mathbf x}.
$$

The amplitude falls to $1/e$ of its surface value at the skin depth $\delta$:

$$
\boxed{\delta=\frac1\alpha=\sqrt{\frac{2}{\mu\sigma\omega}}}.
$$

Since $\mathbf j_f=\sigma\mathbf E$, the current density has the same attenuation depth.
