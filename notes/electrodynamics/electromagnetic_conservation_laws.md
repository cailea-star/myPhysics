# Electromagnetic Conservation Laws

### Energy Conservation

Let $\mathbf j$ be the current density. The power density transferred from the field to charges is

$$
\mathcal P=\mathbf j\cdot\mathbf E.
$$

The magnetic field does no work. The Ampère–Maxwell equation gives

$$
\mathcal P=\frac{1}{\mu_0}\mathbf E\cdot(\nabla\times\mathbf B)-\epsilon_0\mathbf E\cdot\frac{\partial\mathbf E}{\partial t}.
$$

Using Faraday's law and a vector identity,

$$
\mathbf E\cdot(\nabla\times\mathbf B)=\mathbf B\cdot(\nabla\times\mathbf E)-\nabla\cdot(\mathbf E\times\mathbf B)=-\mathbf B\cdot\frac{\partial\mathbf B}{\partial t}-\nabla\cdot(\mathbf E\times\mathbf B).
$$

Define the field energy density $u$ and Poynting vector $\mathbf S_P$ by

$$
\boxed{u=\frac{\epsilon_0}{2}\mathbf E\cdot\mathbf E+\frac{1}{2\mu_0}\mathbf B\cdot\mathbf B,\qquad \mathbf S_P=\frac{1}{\mu_0}\mathbf E\times\mathbf B}.
$$

The local energy balance is

$$
\boxed{\frac{\partial u}{\partial t}+\nabla\cdot\mathbf S_P=-\mathcal P}.
$$

For a fixed volume $V$ with outward-oriented area element $\mathrm{d}\mathbf a$,

$$
\frac{\mathrm{d}}{\mathrm{d}t}\int_V u\,\mathrm{d}V+\oint_{\partial V}\mathbf S_P\cdot\mathrm{d}\mathbf a=-\int_V\mathcal P\,\mathrm{d}V.
$$

### Momentum Conservation

Let $\rho$ and $\mathbf j$ be the charge and current densities. The Lorentz force density $\mathbf f$ exerted on charges is

$$
\mathbf f=\rho\mathbf E+\mathbf j\times\mathbf B.
$$

Eliminating the sources with Maxwell's equations gives

$$
\mathbf f=\epsilon_0(\nabla\cdot\mathbf E)\mathbf E+\frac{1}{\mu_0}(\nabla\times\mathbf B)\times\mathbf B-\epsilon_0\frac{\partial\mathbf E}{\partial t}\times\mathbf B.
$$

Let $\hat{\mathbf I}$ be the identity operator and $\otimes$ denote the outer product. For a vector field $\mathbf F$,

$$
\mathbf F(\nabla\cdot\mathbf F)+(\nabla\times\mathbf F)\times\mathbf F=\nabla\cdot\left(\mathbf F\otimes\mathbf F-\frac{\mathbf F\cdot\mathbf F}{2}\hat{\mathbf I}\right).
$$

Using this identity and Faraday's law, define the field momentum density $\mathbf g$ and momentum flux tensor $\hat{\mathbf T}$. Let $c=1/\sqrt{\epsilon_0\mu_0}$ be the speed of light in vacuum. Then

$$
\boxed{\mathbf g=\epsilon_0\mathbf E\times\mathbf B=\frac{\mathbf S_P}{c^2},\qquad \hat{\mathbf T}=u\hat{\mathbf I}-\epsilon_0\mathbf E\otimes\mathbf E-\frac{1}{\mu_0}\mathbf B\otimes\mathbf B}.
$$

The local momentum balance is

$$
\boxed{\frac{\partial\mathbf g}{\partial t}+\nabla\cdot\hat{\mathbf T}=-\mathbf f}.
$$

For a fixed volume $V$ with outward-oriented area element $\mathrm{d}\mathbf a$,

$$
\frac{\mathrm{d}}{\mathrm{d}t}\int_V\mathbf g\,\mathrm{d}V+\oint_{\partial V}\hat{\mathbf T}\cdot\mathrm{d}\mathbf a=-\int_V\mathbf f\,\mathrm{d}V.
$$

Here $u$ and $\mathbf S_P$ are defined in [Energy Conservation](#energy-conservation). The sign of the tensor called "Maxwell stress" varies by convention; here $\hat{\mathbf T}$ is the outward momentum flux.

### Energy and Momentum Conservation in Matter

Consider a stationary, homogeneous, isotropic, linear, nondispersive medium with constant $\epsilon$ and $\mu$, so that $\mathbf D=\epsilon\mathbf E$ and $\mathbf B=\mu\mathbf H$. Let $\rho_f$ and $\mathbf j_f$ be the free charge and current densities. The power density delivered to free charges is

$$
\mathcal P=\mathbf j_f\cdot\mathbf E=-\nabla\cdot(\mathbf E\times\mathbf H)-\mathbf E\cdot\frac{\partial\mathbf D}{\partial t}-\mathbf H\cdot\frac{\partial\mathbf B}{\partial t}.
$$

The energy density and flux for the combined field–medium system are

$$
\boxed{u=\frac{1}{2}(\mathbf E\cdot\mathbf D+\mathbf B\cdot\mathbf H),\qquad \mathbf S_P=\mathbf E\times\mathbf H}.
$$

Hence,

$$
\boxed{\frac{\partial u}{\partial t}+\nabla\cdot\mathbf S_P=-\mathcal P}.
$$

The force density on free charges is

$$
\mathbf f=\rho_f\mathbf E+\mathbf j_f\times\mathbf B.
$$

The corresponding momentum density and flux tensor are

$$
\boxed{\mathbf g=\mathbf D\times\mathbf B,\qquad \hat{\mathbf T}=u\hat{\mathbf I}-\mathbf D\otimes\mathbf E-\mathbf B\otimes\mathbf H}.
$$

Thus,

$$
\boxed{\frac{\partial\mathbf g}{\partial t}+\nabla\cdot\hat{\mathbf T}=-\mathbf f}.
$$

Here $\mathbf g$ includes momentum assigned to the bound medium; it is not purely field momentum.
