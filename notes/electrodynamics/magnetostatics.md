# Magnetostatics

### Magnetostatic Vector Potential and Interface Conditions

Let $\mathbf A$ be the vector potential and $\mu(\mathbf r)$ the permeability. Then

$$
\mathbf B=\nabla\times\mathbf A,\qquad \mathbf H=\frac{\mathbf B}{\mu(\mathbf r)}=\frac{1}{\mu(\mathbf r)}\nabla\times\mathbf A.
$$

For free current density $\mathbf j_f$, the magnetostatic equations are

$$
\nabla\cdot\mathbf B=0,\qquad \nabla\times\mathbf H=\mathbf j_f+\frac{\partial\mathbf D}{\partial t}=\mathbf j_f.
$$

Substituting the field–potential relations gives

$$
\nabla\times\left[\frac{1}{\mu(\mathbf r)}\nabla\times\mathbf A\right]=\mathbf j_f.
$$

At an interface, $\hat{\mathbf n}$ points from medium 1 to medium 2, $\mathbf K_f$ is the free surface current density, and $\mu_i$ is the permeability on side $i$:

$$
\hat{\mathbf n}\cdot(\mathbf B_2-\mathbf B_1)=0,\qquad \hat{\mathbf n}\times(\mathbf H_2-\mathbf H_1)=\mathbf K_f.
$$

The first condition reads $\hat{\mathbf n}\cdot(\nabla\times\mathbf A_2-\nabla\times\mathbf A_1)=0$. For finite $\mathbf B$ and a nonsingular $\mathbf A$, a narrow surface $S_h$ of width $h$ across the interface gives

$$
\oint_{\partial S_h}\mathbf A\cdot\mathrm d\mathbf l=\int_{S_h}\mathbf B\cdot\mathrm d\mathbf a\longrightarrow0\qquad(h\to0).
$$

Hence $\hat{\mathbf n}\times(\mathbf A_2-\mathbf A_1)=0$. Substituting $\mathbf H_i=\mu_i^{-1}\nabla\times\mathbf A_i$ into the other condition yields

$$
\boxed{\hat{\mathbf n}\times(\mathbf A_2-\mathbf A_1)=0,\qquad \hat{\mathbf n}\times\left(\frac{\nabla\times\mathbf A_2}{\mu_2}-\frac{\nabla\times\mathbf A_1}{\mu_1}\right)=\mathbf K_f}.
$$

### Uniqueness Theorem

Let $V$ be a connected magnetic region with fixed free-current density $\mathbf j_f$ and permeability $\mu(\mathbf r)>0$, and let $\hat{\mathbf n}$ be the outward normal on $\partial V$. Consider two solutions with the same interface conditions and prescribed tangential vector potential on the outer boundary:

$$
\nabla\times\mathbf H=\nabla\times\mathbf H'=\mathbf j_f,\qquad \hat{\mathbf n}\times\mathbf A|_{\partial V}=\hat{\mathbf n}\times\mathbf A'|_{\partial V}.
$$

Define

$$
\delta\mathbf A=\mathbf A-\mathbf A',\qquad \delta\mathbf B=\mathbf B-\mathbf B'=\nabla\times\delta\mathbf A,\qquad \delta\mathbf H=\mathbf H-\mathbf H'=\frac{\delta\mathbf B}{\mu(\mathbf r)}.
$$

Subtracting the field equations and boundary data gives

$$
\nabla\times\delta\mathbf H=0,\qquad \hat{\mathbf n}\times\delta\mathbf A|_{\partial V}=0.
$$

The outer-boundary flux vanishes because

$$
(\delta\mathbf A\times\delta\mathbf H)\cdot\hat{\mathbf n}=(\hat{\mathbf n}\times\delta\mathbf A)\cdot\delta\mathbf H=0.
$$

The common interface conditions cancel the internal surface terms. Applying the divergence theorem in each region gives

$$
0=\oint_{\partial V}(\delta\mathbf A\times\delta\mathbf H)\cdot\mathrm d\mathbf a=\int_V\nabla\cdot(\delta\mathbf A\times\delta\mathbf H)\,\mathrm dV=\int_V\frac{|\delta\mathbf B|^2}{\mu(\mathbf r)}\,\mathrm dV.
$$

Since the permeability is positive,

$$
\boxed{\mathbf B=\mathbf B',\qquad \mathbf H=\mathbf H'}.
$$

The vector potential remains gauge-dependent.

### Two-Dimensional Vector-Potential Method

Assume translation symmetry along $z$. Let the free current and vector potential have only $z$ components:

$$
\mathbf j_f=j_z(x,y)\hat{\mathbf z},\qquad \mathbf A=A_z(x,y)\hat{\mathbf z}.
$$

In a homogeneous region,

$$
\mathbf B=(\partial_y A_z)\hat{\mathbf x}-(\partial_x A_z)\hat{\mathbf y},\qquad \boxed{\nabla_\perp^2 A_z=(\partial_x^2+\partial_y^2)A_z=-\mu j_z}.
$$

##### Magnetic Cylinder in a Uniform Transverse Field

Let a cylinder of radius $R$ and permeability $\mu_{\mathrm{in}}$ lie in a medium of permeability $\mu_{\mathrm{out}}$. There are no free currents at or near the cylinder. For an applied field $\mathbf B_0=B_0\hat{\mathbf x}$, cylindrical coordinates $(\varrho,\phi)$ give

$$
A_{\mathrm{ext},z}=B_0y=B_0\varrho\sin\phi.
$$

Only the $\sin\phi$ mode is needed. Regularity at the axis and the applied-field behavior at infinity suggest

$$
A_{\mathrm{out},z}=\left(B_0\varrho+\frac{a}{\varrho}\right)\sin\phi,\qquad A_{\mathrm{in},z}=b\varrho\sin\phi,
$$

where $a$ and $b$ are coefficients. With no free surface current, the [interface conditions](#magnetostatic-vector-potential-and-interface-conditions) at $\varrho=R$ become

$$
A_{\mathrm{out},z}=A_{\mathrm{in},z},\qquad \frac{1}{\mu_{\mathrm{out}}}\partial_\varrho A_{\mathrm{out},z}=\frac{1}{\mu_{\mathrm{in}}}\partial_\varrho A_{\mathrm{in},z}.
$$

Applying them gives

$$
a=\frac{\mu_{\mathrm{in}}-\mu_{\mathrm{out}}}{\mu_{\mathrm{in}}+\mu_{\mathrm{out}}}B_0R^2,\qquad b=\frac{2\mu_{\mathrm{in}}}{\mu_{\mathrm{in}}+\mu_{\mathrm{out}}}B_0.
$$

Thus,

$$
\boxed{A_{\mathrm{out},z}=B_0\left(\varrho+\frac{\mu_{\mathrm{in}}-\mu_{\mathrm{out}}}{\mu_{\mathrm{in}}+\mu_{\mathrm{out}}}\frac{R^2}{\varrho}\right)\sin\phi,\qquad A_{\mathrm{in},z}=\frac{2\mu_{\mathrm{in}}}{\mu_{\mathrm{in}}+\mu_{\mathrm{out}}}B_0\varrho\sin\phi}.
$$

The field inside the cylinder is uniform:

$$
\mathbf B_{\mathrm{in}}=\frac{2\mu_{\mathrm{in}}}{\mu_{\mathrm{in}}+\mu_{\mathrm{out}}}\mathbf B_0.
$$

### Magnetic Scalar Potential

In a simply connected region without free current, Ampère's law gives

$$
\nabla\times\mathbf H=0\qquad\Longrightarrow\qquad\boxed{\mathbf H=-\nabla\Phi_{\mathrm{mag}}}.
$$

For a linear medium with permeability $\mu(\mathbf r)$, Gauss's law for magnetism gives

$$
0=\nabla\cdot\mathbf B=\nabla\cdot\!\left[\mu(\mathbf r)\mathbf H\right]=-\nabla\cdot\!\left[\mu(\mathbf r)\nabla\Phi_{\mathrm{mag}}(\mathbf r)\right].
$$

Thus $\nabla^2\Phi_{\mathrm{mag}}=0$ within each homogeneous region. At an interface without free surface current, let $\hat{\mathbf n}$ point from medium 1 to medium 2. The interface conditions are

$$
\boxed{\Phi_{\mathrm{mag},1}=\Phi_{\mathrm{mag},2},\qquad \mu_1\partial_{\hat{\mathbf n}}\Phi_{\mathrm{mag},1}=\mu_2\partial_{\hat{\mathbf n}}\Phi_{\mathrm{mag},2}}.
$$

##### Magnetic Sphere in a Uniform Field

Let a sphere of radius $R$ and permeability $\mu_{\mathrm{in}}$ lie in a medium of permeability $\mu_{\mathrm{out}}$, under $\mathbf H_0=H_0\hat{\mathbf z}$. As in the [dielectric sphere in a uniform field](electrostatics_in_matter.md#dielectric-sphere-in-a-uniform-field), spherical coordinates $(r,\theta)$ give the applied-field potential

$$
\Phi_{\mathrm{mag,ext}}=-\mathbf H_0\cdot\mathbf r=-H_0z=-H_0r\cos\theta.
$$

Only the $\ell=1$ mode is needed. With coefficients $c_{\mathrm{out}}$ and $c_{\mathrm{in}}$, regularity at the origin and decay of the induced potential at infinity give

$$
\Phi_{\mathrm{mag,out}}=\left(-H_0r+\frac{c_{\mathrm{out}}}{r^2}\right)\cos\theta,\qquad \Phi_{\mathrm{mag,in}}=c_{\mathrm{in}}r\cos\theta.
$$

Applying the interface conditions at $r=R$ gives

$$
c_{\mathrm{out}}=\frac{\mu_{\mathrm{in}}-\mu_{\mathrm{out}}}{\mu_{\mathrm{in}}+2\mu_{\mathrm{out}}}H_0R^3,\qquad c_{\mathrm{in}}=-\frac{3\mu_{\mathrm{out}}}{\mu_{\mathrm{in}}+2\mu_{\mathrm{out}}}H_0.
$$

Thus,

$$
\boxed{\Phi_{\mathrm{mag,out}}=-H_0r\cos\theta+\frac{\mu_{\mathrm{in}}-\mu_{\mathrm{out}}}{\mu_{\mathrm{in}}+2\mu_{\mathrm{out}}}\frac{H_0R^3}{r^2}\cos\theta,\qquad \Phi_{\mathrm{mag,in}}=-\frac{3\mu_{\mathrm{out}}}{\mu_{\mathrm{in}}+2\mu_{\mathrm{out}}}H_0r\cos\theta}.
$$

With $\mathbf B_0=\mu_{\mathrm{out}}\mathbf H_0$, the magnetic induction inside the sphere is uniform:

$$
\mathbf B_{\mathrm{in}}=\frac{3\mu_{\mathrm{in}}}{\mu_{\mathrm{in}}+2\mu_{\mathrm{out}}}\mathbf B_0.
$$

### Magnetic Multipoles and Dipoles

##### Far-Field Vector Potential

Let $\mathbf j(\mathbf r')$ be a localized steady current density in vacuum, with $r=|\mathbf r|$, $r'=|\mathbf r'|$, and vacuum permeability $\mu_0$. Its vector potential is

$$
\mathbf A(\mathbf r)=\frac{\mu_0}{4\pi}\int\frac{\mathbf j(\mathbf r')}{|\mathbf r-\mathbf r'|}\,\mathrm d^3r'.
$$

For $r\gg r'$, Taylor expansion in $\mathbf r'$ gives

$$
\frac{1}{|\mathbf r-\mathbf r'|}=\frac1r-\mathbf r'\cdot\nabla\left(\frac1r\right)+\frac12(\mathbf r'\otimes\mathbf r'):(\nabla\otimes\nabla)\left(\frac1r\right)+\cdots.
$$

Here $\nabla$ acts on $\mathbf r$; $\otimes$ and $:$ denote outer product and double contraction. The first derivative needed at dipole order is

$$
\nabla\left(\frac1r\right)=-\frac{\mathbf r}{r^3}.
$$

Steady current satisfies $\nabla'\cdot\mathbf j=0$, where $\nabla'$ acts on $\mathbf r'$. Because the current is localized, the zeroth-order term vanishes:

$$
\mathbf A^{(0)}(\mathbf r)=\frac{\mu_0}{4\pi r}\int\mathbf j(\mathbf r')\,\mathrm d^3r'=0.
$$

The first-order term is

$$
\mathbf A^{(1)}(\mathbf r)=\frac{\mu_0}{4\pi r^3}\int(\mathbf r\cdot\mathbf r')\mathbf j(\mathbf r')\,\mathrm d^3r'.
$$

Integration by parts, using $\nabla'\cdot\mathbf j=0$, gives the antisymmetric first-moment identity

$$
\int\left[(\mathbf r\cdot\mathbf r')\mathbf j+(\mathbf r\cdot\mathbf j)\mathbf r'\right]\mathrm d^3r'=0.
$$

Consequently, the first-order term takes the form

$$
\mathbf A^{(1)}(\mathbf r)=\frac{\mu_0}{8\pi r^3}\int\left[(\mathbf r\cdot\mathbf r')\mathbf j-(\mathbf r\cdot\mathbf j)\mathbf r'\right]\mathrm d^3r'=\frac{\mu_0}{8\pi r^3}\int(\mathbf r'\times\mathbf j)\times\mathbf r\,\mathrm d^3r'.
$$

The source integral therefore defines the magnetic dipole moment:

$$
\boxed{\mathbf m=\frac12\int\mathbf r'\times\mathbf j(\mathbf r')\,\mathrm d^3r'}.
$$

For a thin loop with current $I$ and oriented area vector $\mathbf S$, $\mathbf m=I\mathbf S$. Substitution gives the leading far-field vector potential:

$$
\boxed{\mathbf A(\mathbf r)=\frac{\mu_0}{4\pi}\frac{\mathbf m\times\mathbf r}{r^3}+\cdots}.
$$

With $\hat{\mathbf r}=\mathbf r/r$, the corresponding dipole field is

$$
\mathbf B_{\mathrm{dip}}=\frac{\mu_0}{4\pi r^3}\left[3(\mathbf m\cdot\hat{\mathbf r})\hat{\mathbf r}-\mathbf m\right].
$$

Outside the dipole source, it can also be represented by the magnetic scalar potential

$$
\Phi_{\mathrm{mag,dip}}=\frac{\mathbf m\cdot\mathbf r}{4\pi r^3},\qquad \mathbf H_{\mathrm{dip}}=-\nabla\Phi_{\mathrm{mag,dip}}.
$$

##### Interaction with an External Field

Let a distant steady current $\mathbf j_{\mathrm{ext}}$ produce $\mathbf A_{\mathrm{ext}}$ and $\mathbf B_{\mathrm{ext}}=\nabla\times\mathbf A_{\mathrm{ext}}$. The magnetic-field cross-energy is

$$
U_{B,\mathrm{int}}=\frac12\int\left(\mathbf A_{\mathrm{ext}}\cdot\mathbf j+\mathbf A\cdot\mathbf j_{\mathrm{ext}}\right)\mathrm d^3r=\int\mathbf j\cdot\mathbf A_{\mathrm{ext}}\,\mathrm d^3r.
$$

As in the [electric multipole interaction](electrostatics_in_matter.md#interaction-with-an-external-field), expand the external potential about the origin:

$$
\mathbf A_{\mathrm{ext}}(\mathbf r')=\mathbf A_{\mathrm{ext}}(0)+(\mathbf r'\cdot\nabla)\mathbf A_{\mathrm{ext}}(0)+\frac12(\mathbf r'\otimes\mathbf r'):(\nabla\otimes\nabla)\mathbf A_{\mathrm{ext}}(0)+\cdots.
$$

The constant term vanishes because $\int\mathbf j\,\mathrm d^3r'=0$. Using the magnetic moment defined above, the first-order term gives

$$
U_{B,\mathrm{int}}=\mathbf m\cdot\mathbf B_{\mathrm{ext}}(0)+\cdots.
$$

This is the *field* cross-energy. Let $I_{\mathrm{ext}}$ be the external-source current, $\Psi_{\mathrm{ext}}$ the external flux through the dipole loop, and $\Psi$ the dipole flux through the external circuit. If both currents remain fixed, their power supplies do work

$$
\Delta W_{\mathrm{src}}=I\Delta\Psi_{\mathrm{ext}}+I_{\mathrm{ext}}\Delta\Psi=2\Delta U_{B,\mathrm{int}}.
$$

Including the supplies and choosing the zero of energy appropriately, the effective interaction energy is

$$
\boxed{U_{\mathrm{eff}}=-\mathbf m\cdot\mathbf B_{\mathrm{ext}}}.
$$

##### Dipole Force and Torque

Let a rigid dipole with moment $\mathbf m$ be centered at $\mathbf R$. Neglecting higher multipoles, its interaction energy is

$$
U_{\mathrm d}(\mathbf R)=-\mathbf m\cdot\mathbf B_{\mathrm{ext}}(\mathbf R).
$$

Hold $\mathbf m$ fixed when translating the dipole. If the region around it contains no external current, $\nabla_{\mathbf R}\times\mathbf B_{\mathrm{ext}}=0$. With $\boldsymbol\tau$ denoting torque,

$$
\boxed{\mathbf F=-\nabla_{\mathbf R}U_{\mathrm d}=(\mathbf m\cdot\nabla_{\mathbf R})\mathbf B_{\mathrm{ext}}(\mathbf R),\qquad \boldsymbol\tau=\mathbf m\times\mathbf B_{\mathrm{ext}}(\mathbf R)}.
$$

A uniform field can exert a torque, but no net force.
