# Gamma Decay

### Relativistic Kinematics

Let $P$, $D$, and $\gamma$ denote the parent nucleus, daughter nucleus, and photon, with nuclear masses $M_P$ and $M_D$.

The decay process is

$$
P\rightarrow D+\gamma.
$$

Let $E_\gamma$ be the photon energy and $T_D=E_D-M_Dc^2$ the daughter recoil energy. The decay energy is

$$
Q_\gamma=(M_P-M_D)c^2=T_D+E_\gamma.
$$

Four-momentum conservation and the on-shell conditions are

$$
p_P^\mu=p_D^\mu+p_\gamma^\mu,\qquad p^\mu=(E,\mathbf p c),\qquad p^\mu p_\mu=E^2-\mathbf p^2c^2=M^2c^4,\qquad p_\gamma^\mu p_{\gamma\mu}=0.
$$

In the parent-nucleus rest frame, $p_P^\mu=(M_Pc^2,\mathbf0)$.

The photon on-shell condition gives

$$
0=(p_P-p_D)^2=M_P^2c^4+M_D^2c^4-2M_Pc^2E_D.
$$

Therefore,

$$
E_D=\frac{M_P^2+M_D^2}{2M_P}c^2.
$$

The corresponding daughter recoil energy is

$$
T_D=\frac{(M_P-M_D)^2}{2M_P}c^2=\frac{Q_\gamma^2}{2M_Pc^2}.
$$

The photon energy is

$$
\boxed{E_\gamma=M_Pc^2-E_D=\frac{M_P^2-M_D^2}{2M_P}c^2=Q_\gamma\left(1-\frac{Q_\gamma}{2M_Pc^2}\right)}.
$$

For $Q_\gamma\ll M_Pc^2$,

$$
E_\gamma\simeq Q_\gamma,\qquad T_D\ll E_\gamma.
$$

### Electromagnetic Fields

Let $\rho$ and $\mathbf j$ be the charge and current densities. In SI units, Maxwell’s equations are

$$
\nabla\cdot\mathbf E=\frac{\rho}{\epsilon_0},\qquad \nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t},\qquad \nabla\cdot\mathbf B=0,\qquad \nabla\times\mathbf B=\mu_0\mathbf j+\frac{1}{c^2}\frac{\partial\mathbf E}{\partial t},
$$

where $\mu_0\epsilon_0=1/c^2$.

Define the scalar and vector potentials by

$$
\mathbf B=\nabla\times\mathbf A,\qquad \mathbf E=-\nabla\Phi-\frac{\partial\mathbf A}{\partial t}.
$$

The equations for the potentials are

$$
\nabla^2\Phi+\frac{\partial}{\partial t}\left(\nabla\cdot\mathbf A\right)=-\frac{\rho}{\epsilon_0},
$$

$$
\left(\frac{1}{c^2}\frac{\partial^2}{\partial t^2}-\nabla^2\right)\mathbf A+\nabla\left(\nabla\cdot\mathbf A+\frac{1}{c^2}\frac{\partial\Phi}{\partial t}\right)=\mu_0\mathbf j.
$$

For a source-free radiation field,

$$
\rho=0,\qquad \mathbf j=\mathbf0.
$$

In the Coulomb gauge,

$$
\nabla\cdot\mathbf A=0.
$$

With radiation boundary conditions, $\Phi=0$, and the vector potential satisfies

$$
\left(\frac{1}{c^2}\frac{\partial^2}{\partial t^2}-\nabla^2\right)\mathbf A=0,\qquad \nabla\cdot\mathbf A=0.
$$

For a plane-wave mode,

$$
\mathbf A(\mathbf r,t)=\mathbf A_0e^{i(\mathbf k\cdot\mathbf r-\omega t)},
$$

the wave equation and transverse condition give

$$
\mathbf k\cdot\mathbf A_0=0,\qquad \omega=ck,\qquad E_\gamma=\hbar\omega=\hbar ck.
$$

### Vector Spherical Harmonics

The definition, phase convention, and general properties are collected in [Vector Spherical Harmonics](../mathematical_physics/specialfunction_spherical_harmonics.md#vector-spherical-harmonics). The notation used here is $\mathbf Y_{JlM}\equiv\mathbf Y_{JM}^{l}$.

Let $Y_{lm}(\hat{\mathbf r})$ be a scalar spherical harmonic and $\mathbf e_\mu$ a spherical unit vector with $\mu=-1,0,1$. Coupling orbital angular momentum $l$ to the spin-one vector basis gives

$$
\hat{\mathbf J}=\hat{\mathbf L}+\hat{\mathbf S},\qquad \mathbf Y_{JlM}(\hat{\mathbf r})\equiv\sum_{m,\mu}\langle lm,1\mu|JM\rangle Y_{lm}(\hat{\mathbf r})\mathbf e_\mu,
$$

where $J$ and $M$ are the total angular momentum and its projection, with $M=m+\mu$.

The vector spherical harmonics satisfy

$$
\hat{\mathbf J}^2\mathbf Y_{JlM}=\hbar^2J(J+1)\mathbf Y_{JlM},\qquad \hat{\mathbf L}^2\mathbf Y_{JlM}=\hbar^2l(l+1)\mathbf Y_{JlM},
$$

$$
\hat{\mathbf S}^2\mathbf Y_{JlM}=2\hbar^2\mathbf Y_{JlM},\qquad \hat J_z\mathbf Y_{JlM}=\hbar M\mathbf Y_{JlM}.
$$

Their orthonormality relation is

$$
\int d\Omega\,\mathbf Y_{J'l'M'}^*(\hat{\mathbf r})\cdot\mathbf Y_{JlM}(\hat{\mathbf r})=\delta_{JJ'}\delta_{ll'}\delta_{MM'}.
$$
