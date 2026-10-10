# Skyrme effective interaction

## Two-body matrix elements and density matrices

### Single-particle states and coordinate notation

Let $|\alpha\rangle$ denote orthonormal single-particle states with coordinate wave functions

$$
\psi_\alpha(1)=\langle1|\alpha\rangle,\qquad 1=(\mathbf r_1,s_1,q_1),
$$

where $s_1=\uparrow,\downarrow$ and $q_1=n,p$. Define

$$
\int \mathrm d1\equiv\sum_{s_1,q_1}\int \mathrm d^3\mathbf r_1,\qquad \int \mathrm d1\,\psi_\alpha^*(1)\psi_\beta(1)=\delta_{\alpha\beta}.
$$

Here $|\Phi\rangle$ is a normalized Slater determinant constructed from the occupied single-particle states.

### Antisymmetrized matrix elements

For a two-body interaction $\hat v$ symmetric under particle exchange, the many-body interaction is

$$
\hat V=\frac12\sum_{i\ne j}\hat v(i,j)=\sum_{i<j}\hat v(i,j).
$$

Let $P_r,P_\sigma,P_\tau$ denote spatial, spin, and isospin exchange, respectively. Then

$$
P_{12}=P_rP_\sigma P_\tau,\qquad P_{12}|\alpha\beta\rangle=|\beta\alpha\rangle.
$$

Here $|\alpha\beta\rangle=|\alpha\rangle\otimes|\beta\rangle$ is a two-body product state. Define

$$
\bar v=\hat v(1-P_{12}).
$$

The interaction energy is

$$
E_V=\langle\Phi|\hat V|\Phi\rangle=\frac12\sum_{\alpha,\beta\in\mathrm{occ}}\langle\alpha\beta|\bar v|\alpha\beta\rangle=\frac12\sum_{\alpha,\beta\in\mathrm{occ}}\left[\langle\alpha\beta|\hat v|\alpha\beta\rangle-\langle\alpha\beta|\hat v|\beta\alpha\rangle\right].
$$

The two terms in brackets are the direct and exchange contributions.

### Coordinate integrals

Define the general two-body kernel

$$
v(1,2;1',2')=\langle12|\hat v|1'2'\rangle.
$$

The interaction energy is

$$
E_V=\frac12\sum_{\alpha,\beta\in\mathrm{occ}}\int \mathrm d1\,\mathrm d2\,\mathrm d1'\,\mathrm d2'\,\psi_\alpha^*(1)\psi_\beta^*(2)v(1,2;1',2')\left[\psi_\alpha(1')\psi_\beta(2')-\psi_\beta(1')\psi_\alpha(2')\right].
$$

If the interaction can be represented by coordinate, internal-state, and differential operators acting on the wave functions, this becomes

$$
E_V=\frac12\sum_{\alpha,\beta\in\mathrm{occ}}\int \mathrm d1\,\mathrm d2\,\psi_\alpha^*(1)\psi_\beta^*(2)\hat v(1,2)\left[\psi_\alpha(1)\psi_\beta(2)-\psi_\beta(1)\psi_\alpha(2)\right].
$$

Right-acting gradients act on the ket, left-acting gradients on the bra; spin and isospin operators act on the corresponding internal indices. This expression uses the interaction $\hat v$ before antisymmetrization.

### One-body density matrix

Define the one-body density matrix

$$
\rho(1,2)=\sum_{\alpha\in\mathrm{occ}}\psi_\alpha(1)\psi_\alpha^*(2).
$$

Summing over occupied states gives

$$
E_V=\frac12\int \mathrm d1\,\mathrm d2\,\mathrm d1'\,\mathrm d2'\,v(1,2;1',2')\left[\rho(1',1)\rho(2',2)-\rho(2',1)\rho(1',2)\right].
$$

The two density-matrix pairings correspond to the direct and exchange contributions.

Assume no proton–neutron mixing. For nucleon species $q=n,p$,

$$
\rho_{s,s'}(q;\mathbf r,\mathbf r')=\sum_{\alpha\in\mathrm{occ},q}\psi_\alpha(\mathbf r,s,q)\psi_\alpha^*(\mathbf r',s',q).
$$

The sum includes only occupied states of species $q$; $s,s'$ are the row and column indices of the spin matrix.

The total spin density matrix and its local counterpart are defined by

$$
\rho_{s,s'}(\mathbf r,\mathbf r')=\sum_q\rho_{s,s'}(q;\mathbf r,\mathbf r'),\qquad \rho_{s,s'}(q;\mathbf r)\equiv\rho_{s,s'}(q;\mathbf r,\mathbf r).
$$

### Particle and spin densities

Let $\sigma_\nu$ denote the single-particle Pauli matrices, with $\nu=x,y,z$. The scalar particle and spin densities are

$$
\rho(q;\mathbf r,\mathbf r')=\sum_s\rho_{s,s}(q;\mathbf r,\mathbf r'),\qquad
s_\nu(q;\mathbf r,\mathbf r')=\sum_{s,s'}(\sigma_\nu)_{s',s}\rho_{s,s'}(q;\mathbf r,\mathbf r').
$$

The corresponding matrix decomposition is

$$
\rho_{s,s'}(q;\mathbf r,\mathbf r')=\frac12\Bigl[\rho(q;\mathbf r,\mathbf r')\delta_{s,s'}+\sum_\nu s_\nu(q;\mathbf r,\mathbf r')(\sigma_\nu)_{s,s'}\Bigr].
$$

Taking coincident coordinates defines the local densities and their totals:

$$
\rho(q;\mathbf r)=\rho(q;\mathbf r,\mathbf r),\qquad s_\nu(q;\mathbf r)=s_\nu(q;\mathbf r,\mathbf r).
$$

$$
\rho(\mathbf r)=\sum_q\rho(q;\mathbf r),\qquad \mathbf s(\mathbf r)=\sum_q\mathbf s(q;\mathbf r).
$$

Using the Pauli-matrix contraction identities gives

$$
\sum_{s,s'}\rho_{s,s}(q;\mathbf r)\rho_{s',s'}(q';\mathbf r)=\rho(q;\mathbf r)\rho(q';\mathbf r),
$$

$$
\sum_{s,s'}\rho_{s,s'}(q;\mathbf r)\rho_{s',s}(q';\mathbf r)=\frac12\Bigl[\rho(q;\mathbf r)\rho(q';\mathbf r)+\mathbf s(q;\mathbf r)\cdot\mathbf s(q';\mathbf r)\Bigr].
$$

For time-reversal symmetry, $\mathbf s(q;\mathbf r)=0$, so

$$
\rho_{s,s'}(q;\mathbf r)=\frac{\rho(q;\mathbf r)}2\delta_{s,s'},\qquad
\sum_{s,s'}\rho_{s,s'}(q;\mathbf r)\rho_{s',s}(q';\mathbf r)=\frac12\rho(q;\mathbf r)\rho(q';\mathbf r).
$$

This simplification applies only to the local matrix; the nonlocal spin density need not vanish.

## Skyrme interaction

Without explicit tensor terms, the Skyrme interaction is

$$
\boxed{v_{\mathrm{Sk}}(1,2;\rho)=v_{t_0}(1,2)+v_{t_1}(1,2)+v_{t_2}(1,2)+v_{t_3}(1,2;\rho)+v_{\mathrm{SO}}(1,2).}
$$

The $t_0$ term is zero-range; $t_1,t_2$ are momentum-dependent; $t_3$ is density-dependent; and $\mathrm{SO}$ denotes the spin-orbit term. The total particle density is

$$
\rho(\mathbf r)=\sum_{q=n,p}\rho(q;\mathbf r).
$$

The corresponding many-body interaction terms are

$$
\hat V_a=\frac12\sum_{i\ne j}v_a(i,j)=\sum_{i<j}v_a(i,j),\qquad a\in\{t_0,t_1,t_2,t_3,\mathrm{SO}\}.
$$

Here $v_{t_3}$ and $\hat V_{t_3}$ depend on $\rho$. The intrinsic effective Hamiltonian is

$$
\boxed{\hat H[\rho]=\hat T-\hat T_{\mathrm{cm}}+\hat V_{t_0}+\hat V_{t_1}+\hat V_{t_2}+\hat V_{t_3}[\rho]+\hat V_{\mathrm{SO}}+\hat V_{\mathrm C}.}
$$

For $A$ nucleons of common mass $m$, the kinetic and center-of-mass kinetic energies are

$$
\hat T=\sum_{i=1}^{A}\frac{\hat{\mathbf p}_i^2}{2m},\qquad \hat T_{\mathrm{cm}}=\frac{\left(\sum_{i=1}^{A}\hat{\mathbf p}_i\right)^2}{2mA}.
$$

$\hat V_{\mathrm C}$ is the proton–proton Coulomb interaction.

For a Slater determinant $|\Phi\rangle$, each interaction energy is

$$
E_a[\Phi]=\frac12\sum_{\alpha,\beta\in\mathrm{occ}}\langle\alpha\beta|v_a[\rho_\Phi](1-P_{12})|\alpha\beta\rangle.
$$

Only the $t_3$ term explicitly depends on $\rho_\Phi$. The total energy is

$$
E[\Phi]=E_{\mathrm{kin}}-E_{\mathrm{cm}}+E_{t_0}+E_{t_1}+E_{t_2}+E_{t_3}+E_{\mathrm{SO}}+E_{\mathrm C}.
$$

### t0 term

The $t_0$ two-body interaction and its many-body term are

$$
\boxed{\hat v_{t_0}(1,2)=t_0(1+x_0P_\sigma)\delta(\mathbf r_1-\mathbf r_2),\qquad \hat V_{t_0}=\frac12\sum_{i\ne j}\hat v_{t_0}(i,j).}
$$

Here $t_0,x_0$ are the strength and spin-exchange parameters, respectively, and $\delta$ is the three-dimensional delta function.

The spin-exchange operator is

$$
P_\sigma=\frac12\left(I_2\otimes I_2+\sum_{\mu=x,y,z}\sigma_\mu\otimes\sigma_\mu\right),
$$

$$
P_\sigma\Bigl[\psi_\alpha(\mathbf r_1,s_1,q_1)\psi_\beta(\mathbf r_2,s_2,q_2)\Bigr]=\psi_\alpha(\mathbf r_1,s_2,q_1)\psi_\beta(\mathbf r_2,s_1,q_2).
$$

Here $\sigma_\mu$ are the single-particle Pauli matrices; $P_\sigma$ exchanges only the spin indices.

The interaction energy is

$$
E_{t_0}=\frac{t_0}{2}\sum_{\alpha,\beta\in\mathrm{occ}}\int \mathrm d1\,\mathrm d2\,\psi_\alpha^*(1)\psi_\beta^*(2)(1+x_0P_\sigma)\delta(\mathbf r_1-\mathbf r_2)\Bigl[\psi_\alpha(1)\psi_\beta(2)-\psi_\beta(1)\psi_\alpha(2)\Bigr].
$$

Without proton–neutron mixing, the four contractions are written directly in terms of the density matrix defined in the first section. The sums run over occupied states; the expressions omit $t_0/2$, the exchange minus sign, and the factor $x_0$.

**Direct term $t_0$**

$$
\sum_{\alpha,\beta}\psi_\alpha^*(1)\psi_\beta^*(2)\delta(\mathbf r_1-\mathbf r_2)\psi_\alpha(1)\psi_\beta(2)=\delta(\mathbf r_1-\mathbf r_2)\rho_{s_1,s_1}(q_1;\mathbf r_1)\rho_{s_2,s_2}(q_2;\mathbf r_2).
$$

**Direct term $t_0x_0$**

$$
\sum_{\alpha,\beta}\psi_\alpha^*(1)\psi_\beta^*(2)P_\sigma\delta(\mathbf r_1-\mathbf r_2)\psi_\alpha(1)\psi_\beta(2)=\delta(\mathbf r_1-\mathbf r_2)\rho_{s_2,s_1}(q_1;\mathbf r_1)\rho_{s_1,s_2}(q_2;\mathbf r_2).
$$

**Exchange term $t_0$**

$$
\sum_{\alpha,\beta}\psi_\alpha^*(1)\psi_\beta^*(2)\delta(\mathbf r_1-\mathbf r_2)\psi_\beta(1)\psi_\alpha(2)=\delta(\mathbf r_1-\mathbf r_2)\delta_{q_1,q_2}\rho_{s_2,s_1}(q_1;\mathbf r_2,\mathbf r_1)\rho_{s_1,s_2}(q_1;\mathbf r_1,\mathbf r_2).
$$

**Exchange term $t_0x_0$**

$$
\sum_{\alpha,\beta}\psi_\alpha^*(1)\psi_\beta^*(2)P_\sigma\delta(\mathbf r_1-\mathbf r_2)\psi_\beta(1)\psi_\alpha(2)=\delta(\mathbf r_1-\mathbf r_2)\delta_{q_1,q_2}\rho_{s_1,s_1}(q_1;\mathbf r_2,\mathbf r_1)\rho_{s_2,s_2}(q_1;\mathbf r_1,\mathbf r_2).
$$

Use the delta function to set $\mathbf r_1=\mathbf r_2=\mathbf r$, and perform the spin sums using the particle and spin densities of the first section. Define

$$
\rho(\mathbf r)=\sum_q\rho(q;\mathbf r),\qquad \mathbf s(\mathbf r)=\sum_q\mathbf s(q;\mathbf r).
$$

Adding the two direct terms and summing over internal indices gives the Hartree energy:

$$
E_{t_0}^{\mathrm{Hartree}}=\frac{t_0}{2}\int \mathrm d^3\mathbf r\,\sum_{q,q'}\Bigl[\rho(q;\mathbf r)\rho(q';\mathbf r)+\frac{x_0}{2}\bigl(\rho(q;\mathbf r)\rho(q';\mathbf r)+\mathbf s(q;\mathbf r)\cdot\mathbf s(q';\mathbf r)\bigr)\Bigr]=\frac{t_0}{4}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_0)\rho(\mathbf r)^2+x_0|\mathbf s(\mathbf r)|^2\Bigr].
$$

The exchange factor $\delta_{q_1,q_2}$ restricts the contractions to the same species, giving

$$
E_{t_0}^{\mathrm{Fock}}=-\frac{t_0}{2}\int \mathrm d^3\mathbf r\,\sum_q\Bigl[\frac12\bigl(\rho(q;\mathbf r)^2+|\mathbf s(q;\mathbf r)|^2\bigr)+x_0\rho(q;\mathbf r)^2\Bigr]=-\frac{t_0}{4}\int \mathrm d^3\mathbf r\,\sum_q\Bigl[(1+2x_0)\rho(q;\mathbf r)^2+|\mathbf s(q;\mathbf r)|^2\Bigr].
$$

The total $t_0$ energy is

$$
E_{t_0}=E_{t_0}^{\mathrm{Hartree}}+E_{t_0}^{\mathrm{Fock}}=\frac{t_0}{4}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_0)\rho(\mathbf r)^2-(1+2x_0)\sum_q\rho(q;\mathbf r)^2+x_0|\mathbf s(\mathbf r)|^2-\sum_q|\mathbf s(q;\mathbf r)|^2\Bigr].
$$

Varying the particle and spin densities separately while holding the remaining densities fixed gives

$$
\frac{\delta E_{t_0}}{\delta\rho(q;\mathbf r)}=\frac{t_0}{2}\Bigl[(2+x_0)\rho(\mathbf r)-(1+2x_0)\rho(q;\mathbf r)\Bigr],
$$

$$
\frac{\delta E_{t_0}}{\delta\mathbf s(q;\mathbf r)}=\frac{t_0}{2}\Bigl[x_0\mathbf s(\mathbf r)-\mathbf s(q;\mathbf r)\Bigr].
$$

For **time-reversal symmetry**,

$$
\mathbf s(q;\mathbf r)=0.
$$

The energy reduces to

$$
\boxed{E_{t_0}=\frac{t_0}{4}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_0)\rho(\mathbf r)^2-(1+2x_0)\sum_q\rho(q;\mathbf r)^2\Bigr].}
$$

The corresponding variations are

$$
\frac{\delta E_{t_0}}{\delta\rho(q;\mathbf r)}=\frac{t_0}{2}\Bigl[(2+x_0)\rho(\mathbf r)-(1+2x_0)\rho(q;\mathbf r)\Bigr],\qquad \frac{\delta E_{t_0}}{\delta\mathbf s(q;\mathbf r)}=0.
$$

### t1 term

The $t_1$ two-body interaction and its many-body term are

$$
\boxed{v_{t_1}(1,2)=\frac{t_1}{2}(1+x_1P_\sigma)\Bigl[\overleftarrow{\mathbf k}^{\,2}\delta(\mathbf r_1-\mathbf r_2)+\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}^{\,2}\Bigr],\qquad \hat V_{t_1}=\frac12\sum_{i\ne j}v_{t_1}(i,j).}
$$

The relative wave-number operators are

$$
\overrightarrow{\mathbf k}=\frac{\overrightarrow{\boldsymbol\nabla}_1-\overrightarrow{\boldsymbol\nabla}_2}{2i},\qquad \overleftarrow{\mathbf k}=-\frac{\overleftarrow{\boldsymbol\nabla}_1-\overleftarrow{\boldsymbol\nabla}_2}{2i}.
$$

The left and right derivatives act on the corresponding wave functions, not on the intervening $\delta(\mathbf r_1-\mathbf r_2)$.

The corresponding energy is

$$
E_{t_1}=\frac{t_1}{4}\sum_{\alpha,\beta\in\mathrm{occ}}\int \mathrm d1\,\mathrm d2\,\psi_\alpha^*(1)\psi_\beta^*(2)(1+x_1P_\sigma)\Bigl[\overleftarrow{\mathbf k}^{\,2}\delta(\mathbf r_1-\mathbf r_2)+\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}^{\,2}\Bigr]\Bigl[\psi_\alpha(1)\psi_\beta(2)-\psi_\beta(1)\psi_\alpha(2)\Bigr],
$$

The following expansions omit external coefficients, occupied-state sums, and integrals. The exchange minus sign remains in the energy expression above.

**Direct term $t_1$**

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\overleftarrow{\mathbf k}^{\,2}\delta(\mathbf r_1-\mathbf r_2)\psi_\alpha(1)\psi_\beta(2)=-\frac14\delta(\mathbf r_1-\mathbf r_2)\Bigl[(\Delta_1\psi_\alpha^*(1))\psi_\beta^*(2)+\psi_\alpha^*(1)(\Delta_2\psi_\beta^*(2))-2\nabla_1\psi_\alpha^*(1)\cdot\nabla_2\psi_\beta^*(2)\Bigr]\psi_\alpha(1)\psi_\beta(2).
$$

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}^{\,2}\psi_\alpha(1)\psi_\beta(2)=-\frac14\delta(\mathbf r_1-\mathbf r_2)\psi_\alpha^*(1)\psi_\beta^*(2)\Bigl[(\Delta_1\psi_\alpha(1))\psi_\beta(2)+\psi_\alpha(1)(\Delta_2\psi_\beta(2))-2\nabla_1\psi_\alpha(1)\cdot\nabla_2\psi_\beta(2)\Bigr].
$$

**Direct term $t_1x_1$**

$P_\sigma$ exchanges the spin indices of the two ket wave functions.

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\overleftarrow{\mathbf k}^{\,2}\delta(\mathbf r_1-\mathbf r_2)P_\sigma\psi_\alpha(1)\psi_\beta(2)=-\frac14\delta(\mathbf r_1-\mathbf r_2)\Bigl[(\Delta_1\psi_\alpha^*(1))\psi_\beta^*(2)+\psi_\alpha^*(1)(\Delta_2\psi_\beta^*(2))-2\nabla_1\psi_\alpha^*(1)\cdot\nabla_2\psi_\beta^*(2)\Bigr]\psi_\alpha(\mathbf r_1,s_2,q_1)\psi_\beta(\mathbf r_2,s_1,q_2).
$$

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}^{\,2}P_\sigma\psi_\alpha(1)\psi_\beta(2)=-\frac14\delta(\mathbf r_1-\mathbf r_2)\psi_\alpha^*(1)\psi_\beta^*(2)\Bigl[(\Delta_1\psi_\alpha(\mathbf r_1,s_2,q_1))\psi_\beta(\mathbf r_2,s_1,q_2)+\psi_\alpha(\mathbf r_1,s_2,q_1)(\Delta_2\psi_\beta(\mathbf r_2,s_1,q_2))-2\nabla_1\psi_\alpha(\mathbf r_1,s_2,q_1)\cdot\nabla_2\psi_\beta(\mathbf r_2,s_1,q_2)\Bigr].
$$

**Exchange term $t_1$**

The ket states are exchanged to $\psi_\beta(1)\psi_\alpha(2)$.

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\overleftarrow{\mathbf k}^{\,2}\delta(\mathbf r_1-\mathbf r_2)\psi_\beta(1)\psi_\alpha(2)=-\frac14\delta(\mathbf r_1-\mathbf r_2)\Bigl[(\Delta_1\psi_\alpha^*(1))\psi_\beta^*(2)+\psi_\alpha^*(1)(\Delta_2\psi_\beta^*(2))-2\nabla_1\psi_\alpha^*(1)\cdot\nabla_2\psi_\beta^*(2)\Bigr]\psi_\beta(1)\psi_\alpha(2).
$$

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}^{\,2}\psi_\beta(1)\psi_\alpha(2)=-\frac14\delta(\mathbf r_1-\mathbf r_2)\psi_\alpha^*(1)\psi_\beta^*(2)\Bigl[(\Delta_1\psi_\beta(1))\psi_\alpha(2)+\psi_\beta(1)(\Delta_2\psi_\alpha(2))-2\nabla_1\psi_\beta(1)\cdot\nabla_2\psi_\alpha(2)\Bigr].
$$

**Exchange term $t_1x_1$**

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\overleftarrow{\mathbf k}^{\,2}\delta(\mathbf r_1-\mathbf r_2)P_\sigma\psi_\beta(1)\psi_\alpha(2)=-\frac14\delta(\mathbf r_1-\mathbf r_2)\Bigl[(\Delta_1\psi_\alpha^*(1))\psi_\beta^*(2)+\psi_\alpha^*(1)(\Delta_2\psi_\beta^*(2))-2\nabla_1\psi_\alpha^*(1)\cdot\nabla_2\psi_\beta^*(2)\Bigr]\psi_\beta(\mathbf r_1,s_2,q_1)\psi_\alpha(\mathbf r_2,s_1,q_2).
$$

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}^{\,2}P_\sigma\psi_\beta(1)\psi_\alpha(2)=-\frac14\delta(\mathbf r_1-\mathbf r_2)\psi_\alpha^*(1)\psi_\beta^*(2)\Bigl[(\Delta_1\psi_\beta(\mathbf r_1,s_2,q_1))\psi_\alpha(\mathbf r_2,s_1,q_2)+\psi_\beta(\mathbf r_1,s_2,q_1)(\Delta_2\psi_\alpha(\mathbf r_2,s_1,q_2))-2\nabla_1\psi_\beta(\mathbf r_1,s_2,q_1)\cdot\nabla_2\psi_\alpha(\mathbf r_2,s_1,q_2)\Bigr].
$$

**Local densities and derivative contractions**

Using the nonlocal particle and spin densities from the first section, define

$$
\tau(q;\mathbf r)=\nabla\cdot\nabla'\rho(q;\mathbf r,\mathbf r'),\qquad \mathbf j(q;\mathbf r)=\frac{\nabla-\nabla'}{2i}\rho(q;\mathbf r,\mathbf r'),
$$

$$
T_\nu(q;\mathbf r)=\nabla\cdot\nabla's_\nu(q;\mathbf r,\mathbf r'),\qquad J_{\mu\nu}(q;\mathbf r)=\frac{\partial_\mu-\partial_\mu'}{2i}s_\nu(q;\mathbf r,\mathbf r').
$$

On each right-hand side, take the derivatives before setting $\mathbf r'=\mathbf r$; $\mu,\nu=x,y,z$. Total densities are sums over nucleon species.

The particle-density derivative contractions are

$$
\sum_{\alpha\in\mathrm{occ},q}\sum_s\nabla\psi_\alpha^*(\mathbf r,s,q)\cdot\nabla\psi_\alpha(\mathbf r,s,q)=\tau(q;\mathbf r),
$$

$$
\sum_{\alpha\in\mathrm{occ},q}\sum_s\psi_\alpha^*(\mathbf r,s,q)\nabla\psi_\alpha(\mathbf r,s,q)=\frac12\nabla\rho(q;\mathbf r)+i\mathbf j(q;\mathbf r),
$$

$$
\sum_{\alpha\in\mathrm{occ},q}\sum_s\Bigl[\psi_\alpha^*(\mathbf r,s,q)\Delta\psi_\alpha(\mathbf r,s,q)+(\Delta\psi_\alpha^*(\mathbf r,s,q))\psi_\alpha(\mathbf r,s,q)\Bigr]=\Delta\rho(q;\mathbf r)-2\tau(q;\mathbf r).
$$

The spin-density derivative contractions are

$$
\sum_{\alpha\in\mathrm{occ},q}\sum_{s,s'}\nabla\psi_\alpha^*(\mathbf r,s',q)\cdot\nabla\psi_\alpha(\mathbf r,s,q)(\sigma_\nu)_{s',s}=T_\nu(q;\mathbf r),
$$

$$
\sum_{\alpha\in\mathrm{occ},q}\sum_{s,s'}\psi_\alpha^*(\mathbf r,s',q)\partial_\mu\psi_\alpha(\mathbf r,s,q)(\sigma_\nu)_{s',s}=\frac12\partial_\mu s_\nu(q;\mathbf r)+iJ_{\mu\nu}(q;\mathbf r),
$$

$$
\sum_{\alpha\in\mathrm{occ},q}\sum_{s,s'}\Bigl[\psi_\alpha^*(\mathbf r,s',q)\Delta\psi_\alpha(\mathbf r,s,q)+(\Delta\psi_\alpha^*(\mathbf r,s',q))\psi_\alpha(\mathbf r,s,q)\Bigr](\sigma_\nu)_{s',s}=\Delta s_\nu(q;\mathbf r)-2T_\nu(q;\mathbf r).
$$

**Energy contributions**

Use these contractions, the spin-density matrix decomposition from the first section, and integration by parts:

$$
\int \mathrm d^3\mathbf r\,\rho(\mathbf r)\Delta\rho(\mathbf r)=-\int \mathrm d^3\mathbf r\,|\nabla\rho(\mathbf r)|^2.
$$

Boundary terms vanish. Below, all densities are evaluated at $\mathbf r$, and $|\nabla\mathbf s|^2=\sum_{\mu,\nu}(\partial_\mu s_\nu)^2$.

$$
E_{t_1}^{\mathrm{Hartree}}=\frac{t_1}{8}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_1)(\rho\tau-|\mathbf j|^2)+x_1\Bigl(\mathbf s\cdot\mathbf T-\sum_{\mu,\nu}J_{\mu\nu}^2\Bigr)\Bigr]+\frac{3t_1}{32}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_1)|\nabla\rho|^2+x_1|\nabla\mathbf s|^2\Bigr].
$$

$$
E_{t_1}^{\mathrm{Fock}}=-\frac{t_1}{8}\int \mathrm d^3\mathbf r\,\sum_q\Bigl[(1+2x_1)\bigl(\rho(q)\tau(q)-|\mathbf j(q)|^2\bigr)+\mathbf s(q)\cdot\mathbf T(q)-\sum_{\mu,\nu}J_{\mu\nu}(q)^2\Bigr]-\frac{3t_1}{32}\int \mathrm d^3\mathbf r\,\sum_q\Bigl[(1+2x_1)|\nabla\rho(q)|^2+|\nabla\mathbf s(q)|^2\Bigr].
$$

These expressions do not assume time-reversal symmetry.

$$
\boxed{E_{t_1}=E_{t_1}^{\mathrm{Hartree}}+E_{t_1}^{\mathrm{Fock}}.}
$$

Varying each local density while holding the remaining densities fixed gives

$$
\frac{\delta E_{t_1}}{\delta\rho(q;\mathbf r)}=\frac{t_1}{8}\Bigl[(2+x_1)\tau-(1+2x_1)\tau(q)\Bigr]-\frac{3t_1}{16}\Bigl[(2+x_1)\Delta\rho-(1+2x_1)\Delta\rho(q)\Bigr],
$$

$$
\frac{\delta E_{t_1}}{\delta\mathbf s(q;\mathbf r)}=\frac{t_1}{8}\Bigl[x_1\mathbf T-\mathbf T(q)\Bigr]-\frac{3t_1}{16}\Bigl[x_1\Delta\mathbf s-\Delta\mathbf s(q)\Bigr].
$$

$$
\frac{\delta E_{t_1}}{\delta\tau(q;\mathbf r)}=\frac{t_1}{8}\Bigl[(2+x_1)\rho-(1+2x_1)\rho(q)\Bigr],
$$

$$
\frac{\delta E_{t_1}}{\delta\mathbf j(q;\mathbf r)}=-\frac{t_1}{4}\Bigl[(2+x_1)\mathbf j-(1+2x_1)\mathbf j(q)\Bigr],
$$

$$
\frac{\delta E_{t_1}}{\delta\mathbf T(q;\mathbf r)}=\frac{t_1}{8}\Bigl[x_1\mathbf s-\mathbf s(q)\Bigr],
$$

$$
\frac{\delta E_{t_1}}{\delta J_{\mu\nu}(q;\mathbf r)}=-\frac{t_1}{4}\Bigl[x_1J_{\mu\nu}-J_{\mu\nu}(q)\Bigr].
$$

All densities are evaluated at $\mathbf r$. Gradient contributions are included in the $\rho$ and $\mathbf s$ variations by integration by parts.

For **time-reversal symmetry**,

$$
\mathbf s(q;\mathbf r)=\mathbf j(q;\mathbf r)=\mathbf T(q;\mathbf r)=0.
$$

The spin-current tensor $J_{\mu\nu}(q;\mathbf r)$ need not vanish. The energy reduces to

$$
\boxed{E_{t_1}=E_{t_1}^{\mathrm{Hartree}}+E_{t_1}^{\mathrm{Fock}}=\frac{t_1}{8}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_1)\rho\tau-(1+2x_1)\sum_q\rho(q)\tau(q)-x_1\sum_{\mu,\nu}J_{\mu\nu}^2+\sum_q\sum_{\mu,\nu}J_{\mu\nu}(q)^2\Bigr]+\frac{3t_1}{32}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_1)|\nabla\rho|^2-(1+2x_1)\sum_q|\nabla\rho(q)|^2\Bigr].}
$$

The corresponding variations are

$$
\frac{\delta E_{t_1}}{\delta\rho(q;\mathbf r)}=\frac{t_1}{8}\Bigl[(2+x_1)\tau-(1+2x_1)\tau(q)\Bigr]-\frac{3t_1}{16}\Bigl[(2+x_1)\Delta\rho-(1+2x_1)\Delta\rho(q)\Bigr],
$$

$$
\frac{\delta E_{t_1}}{\delta\tau(q;\mathbf r)}=\frac{t_1}{8}\Bigl[(2+x_1)\rho-(1+2x_1)\rho(q)\Bigr],
$$

$$
\frac{\delta E_{t_1}}{\delta J_{\mu\nu}(q;\mathbf r)}=-\frac{t_1}{4}\Bigl[x_1J_{\mu\nu}-J_{\mu\nu}(q)\Bigr],
$$

$$
\frac{\delta E_{t_1}}{\delta\mathbf s(q;\mathbf r)}=\frac{\delta E_{t_1}}{\delta\mathbf j(q;\mathbf r)}=\frac{\delta E_{t_1}}{\delta\mathbf T(q;\mathbf r)}=0.
$$

All densities are evaluated at $\mathbf r$.

### t2 term

The $t_2$ two-body interaction and its many-body term are

$$
\boxed{v_{t_2}(1,2)=t_2(1+x_2P_\sigma)\overleftarrow{\mathbf k}\cdot\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k},\qquad \hat V_{t_2}=\frac12\sum_{i\ne j}v_{t_2}(i,j).}
$$

Use the relative wave-number operators defined in the $t_1$ term.

The corresponding energy is

$$
E_{t_2}=\frac{t_2}{2}\sum_{\alpha,\beta\in\mathrm{occ}}\int \mathrm d1\,\mathrm d2\,\psi_\alpha^*(1)\psi_\beta^*(2)(1+x_2P_\sigma)\overleftarrow{\mathbf k}\cdot\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}\Bigl[\psi_\alpha(1)\psi_\beta(2)-\psi_\beta(1)\psi_\alpha(2)\Bigr].
$$

The following expansions omit external coefficients, occupied-state sums, and integrals. The exchange minus sign remains in the energy expression above.

**Direct term $t_2$**

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\overleftarrow{\mathbf k}\cdot\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}\psi_\alpha(1)\psi_\beta(2)=\frac14\delta(\mathbf r_1-\mathbf r_2)\Bigl[(\nabla_1\psi_\alpha^*(1))\psi_\beta^*(2)-\psi_\alpha^*(1)\nabla_2\psi_\beta^*(2)\Bigr]\cdot\Bigl[(\nabla_1\psi_\alpha(1))\psi_\beta(2)-\psi_\alpha(1)\nabla_2\psi_\beta(2)\Bigr].
$$

**Direct term $t_2x_2$**

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\overleftarrow{\mathbf k}\cdot\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}P_\sigma\psi_\alpha(1)\psi_\beta(2)=\frac14\delta(\mathbf r_1-\mathbf r_2)\Bigl[(\nabla_1\psi_\alpha^*(1))\psi_\beta^*(2)-\psi_\alpha^*(1)\nabla_2\psi_\beta^*(2)\Bigr]\cdot\Bigl[(\nabla_1\psi_\alpha(\mathbf r_1,s_2,q_1))\psi_\beta(\mathbf r_2,s_1,q_2)-\psi_\alpha(\mathbf r_1,s_2,q_1)\nabla_2\psi_\beta(\mathbf r_2,s_1,q_2)\Bigr].
$$

**Exchange term $t_2$**

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\overleftarrow{\mathbf k}\cdot\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}\psi_\beta(1)\psi_\alpha(2)=\frac14\delta(\mathbf r_1-\mathbf r_2)\Bigl[(\nabla_1\psi_\alpha^*(1))\psi_\beta^*(2)-\psi_\alpha^*(1)\nabla_2\psi_\beta^*(2)\Bigr]\cdot\Bigl[(\nabla_1\psi_\beta(1))\psi_\alpha(2)-\psi_\beta(1)\nabla_2\psi_\alpha(2)\Bigr].
$$

**Exchange term $t_2x_2$**

$$
\psi_\alpha^*(1)\psi_\beta^*(2)\overleftarrow{\mathbf k}\cdot\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}P_\sigma\psi_\beta(1)\psi_\alpha(2)=\frac14\delta(\mathbf r_1-\mathbf r_2)\Bigl[(\nabla_1\psi_\alpha^*(1))\psi_\beta^*(2)-\psi_\alpha^*(1)\nabla_2\psi_\beta^*(2)\Bigr]\cdot\Bigl[(\nabla_1\psi_\beta(\mathbf r_1,s_2,q_1))\psi_\alpha(\mathbf r_2,s_1,q_2)-\psi_\beta(\mathbf r_1,s_2,q_1)\nabla_2\psi_\alpha(\mathbf r_2,s_1,q_2)\Bigr].
$$

**Energy contributions**

Using the local densities and derivative contractions defined in the $t_1$ term, with vanishing boundary terms, gives

$$
E_{t_2}^{\mathrm{Hartree}}=\frac{t_2}{8}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_2)(\rho\tau-|\mathbf j|^2)+x_2\Bigl(\mathbf s\cdot\mathbf T-\sum_{\mu,\nu}J_{\mu\nu}^2\Bigr)\Bigr]-\frac{t_2}{32}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_2)|\nabla\rho|^2+x_2|\nabla\mathbf s|^2\Bigr].
$$

$$
E_{t_2}^{\mathrm{Fock}}=\frac{t_2}{8}\int \mathrm d^3\mathbf r\,\sum_q\Bigl[(1+2x_2)\bigl(\rho(q)\tau(q)-|\mathbf j(q)|^2\bigr)+\mathbf s(q)\cdot\mathbf T(q)-\sum_{\mu,\nu}J_{\mu\nu}(q)^2\Bigr]-\frac{t_2}{32}\int \mathrm d^3\mathbf r\,\sum_q\Bigl[(1+2x_2)|\nabla\rho(q)|^2+|\nabla\mathbf s(q)|^2\Bigr].
$$

These expressions do not assume time-reversal symmetry.

$$
\boxed{E_{t_2}=E_{t_2}^{\mathrm{Hartree}}+E_{t_2}^{\mathrm{Fock}}.}
$$

Varying each local density while holding the remaining densities fixed gives

$$
\frac{\delta E_{t_2}}{\delta\rho(q;\mathbf r)}=\frac{t_2}{8}\Bigl[(2+x_2)\tau+(1+2x_2)\tau(q)\Bigr]+\frac{t_2}{16}\Bigl[(2+x_2)\Delta\rho+(1+2x_2)\Delta\rho(q)\Bigr],
$$

$$
\frac{\delta E_{t_2}}{\delta\mathbf s(q;\mathbf r)}=\frac{t_2}{8}\Bigl[x_2\mathbf T+\mathbf T(q)\Bigr]+\frac{t_2}{16}\Bigl[x_2\Delta\mathbf s+\Delta\mathbf s(q)\Bigr].
$$

$$
\frac{\delta E_{t_2}}{\delta\tau(q;\mathbf r)}=\frac{t_2}{8}\Bigl[(2+x_2)\rho+(1+2x_2)\rho(q)\Bigr],
$$

$$
\frac{\delta E_{t_2}}{\delta\mathbf j(q;\mathbf r)}=-\frac{t_2}{4}\Bigl[(2+x_2)\mathbf j+(1+2x_2)\mathbf j(q)\Bigr],
$$

$$
\frac{\delta E_{t_2}}{\delta\mathbf T(q;\mathbf r)}=\frac{t_2}{8}\Bigl[x_2\mathbf s+\mathbf s(q)\Bigr],
$$

$$
\frac{\delta E_{t_2}}{\delta J_{\mu\nu}(q;\mathbf r)}=-\frac{t_2}{4}\Bigl[x_2J_{\mu\nu}+J_{\mu\nu}(q)\Bigr].
$$

All densities are evaluated at $\mathbf r$. Gradient contributions are included in the $\rho$ and $\mathbf s$ variations by integration by parts.

For **time-reversal symmetry**,

$$
\mathbf s(q;\mathbf r)=\mathbf j(q;\mathbf r)=\mathbf T(q;\mathbf r)=0.
$$

The spin-current tensor $J_{\mu\nu}(q;\mathbf r)$ need not vanish. The energy reduces to

$$
\boxed{E_{t_2}=E_{t_2}^{\mathrm{Hartree}}+E_{t_2}^{\mathrm{Fock}}=\frac{t_2}{8}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_2)\rho\tau+(1+2x_2)\sum_q\rho(q)\tau(q)-x_2\sum_{\mu,\nu}J_{\mu\nu}^2-\sum_q\sum_{\mu,\nu}J_{\mu\nu}(q)^2\Bigr]-\frac{t_2}{32}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_2)|\nabla\rho|^2+(1+2x_2)\sum_q|\nabla\rho(q)|^2\Bigr].}
$$

The corresponding variations are

$$
\frac{\delta E_{t_2}}{\delta\rho(q;\mathbf r)}=\frac{t_2}{8}\Bigl[(2+x_2)\tau+(1+2x_2)\tau(q)\Bigr]+\frac{t_2}{16}\Bigl[(2+x_2)\Delta\rho+(1+2x_2)\Delta\rho(q)\Bigr],
$$

$$
\frac{\delta E_{t_2}}{\delta\tau(q;\mathbf r)}=\frac{t_2}{8}\Bigl[(2+x_2)\rho+(1+2x_2)\rho(q)\Bigr],
$$

$$
\frac{\delta E_{t_2}}{\delta J_{\mu\nu}(q;\mathbf r)}=-\frac{t_2}{4}\Bigl[x_2J_{\mu\nu}+J_{\mu\nu}(q)\Bigr],
$$

$$
\frac{\delta E_{t_2}}{\delta\mathbf s(q;\mathbf r)}=\frac{\delta E_{t_2}}{\delta\mathbf j(q;\mathbf r)}=\frac{\delta E_{t_2}}{\delta\mathbf T(q;\mathbf r)}=0.
$$

All densities are evaluated at $\mathbf r$.

### t3 term

The density-dependent two-body interaction and its many-body term are

$$
\boxed{v_{t_3}(1,2;\rho)=\frac{t_3}{6}(1+x_3P_\sigma)\rho^\gamma\!\left(\frac{\mathbf r_1+\mathbf r_2}{2}\right)\delta(\mathbf r_1-\mathbf r_2),\qquad \hat V_{t_3}[\rho]=\frac12\sum_{i\ne j}v_{t_3}(i,j;\rho).}
$$

Here $t_3,x_3,\gamma$ are the strength, spin-exchange parameter, and density exponent, respectively, with $\rho(\mathbf r)=\sum_q\rho(q;\mathbf r)$. Use the definition of $P_\sigma$ from the $t_0$ term.

The corresponding energy is

$$
E_{t_3}=\frac{t_3}{12}\sum_{\alpha,\beta\in\mathrm{occ}}\int \mathrm d1\,\mathrm d2\,\rho^\gamma\!\left(\frac{\mathbf r_1+\mathbf r_2}{2}\right)\delta(\mathbf r_1-\mathbf r_2)\psi_\alpha^*(1)\psi_\beta^*(2)(1+x_3P_\sigma)\Bigl[\psi_\alpha(1)\psi_\beta(2)-\psi_\beta(1)\psi_\alpha(2)\Bigr].
$$

Without proton–neutron mixing, the four contractions are the same as in the $t_0$ term. The following expressions omit $t_3/12$, the density factor $\rho^\gamma$, the exchange minus sign, and $x_3$.

**Direct term $t_3$**

$$
\delta(\mathbf r_1-\mathbf r_2)\sum_{\alpha,\beta\in\mathrm{occ}}\psi_\alpha^*(1)\psi_\beta^*(2)\psi_\alpha(1)\psi_\beta(2)=\delta(\mathbf r_1-\mathbf r_2)\rho_{s_1,s_1}(q_1;\mathbf r_1)\rho_{s_2,s_2}(q_2;\mathbf r_2).
$$

**Direct term $t_3x_3$**

$$
\delta(\mathbf r_1-\mathbf r_2)\sum_{\alpha,\beta\in\mathrm{occ}}\psi_\alpha^*(1)\psi_\beta^*(2)\psi_\alpha(\mathbf r_1,s_2,q_1)\psi_\beta(\mathbf r_2,s_1,q_2)=\delta(\mathbf r_1-\mathbf r_2)\rho_{s_2,s_1}(q_1;\mathbf r_1)\rho_{s_1,s_2}(q_2;\mathbf r_2).
$$

**Exchange term $t_3$**

$$
\delta(\mathbf r_1-\mathbf r_2)\sum_{\alpha,\beta\in\mathrm{occ}}\psi_\alpha^*(1)\psi_\beta^*(2)\psi_\beta(1)\psi_\alpha(2)=\delta(\mathbf r_1-\mathbf r_2)\delta_{q_1,q_2}\rho_{s_2,s_1}(q_1;\mathbf r_2,\mathbf r_1)\rho_{s_1,s_2}(q_1;\mathbf r_1,\mathbf r_2).
$$

**Exchange term $t_3x_3$**

$$
\delta(\mathbf r_1-\mathbf r_2)\sum_{\alpha,\beta\in\mathrm{occ}}\psi_\alpha^*(1)\psi_\beta^*(2)\psi_\beta(\mathbf r_1,s_2,q_1)\psi_\alpha(\mathbf r_2,s_1,q_2)=\delta(\mathbf r_1-\mathbf r_2)\delta_{q_1,q_2}\rho_{s_1,s_1}(q_1;\mathbf r_2,\mathbf r_1)\rho_{s_2,s_2}(q_1;\mathbf r_1,\mathbf r_2).
$$

**Energy contributions**

Integrate using $\delta(\mathbf r_1-\mathbf r_2)$ and sum over spin and nucleon species. Below, all densities are evaluated at $\mathbf r$, with $\mathbf s=\sum_q\mathbf s(q)$.

$$
E_{t_3}^{\mathrm{Hartree}}=\frac{t_3}{12}\int \mathrm d^3\mathbf r\,\rho^\gamma\sum_{q,q'}\Bigl[\rho(q)\rho(q')+\frac{x_3}{2}\bigl(\rho(q)\rho(q')+\mathbf s(q)\cdot\mathbf s(q')\bigr)\Bigr]=\frac{t_3}{24}\int \mathrm d^3\mathbf r\,\rho^\gamma\Bigl[(2+x_3)\rho^2+x_3|\mathbf s|^2\Bigr].
$$

$$
E_{t_3}^{\mathrm{Fock}}=-\frac{t_3}{12}\int \mathrm d^3\mathbf r\,\rho^\gamma\sum_q\Bigl[\frac12\bigl(\rho(q)^2+|\mathbf s(q)|^2\bigr)+x_3\rho(q)^2\Bigr]=-\frac{t_3}{24}\int \mathrm d^3\mathbf r\,\rho^\gamma\sum_q\Bigl[(1+2x_3)\rho(q)^2+|\mathbf s(q)|^2\Bigr].
$$

The total energy is

$$
\boxed{E_{t_3}=E_{t_3}^{\mathrm{Hartree}}+E_{t_3}^{\mathrm{Fock}}=\frac{t_3}{24}\int \mathrm d^3\mathbf r\,\rho^\gamma\Bigl[(2+x_3)\rho^2-(1+2x_3)\sum_q\rho(q)^2+x_3|\mathbf s|^2-\sum_q|\mathbf s(q)|^2\Bigr].}
$$

Varying the particle and spin densities separately while holding the remaining densities fixed gives

$$
\frac{\delta E_{t_3}}{\delta\rho(q;\mathbf r)}=\frac{t_3}{12}\rho^\gamma\Bigl[(2+x_3)\rho-(1+2x_3)\rho(q)\Bigr]+\frac{\gamma t_3}{24}\rho^{\gamma-1}\Bigl[(2+x_3)\rho^2-(1+2x_3)\sum_{q'}\rho(q')^2+x_3|\mathbf s|^2-\sum_{q'}|\mathbf s(q')|^2\Bigr],
$$

$$
\frac{\delta E_{t_3}}{\delta\mathbf s(q;\mathbf r)}=\frac{t_3}{12}\rho^\gamma\Bigl[x_3\mathbf s-\mathbf s(q)\Bigr].
$$

The second term in the particle-density variation comes from differentiating $\rho^\gamma$ and is the rearrangement term.

For **time-reversal symmetry**,

$$
\mathbf s(q;\mathbf r)=0.
$$

The energy reduces to

$$
\boxed{E_{t_3}=\frac{t_3}{24}\int \mathrm d^3\mathbf r\,\rho^\gamma\Bigl[(2+x_3)\rho^2-(1+2x_3)\sum_q\rho(q)^2\Bigr].}
$$

The corresponding variations are

$$
\frac{\delta E_{t_3}}{\delta\rho(q;\mathbf r)}=\frac{t_3}{12}\rho^\gamma\Bigl[(2+x_3)\rho-(1+2x_3)\rho(q)\Bigr]+\frac{\gamma t_3}{24}\rho^{\gamma-1}\Bigl[(2+x_3)\rho^2-(1+2x_3)\sum_{q'}\rho(q')^2\Bigr],
$$

$$
\frac{\delta E_{t_3}}{\delta\mathbf s(q;\mathbf r)}=0.
$$

### SO term

The spin-orbit two-body interaction and its many-body term are

$$
\boxed{v_{\mathrm{SO}}(1,2)=iW_0(\boldsymbol\sigma_1+\boldsymbol\sigma_2)\cdot\Bigl[\overleftarrow{\mathbf k}\times\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}\Bigr],\qquad \hat V_{\mathrm{SO}}=\frac12\sum_{i\ne j}v_{\mathrm{SO}}(i,j).}
$$

Here $W_0$ is the spin-orbit strength. Use the left- and right-acting wave-number operators defined in the $t_1$ term. The two-nucleon spin operators are

$$
\sigma_{1\nu}=\sigma_\nu\otimes I_2,\qquad \sigma_{2\nu}=I_2\otimes\sigma_\nu.
$$

The corresponding energy is

$$
E_{\mathrm{SO}}=\frac{iW_0}{2}\sum_{\alpha,\beta\in\mathrm{occ}}\int \mathrm d1\,\mathrm d2\,\psi_\alpha^*(1)\psi_\beta^*(2)(\boldsymbol\sigma_1+\boldsymbol\sigma_2)\cdot\Bigl[\overleftarrow{\mathbf k}\times\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}\Bigr]\Bigl[\psi_\alpha(1)\psi_\beta(2)-\psi_\beta(1)\psi_\alpha(2)\Bigr].
$$

The following expansions omit $iW_0/2$, occupied-state sums, and integrals. The exchange minus sign remains in the energy expression above.

**Direct term $W_0$**

$$
\psi_\alpha^*(1)\psi_\beta^*(2)(\boldsymbol\sigma_1+\boldsymbol\sigma_2)\cdot\Bigl[\overleftarrow{\mathbf k}\times\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}\Bigr]\psi_\alpha(1)\psi_\beta(2)=\frac14\delta(\mathbf r_1-\mathbf r_2)\sum_{\nu,\mu,\lambda}\epsilon_{\nu\mu\lambda}\Bigl[(\partial_{1\mu}\psi_\alpha^*(1))\psi_\beta^*(2)-\psi_\alpha^*(1)\partial_{2\mu}\psi_\beta^*(2)\Bigr](\sigma_{1\nu}+\sigma_{2\nu})\Bigl[(\partial_{1\lambda}\psi_\alpha(1))\psi_\beta(2)-\psi_\alpha(1)\partial_{2\lambda}\psi_\beta(2)\Bigr].
$$

**Exchange term $W_0$**

$$
\psi_\alpha^*(1)\psi_\beta^*(2)(\boldsymbol\sigma_1+\boldsymbol\sigma_2)\cdot\Bigl[\overleftarrow{\mathbf k}\times\delta(\mathbf r_1-\mathbf r_2)\overrightarrow{\mathbf k}\Bigr]\psi_\beta(1)\psi_\alpha(2)=\frac14\delta(\mathbf r_1-\mathbf r_2)\sum_{\nu,\mu,\lambda}\epsilon_{\nu\mu\lambda}\Bigl[(\partial_{1\mu}\psi_\alpha^*(1))\psi_\beta^*(2)-\psi_\alpha^*(1)\partial_{2\mu}\psi_\beta^*(2)\Bigr](\sigma_{1\nu}+\sigma_{2\nu})\Bigl[(\partial_{1\lambda}\psi_\beta(1))\psi_\alpha(2)-\psi_\beta(1)\partial_{2\lambda}\psi_\alpha(2)\Bigr].
$$

Here $\epsilon_{\nu\mu\lambda}$ is the Levi-Civita symbol, with $\epsilon_{xyz}=1$; Cartesian indices run over $x,y,z$. The spin matrices act on the spin indices of the wave functions to their right.

**Spin-current vector**

Using the spin-current tensor from the $t_1$ term, define

$$
J_\kappa(q;\mathbf r)=\sum_{\mu,\nu}\epsilon_{\kappa\mu\nu}J_{\mu\nu}(q;\mathbf r),\qquad \mathbf J(\mathbf r)=\sum_q\mathbf J(q;\mathbf r).
$$

**Energy contributions**

Integrate using $\delta(\mathbf r_1-\mathbf r_2)$, perform the spin contractions, and integrate by parts with vanishing boundary terms. Below, all densities are evaluated at $\mathbf r$.

$$
E_{\mathrm{SO}}^{\mathrm{Hartree}}=-\frac{W_0}{2}\int \mathrm d^3\mathbf r\,\Bigl[\rho\,\nabla\cdot\mathbf J+\mathbf s\cdot(\nabla\times\mathbf j)\Bigr],
$$

$$
E_{\mathrm{SO}}^{\mathrm{Fock}}=-\frac{W_0}{2}\int \mathrm d^3\mathbf r\,\sum_q\Bigl[\rho(q)\,\nabla\cdot\mathbf J(q)+\mathbf s(q)\cdot(\nabla\times\mathbf j(q))\Bigr].
$$

The total energy is

$$
\boxed{E_{\mathrm{SO}}=E_{\mathrm{SO}}^{\mathrm{Hartree}}+E_{\mathrm{SO}}^{\mathrm{Fock}}.}
$$

Varying each density while holding the remaining densities fixed gives

$$
\frac{\delta E_{\mathrm{SO}}}{\delta\rho(q;\mathbf r)}=-\frac{W_0}{2}\Bigl[\nabla\cdot\mathbf J+\nabla\cdot\mathbf J(q)\Bigr],
$$

$$
\frac{\delta E_{\mathrm{SO}}}{\delta\mathbf s(q;\mathbf r)}=-\frac{W_0}{2}\Bigl[\nabla\times\mathbf j+\nabla\times\mathbf j(q)\Bigr],
$$

$$
\frac{\delta E_{\mathrm{SO}}}{\delta\mathbf j(q;\mathbf r)}=-\frac{W_0}{2}\Bigl[\nabla\times\mathbf s+\nabla\times\mathbf s(q)\Bigr],
$$

$$
\frac{\delta E_{\mathrm{SO}}}{\delta J_{\mu\nu}(q;\mathbf r)}=\frac{W_0}{2}\sum_\kappa\epsilon_{\kappa\mu\nu}\Bigl[\partial_\kappa\rho+\partial_\kappa\rho(q)\Bigr].
$$

For **time-reversal symmetry**,

$$
\mathbf s(q;\mathbf r)=\mathbf j(q;\mathbf r)=0.
$$

The spin-current tensor need not vanish. The energy reduces to

$$
\boxed{E_{\mathrm{SO}}=-\frac{W_0}{2}\int \mathrm d^3\mathbf r\,\Bigl[\rho\,\nabla\cdot\mathbf J+\sum_q\rho(q)\,\nabla\cdot\mathbf J(q)\Bigr].}
$$

The corresponding variations are

$$
\frac{\delta E_{\mathrm{SO}}}{\delta\rho(q;\mathbf r)}=-\frac{W_0}{2}\Bigl[\nabla\cdot\mathbf J+\nabla\cdot\mathbf J(q)\Bigr],
$$

$$
\frac{\delta E_{\mathrm{SO}}}{\delta J_{\mu\nu}(q;\mathbf r)}=\frac{W_0}{2}\sum_\kappa\epsilon_{\kappa\mu\nu}\Bigl[\partial_\kappa\rho+\partial_\kappa\rho(q)\Bigr],
$$

$$
\frac{\delta E_{\mathrm{SO}}}{\delta\mathbf s(q;\mathbf r)}=\frac{\delta E_{\mathrm{SO}}}{\delta\mathbf j(q;\mathbf r)}=0.
$$

The interaction and energy coefficients follow the [standard Skyrme convention](https://users.fuw.edu.pl/~dobaczew/landau19w/node18.html).

### Coulomb term

The Coulomb two-body interaction and its many-body term are

$$
\boxed{v_{\mathrm C}(1,2)=\frac{e^2}{|\mathbf r_1-\mathbf r_2|}\delta_{q_1,p}\delta_{q_2,p},\qquad \hat V_{\mathrm C}=\frac12\sum_{i\ne j}v_{\mathrm C}(i,j).}
$$

Here $e^2$ includes the Coulomb constant $1/(4\pi\epsilon_0)$; $\delta_{q,p}$ selects protons.

The corresponding energy is

$$
E_{\mathrm C}=\frac{e^2}{2}\sum_{\alpha,\beta\in\mathrm{occ}}\int \mathrm d1\,\mathrm d2\,\frac{\delta_{q_1,p}\delta_{q_2,p}}{|\mathbf r_1-\mathbf r_2|}\psi_\alpha^*(1)\psi_\beta^*(2)\Bigl[\psi_\alpha(1)\psi_\beta(2)-\psi_\beta(1)\psi_\alpha(2)\Bigr].
$$

The following expressions omit $e^2/2$, the Coulomb kernel, the proton selectors, and the exchange minus sign.

**Direct term**

$$
\sum_{\alpha,\beta\in\mathrm{occ},p}\psi_\alpha^*(\mathbf r_1,s_1,p)\psi_\beta^*(\mathbf r_2,s_2,p)\psi_\alpha(\mathbf r_1,s_1,p)\psi_\beta(\mathbf r_2,s_2,p)=\rho_{s_1,s_1}(p;\mathbf r_1)\rho_{s_2,s_2}(p;\mathbf r_2).
$$

**Exchange term**

$$
\sum_{\alpha,\beta\in\mathrm{occ},p}\psi_\alpha^*(\mathbf r_1,s_1,p)\psi_\beta^*(\mathbf r_2,s_2,p)\psi_\beta(\mathbf r_1,s_1,p)\psi_\alpha(\mathbf r_2,s_2,p)=\rho_{s_2,s_1}(p;\mathbf r_2,\mathbf r_1)\rho_{s_1,s_2}(p;\mathbf r_1,\mathbf r_2).
$$

**Energy contributions**

Summing over spin gives

$$
E_{\mathrm C}^{\mathrm{Hartree}}=\frac{e^2}{2}\int \mathrm d^3\mathbf r\,\mathrm d^3\mathbf r'\,\frac{\rho(p;\mathbf r)\rho(p;\mathbf r')}{|\mathbf r-\mathbf r'|},
$$

$$
E_{\mathrm C}^{\mathrm{Fock}}=-\frac{e^2}{2}\int \mathrm d^3\mathbf r\,\mathrm d^3\mathbf r'\,\frac{\sum_{s,s'}\rho_{s,s'}(p;\mathbf r,\mathbf r')\rho_{s',s}(p;\mathbf r',\mathbf r)}{|\mathbf r-\mathbf r'|}.
$$

The total energy is

$$
\boxed{E_{\mathrm C}=E_{\mathrm C}^{\mathrm{Hartree}}+E_{\mathrm C}^{\mathrm{Fock}}.}
$$

The particle-density variation of the direct term is

$$
\frac{\delta E_{\mathrm C}^{\mathrm{Hartree}}}{\delta\rho(q;\mathbf r)}=\delta_{q,p}\,e^2\int \mathrm d^3\mathbf r'\,\frac{\rho(p;\mathbf r')}{|\mathbf r-\mathbf r'|}.
$$

The exact exchange term depends on the nonlocal density matrix; its variation generates a nonlocal exchange operator.

**Time-reversal symmetry and the Slater approximation**

For **time-reversal symmetry**,

$$
\mathbf s(p;\mathbf r)=0.
$$

Applying the Slater local-density approximation to the exchange term gives

$$
E_{\mathrm C}^{\mathrm{Fock}}\approx E_{\mathrm C}^{\mathrm{Slater}}=-\frac{3e^2}{4}\left(\frac{3}{\pi}\right)^{1/3}\int \mathrm d^3\mathbf r\,\rho(p;\mathbf r)^{4/3}.
$$

This is an additional approximation to exchange, not a consequence of time-reversal symmetry alone. [Slater approximation reference](https://arxiv.org/abs/1210.3162)

The total energy becomes

$$
\boxed{E_{\mathrm C}\approx E_{\mathrm C}^{\mathrm{Hartree}}+E_{\mathrm C}^{\mathrm{Slater}}=\frac{e^2}{2}\int \mathrm d^3\mathbf r\,\mathrm d^3\mathbf r'\,\frac{\rho(p;\mathbf r)\rho(p;\mathbf r')}{|\mathbf r-\mathbf r'|}-\frac{3e^2}{4}\left(\frac{3}{\pi}\right)^{1/3}\int \mathrm d^3\mathbf r\,\rho(p;\mathbf r)^{4/3}.}
$$

The corresponding particle-density variation is

$$
\frac{\delta\bigl(E_{\mathrm C}^{\mathrm{Hartree}}+E_{\mathrm C}^{\mathrm{Slater}}\bigr)}{\delta\rho(q;\mathbf r)}=\delta_{q,p}\Biggl[e^2\int \mathrm d^3\mathbf r'\,\frac{\rho(p;\mathbf r')}{|\mathbf r-\mathbf r'|}-e^2\left(\frac{3}{\pi}\right)^{1/3}\rho(p;\mathbf r)^{1/3}\Biggr].
$$

## Time-reversal symmetry

### Time reversal and local densities

The time-reversal operator for a spin-$1/2$ nucleon is

$$
\Theta=i\sigma_yK,\qquad \Theta^2=-1,
$$

where $K$ denotes complex conjugation. Time reversal does not change the nucleon species $q$. In a time-reversal-invariant Slater determinant, occupied states occur in pairs:

$$
\psi_{\bar\alpha}(\mathbf r,s,q)=\sum_{s'}(i\sigma_y)_{s,s'}\psi_\alpha^*(\mathbf r,s',q).
$$

The local densities have the following time-reversal properties:

| Time-even | Time-odd |
|---|---|
| $\rho,\tau,J_{\mu\nu}$ | $\mathbf s,\mathbf j,\mathbf T$ |

Time-odd densities cancel between paired occupied states, giving

$$
\boxed{\mathbf s(q;\mathbf r)=\mathbf j(q;\mathbf r)=\mathbf T(q;\mathbf r)=0.}
$$

Time-even densities may remain nonzero. In particular, time-reversal symmetry does not require $J_{\mu\nu}=0$.

### Energy functional

For time-reversal symmetry, the Skyrme energy depends on $\rho,\tau,J_{\mu\nu}$. With the Slater approximation for Coulomb exchange, the total energy is

$$
E=E_{\mathrm{kin}}-E_{\mathrm{cm}}+E_{t_0}+E_{t_1}+E_{t_2}+E_{t_3}+E_{\mathrm{SO}}+E_{\mathrm C}^{\mathrm{Hartree}}+E_{\mathrm C}^{\mathrm{Slater}}.
$$

Here

$$
E_{\mathrm{kin}}=\frac{\hbar^2}{2m}\int \mathrm d^3\mathbf r\,\tau(\mathbf r).
$$

Use the time-reversal-symmetric expressions given above for each interaction term. The Skyrme terms retain the full spin-current tensor; the SO term uses only

$$
J_\kappa(q;\mathbf r)=\sum_{\mu,\nu}\epsilon_{\kappa\mu\nu}J_{\mu\nu}(q;\mathbf r).
$$

### Remaining density variations

First vary the general functional, then impose time-reversal symmetry. Without the center-of-mass correction, the three time-odd density variations vanish:

$$
\frac{\delta(E+E_{\mathrm{cm}})}{\delta\mathbf s(q;\mathbf r)}=\frac{\delta(E+E_{\mathrm{cm}})}{\delta\mathbf j(q;\mathbf r)}=\frac{\delta(E+E_{\mathrm{cm}})}{\delta\mathbf T(q;\mathbf r)}=0.
$$

Below, the variations with respect to $\rho,\tau,J_{\mu\nu}$ are listed in the order t0, t1, t2, t3, SO, and Coulomb. All densities are evaluated at $\mathbf r$.

**t0 term**

$$
E_{t_0}=\frac{t_0}{4}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_0)\rho(\mathbf r)^2-(1+2x_0)\sum_q\rho(q;\mathbf r)^2\Bigr].
$$

$$
\frac{\delta E_{t_0}}{\delta\rho(q;\mathbf r)}=\frac{t_0}{2}\Bigl[(2+x_0)\rho(\mathbf r)-(1+2x_0)\rho(q;\mathbf r)\Bigr],
$$

$$
\frac{\delta E_{t_0}}{\delta\tau(q;\mathbf r)}=\frac{\delta E_{t_0}}{\delta J_{\mu\nu}(q;\mathbf r)}=0.
$$

**t1 term**

$$
E_{t_1}=\frac{t_1}{8}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_1)\rho(\mathbf r)\tau(\mathbf r)-(1+2x_1)\sum_q\rho(q;\mathbf r)\tau(q;\mathbf r)-x_1\sum_{\mu,\nu}J_{\mu\nu}(\mathbf r)^2+\sum_q\sum_{\mu,\nu}J_{\mu\nu}(q;\mathbf r)^2\Bigr]+\frac{3t_1}{32}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_1)|\nabla\rho(\mathbf r)|^2-(1+2x_1)\sum_q|\nabla\rho(q;\mathbf r)|^2\Bigr].
$$

$$
\frac{\delta E_{t_1}}{\delta\rho(q;\mathbf r)}=\frac{t_1}{8}\Bigl[(2+x_1)\tau(\mathbf r)-(1+2x_1)\tau(q;\mathbf r)\Bigr]-\frac{3t_1}{16}\Bigl[(2+x_1)\Delta\rho(\mathbf r)-(1+2x_1)\Delta\rho(q;\mathbf r)\Bigr],
$$

$$
\frac{\delta E_{t_1}}{\delta\tau(q;\mathbf r)}=\frac{t_1}{8}\Bigl[(2+x_1)\rho(\mathbf r)-(1+2x_1)\rho(q;\mathbf r)\Bigr],
$$

$$
\frac{\delta E_{t_1}}{\delta J_{\mu\nu}(q;\mathbf r)}=-\frac{t_1}{4}\Bigl[x_1J_{\mu\nu}(\mathbf r)-J_{\mu\nu}(q;\mathbf r)\Bigr].
$$

**t2 term**

$$
E_{t_2}=\frac{t_2}{8}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_2)\rho(\mathbf r)\tau(\mathbf r)+(1+2x_2)\sum_q\rho(q;\mathbf r)\tau(q;\mathbf r)-x_2\sum_{\mu,\nu}J_{\mu\nu}(\mathbf r)^2-\sum_q\sum_{\mu,\nu}J_{\mu\nu}(q;\mathbf r)^2\Bigr]-\frac{t_2}{32}\int \mathrm d^3\mathbf r\,\Bigl[(2+x_2)|\nabla\rho(\mathbf r)|^2+(1+2x_2)\sum_q|\nabla\rho(q;\mathbf r)|^2\Bigr].
$$

$$
\frac{\delta E_{t_2}}{\delta\rho(q;\mathbf r)}=\frac{t_2}{8}\Bigl[(2+x_2)\tau(\mathbf r)+(1+2x_2)\tau(q;\mathbf r)\Bigr]+\frac{t_2}{16}\Bigl[(2+x_2)\Delta\rho(\mathbf r)+(1+2x_2)\Delta\rho(q;\mathbf r)\Bigr],
$$

$$
\frac{\delta E_{t_2}}{\delta\tau(q;\mathbf r)}=\frac{t_2}{8}\Bigl[(2+x_2)\rho(\mathbf r)+(1+2x_2)\rho(q;\mathbf r)\Bigr],
$$

$$
\frac{\delta E_{t_2}}{\delta J_{\mu\nu}(q;\mathbf r)}=-\frac{t_2}{4}\Bigl[x_2J_{\mu\nu}(\mathbf r)+J_{\mu\nu}(q;\mathbf r)\Bigr].
$$

**t3 term**

$$
E_{t_3}=\frac{t_3}{24}\int \mathrm d^3\mathbf r\,\rho(\mathbf r)^\gamma\Bigl[(2+x_3)\rho(\mathbf r)^2-(1+2x_3)\sum_q\rho(q;\mathbf r)^2\Bigr].
$$

$$
\frac{\delta E_{t_3}}{\delta\rho(q;\mathbf r)}=\frac{t_3}{12}\rho(\mathbf r)^\gamma\Bigl[(2+x_3)\rho(\mathbf r)-(1+2x_3)\rho(q;\mathbf r)\Bigr]+\frac{\gamma t_3}{24}\rho(\mathbf r)^{\gamma-1}\Bigl[(2+x_3)\rho(\mathbf r)^2-(1+2x_3)\sum_{q'}\rho(q';\mathbf r)^2\Bigr],
$$

$$
\frac{\delta E_{t_3}}{\delta\tau(q;\mathbf r)}=\frac{\delta E_{t_3}}{\delta J_{\mu\nu}(q;\mathbf r)}=0.
$$

**SO term**

$$
E_{\mathrm{SO}}=-\frac{W_0}{2}\int \mathrm d^3\mathbf r\,\Bigl[\rho(\mathbf r)\nabla\cdot\mathbf J(\mathbf r)+\sum_q\rho(q;\mathbf r)\nabla\cdot\mathbf J(q;\mathbf r)\Bigr].
$$

$$
\frac{\delta E_{\mathrm{SO}}}{\delta\rho(q;\mathbf r)}=-\frac{W_0}{2}\Bigl[\nabla\cdot\mathbf J(\mathbf r)+\nabla\cdot\mathbf J(q;\mathbf r)\Bigr],
$$

$$
\frac{\delta E_{\mathrm{SO}}}{\delta\tau(q;\mathbf r)}=0,
$$

$$
\frac{\delta E_{\mathrm{SO}}}{\delta J_{\mu\nu}(q;\mathbf r)}=\frac{W_0}{2}\sum_\kappa\epsilon_{\kappa\mu\nu}\Bigl[\partial_\kappa\rho(\mathbf r)+\partial_\kappa\rho(q;\mathbf r)\Bigr].
$$

**Coulomb term**

With the Slater approximation for exchange,

$$
E_{\mathrm C}\approx E_{\mathrm C}^{\mathrm{Hartree}}+E_{\mathrm C}^{\mathrm{Slater}}=\frac{e^2}{2}\int \mathrm d^3\mathbf r\,\mathrm d^3\mathbf r'\,\frac{\rho(p;\mathbf r)\rho(p;\mathbf r')}{|\mathbf r-\mathbf r'|}-\frac{3e^2}{4}\left(\frac{3}{\pi}\right)^{1/3}\int \mathrm d^3\mathbf r\,\rho(p;\mathbf r)^{4/3}.
$$

$$
\frac{\delta(E_{\mathrm C}^{\mathrm{Hartree}}+E_{\mathrm C}^{\mathrm{Slater}})}{\delta\rho(q;\mathbf r)}=\delta_{q,p}\Biggl[e^2\int \mathrm d^3\mathbf r'\,\frac{\rho(p;\mathbf r')}{|\mathbf r-\mathbf r'|}-e^2\left(\frac{3}{\pi}\right)^{1/3}\rho(p;\mathbf r)^{1/3}\Biggr],
$$

$$
\frac{\delta(E_{\mathrm C}^{\mathrm{Hartree}}+E_{\mathrm C}^{\mathrm{Slater}})}{\delta\tau(q;\mathbf r)}=\frac{\delta(E_{\mathrm C}^{\mathrm{Hartree}}+E_{\mathrm C}^{\mathrm{Slater}})}{\delta J_{\mu\nu}(q;\mathbf r)}=0.
$$

Sum the contributions to obtain the total variations.

**Kinetic term**

$$
E_{\mathrm{kin}}=\frac{\hbar^2}{2m}\int \mathrm d^3\mathbf r\,\tau(\mathbf r).
$$

$$
\frac{\delta E_{\mathrm{kin}}}{\delta\rho(q;\mathbf r)}=0,\qquad \frac{\delta E_{\mathrm{kin}}}{\delta\tau(q;\mathbf r)}=\frac{\hbar^2}{2m},\qquad \frac{\delta E_{\mathrm{kin}}}{\delta J_{\mu\nu}(q;\mathbf r)}=0.
$$

## UNEDF energy density functional

The energy density functional is specified directly through density couplings. The general form below follows the code's EDF parameterization; UNEDF1 selects a particular set of these couplings.

### Densities and coupling constants

For each density $X$, define the isoscalar and isovector combinations

$$
X_0(\mathbf r)=X_n(\mathbf r)+X_p(\mathbf r),\qquad X_1(\mathbf r)=X_n(\mathbf r)-X_p(\mathbf r),\qquad X\in\{\rho,\mathbf s,\tau,\mathbf T,\mathbf j,J_{\mu\nu}\}.
$$

Here $t=0,1$ labels the isospin channel. The spin-current vector is

$$
J_{t,\kappa}(\mathbf r)=\sum_{\mu,\nu}\epsilon_{\kappa\mu\nu}J_{t,\mu\nu}(\mathbf r).
$$

Coupling superscripts identify the associated density products. The density-dependent couplings are

$$
C_t^{\rho\rho}[\rho_0(\mathbf r)]=C_{t,0}^{\rho\rho}+C_{t,D}^{\rho\rho}\rho_0^\alpha(\mathbf r),\qquad C_t^{ss}[\rho_0(\mathbf r)]=C_{t,0}^{ss}+C_{t,D}^{ss}\rho_0^\alpha(\mathbf r).
$$

Here $\alpha$ is the density exponent, corresponding to $\gamma$ in the preceding Skyrme interaction.

### Time-even energy density

All densities below are evaluated at $\mathbf r$. The time-even interaction energy density contains eight density combinations:

$$
\boxed{\mathcal H_{\mathrm{even}}(\mathbf r)=\sum_{t=0,1}\Bigl[C_t^{\rho\rho}[\rho_0(\mathbf r)]\rho_t^2(\mathbf r)+C_t^{\rho\tau}\rho_t(\mathbf r)\tau_t(\mathbf r)+C_t^{\rho\Delta\rho}\rho_t(\mathbf r)\Delta\rho_t(\mathbf r)+C_t^{\rho\nabla J}\rho_t(\mathbf r)\nabla\cdot\mathbf J_t(\mathbf r)+C_t^{JJ}\sum_{\mu,\nu}J_{t,\mu\nu}^2(\mathbf r)+C_t^{JJ^{\mathsf T}}\sum_{\mu,\nu}J_{t,\mu\nu}(\mathbf r)J_{t,\nu\mu}(\mathbf r)+C_t^{\nabla\rho\nabla\rho}|\nabla\rho_t(\mathbf r)|^2+C_t^{J\nabla\rho}\mathbf J_t(\mathbf r)\cdot\nabla\rho_t(\mathbf r)\Bigr].}
$$

The two tensor contractions are distinct. The vector $\mathbf J_t$ contains only the antisymmetric part of $J_{t,\mu\nu}$.

For constant gradient couplings and vanishing boundary terms,

$$
\int \mathrm d^3\mathbf r\,\rho_t(\mathbf r)\Delta\rho_t(\mathbf r)=-\int \mathrm d^3\mathbf r\,|\nabla\rho_t(\mathbf r)|^2,
$$

$$
\int \mathrm d^3\mathbf r\,\rho_t(\mathbf r)\nabla\cdot\mathbf J_t(\mathbf r)=-\int \mathrm d^3\mathbf r\,\mathbf J_t(\mathbf r)\cdot\nabla\rho_t(\mathbf r).
$$

Thus the integrated energy depends on the combinations

$$
C_t^{\rho\Delta\rho}-C_t^{\nabla\rho\nabla\rho},\qquad C_t^{\rho\nabla J}-C_t^{J\nabla\rho}.
$$

The code retains both representations explicitly.

### Time-odd energy density

The general parameterization also includes the tensor-kinetic density $\mathbf F_t$. In terms of the nonlocal spin density,

$$
F_{t,\mu}(\mathbf r)=\frac12\sum_\nu\left[(\partial_\mu\partial_\nu'+\partial_\mu'\partial_\nu)s_{t,\nu}(\mathbf r,\mathbf r')\right]_{\mathbf r'=\mathbf r}.
$$

The time-odd interaction energy density is

$$
\boxed{\mathcal H_{\mathrm{odd}}(\mathbf r)=\sum_{t=0,1}\Bigl[C_t^{ss}[\rho_0(\mathbf r)]|\mathbf s_t(\mathbf r)|^2+C_t^{jj}|\mathbf j_t(\mathbf r)|^2+C_t^{s\Delta s}\mathbf s_t(\mathbf r)\cdot\Delta\mathbf s_t(\mathbf r)+C_t^{s\nabla j}\mathbf s_t(\mathbf r)\cdot(\nabla\times\mathbf j_t(\mathbf r))+C_t^{sT}\mathbf s_t(\mathbf r)\cdot\mathbf T_t(\mathbf r)+C_t^{\nabla s\nabla s}(\nabla\cdot\mathbf s_t(\mathbf r))^2+C_t^{sF}\mathbf s_t(\mathbf r)\cdot\mathbf F_t(\mathbf r)\Bigr].}
$$

These couplings are declared in the general parameter container; their presence does not establish a full time-reversal-breaking solver implementation.

For time-reversal symmetry,

$$
\mathbf s_t(\mathbf r)=\mathbf j_t(\mathbf r)=\mathbf T_t(\mathbf r)=\mathbf F_t(\mathbf r)=0,\qquad \mathcal H_{\mathrm{odd}}(\mathbf r)=0.
$$

### Conversion between Skyrme parameters and EDF couplings

For the standard Skyrme interaction without explicit tensor terms, use $\alpha=\gamma$ and the $\rho_t\Delta\rho_t$, $\rho_t\nabla\cdot\mathbf J_t$ representation.

**Time-even couplings**

The particle-density couplings are

$$
C_{0,0}^{\rho\rho}=\frac38t_0,\qquad C_{1,0}^{\rho\rho}=-\frac18t_0(1+2x_0),
$$

$$
C_{0,D}^{\rho\rho}=\frac1{16}t_3,\qquad C_{1,D}^{\rho\rho}=-\frac1{48}t_3(1+2x_3).
$$

The kinetic-density and surface couplings are

$$
C_0^{\rho\tau}=\frac1{16}\bigl[3t_1+(5+4x_2)t_2\bigr],\qquad C_1^{\rho\tau}=\frac1{16}\bigl[-(1+2x_1)t_1+(1+2x_2)t_2\bigr],
$$

$$
C_0^{\rho\Delta\rho}=\frac1{64}\bigl[-9t_1+(5+4x_2)t_2\bigr],\qquad C_1^{\rho\Delta\rho}=\frac1{64}\bigl[3(1+2x_1)t_1+(1+2x_2)t_2\bigr].
$$

The spin-current and spin-orbit couplings are

$$
C_0^{JJ}=\frac1{16}\bigl[(1-2x_1)t_1-(1+2x_2)t_2\bigr],\qquad C_1^{JJ}=\frac1{16}(t_1-t_2),
$$

$$
C_0^{\rho\nabla J}=-\frac34W_0,\qquad C_1^{\rho\nabla J}=-\frac14W_0.
$$

The remaining couplings vanish:

$$
C_t^{JJ^{\mathsf T}}=C_t^{\nabla\rho\nabla\rho}=C_t^{J\nabla\rho}=0.
$$

The $JJ$ term arises from the momentum-dependent central interaction even without explicit tensor terms.

**Time-odd couplings**

The spin-density couplings are

$$
C_{0,0}^{ss}=-\frac18t_0(1-2x_0),\qquad C_{1,0}^{ss}=-\frac18t_0,
$$

$$
C_{0,D}^{ss}=-\frac1{48}t_3(1-2x_3),\qquad C_{1,D}^{ss}=-\frac1{48}t_3.
$$

The spin-gradient couplings are

$$
C_0^{s\Delta s}=\frac1{64}\bigl[3(1-2x_1)t_1+(1+2x_2)t_2\bigr],\qquad C_1^{s\Delta s}=\frac1{64}(3t_1+t_2).
$$

The remaining nonzero couplings satisfy

$$
\boxed{C_t^{jj}=-C_t^{\rho\tau},\qquad C_t^{sT}=-C_t^{JJ},\qquad C_t^{s\nabla j}=C_t^{\rho\nabla J}.}
$$

Without explicit tensor terms,

$$
C_t^{\nabla s\nabla s}=C_t^{sF}=0.
$$

These relations apply to an EDF generated by the specified Skyrme interaction. Independently fitted EDF couplings, including UNEDF parameter sets, need not satisfy all of them; the reverse conversion reproduces the full functional only when the corresponding constraints hold.

### Pairing term

For each nucleon species $q=n,p$, introduce the local pairing density $\kappa_q(\mathbf r)$ and the density-dependent pairing coupling $C_q^{\mathrm{pair}}[\rho_0(\mathbf r)]$. Here $C_q^{V0}$ is the pairing strength, $C_q^{V1}$ is the dimensionless density factor, and $\rho_c=0.16\,\mathrm{fm}^{-3}$ is the reference density used in the code.

$$
C_q^{\mathrm{pair}}[\rho_0(\mathbf r)]=C_q^{V0}\left[1-C_q^{V1}\frac{\rho_0(\mathbf r)}{\rho_c}\right].
$$

With the code's normalization convention, the local pairing field $\Delta_q(\mathbf r)$ and pairing energy are

$$
\boxed{\Delta_q(\mathbf r)=C_q^{\mathrm{pair}}[\rho_0(\mathbf r)]\,\kappa_q(\mathbf r),\qquad E_{\mathrm{pair}}=\sum_{q=n,p}\int\mathrm d^3\mathbf r\,\kappa_q^*(\mathbf r)\Delta_q(\mathbf r).}
$$

The pairing density $\kappa_q$ supplements the normal densities; $C_q^{\mathrm{pair}}$ is a coupling, not a density. The pairing field $\Delta_q$ is distinct from the Laplacian operator $\Delta$.
