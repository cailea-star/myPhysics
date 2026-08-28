# Hartree–Fock–Bogoliubov Method

### Many-Body Hamiltonian

Use the fermionic conventions defined in [Second Quantization](../quantum_mechanics/second_quantization.md#fermionic-second-quantization). Let $\{|\phi_\alpha\rangle\}_{\alpha=1}^{M}$ be an orthonormal basis of $M$ single-particle modes, and let $|0\rangle$ be the particle vacuum. The operators $\hat c_\alpha^\dagger$ and $\hat c_\alpha$ create and annihilate a fermion in $|\phi_\alpha\rangle$:

$$
\hat c_\alpha^\dagger|0\rangle=|\phi_\alpha\rangle,\qquad \hat c_\alpha|0\rangle=0.
$$

Their canonical anticommutation relations are

$$
\{\hat c_\alpha,\hat c_\beta^\dagger\}=\delta_{\alpha\beta},\qquad \{\hat c_\alpha,\hat c_\beta\}=\{\hat c_\alpha^\dagger,\hat c_\beta^\dagger\}=0,
$$

where $\delta_{\alpha\beta}$ is the Kronecker delta.

Let $\hat t$ be the one-body kinetic operator and $\hat v$ a symmetric two-particle interaction. Their matrix elements are

$$
t_{\alpha\beta}\equiv\langle\phi_\alpha|\hat t|\phi_\beta\rangle,\qquad v_{\alpha\beta;\gamma\delta}\equiv\langle\phi_\alpha\phi_\beta|\hat v|\phi_\gamma\phi_\delta\rangle.
$$

The antisymmetrized two-body matrix element is

$$
\bar v_{\alpha\beta;\gamma\delta}\equiv v_{\alpha\beta;\gamma\delta}-v_{\alpha\beta;\delta\gamma}.
$$

For a Hermitian two-particle interaction symmetric under particle exchange,

$$
\bar v_{\alpha\beta;\gamma\delta}=-\bar v_{\beta\alpha;\gamma\delta}=-\bar v_{\alpha\beta;\delta\gamma}=\bar v_{\beta\alpha;\delta\gamma},\qquad \bar v_{\alpha\beta;\gamma\delta}^*=\bar v_{\gamma\delta;\alpha\beta}.
$$

The second-quantized many-body Hamiltonian is

$$
\boxed{\hat H=\sum_{\alpha,\beta=1}^{M}t_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta+\frac14\sum_{\alpha,\beta,\gamma,\delta=1}^{M}\bar v_{\alpha\beta;\gamma\delta}\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma}.
$$

The particle-number operator is

$$
\hat N=\sum_{\alpha=1}^{M}\hat c_\alpha^\dagger\hat c_\alpha.
$$

### Bogoliubov Transformation

Let $\mu,\nu=1,\ldots,M$ label quasiparticle modes. Complex conjugation, transpose, and Hermitian conjugation are denoted by $*$, $T$, and $\dagger$, respectively. Define the particle and quasiparticle Nambu columns by

$$
\hat{\mathcal C}\equiv\begin{pmatrix}\hat c_1&\cdots&\hat c_M&\hat c_1^\dagger&\cdots&\hat c_M^\dagger\end{pmatrix}^{T},\qquad \hat{\mathcal B}\equiv\begin{pmatrix}\hat\beta_1&\cdots&\hat\beta_M&\hat\beta_1^\dagger&\cdots&\hat\beta_M^\dagger\end{pmatrix}^{T}.
$$

Let $U,V\in\mathbb C^{M\times M}$ be the particle and hole amplitude matrices. The Bogoliubov transformation and its inverse are

$$
\boxed{\hat{\mathcal B}=\mathcal W^\dagger\hat{\mathcal C},\qquad \hat{\mathcal C}=\mathcal W\hat{\mathcal B},\qquad \mathcal W\equiv\begin{pmatrix}U&V^*\\V&U^*\end{pmatrix}}.
$$

Let $I_n$ denote the $n\times n$ identity matrix, and define the particle–hole exchange matrix by

$$
\tau_x\equiv\begin{pmatrix}0&I_M\\I_M&0\end{pmatrix}.
$$

The particle and hole components of a Nambu column are not independent, and the Bogoliubov transformation preserves this conjugation structure:

$$
\hat{\mathcal C}=\tau_x\hat{\mathcal C}^*,\qquad \hat{\mathcal B}=\tau_x\hat{\mathcal B}^*,\qquad \tau_x\mathcal W^*\tau_x=\mathcal W.
$$

This constraint expresses the redundancy introduced by Nambu doubling and is not an additional physical symmetry.

In components,

$$
\hat\beta_\mu=\sum_{\alpha=1}^{M}\left(U_{\alpha\mu}^*\hat c_\alpha+V_{\alpha\mu}^*\hat c_\alpha^\dagger\right).
$$

$$
\hat c_\alpha=\sum_{\mu=1}^{M}\left(U_{\alpha\mu}\hat\beta_\mu+V_{\alpha\mu}^*\hat\beta_\mu^\dagger\right).
$$

The quasiparticle operators must satisfy the fermionic canonical anticommutation relations

$$
\{\hat\beta_\mu,\hat\beta_\nu^\dagger\}=\delta_{\mu\nu},\qquad \{\hat\beta_\mu,\hat\beta_\nu\}=\{\hat\beta_\mu^\dagger,\hat\beta_\nu^\dagger\}=0.
$$

Preservation of the anticommutation relations requires

$$
\mathcal W^\dagger\mathcal W=\mathcal W\mathcal W^\dagger=I_{2M}.
$$

Equivalently,

$$
U^\dagger U+V^\dagger V=I_M,\qquad U^TV+V^TU=0,
$$

and

$$
UU^\dagger+V^*V^T=I_M,\qquad UV^\dagger+V^*U^T=0.
$$

The quasiparticle vacuum $|\Phi\rangle$ is defined by

$$
\boxed{\hat\beta_\mu|\Phi\rangle=0,\qquad \mu=1,\ldots,M}.
$$

Because each $\hat\beta_\mu$ mixes particle creation and annihilation operators, $|\Phi\rangle$ is generally a superposition of different particle-number sectors and is not an eigenstate of $\hat N$.

### Generalized Density Matrix

For the normalized quasiparticle vacuum $|\Phi\rangle$, define the normal density matrix $\rho$ and pairing tensor $\kappa$ by

$$
\rho_{\alpha\beta}\equiv\langle\Phi|\hat c_\beta^\dagger\hat c_\alpha|\Phi\rangle,\qquad \kappa_{\alpha\beta}\equiv\langle\Phi|\hat c_\beta\hat c_\alpha|\Phi\rangle.
$$

Using the inverse Bogoliubov transformation,

$$
\rho_{\alpha\beta}=\sum_{\mu=1}^{M}V_{\alpha\mu}^*V_{\beta\mu},\qquad \kappa_{\alpha\beta}=\sum_{\mu=1}^{M}V_{\alpha\mu}^*U_{\beta\mu}.
$$

Hence,

$$
\rho=V^*V^T,\qquad \kappa=V^*U^T,\qquad \rho^\dagger=\rho,\qquad \kappa^T=-\kappa.
$$

The generalized density matrix is

$$
\mathcal R\equiv\begin{pmatrix}\langle\Phi|\hat c_\beta^\dagger\hat c_\alpha|\Phi\rangle&\langle\Phi|\hat c_\beta\hat c_\alpha|\Phi\rangle\\\langle\Phi|\hat c_\beta^\dagger\hat c_\alpha^\dagger|\Phi\rangle&\langle\Phi|\hat c_\beta\hat c_\alpha^\dagger|\Phi\rangle\end{pmatrix}=\begin{pmatrix}V^*\\U^*\end{pmatrix}\begin{pmatrix}V^T&U^T\end{pmatrix}=\begin{pmatrix}\rho&\kappa\\-\kappa^*&I_M-\rho^*\end{pmatrix}.
$$

It inherits the Nambu constraint

$$
\tau_x\mathcal R^*\tau_x=I_{2M}-\mathcal R.
$$

In the quasiparticle basis,

$$
\mathcal W^\dagger\mathcal R\mathcal W=\begin{pmatrix}\langle\Phi|\hat\beta_\nu^\dagger\hat\beta_\mu|\Phi\rangle&\langle\Phi|\hat\beta_\nu\hat\beta_\mu|\Phi\rangle\\\langle\Phi|\hat\beta_\nu^\dagger\hat\beta_\mu^\dagger|\Phi\rangle&\langle\Phi|\hat\beta_\nu\hat\beta_\mu^\dagger|\Phi\rangle\end{pmatrix}=\begin{pmatrix}0\\I_M\end{pmatrix}\begin{pmatrix}0&I_M\end{pmatrix}=\begin{pmatrix}0&0\\0&I_M\end{pmatrix}.
$$

The generalized density matrix is therefore Hermitian and idempotent:

$$
\boxed{\mathcal R^\dagger=\mathcal R,\qquad \mathcal R^2=\mathcal R}.
$$

The block form of $\mathcal R^2=\mathcal R$ gives

$$
\rho-\rho^2=\kappa\kappa^\dagger,\qquad \rho\kappa=\kappa\rho^*.
$$

The average particle number is

$$
\langle\Phi|\hat N|\Phi\rangle=\operatorname{Tr}\rho=\operatorname{Tr}(V^\dagger V).
$$

### $\Gamma$ and $\Delta$ Fields

The generalized Wick theorem gives

$$
\langle\Phi|\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma|\Phi\rangle=\rho_{\gamma\alpha}\rho_{\delta\beta}-\rho_{\delta\alpha}\rho_{\gamma\beta}+\kappa_{\alpha\beta}^*\kappa_{\gamma\delta}.
$$

The expectation value of the many-body Hamiltonian defines the energy functional.

The allowed variations satisfy $\delta\rho^\dagger=\delta\rho$ and $\delta\kappa^T=-\delta\kappa$, while $\kappa$ and $\kappa^*$ are treated as formally independent variables in the functional variation.

$$
E[\rho,\kappa,\kappa^*]\equiv\langle\Phi|\hat H|\Phi\rangle=\sum_{\alpha,\beta=1}^{M}t_{\alpha\beta}\rho_{\beta\alpha}+\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{M}\bar v_{\alpha\beta;\gamma\delta}\rho_{\gamma\alpha}\rho_{\delta\beta}+\frac14\sum_{\alpha,\beta,\gamma,\delta=1}^{M}\bar v_{\alpha\beta;\gamma\delta}\kappa_{\alpha\beta}^*\kappa_{\gamma\delta}.
$$

Define the Hartree–Fock field $\Gamma$, single-particle field $h$, and pairing field $\Delta$ by

$$
\boxed{\Gamma_{\alpha\beta}\equiv\sum_{\gamma,\delta=1}^{M}\bar v_{\alpha\gamma;\beta\delta}\rho_{\delta\gamma},\qquad h_{\alpha\beta}\equiv t_{\alpha\beta}+\Gamma_{\alpha\beta},\qquad \Delta_{\alpha\beta}\equiv\frac12\sum_{\gamma,\delta=1}^{M}\bar v_{\alpha\beta;\gamma\delta}\kappa_{\gamma\delta}}.
$$

Their symmetry properties are

$$
\Gamma^\dagger=\Gamma,\qquad h^\dagger=h,\qquad \Delta^T=-\Delta.
$$

The energy functional becomes

$$
E[\rho,\kappa,\kappa^*]=\operatorname{Tr}(t\rho)+\frac12\operatorname{Tr}(\Gamma\rho)+\frac12\operatorname{Tr}(\Delta\kappa^\dagger).
$$

Its first variation is

$$
\boxed{\delta E=\operatorname{Tr}(h\,\delta\rho)+\frac12\operatorname{Tr}(\Delta\,\delta\kappa^\dagger)+\frac12\operatorname{Tr}(\Delta^\dagger\delta\kappa)}.
$$

For a general energy density functional, $h$ and $\Delta$ are defined by the same first-variation relation.

### Generalized Single-Particle Hamiltonian

Let $N_0$ be the target particle number and $\lambda\in\mathbb R$ its Lagrange multiplier. Define

$$
\mathcal L[\rho,\kappa,\kappa^*;\lambda]\equiv E[\rho,\kappa,\kappa^*]-\lambda\left(\langle\Phi|\hat N|\Phi\rangle-N_0\right)=E[\rho,\kappa,\kappa^*]-\lambda\left(\operatorname{Tr}_M\rho-N_0\right).
$$

Here $\operatorname{Tr}_M$ denotes the trace over the $M$ single-particle modes. Stationarity with respect to $\lambda$ imposes

$$
\frac{\partial\mathcal L}{\partial\lambda}=-\left(\operatorname{Tr}_M\rho-N_0\right)=0.
$$

Its first variation with respect to the densities is

$$
\delta\mathcal L=\operatorname{Tr}_M\!\left[(h-\lambda I_M)\delta\rho\right]+\frac12\operatorname{Tr}_M(\Delta\,\delta\kappa^\dagger)+\frac12\operatorname{Tr}_M(\Delta^\dagger\delta\kappa).
$$

The variation of the generalized density matrix is

$$
\delta\mathcal R=\begin{pmatrix}\delta\rho&\delta\kappa\\-\delta\kappa^*&-\delta\rho^*\end{pmatrix}.
$$

Define the generalized single-particle Hamiltonian by

$$
\boxed{\mathcal H\equiv\begin{pmatrix}h-\lambda I_M&\Delta\\-\Delta^*&-h^*+\lambda I_M\end{pmatrix}}.
$$

It is Hermitian and inherits the Nambu constraint

$$
\mathcal H^\dagger=\mathcal H,\qquad \tau_x\mathcal H^*\tau_x=-\mathcal H.
$$

Let $\operatorname{Tr}_{2M}$ denote the trace over Nambu space. Expanding the block-matrix trace gives

$$
\operatorname{Tr}_{2M}(\mathcal H\,\delta\mathcal R)=\operatorname{Tr}_M\!\left[(h-\lambda I_M)\delta\rho-\Delta\,\delta\kappa^*-\Delta^*\delta\kappa+(h^*-\lambda I_M)\delta\rho^*\right].
$$

The fields and allowed variations satisfy

$$
h^\dagger=h,\qquad \Delta^T=-\Delta,\qquad \delta\rho^\dagger=\delta\rho,\qquad \delta\kappa^T=-\delta\kappa.
$$

Hence,

$$
\operatorname{Tr}_{2M}(\mathcal H\,\delta\mathcal R)=2\operatorname{Tr}_M\!\left[(h-\lambda I_M)\delta\rho\right]+\operatorname{Tr}_M(\Delta\,\delta\kappa^\dagger)+\operatorname{Tr}_M(\Delta^\dagger\delta\kappa)=2\delta\mathcal L.
$$

Thus,

$$
\boxed{\delta\mathcal L=\frac12\operatorname{Tr}_{2M}(\mathcal H\,\delta\mathcal R)}.
$$

The factor $1/2$ removes the Nambu-space double counting.

### Variational HFB Equation

Since $\mathcal R$ is a Hermitian projector, its variation is restricted to a unitary rotation that preserves the Nambu structure. Let $\mathcal G$ be an admissible Hermitian generator:

$$
\mathcal G^\dagger=\mathcal G,\qquad \tau_x\mathcal G^*\tau_x=-\mathcal G.
$$

For each variation direction $\mathcal G$ and real parameter $\epsilon$, define

$$
\mathcal R(\epsilon)=e^{i\epsilon\mathcal G}\mathcal R e^{-i\epsilon\mathcal G}.
$$

The corresponding infinitesimal variation is

$$
\delta\mathcal R=\left.\frac{d\mathcal R(\epsilon)}{d\epsilon}\right|_{\epsilon=0}=i[\mathcal G,\mathcal R].
$$

The stationary condition is

$$
0=\delta\mathcal L=\frac{i}{2}\operatorname{Tr}_{2M}\!\left(\mathcal H[\mathcal G,\mathcal R]\right)=\frac{i}{2}\operatorname{Tr}_{2M}\!\left([\mathcal R,\mathcal H]\mathcal G\right).
$$

Stationarity with respect to every admissible $\mathcal G$ gives the variational HFB equation

$$
\boxed{[\mathcal H,\mathcal R]=0}.
$$

Since $\mathcal H$ and $\mathcal R$ are Hermitian and commute, their eigenvectors can be chosen in common. Choose the Bogoliubov matrix such that

$$
\mathcal W^\dagger\mathcal R\mathcal W=\begin{pmatrix}0&0\\0&I_M\end{pmatrix}.
$$

Let

$$
\Psi_\mu\equiv\begin{pmatrix}U_\mu\\V_\mu\end{pmatrix},\qquad \widetilde\Psi_\mu\equiv\tau_x\Psi_\mu^*=\begin{pmatrix}V_\mu^*\\U_\mu^*\end{pmatrix},\qquad \mathcal R=\sum_{\nu=1}^{M}\widetilde\Psi_\nu\widetilde\Psi_\nu^\dagger.
$$

The Nambu constraint pairs the two eigenspaces:

$$
\mathcal R\Psi_\mu=\left(\sum_{\nu=1}^{M}\widetilde\Psi_\nu\widetilde\Psi_\nu^\dagger\right)\Psi_\mu=0,\qquad \mathcal R\widetilde\Psi_\mu=\left(\sum_{\nu=1}^{M}\widetilde\Psi_\nu\widetilde\Psi_\nu^\dagger\right)\widetilde\Psi_\mu=\widetilde\Psi_\mu.
$$

For a stable solution without zero modes, let $E_\mu>0$ denote the positive member of each Nambu pair. The constraint on $\mathcal H$ gives

$$
\mathcal H\Psi_\mu=E_\mu\Psi_\mu,\qquad \mathcal H\widetilde\Psi_\mu=-E_\mu\widetilde\Psi_\mu.
$$

With $\Psi_\mu$ and $\widetilde\Psi_\mu$ arranged as the columns of $\mathcal W$,

$$
\mathcal W^\dagger\mathcal H\mathcal W=\begin{pmatrix}E&0\\0&-E\end{pmatrix},\qquad E\equiv\operatorname{diag}(E_1,\ldots,E_M).
$$

The positive-energy HFB equation is

$$
\begin{pmatrix}h-\lambda I_M&\Delta\\-\Delta^*&-h^*+\lambda I_M\end{pmatrix}\begin{pmatrix}U_\mu\\V_\mu\end{pmatrix}=E_\mu\begin{pmatrix}U_\mu\\V_\mu\end{pmatrix}.
$$

### Quasiparticle Vacuum

The paired Nambu vectors define the quasiparticle operators, vacuum, and one-quasiparticle states:

$$
\hat\beta_\mu=\Psi_\mu^\dagger\hat{\mathcal C},\qquad \hat\beta_\mu^\dagger=\widetilde\Psi_\mu^\dagger\hat{\mathcal C},\qquad \hat\beta_\mu|\Phi\rangle=0,\qquad |\mu\rangle\equiv\hat\beta_\mu^\dagger|\Phi\rangle.
$$

The quadratic mean-field Hamiltonian associated with $\mathcal H$, distinct from the original interacting Hamiltonian $\hat H$, is

$$
\hat H_{\mathrm{qp}}=E_{\mathrm{vac}}+\sum_{\mu=1}^{M}E_\mu\hat\beta_\mu^\dagger\hat\beta_\mu,\qquad E_\mu>0.
$$

Hence,

$$
\hat H_{\mathrm{qp}}|\Phi\rangle=E_{\mathrm{vac}}|\Phi\rangle,\qquad \hat H_{\mathrm{qp}}|\mu\rangle=(E_{\mathrm{vac}}+E_\mu)|\mu\rangle.
$$

Although $\widetilde\Psi_\mu$ carries the Nambu eigenvalue $-E_\mu$, it represents the creation operator $\hat\beta_\mu^\dagger$, which generates the many-body excitation $\hat\beta_\mu^\dagger|\Phi\rangle$ with energy $+E_\mu$.

The zero-eigenvalue subspace represents unoccupied quasiparticle modes.

$$
\mathcal R\Psi_\mu=0\quad\Longleftrightarrow\quad\langle\Phi|\hat\beta_\mu^\dagger\hat\beta_\mu|\Phi\rangle=0.
$$

The fermionic anticommutation relation gives the complementary Nambu contraction:

$$
\mathcal R\widetilde\Psi_\mu=\widetilde\Psi_\mu\quad\Longleftrightarrow\quad\langle\Phi|\hat\beta_\mu\hat\beta_\mu^\dagger|\Phi\rangle=1.
$$

Only the positive-energy solutions are independent; the negative-energy branch follows from Nambu conjugation.

### Self-Consistent Iteration

The HFB equation is solved through nested outer and inner iterations:

1. Construct the fields from the densities at iteration $n$:

$$
h^{(n)}=t+\Gamma[\rho^{(n)}],\qquad \Delta^{(n)}=\Delta[\kappa^{(n)}].
$$

2. At fixed $h^{(n)}$ and $\Delta^{(n)}$, choose $\lambda$, construct and diagonalize $\mathcal H^{(n)}(\lambda)$, and retain only the independent solutions with $E_{\mu,\lambda}^{(n)}>0$:

$$
\mathcal H^{(n)}(\lambda)=\begin{pmatrix}h^{(n)}-\lambda I_M&\Delta^{(n)}\\-\Delta^{(n)*}&-h^{(n)*}+\lambda I_M\end{pmatrix}.
$$

$$
\mathcal H^{(n)}(\lambda)\begin{pmatrix}U_{\mu,\lambda}^{(n)}\\V_{\mu,\lambda}^{(n)}\end{pmatrix}=E_{\mu,\lambda}^{(n)}\begin{pmatrix}U_{\mu,\lambda}^{(n)}\\V_{\mu,\lambda}^{(n)}\end{pmatrix}.
$$

3. Construct the densities and particle number from the positive-energy solutions:

$$
\rho_\lambda^{(n)}=V_\lambda^{(n)*}V_\lambda^{(n)T},\qquad \kappa_\lambda^{(n)}=V_\lambda^{(n)*}U_\lambda^{(n)T},\qquad N_\lambda^{(n)}=\operatorname{Tr}_M\rho_\lambda^{(n)}.
$$

4. Repeat steps 2–3 and adjust $\lambda$ until the root $\lambda_*^{(n)}$ satisfies

$$
N_{\lambda_*^{(n)}}^{(n)}-N_0=0.
$$

5. Update the outer iteration with the particle-number-converged densities and return to step 1.

6. The iteration terminates when

$$
\rho^{(n+1)}-\rho^{(n)}\rightarrow0,\qquad \kappa^{(n+1)}-\kappa^{(n)}\rightarrow0,\qquad N_{\lambda_*^{(n)}}^{(n)}-N_0\rightarrow0.
$$
