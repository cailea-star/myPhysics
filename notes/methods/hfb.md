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
\hat{\mathcal C}\equiv\begin{pmatrix}\hat c_1&\cdots&\hat c_M&\hat c_1^\dagger&\cdots&\hat c_M^\dagger\end{pmatrix}^T,\qquad \hat{\mathcal B}\equiv\begin{pmatrix}\hat\beta_1&\cdots&\hat\beta_M&\hat\beta_1^\dagger&\cdots&\hat\beta_M^\dagger\end{pmatrix}^T.
$$

Let $U,V\in\mathbb C^{M\times M}$ be the particle and hole amplitude matrices. The Bogoliubov transformation and its inverse are

$$
\boxed{\hat{\mathcal B}=\mathcal W^\dagger \hat{\mathcal C},\qquad \hat{\mathcal C}=\mathcal W \hat{\mathcal B},\qquad \mathcal W\equiv\begin{pmatrix}U&V^*\\V&U^*\end{pmatrix}}.
$$

We can also write the Bogoliubov transformation in block form:

$$
\hat{\boldsymbol c}\equiv\begin{pmatrix}\hat c_1&\cdots&\hat c_M\end{pmatrix}^T,\qquad \hat{\boldsymbol\beta}\equiv\begin{pmatrix}\hat\beta_1&\cdots&\hat\beta_M\end{pmatrix}^T.
$$

$$
\begin{pmatrix}\hat{\boldsymbol\beta}\\\hat{\boldsymbol\beta}^{\dagger T}\end{pmatrix}=\begin{pmatrix}U^\dagger&V^\dagger\\V^T&U^T\end{pmatrix}\begin{pmatrix}\hat{\boldsymbol c}\\\hat{\boldsymbol c}^{\dagger T}\end{pmatrix},\qquad \begin{pmatrix}\hat{\boldsymbol c}\\\hat{\boldsymbol c}^{\dagger T}\end{pmatrix}=\begin{pmatrix}U&V^*\\V&U^*\end{pmatrix}\begin{pmatrix}\hat{\boldsymbol\beta}\\\hat{\boldsymbol\beta}^{\dagger T}\end{pmatrix}.
$$

In components,

$$
\hat\beta_\mu=\sum_{\alpha=1}^{M}\left(U_{\alpha\mu}^*\hat c_\alpha+V_{\alpha\mu}^*\hat c_\alpha^\dagger\right),\qquad \hat c_\alpha=\sum_{\mu=1}^{M}\left(U_{\alpha\mu}\hat\beta_\mu+V_{\alpha\mu}^*\hat\beta_\mu^\dagger\right).
$$


Let $I_n$ denote the $n\times n$ identity matrix, and define the particle–hole exchange matrix by

$$
\tau_x\equiv\begin{pmatrix}0&I_M\\I_M&0\end{pmatrix}.
$$

The particle and hole components of a Nambu column are not independent, and the Bogoliubov transformation preserves this conjugation structure:

$$
\hat{\mathcal C}^*=\tau_x\hat{\mathcal C},\qquad \hat{\mathcal B}^*=\tau_x\hat{\mathcal B},\qquad \mathcal W^*=\tau_x\mathcal W\tau_x.
$$

This constraint expresses the redundancy introduced by Nambu doubling and is not an additional physical symmetry.

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
\boxed{E[\rho,\kappa,\kappa^*]=\operatorname{Tr}(t\rho)+\frac12\operatorname{Tr}(\Gamma\rho)+\frac12\operatorname{Tr}(\Delta\kappa^\dagger)}.
$$

Its first variation is

$$
\boxed{\delta E=\operatorname{Tr}(h\,\delta\rho)+\frac12\operatorname{Tr}(\Delta\,\delta\kappa^\dagger)+\frac12\operatorname{Tr}(\Delta^\dagger\delta\kappa)}.
$$

For a general energy density functional, $h$ and $\Delta$ are defined by the same first-variation relation.

### Generalized Single-Particle Hamiltonian

Let $N_0$ be the target particle number and $\lambda\in\mathbb R$ its Lagrange multiplier. Define

$$
\mathcal L[\rho,\kappa,\kappa^*;\lambda]\equiv E[\rho,\kappa,\kappa^*]-\lambda\big(\langle\Phi|\hat N|\Phi\rangle-N_0\big).
$$

We can expand the Lagrangian in terms of the fields:

$$
\boxed{\mathcal L[\rho,\kappa,\kappa^*;\lambda]=\operatorname{Tr}(t\rho)+\frac12\operatorname{Tr}(\Gamma\rho)+\frac12\operatorname{Tr}(\Delta\kappa^\dagger)-\lambda\big(\operatorname{Tr}\rho-N_0\big)}.
$$

Stationarity with respect to $\lambda$ imposes

$$
\frac{\partial\mathcal L}{\partial\lambda}=-\left(\operatorname{Tr}\rho-N_0\right)=0.
$$

Its first variation with respect to the densities is

$$
\boxed{\delta\mathcal L=\operatorname{Tr}\!\left[(h-\lambda I)\delta\rho\right]+\frac12\operatorname{Tr}_M(\Delta\,\delta\kappa^\dagger)+\frac12\operatorname{Tr}_M(\Delta^\dagger\delta\kappa)}.
$$

The variation of the generalized density matrix is

$$
\delta\mathcal R=\begin{pmatrix}\delta\rho&\delta\kappa\\-\delta\kappa^*&-\delta\rho^*\end{pmatrix}.
$$

Define the generalized single-particle Hamiltonian by

$$
\boxed{\mathcal H\equiv\begin{pmatrix}h-\lambda I&\Delta\\-\Delta^*&-h^*+\lambda I\end{pmatrix}}.
$$

It is Hermitian and inherits the Nambu constraint

$$
\mathcal H^\dagger=\mathcal H,\qquad \tau_x\mathcal H^*\tau_x=-\mathcal H.
$$

Expanding the block-matrix trace gives

$$
\operatorname{Tr}(\mathcal H\,\delta\mathcal R)=\operatorname{Tr}\!\left[(h-\lambda I)\delta\rho-\Delta\,\delta\kappa^*-\Delta^*\delta\kappa+(h^*-\lambda I)\delta\rho^*\right].
$$

The fields and allowed variations satisfy

$$
h^\dagger=h,\qquad \Delta^T=-\Delta,\qquad \delta\rho^\dagger=\delta\rho,\qquad \delta\kappa^T=-\delta\kappa.
$$

Hence,

$$
\operatorname{Tr}(\mathcal H\,\delta\mathcal R)=2\operatorname{Tr}\!\left[(h-\lambda I)\delta\rho\right]+\operatorname{Tr}(\Delta\,\delta\kappa^\dagger)+\operatorname{Tr}(\Delta^\dagger\delta\kappa)=2\delta\mathcal L.
$$

Thus,

$$
\boxed{\delta\mathcal L=\frac12\operatorname{Tr}(\mathcal H\,\delta\mathcal R)}.
$$

The factor $1/2$ removes the Nambu-space double counting.

### Solution of the Generalized Single-Particle Eigenproblem

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

Since $[\mathcal H,\mathcal R]=0$, choose $\mathcal W$ to diagonalize both matrices. For a stable solution without zero modes, let

$$
E\equiv\operatorname{diag}(E_1,\ldots,E_M),\qquad E_\mu>0.
$$

Then

$$
\mathcal W^\dagger\mathcal H\mathcal W=\begin{pmatrix}E&0\\0&-E\end{pmatrix},\qquad \mathcal W^\dagger\mathcal R\mathcal W=\begin{pmatrix}0&0\\0&I\end{pmatrix}.
$$

The Bogoliubov matrix and its paired columns are

$$
\mathcal W=\begin{pmatrix}U&V^*\\V&U^*\end{pmatrix},\qquad \Psi_\mu(E_\mu)\equiv\begin{pmatrix}U_\mu\\V_\mu\end{pmatrix},\qquad \Psi_\mu(-E_\mu)\equiv\tau_x\Psi_\mu^*(E_\mu)=\begin{pmatrix}V_\mu^*\\U_\mu^*\end{pmatrix}.
$$

The positive-energy columns satisfy

$$
\mathcal H\Psi_\mu(E_\mu)=E_\mu\Psi_\mu(E_\mu),\qquad \mathcal R\Psi_\mu(E_\mu)=0.
$$

The negative-energy columns satisfy

$$
\mathcal H\Psi_\mu(-E_\mu)=-E_\mu\Psi_\mu(-E_\mu),\qquad \mathcal R\Psi_\mu(-E_\mu)=\Psi_\mu(-E_\mu).
$$

### Quasiparticle Diagonalization

At the stationary fields, define the quadratic HFB operator by

$$
\boxed{\hat{\mathcal L}_{\mathrm{HFB}}\equiv\sum_{\alpha,\beta=1}^{M}\left(h_{\alpha\beta}-\lambda\delta_{\alpha\beta}\right)\hat c_\alpha^\dagger\hat c_\beta+\frac12\sum_{\alpha,\beta=1}^{M}\left(\Delta_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta^\dagger+\Delta_{\alpha\beta}^*\hat c_\beta\hat c_\alpha\right)}.
$$

Its Nambu representation is

$$
\hat{\mathcal L}_{\mathrm{HFB}}=\frac12\hat{\mathcal C}^\dagger\mathcal H\hat{\mathcal C}+\frac12\operatorname{Tr}(h-\lambda I).
$$

The trace term cancels the constant generated by normal ordering the Nambu quadratic form.

The quasiparticle transformation and diagonalization give

$$
\hat{\mathcal C}=\mathcal W\hat{\mathcal B},\qquad \mathcal W^\dagger\mathcal H\mathcal W=\begin{pmatrix}E&0\\0&-E\end{pmatrix}.
$$

Expanding the Nambu components gives

$$
\hat{\mathcal L}_{\mathrm{HFB}}=\frac12\operatorname{Tr}(h-\lambda I)+\frac12\sum_{\mu=1}^{M}E_\mu\left(\hat\beta_\mu^\dagger\hat\beta_\mu-\hat\beta_\mu\hat\beta_\mu^\dagger\right).
$$

Using the anticommutation relation $\{\hat\beta_\mu,\hat\beta_\mu^\dagger\}=1$ gives

$$
\hat{\mathcal L}_{\mathrm{HFB}}=\frac12\operatorname{Tr}(h-\lambda I)-\frac12\sum_{\mu=1}^{M}E_\mu+\sum_{\mu=1}^{M}E_\mu\hat\beta_\mu^\dagger\hat\beta_\mu.
$$

Since $\hat\beta_\mu|\Phi\rangle=0$, define

$$
\boxed{E_{\mathrm{vac}}\equiv\langle\Phi|\hat{\mathcal L}_{\mathrm{HFB}}|\Phi\rangle=\frac12\operatorname{Tr}(h-\lambda I)-\frac12\sum_{\mu=1}^{M}E_\mu}.
$$

The quadratic HFB operator is therefore

$$
\boxed{\hat{\mathcal L}_{\mathrm{HFB}}=E_{\mathrm{vac}}+\sum_{\mu=1}^{M}E_\mu\hat\beta_\mu^\dagger\hat\beta_\mu}.
$$

Define the one-quasiparticle state by

$$
|\mu\rangle\equiv\hat\beta_\mu^\dagger|\Phi\rangle.
$$

Then,

$$
\hat{\mathcal L}_{\mathrm{HFB}}|\Phi\rangle=E_{\mathrm{vac}}|\Phi\rangle,\qquad \hat{\mathcal L}_{\mathrm{HFB}}|\mu\rangle=(E_{\mathrm{vac}}+E_\mu)|\mu\rangle.
$$

Thus, $E_\mu$ is the excitation energy relative to the quasiparticle vacuum. Although $\Psi_\mu(-E_\mu)$ has Nambu eigenvalue $-E_\mu$, its Hermitian conjugate defines $\hat\beta_\mu^\dagger$ and creates an excitation of energy $E_\mu$.

The physical HFB energy is evaluated from

$$
E_{\mathrm{HFB}}\equiv E[\rho,\kappa,\kappa^*],
$$

and is not generally equal to $E_{\mathrm{vac}}$ or $\sum_\mu E_\mu$.

### Quasiparticle Vacuum

The positive-energy Nambu covector $\Psi_\mu^\dagger(E_\mu)$ defines the quasiparticle annihilation operator

$$
\hat\beta_\mu=\Psi_\mu^\dagger(E_\mu)\hat{\mathcal C}=\sum_{\alpha=1}^{M}\left(U_{\alpha\mu}^*\hat c_\alpha+V_{\alpha\mu}^*\hat c_\alpha^\dagger\right).
$$

The negative-energy Nambu covector $\Psi_\mu^\dagger(-E_\mu)$ defines the quasiparticle creation operator

$$
\hat\beta_\mu^\dagger=\Psi_\mu^\dagger(-E_\mu)\hat{\mathcal C}=\sum_{\alpha=1}^{M}\left(V_{\alpha\mu}\hat c_\alpha+U_{\alpha\mu}\hat c_\alpha^\dagger\right).
$$

The normalized quasiparticle vacuum $|\Phi\rangle$ is defined by

$$
\boxed{\hat\beta_\mu|\Phi\rangle=0,\qquad \langle\Phi|\Phi\rangle=1,\qquad \mu=1,\ldots,M}.
$$

Assume that $U$ is invertible. The Bogoliubov relation gives

$$
U^TV+V^TU=0\;\to\;U^{-T}V^T=-VU^{-1}.
$$

Define the Thouless matrix by

$$
Z\equiv(VU^{-1})^*=V^*U^{-*}.
$$

It is antisymmetric:

$$
Z^T=U^{-\dagger}V^\dagger=-V^*U^{-*}=-Z.
$$

Let $|0\rangle$ denote the particle vacuum, and define its amplitude in $|\Phi\rangle$ by

$$
\nu\equiv\langle0|\Phi\rangle.
$$

The quasiparticle vacuum has the Thouless form

$$
\boxed{|\Phi\rangle=\nu\exp\left[\frac12\sum_{\alpha,\beta=1}^{M}Z_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta^\dagger\right]|0\rangle}.
$$

Define the pair-creation operator

$$
\hat A\equiv\frac12\sum_{\gamma,\delta=1}^{M}Z_{\gamma\delta}\hat c_\gamma^\dagger\hat c_\delta^\dagger.
$$

Its commutator with an annihilation operator is

$$
[\hat c_\alpha,\hat A]=\frac12\sum_{\gamma,\delta=1}^{M}Z_{\gamma\delta}\left(\delta_{\alpha\gamma}\hat c_\delta^\dagger-\delta_{\alpha\delta}\hat c_\gamma^\dagger\right)=\sum_{\beta=1}^{M}Z_{\alpha\beta}\hat c_\beta^\dagger.
$$

Its commutator with a creation operator is

$$
[\hat c_\alpha^\dagger,\hat A]=\frac12\sum_{\gamma,\delta=1}^{M}Z_{\gamma\delta}\left(\hat c_\alpha^\dagger\hat c_\gamma^\dagger\hat c_\delta^\dagger-\hat c_\gamma^\dagger\hat c_\delta^\dagger\hat c_\alpha^\dagger\right)=0.
$$

The Baker–Campbell–Hausdorff series therefore terminates at first order:

$$
[\hat c_\alpha^\dagger,\hat A]=0\;\to\;[[\hat c_\alpha,\hat A],\hat A]=0\;\to\;\hat c_\alpha e^{\hat A}=e^{\hat A}\left(\hat c_\alpha+[\hat c_\alpha,\hat A]\right).
$$

Using $\hat c_\alpha|0\rangle=0$ gives

$$
\hat c_\alpha|\Phi\rangle=\sum_{\beta=1}^{M}Z_{\alpha\beta}\hat c_\beta^\dagger|\Phi\rangle,\qquad \hat{\boldsymbol c}|\Phi\rangle=Z\hat{\boldsymbol c}^\dagger|\Phi\rangle.
$$

The required operator and matrix relations are

$$
\hat{\boldsymbol\beta}=U^\dagger\hat{\boldsymbol c}+V^\dagger\hat{\boldsymbol c}^\dagger,\qquad U^\dagger V^*+V^\dagger U^*=0,\qquad Z=V^*U^{-*}.
$$

Hence,

$$
\hat{\boldsymbol\beta}|\Phi\rangle=(U^\dagger Z+V^\dagger)\hat{\boldsymbol c}^\dagger|\Phi\rangle=(U^\dagger V^*U^{-*}+V^\dagger)\hat{\boldsymbol c}^\dagger|\Phi\rangle=0.
$$

The Thouless state therefore satisfies the quasiparticle-vacuum condition.

Expanding the exponential gives

$$
|\Phi\rangle=\nu\sum_{p=0}^{\lfloor M/2\rfloor}\frac{\hat A^p}{p!}|0\rangle=\sum_{n=0}^{\lfloor M/2\rfloor}|\Phi^{(2n)}\rangle.
$$

Since $\hat A$ creates two particles, every odd-particle component vanishes:

$$
\boxed{|\Phi^{(2n+1)}\rangle=0,\qquad 0\leq2n+1\leq M}.
$$

Fermionic creation operators obey

$$
\{\hat c_\alpha^\dagger,\hat c_\beta^\dagger\}=0,\qquad (\hat c_\alpha^\dagger)^2=0,
$$

so a single-particle orbital cannot occur more than once. The zero-particle component is

$$
|\Phi^{(0)}\rangle=\nu|0\rangle.
$$

The two-particle component is

$$
|\Phi^{(2)}\rangle=\nu\sum_{i<j}Z_{ij}\hat c_i^\dagger\hat c_j^\dagger|0\rangle.
$$

The four-particle component is

$$
|\Phi^{(4)}\rangle=\nu\sum_{i<j<k<l}\left(Z_{ij}Z_{kl}-Z_{ik}Z_{jl}+Z_{il}Z_{jk}\right)\hat c_i^\dagger\hat c_j^\dagger\hat c_k^\dagger\hat c_l^\dagger|0\rangle.
$$

The three terms are the complete pairings of four orbitals, with relative signs fixed by fermionic exchange.

For an ordered orbital set $\mathcal I=(i_1<\cdots<i_{2n})$ with $n\geq1$, let $Z_{\mathcal I}$ be the corresponding principal submatrix of $Z$. The general $2n$-particle component is

$$
|\Phi^{(2n)}\rangle=\nu\sum_{|\mathcal I|=2n}\operatorname{pf}(Z_{\mathcal I})\hat c_{i_1}^\dagger\cdots\hat c_{i_{2n}}^\dagger|0\rangle.
$$

The lowest-order Pfaffians, including the empty-set convention, are

$$
\operatorname{pf}(Z_\varnothing)=1,
$$

$$
\operatorname{pf}(Z_{\{i,j\}})=Z_{ij},
$$

$$
\operatorname{pf}(Z_{\{i,j,k,l\}})=Z_{ij}Z_{kl}-Z_{ik}Z_{jl}+Z_{il}Z_{jk}.
$$

Thus, the complete quasiparticle vacuum is

$$
\boxed{|\Phi\rangle=\nu\sum_{n=0}^{\lfloor M/2\rfloor}\sum_{|\mathcal I|=2n}\operatorname{pf}(Z_{\mathcal I})\hat c_{i_1}^\dagger\cdots\hat c_{i_{2n}}^\dagger|0\rangle}.
$$

The Pfaffian is therefore the natural fixed-particle-number coefficient of the fermionic pair amplitudes.

The Bogoliubov relations give the norm

$$
\langle\Phi|\Phi\rangle=|\nu|^2\sqrt{\det(I_M+Z^\dagger Z)}=|\nu|^2\sqrt{\det(U^{-T}U^{-*})}=\frac{|\nu|^2}{|\det U|}.
$$

Normalization requires

$$
\boxed{|\nu|^2=|\det U|}.
$$

This condition fixes only the magnitude of $\nu$; its phase is the overall phase of $|\Phi\rangle$.

If $\det U=0$, then

$$
\nu=\langle0|\Phi\rangle=0.
$$

The Thouless matrix relative to $|0\rangle$ is then undefined, although the quasiparticle-vacuum condition remains valid. One must choose another reference vacuum or use a Pfaffian expression constructed directly from $U$ and $V$.

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
