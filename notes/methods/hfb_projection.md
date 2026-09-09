# HFB Quantum-Number Projection

### Configurations and Projected States

Use the normalized vacuum and Bogoliubov convention in [HFB Pfaffian Algorithm](hfb_pfaffian.md). Let $\kappa=(\mu_1,\ldots,\mu_{n_\kappa})$ label an ordered list of distinct quasiparticle modes. Define

$$
|\Phi_\kappa\rangle=\hat\beta_{\mu_1}^\dagger\cdots\hat\beta_{\mu_{n_\kappa}}^\dagger|\Phi\rangle,\qquad \langle\Phi_\kappa|=\langle\Phi|\hat\beta_{\mu_{n_\kappa}}\cdots\hat\beta_{\mu_1}.
$$

All configurations share the same vacuum; the zero-quasiparticle configuration is $|\Phi_0\rangle=|\Phi\rangle$.

**Particle-number projection.** For a nonnegative integer $N$,

$$
\hat N=\sum_{\alpha=1}^{N_{\mathrm{sp}}}\hat c_\alpha^\dagger\hat c_\alpha,\qquad \hat P^N=\frac{1}{2\pi}\int_0^{2\pi}e^{i\varphi(\hat N-N)}\,\mathrm d\varphi.
$$

$$
(\hat P^N)^\dagger=\hat P^N,\qquad \hat P^{N_1}\hat P^{N_2}=\delta_{N_1N_2}\hat P^{N_1}.
$$

**Angular-momentum projection.** For Euler angles $\Omega=(\alpha,\beta,\gamma)$,

$$
\hat R(\Omega)=e^{-i\alpha\hat I_z/\hbar}e^{-i\beta\hat I_y/\hbar}e^{-i\gamma\hat I_z/\hbar},\qquad D^I_{MK}(\Omega)=\langle IM|\hat R(\Omega)|IK\rangle.
$$

Use the integration convention in [SO(3) and SU(2)](../group_theory/so3_su2.md#wigner-d-functions):

$$
\mathrm d\Omega=\sin\beta\,\mathrm d\alpha\,\mathrm d\beta\,\mathrm d\gamma,\qquad \alpha\in[0,2\pi),\quad\beta\in[0,\pi],\quad\gamma\in[0,\gamma_{\max}),\qquad \mathcal V_\Omega=4\pi\gamma_{\max}.
$$

Here $\gamma_{\max}=2\pi$ for integer angular momentum and $\gamma_{\max}=4\pi$ for full $SU(2)$ integration.

$$
\hat P^I_{MK}=\frac{2I+1}{\mathcal V_\Omega}\int\mathrm d\Omega\,D^{I*}_{MK}(\Omega)\hat R(\Omega).
$$

$$
(\hat P^I_{MK})^\dagger=\hat P^I_{KM},\qquad (\hat P^{I_1}_{M_1K_1})^\dagger\hat P^{I_2}_{M_2K_2}=\delta_{I_1I_2}\delta_{M_1M_2}\hat P^{I_1}_{K_1K_2}.
$$

**Projected states.** Since $[\hat N,\hat R(\Omega)]=0$, the two projections commute:

$$
[\hat P^N,\hat P^I_{MK}]=0.
$$

Define the unnormalized projected state by

$$
\boxed{|\Psi^{I;N}_{M,K;\kappa}\rangle=\hat P^I_{MK}\hat P^N|\Phi_\kappa\rangle=\frac{2I+1}{2\pi\mathcal V_\Omega}\int\mathrm d\Omega\int_0^{2\pi}\mathrm d\varphi\,D^{I*}_{MK}(\Omega)e^{-iN\varphi}\hat R(\Omega)e^{i\varphi\hat N}|\Phi_\kappa\rangle}.
$$

Here $I$ is total angular momentum, $M$ its laboratory-frame projection, and $K$ the projector’s second index. No axial symmetry is assumed.

### Hamiltonian and Projected Matrix Elements

For unsymmetrized two-body matrix elements, write

$$
\hat H=\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}\langle\alpha|\hat h|\beta\rangle\hat c_\alpha^\dagger\hat c_\beta+\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}\langle\alpha\otimes\beta|\hat v|\gamma\otimes\delta\rangle\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma.
$$

Assume rotational invariance and particle-number conservation:

$$
[\hat H,\hat R(\Omega)]=0,\qquad [\hat H,\hat N]=0.
$$

For a common laboratory projection $M$, the projector product relations give

$$
\langle\Psi^{I_1;N_1}_{M,K_1;\kappa_1}|\hat H|\Psi^{I_2;N_2}_{M,K_2;\kappa_2}\rangle=\langle\Phi_{\kappa_1}|\hat P^{N_1}(\hat P^{I_1}_{MK_1})^\dagger\hat H\hat P^{I_2}_{MK_2}\hat P^{N_2}|\Phi_{\kappa_2}\rangle=\delta_{N_1N_2}\delta_{I_1I_2}\langle\Phi_{\kappa_1}|\hat H\hat P^{N_1}\hat P^{I_1}_{K_1K_2}|\Phi_{\kappa_2}\rangle.
$$

For fixed $I,N$, define the Hamiltonian matrix

$$
\mathcal H^{I;N}_{K_1,\kappa_1;K_2,\kappa_2}\equiv\langle\Phi_{\kappa_1}|\hat H\hat P^N\hat P^I_{K_1K_2}|\Phi_{\kappa_2}\rangle=\frac{2I+1}{2\pi\mathcal V_\Omega}\int\mathrm d\Omega\int_0^{2\pi}\mathrm d\varphi\,D^{I*}_{K_1K_2}(\Omega)e^{-iN\varphi}\langle\Phi_{\kappa_1}|\hat H\hat R(\Omega)e^{i\varphi\hat N}|\Phi_{\kappa_2}\rangle.
$$

The norm matrix is

$$
\mathcal N^{I;N}_{K_1,\kappa_1;K_2,\kappa_2}\equiv\langle\Phi_{\kappa_1}|\hat P^N\hat P^I_{K_1K_2}|\Phi_{\kappa_2}\rangle=\frac{2I+1}{2\pi\mathcal V_\Omega}\int\mathrm d\Omega\int_0^{2\pi}\mathrm d\varphi\,D^{I*}_{K_1K_2}(\Omega)e^{-iN\varphi}\langle\Phi_{\kappa_1}|\hat R(\Omega)e^{i\varphi\hat N}|\Phi_{\kappa_2}\rangle.
$$

Define the mixed state with coefficients $f_{K,\kappa}$ by

$$
|\Psi^{I;N}_M\rangle=\sum_{K,\kappa}f_{K,\kappa}|\Psi^{I;N}_{M,K;\kappa}\rangle,\qquad f^\dagger\mathcal N^{I;N}f=1.
$$

Stationarity of the energy gives

$$
\boxed{\sum_{K_2,\kappa_2}\left(\mathcal H^{I;N}_{K_1,\kappa_1;K_2,\kappa_2}-E\,\mathcal N^{I;N}_{K_1,\kappa_1;K_2,\kappa_2}\right)f_{K_2,\kappa_2}=0}.
$$

### Rotated Quasiparticle Operators and Vacua

Define the single-particle rotation matrix by

$$
D_{\alpha\beta}(\Omega)\equiv\langle\alpha|\hat R(\Omega)|\beta\rangle.
$$

Assume the single-particle space is closed under rotations. Then

$$
\hat R(\Omega)\hat c_\alpha^\dagger\hat R^\dagger(\Omega)=\sum_\beta D_{\beta\alpha}(\Omega)\hat c_\beta^\dagger,\qquad \hat R(\Omega)\hat c_\alpha\hat R^\dagger(\Omega)=\sum_\beta D_{\beta\alpha}^*(\Omega)\hat c_\beta.
$$

The rotated quasiparticle annihilation operators are

$$
\hat{\boldsymbol\beta}(\Omega)\equiv\hat R(\Omega)\hat{\boldsymbol\beta}\hat R^\dagger(\Omega)=U^\dagger(\Omega)\hat{\boldsymbol c}+V^\dagger(\Omega)\hat{\boldsymbol c}^*.
$$

Taking the Hermitian adjoint gives

$$
\hat{\boldsymbol\beta}^\dagger(\Omega)\equiv\hat R(\Omega)\hat{\boldsymbol\beta}^\dagger\hat R^\dagger(\Omega)=\hat{\boldsymbol c}^\dagger U(\Omega)+\hat{\boldsymbol c}^TV(\Omega).
$$

The transformed Bogoliubov matrices follow from the particle rotation:

$$
\boxed{U(\Omega)=D(\Omega)U,\qquad V(\Omega)=D^*(\Omega)V}.
$$

The rotated vacuum satisfies

$$
|\Phi(\Omega)\rangle\equiv\hat R(\Omega)|\Phi\rangle,\qquad \hat{\boldsymbol\beta}(\Omega)|\Phi(\Omega)\rangle=\hat R(\Omega)\hat{\boldsymbol\beta}|\Phi\rangle=0.
$$

For a multiquasiparticle configuration,

$$
|\Phi_\kappa(\Omega)\rangle\equiv\hat R(\Omega)|\Phi_\kappa\rangle=\hat\beta_{\mu_1}^\dagger(\Omega)\cdots\hat\beta_{\mu_{n_\kappa}}^\dagger(\Omega)|\Phi(\Omega)\rangle.
$$

### One-Body Matrix Elements under Rotation

For single-particle matrix elements $\langle\alpha|\hat O|\beta\rangle$, define

$$
\hat O^{(1)}=\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}\langle\alpha|\hat O|\beta\rangle\hat c_\alpha^\dagger\hat c_\beta.
$$

Let $\kappa_1=(\mu_1,\ldots,\mu_{n_{\kappa_1}})$ and $\kappa_2=(\nu_1,\ldots,\nu_{n_{\kappa_2}})$ label the two configurations. Using the rotated quasiparticle operators,

$$
\langle\Phi_{\kappa_1}|\hat O^{(1)}\hat R(\Omega)|\Phi_{\kappa_2}\rangle=\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}\langle\alpha|\hat O|\beta\rangle\langle\Phi|\hat\beta_{\mu_{n_{\kappa_1}}}\cdots\hat\beta_{\mu_1}\hat c_\alpha^\dagger\hat c_\beta\hat\beta_{\nu_1}^\dagger(\Omega)\cdots\hat\beta_{\nu_{n_{\kappa_2}}}^\dagger(\Omega)|\Phi(\Omega)\rangle.
$$

For each $(\alpha,\beta,\Omega)$, use the ordered operator column

$$
\hat{\boldsymbol a}(\alpha,\beta)=\begin{pmatrix}\hat\beta_{\mu_{n_{\kappa_1}}}&\cdots&\hat\beta_{\mu_1}&\hat c_\alpha^\dagger&\hat c_\beta&\hat\beta_{\nu_1}^\dagger(\Omega)&\cdots&\hat\beta_{\nu_{n_{\kappa_2}}}^\dagger(\Omega)\end{pmatrix}^T,\qquad L=n_{\kappa_1}+n_{\kappa_2}+2.
$$

Assume $n(\Omega)\equiv\langle\Phi|\Phi(\Omega)\rangle\ne0$. Define $S(\alpha,\beta;\Omega)\in\mathbb C^{L\times L}$ by

$$
S_{ij}(\alpha,\beta;\Omega)\equiv\frac{\langle\Phi|\hat a_i(\alpha,\beta)\hat a_j(\alpha,\beta)|\Phi(\Omega)\rangle}{n(\Omega)}\quad(i<j),\qquad S(\alpha,\beta;\Omega)=-S^T(\alpha,\beta;\Omega).
$$

Evaluate its entries using [Transition Densities](hfb_pfaffian.md#transition-densities), with left matrices $(U,V)$ and right matrices $(U(\Omega),V(\Omega))$.

For even $n_{\kappa_1}+n_{\kappa_2}$, the generalized Wick theorem gives

$$
\boxed{\langle\Phi_{\kappa_1}|\hat O^{(1)}\hat R(\Omega)|\Phi_{\kappa_2}\rangle=n(\Omega)\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}\langle\alpha|\hat O|\beta\rangle\operatorname{pf}[S(\alpha,\beta;\Omega)]}.
$$

For odd $n_{\kappa_1}+n_{\kappa_2}$, the matrix element vanishes.

For $\kappa_1=\kappa_2=0$, $S_{12}(\alpha,\beta;\Omega)=\rho_{\beta\alpha}(\Omega)$, recovering

$$
\langle\Phi|\hat O^{(1)}\hat R(\Omega)|\Phi\rangle=n(\Omega)\operatorname{Tr}[O\rho(\Omega)].
$$

### Two-Body Matrix Elements under Rotation

For unsymmetrized two-body matrix elements, define

$$
\hat O^{(2)}=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}\langle\alpha\otimes\beta|\hat O|\gamma\otimes\delta\rangle\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma.
$$

Use the configurations and operator ordering defined in [One-Body Matrix Elements under Rotation](#one-body-matrix-elements-under-rotation). Then

$$
\langle\Phi_{\kappa_1}|\hat O^{(2)}\hat R(\Omega)|\Phi_{\kappa_2}\rangle=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}\langle\alpha\otimes\beta|\hat O|\gamma\otimes\delta\rangle\langle\Phi|\hat\beta_{\mu_{n_{\kappa_1}}}\cdots\hat\beta_{\mu_1}\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma\hat\beta_{\nu_1}^\dagger(\Omega)\cdots\hat\beta_{\nu_{n_{\kappa_2}}}^\dagger(\Omega)|\Phi(\Omega)\rangle.
$$

For each $(\alpha,\beta,\gamma,\delta,\Omega)$, use the ordered operator column

$$
\hat{\boldsymbol a}(\alpha,\beta,\gamma,\delta)=\begin{pmatrix}\hat\beta_{\mu_{n_{\kappa_1}}}&\cdots&\hat\beta_{\mu_1}&\hat c_\alpha^\dagger&\hat c_\beta^\dagger&\hat c_\delta&\hat c_\gamma&\hat\beta_{\nu_1}^\dagger(\Omega)&\cdots&\hat\beta_{\nu_{n_{\kappa_2}}}^\dagger(\Omega)\end{pmatrix}^T,\qquad L=n_{\kappa_1}+n_{\kappa_2}+4.
$$

Assume $n(\Omega)\ne0$. Define $S(\alpha,\beta,\gamma,\delta;\Omega)\in\mathbb C^{L\times L}$ by

$$
S_{ij}(\alpha,\beta,\gamma,\delta;\Omega)\equiv\frac{\langle\Phi|\hat a_i(\alpha,\beta,\gamma,\delta)\hat a_j(\alpha,\beta,\gamma,\delta)|\Phi(\Omega)\rangle}{n(\Omega)}\quad(i<j),\qquad S(\alpha,\beta,\gamma,\delta;\Omega)=-S^T(\alpha,\beta,\gamma,\delta;\Omega).
$$

Evaluate its entries using [Transition Densities](hfb_pfaffian.md#transition-densities), with left matrices $(U,V)$ and right matrices $(U(\Omega),V(\Omega))$.

For even $n_{\kappa_1}+n_{\kappa_2}$, the generalized Wick theorem gives

$$
\boxed{\langle\Phi_{\kappa_1}|\hat O^{(2)}\hat R(\Omega)|\Phi_{\kappa_2}\rangle=\frac{n(\Omega)}2\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}\langle\alpha\otimes\beta|\hat O|\gamma\otimes\delta\rangle\operatorname{pf}[S(\alpha,\beta,\gamma,\delta;\Omega)]}.
$$

For odd $n_{\kappa_1}+n_{\kappa_2}$, the matrix element vanishes.

For $\kappa_1=\kappa_2=0$,

$$
\operatorname{pf}[S(\alpha,\beta,\gamma,\delta;\Omega)]=\rho_{\gamma\alpha}(\Omega)\rho_{\delta\beta}(\Omega)-\rho_{\delta\alpha}(\Omega)\rho_{\gamma\beta}(\Omega)-\bar\kappa_{\alpha\beta}(\Omega)\kappa_{\gamma\delta}(\Omega),
$$

recovering

$$
\langle\Phi|\hat O^{(2)}\hat R(\Omega)|\Phi\rangle=\frac{n(\Omega)}2\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}\langle\alpha\otimes\beta|\hat O|\gamma\otimes\delta\rangle\left[\rho_{\gamma\alpha}(\Omega)\rho_{\delta\beta}(\Omega)-\rho_{\delta\alpha}(\Omega)\rho_{\gamma\beta}(\Omega)-\bar\kappa_{\alpha\beta}(\Omega)\kappa_{\gamma\delta}(\Omega)\right].
$$

### Irreducible Tensor Operators and Projected Matrix Elements

Let $\hat T_{\lambda\mu}$ be an irreducible spherical tensor of rank $\lambda$ and component $\mu$, conserving particle number:

$$
\hat R(\Omega)\hat T_{\lambda\mu}\hat R^\dagger(\Omega)=\sum_\nu D^\lambda_{\nu\mu}(\Omega)\hat T_{\lambda\nu},\qquad [\hat N,\hat T_{\lambda\mu}]=0.
$$

Using the angular-momentum projector definition,

$$
(\hat P^{I_1}_{M_1K_1})^\dagger\hat T_{\lambda\mu}=\frac{2I_1+1}{\mathcal V_\Omega}\sum_\nu\int\mathrm d\Omega\,D^{I_1*}_{K_1M_1}(\Omega)D^\lambda_{\nu\mu}(\Omega)\hat T_{\lambda\nu}\hat R(\Omega).
$$

The Wigner-$D$ product expansion gives

$$
(\hat P^{I_1}_{M_1K_1})^\dagger\hat T_{\lambda\mu}=\sum_{Jmk\nu}\langle Jm;\lambda\mu|I_1M_1\rangle\langle Jk;\lambda\nu|I_1K_1\rangle\hat T_{\lambda\nu}\hat P^J_{km},
$$

where the bracketed coefficients are Clebsch–Gordan coefficients. Multiplying by $\hat P^{I_2}_{M_2K_2}$ and using $k+\nu=K_1$ gives

$$
\boxed{(\hat P^{I_1}_{M_1K_1})^\dagger\hat T_{\lambda\mu}\hat P^{I_2}_{M_2K_2}=\langle I_2M_2;\lambda\mu|I_1M_1\rangle\sum_\nu\langle I_2,K_1-\nu;\lambda\nu|I_1K_1\rangle\hat T_{\lambda\nu}\hat P^{I_2}_{K_1-\nu,K_2}}.
$$

Including particle-number projection,

$$
\langle\Psi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat T_{\lambda\mu}|\Psi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle=\delta_{N_1N_2}\langle I_2M_2;\lambda\mu|I_1M_1\rangle\sum_\nu\langle I_2,K_1-\nu;\lambda\nu|I_1K_1\rangle\langle\Phi_{\kappa_1}|\hat T_{\lambda\nu}\hat P^{I_2}_{K_1-\nu,K_2}\hat P^{N_2}|\Phi_{\kappa_2}\rangle.
$$

Use the Wigner–Eckart convention

$$
\langle\Psi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat T_{\lambda\mu}|\Psi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle=\frac{\langle I_2M_2;\lambda\mu|I_1M_1\rangle}{\sqrt{2I_1+1}}\langle\Psi^{I_1;N_1}_{K_1;\kappa_1}\Vert\hat T_\lambda\Vert\Psi^{I_2;N_2}_{K_2;\kappa_2}\rangle.
$$

Comparison yields

$$
\boxed{\langle\Psi^{I_1;N_1}_{K_1;\kappa_1}\Vert\hat T_\lambda\Vert\Psi^{I_2;N_2}_{K_2;\kappa_2}\rangle=\delta_{N_1N_2}\sqrt{2I_1+1}\sum_\nu\langle I_2,K_1-\nu;\lambda\nu|I_1K_1\rangle\langle\Phi_{\kappa_1}|\hat T_{\lambda\nu}\hat P^{I_2}_{K_1-\nu,K_2}\hat P^{N_2}|\Phi_{\kappa_2}\rangle}.
$$

The remaining matrix element is

$$
\langle\Phi_{\kappa_1}|\hat T_{\lambda\nu}\hat P^{I_2}_{K_1-\nu,K_2}\hat P^{N_2}|\Phi_{\kappa_2}\rangle=\frac{2I_2+1}{2\pi\mathcal V_\Omega}\int\mathrm d\Omega\int_0^{2\pi}\mathrm d\varphi\,D^{I_2*}_{K_1-\nu,K_2}(\Omega)e^{-iN_2\varphi}\langle\Phi_{\kappa_1}|\hat T_{\lambda\nu}\hat R(\Omega)e^{i\varphi\hat N}|\Phi_{\kappa_2}\rangle.
$$
