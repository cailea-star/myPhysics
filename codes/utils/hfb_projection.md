# HFB Quantum-Number Projection

### Configurations and Projected States

Use the normalized vacuum and Bogoliubov convention in [HFB Pfaffian Algorithm](hfb_pfaffian.md). Let $\kappa=(\mu_1,\ldots,\mu_{n_\kappa})$ label an ordered list of distinct quasiparticle modes. Define

$$
|\Phi_\kappa\rangle=\hat\beta_{\mu_1}^*\cdots\hat\beta_{\mu_{n_\kappa}}^*|\Phi\rangle,\qquad \langle\Phi_\kappa|=\langle\Phi|\hat\beta_{\mu_{n_\kappa}}\cdots\hat\beta_{\mu_1}.
$$

All configurations share the same vacuum; the zero-quasiparticle configuration is $|\Phi_0\rangle=|\Phi\rangle$.

**Particle-number projection.** For a nonnegative integer $N$,

$$
\hat N=\sum_{\alpha=1}^{N_{\mathrm{sp}}}\hat c_\alpha^*\hat c_\alpha,\qquad \hat P^N=\frac{1}{2\pi}\int_0^{2\pi}e^{-i\varphi(\hat N-N)}\,\mathrm d\varphi.
$$

$$
(\hat P^N)^\dagger=\hat P^N,\qquad \hat P^{N_1}\hat P^{N_2}=\delta_{N_1N_2}\hat P^{N_1}.
$$

**Angular-momentum projection.** For Euler angles $\Omega=(\alpha,\beta,\gamma)$,

$$
\hat R(\Omega)=e^{-i\alpha\hat I_z/\hbar}e^{-i\beta\hat I_y/\hbar}e^{-i\gamma\hat I_z/\hbar},\qquad D^I_{MK}(\Omega)=\langle IM|\hat R(\Omega)|IK\rangle.
$$

Use the integration convention in [SO(3) and SU(2)](../../notes/group_theory/so3_su2.md#wigner-d-functions):

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
\boxed{|\Phi^{I;N}_{M,K;\kappa}\rangle=\hat P^I_{MK}\hat P^N|\Phi_\kappa\rangle=\frac{2I+1}{2\pi\mathcal V_\Omega}\int\mathrm d\Omega\int_0^{2\pi}\mathrm d\varphi\,D^{I*}_{MK}(\Omega)e^{iN\varphi}\hat R(\Omega)e^{-i\varphi\hat N}|\Phi_\kappa\rangle}.
$$

Here $I$ is total angular momentum, $M$ its laboratory-frame projection, and $K$ the projector’s second index. No axial symmetry is assumed.

### Hamiltonian and Projected Matrix Elements

For unsymmetrized two-body matrix elements, write

$$
\hat H=\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}\langle\alpha|\hat h|\beta\rangle\hat c_\alpha^*\hat c_\beta+\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}\langle\alpha\otimes\beta|\hat v|\gamma\otimes\delta\rangle\hat c_\alpha^*\hat c_\beta^*\hat c_\delta\hat c_\gamma.
$$

Assume rotational invariance and particle-number conservation:

$$
[\hat H,\hat R(\Omega)]=0,\qquad [\hat H,\hat N]=0.
$$

For a common laboratory projection $M$, the projector product relations give

$$
\langle\Phi^{I_1;N_1}_{M,K_1;\kappa_1}|\hat H|\Phi^{I_2;N_2}_{M,K_2;\kappa_2}\rangle=\langle\Phi_{\kappa_1}|\hat P^{N_1}(\hat P^{I_1}_{MK_1})^\dagger\hat H\hat P^{I_2}_{MK_2}\hat P^{N_2}|\Phi_{\kappa_2}\rangle=\delta_{N_1N_2}\delta_{I_1I_2}\langle\Phi_{\kappa_1}|\hat H\hat P^{N_1}\hat P^{I_1}_{K_1K_2}|\Phi_{\kappa_2}\rangle.
$$

For fixed $I,N$, define the Hamiltonian matrix

$$
\mathcal H^{I;N}_{K_1,\kappa_1;K_2,\kappa_2}\equiv\langle\Phi_{\kappa_1}|\hat H\hat P^N\hat P^I_{K_1K_2}|\Phi_{\kappa_2}\rangle=\frac{2I+1}{2\pi\mathcal V_\Omega}\int\mathrm d\Omega\int_0^{2\pi}\mathrm d\varphi\,D^{I*}_{K_1K_2}(\Omega)e^{iN\varphi}\langle\Phi_{\kappa_1}|\hat H\hat R(\Omega)e^{-i\varphi\hat N}|\Phi_{\kappa_2}\rangle.
$$

The norm matrix is

$$
\mathcal N^{I;N}_{K_1,\kappa_1;K_2,\kappa_2}\equiv\langle\Phi_{\kappa_1}|\hat P^N\hat P^I_{K_1K_2}|\Phi_{\kappa_2}\rangle=\frac{2I+1}{2\pi\mathcal V_\Omega}\int\mathrm d\Omega\int_0^{2\pi}\mathrm d\varphi\,D^{I*}_{K_1K_2}(\Omega)e^{iN\varphi}\langle\Phi_{\kappa_1}|\hat R(\Omega)e^{-i\varphi\hat N}|\Phi_{\kappa_2}\rangle.
$$

Define the mixed state with coefficients $f_{K,\kappa}$ by

$$
|\Phi^{I;N}_M\rangle=\sum_{K,\kappa}f_{K,\kappa}|\Phi^{I;N}_{M,K;\kappa}\rangle,\qquad f^\dagger\mathcal N^{I;N}f=1.
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
\hat R(\Omega)\hat c_\alpha^*\hat R^\dagger(\Omega)=\sum_\beta D_{\beta\alpha}(\Omega)\hat c_\beta^*,\qquad \hat R(\Omega)\hat c_\alpha\hat R^\dagger(\Omega)=\sum_\beta D_{\beta\alpha}^*(\Omega)\hat c_\beta.
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
|\Phi_\kappa(\Omega)\rangle\equiv\hat R(\Omega)|\Phi_\kappa\rangle=\hat\beta_{\mu_1}^*(\Omega)\cdots\hat\beta_{\mu_{n_\kappa}}^*(\Omega)|\Phi(\Omega)\rangle.
$$

Apply the [Generalized Wick Theorem](hfb_pfaffian.md#generalized-wick-theorem) with left matrices $(U,V)$ and right matrices $(U(\Omega),V(\Omega))$ to calculate one-body and two-body matrix elements, preserving the operator order. Normalized contractions require $\langle\Phi|\Phi(\Omega)\rangle\ne0$.

### Irreducible Tensor Operators and Projected Matrix Elements

Let $\hat T_{\lambda\mu}$ be an irreducible spherical tensor conserving particle number:

$$
\hat R(\Omega)\hat T_{\lambda\mu}\hat R^\dagger(\Omega)=\sum_{\mu'} D^\lambda_{\mu'\mu}(\Omega)\hat T_{\lambda\mu'},\qquad [\hat N,\hat T_{\lambda\mu}]=0.
$$

Expanding both angular-momentum projectors gives

$$
\hat P^{I_1}_{K_1M_1}\hat T_{\lambda\mu}\hat P^{I_2}_{M_2K_2}=\frac{(2I_1+1)(2I_2+1)}{\mathcal V_\Omega^2}\int\mathrm d\Omega_1\,\mathrm d\Omega_2\,D^{I_1*}_{K_1M_1}(\Omega_1)D^{I_2*}_{M_2K_2}(\Omega_2)\hat R(\Omega_1)\hat T_{\lambda\mu}\hat R(\Omega_2).
$$

Using the tensor transformation,

$$
\hat P^{I_1}_{K_1M_1}\hat T_{\lambda\mu}\hat P^{I_2}_{M_2K_2}=\frac{(2I_1+1)(2I_2+1)}{\mathcal V_\Omega^2}\sum_{\mu'}\int\mathrm d\Omega_1\,\mathrm d\Omega_2\,D^{I_1*}_{K_1M_1}(\Omega_1)D^\lambda_{\mu'\mu}(\Omega_1)D^{I_2*}_{M_2K_2}(\Omega_2)\hat T_{\lambda\mu'}\hat R(\Omega_1)\hat R(\Omega_2).
$$

Introduce the composite rotation. For fixed $\Omega_1$, invariance of the Haar measure implies

$$
\hat R(\Omega)=\hat R(\Omega_1)\hat R(\Omega_2),\qquad \Omega_2=\Omega_1^{-1}\Omega,\qquad \mathrm d\Omega_2=\mathrm d\Omega.
$$

The representation matrices satisfy

$$
D^{I_2*}_{M_2K_2}(\Omega_2)=\sum_k D^{I_2}_{kM_2}(\Omega_1)D^{I_2*}_{kK_2}(\Omega).
$$

Separating the integrals gives

$$
\hat P^{I_1}_{K_1M_1}\hat T_{\lambda\mu}\hat P^{I_2}_{M_2K_2}=\frac{(2I_1+1)(2I_2+1)}{\mathcal V_\Omega^2}\sum_{\mu',k}\left[\int\mathrm d\Omega_1\,D^{I_1*}_{K_1M_1}(\Omega_1)D^{I_2}_{kM_2}(\Omega_1)D^\lambda_{\mu'\mu}(\Omega_1)\right]\left[\int\mathrm d\Omega\,D^{I_2*}_{kK_2}(\Omega)\hat T_{\lambda\mu'}\hat R(\Omega)\right].
$$

The Clebsch–Gordan decomposition is

$$
D^{I_2}_{kM_2}(\Omega_1)D^\lambda_{\mu'\mu}(\Omega_1)=\sum_{IMK}\langle I_2k;\lambda\mu'|IK\rangle\langle I_2M_2;\lambda\mu|IM\rangle D^I_{KM}(\Omega_1).
$$

Orthogonality selects $I=I_1$, $K=K_1$, and $M=M_1$:

$$
\int\mathrm d\Omega_1\,D^{I_1*}_{K_1M_1}(\Omega_1)D^{I_2}_{kM_2}(\Omega_1)D^\lambda_{\mu'\mu}(\Omega_1)=\frac{\mathcal V_\Omega}{2I_1+1}\langle I_2k;\lambda\mu'|I_1K_1\rangle\langle I_2M_2;\lambda\mu|I_1M_1\rangle.
$$

The remaining integral defines $\hat P^{I_2}_{kK_2}$. Using $k+\mu'=K_1$ gives

$$
\boxed{\hat P^{I_1}_{K_1M_1}\hat T_{\lambda\mu}\hat P^{I_2}_{M_2K_2}=\langle I_2M_2;\lambda\mu|I_1M_1\rangle\sum_{\mu'}\langle I_2,K_1-\mu';\lambda\mu'|I_1K_1\rangle\hat T_{\lambda\mu'}\hat P^{I_2}_{K_1-\mu',K_2}}.
$$

Including particle-number projection,

$$
\langle\Phi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat T_{\lambda\mu}|\Phi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle=\delta_{N_1N_2}\langle I_2M_2;\lambda\mu|I_1M_1\rangle\sum_{\mu'}\langle I_2,K_1-\mu';\lambda\mu'|I_1K_1\rangle\langle\Phi_{\kappa_1}|\hat T_{\lambda\mu'}\hat P^{I_2}_{K_1-\mu',K_2}\hat P^{N_2}|\Phi_{\kappa_2}\rangle.
$$

Use the Wigner–Eckart convention

$$
\langle\Phi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat T_{\lambda\mu}|\Phi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle=\frac{\langle I_2M_2;\lambda\mu|I_1M_1\rangle}{\sqrt{2I_1+1}}\langle\Phi^{I_1;N_1}_{K_1;\kappa_1}\Vert\hat T_\lambda\Vert\Phi^{I_2;N_2}_{K_2;\kappa_2}\rangle.
$$

Comparison yields

$$
\boxed{\langle\Phi^{I_1;N_1}_{K_1;\kappa_1}\Vert\hat T_\lambda\Vert\Phi^{I_2;N_2}_{K_2;\kappa_2}\rangle=\delta_{N_1N_2}\sqrt{2I_1+1}\sum_{\mu'}\langle I_2,K_1-\mu';\lambda\mu'|I_1K_1\rangle\langle\Phi_{\kappa_1}|\hat T_{\lambda\mu'}\hat P^{I_2}_{K_1-\mu',K_2}\hat P^{N_2}|\Phi_{\kappa_2}\rangle}.
$$

The remaining matrix element is

$$
\langle\Phi_{\kappa_1}|\hat T_{\lambda\mu'}\hat P^{I_2}_{K_1-\mu',K_2}\hat P^{N_2}|\Phi_{\kappa_2}\rangle=\frac{2I_2+1}{2\pi\mathcal V_\Omega}\int\mathrm d\Omega\int_0^{2\pi}\mathrm d\varphi\,D^{I_2*}_{K_1-\mu',K_2}(\Omega)e^{iN_2\varphi}\langle\Phi_{\kappa_1}|\hat T_{\lambda\mu'}\hat R(\Omega)e^{-i\varphi\hat N}|\Phi_{\kappa_2}\rangle.
$$

### Projected Matrix Elements of Multipole Interactions

Use the multipole interaction and angular-momentum conventions in [Multipole Interactions and Matrix Elements](../../notes/quantum_mechanics/second_quantization_multipole.md#multipole-interactions-and-matrix-elements). Define

$$
\hat H_\lambda=\frac{\chi_\lambda}{2}\sum_{\mu=-\lambda}^{\lambda}\hat Q_{\lambda\mu}^\dagger\hat Q_{\lambda\mu},\qquad \hat Q_{\lambda\mu}=\sum_{\alpha,\beta}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\hat c_\alpha^*\hat c_\beta,\qquad \hat Q_{\lambda\mu}^\dagger=(-1)^\mu\hat Q_{\lambda,-\mu}.
$$

The multipole operator is

$$
\boxed{\hat Q_{\lambda\mu}=\sum_{\alpha,\beta}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\hat c_\alpha^*\hat c_\beta=\sum_{a,b}\frac{\langle a\Vert\hat q_\lambda\Vert b\rangle}{\sqrt{2\lambda+1}}[\hat c_a^*\otimes\tilde c_b]_{\lambda\mu}}.
$$

The matrix element between projected states is

$$
\langle\Phi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat Q_{\lambda\mu}|\Phi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle=\sum_{\alpha,\beta}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\langle\Phi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat c_\alpha^*\hat c_\beta|\Phi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle.
$$

The reduced matrix element between projected states is

$$
\boxed{\langle\Phi^{I_1;N_1}_{K_1;\kappa_1}\Vert\hat Q_\lambda\Vert\Phi^{I_2;N_2}_{K_2;\kappa_2}\rangle=\sum_{a,b}\frac{\langle a\Vert\hat q_\lambda\Vert b\rangle}{\sqrt{2\lambda+1}}\langle\Phi^{I_1;N_1}_{K_1;\kappa_1}\Vert[\hat c_a^*\otimes\tilde c_b]_\lambda\Vert\Phi^{I_2;N_2}_{K_2;\kappa_2}\rangle}.
$$

The interaction matrix element between projected states is

$$
\langle\Phi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat H_\lambda|\Phi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle=\frac{\chi_\lambda}{2}\sum_{\mu,\alpha,\beta,\gamma,\delta}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle^*\langle\gamma|\hat q_{\lambda\mu}|\delta\rangle\langle\Phi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat c_\beta^*\hat c_\alpha\hat c_\gamma^*\hat c_\delta|\Phi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle.
$$

For its reduced form, use the scalar operator

$$
\hat X^{(0)}_{\lambda;ab,cd}\equiv\sum_{\mu=-\lambda}^{\lambda}\left([\hat c_a^*\otimes\tilde c_b]_{\lambda\mu}\right)^\dagger[\hat c_c^*\otimes\tilde c_d]_{\lambda\mu}=\sum_{\mu,m_a,m_b,m_c,m_d}(-1)^{j_b-m_b+j_d-m_d}\langle j_a m_a;j_b,-m_b|\lambda\mu\rangle\langle j_c m_c;j_d,-m_d|\lambda\mu\rangle\hat c_{b m_b}^*\hat c_{a m_a}\hat c_{c m_c}^*\hat c_{d m_d}.
$$

Then

$$
\boxed{\langle\Phi^{I_1;N_1}_{K_1;\kappa_1}\Vert\hat H_\lambda\Vert\Phi^{I_2;N_2}_{K_2;\kappa_2}\rangle=\frac{\chi_\lambda}{2(2\lambda+1)}\sum_{a,b,c,d}\langle a\Vert\hat q_\lambda\Vert b\rangle^*\langle c\Vert\hat q_\lambda\Vert d\rangle\langle\Phi^{I_1;N_1}_{K_1;\kappa_1}\Vert\hat X^{(0)}_{\lambda;ab,cd}\Vert\Phi^{I_2;N_2}_{K_2;\kappa_2}\rangle}.
$$

### Projected Matrix Elements of Multipole Pairing Interactions

Use the pairing interaction and angular-momentum conventions in [Multipole Pairing Interactions and Matrix Elements](../../notes/quantum_mechanics/second_quantization_multipole.md#multipole-pairing-interactions-and-matrix-elements). Define

$$
\hat H_\lambda^{\mathrm{pair}}=-G_\lambda\sum_{\mu=-\lambda}^{\lambda}\hat P_{\lambda\mu}^\dagger\hat P_{\lambda\mu},\qquad \hat P_{\lambda\mu}^\dagger=\frac12\sum_{\alpha,\beta}(-1)^{j_b-m_b}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\hat c_{a m_a}^*\hat c_{b,-m_b}^*.
$$

The pair-creation operator is

$$
\boxed{\hat P_{\lambda\mu}^\dagger=\frac12\sum_{\alpha,\beta}(-1)^{j_b-m_b}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\hat c_{a m_a}^*\hat c_{b,-m_b}^*=\frac12\sum_{a,b}\frac{\langle a\Vert\hat q_\lambda\Vert b\rangle}{\sqrt{2\lambda+1}}[\hat c_a^*\otimes\hat c_b^*]_{\lambda\mu}}.
$$

The matrix element between projected states is

$$
\langle\Phi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat P_{\lambda\mu}^\dagger|\Phi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle=\frac12\sum_{\alpha,\beta}(-1)^{j_b-m_b}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\langle\Phi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat c_{a m_a}^*\hat c_{b,-m_b}^*|\Phi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle.
$$

The reduced matrix element between projected states is

$$
\boxed{\langle\Phi^{I_1;N_1}_{K_1;\kappa_1}\Vert\hat P_\lambda^\dagger\Vert\Phi^{I_2;N_2}_{K_2;\kappa_2}\rangle=\frac12\sum_{a,b}\frac{\langle a\Vert\hat q_\lambda\Vert b\rangle}{\sqrt{2\lambda+1}}\langle\Phi^{I_1;N_1}_{K_1;\kappa_1}\Vert[\hat c_a^*\otimes\hat c_b^*]_\lambda\Vert\Phi^{I_2;N_2}_{K_2;\kappa_2}\rangle}.
$$

These pair-addition matrix elements vanish unless $N_1=N_2+2$.

The interaction matrix element between projected states is

$$
\langle\Phi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat H_\lambda^{\mathrm{pair}}|\Phi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle=-\frac{G_\lambda}{4}\sum_{\mu,\alpha,\beta,\gamma,\delta}(-1)^{j_b-m_b+j_d-m_d}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\langle\gamma|\hat q_{\lambda\mu}|\delta\rangle^*\langle\Phi^{I_1;N_1}_{M_1,K_1;\kappa_1}|\hat c_{a m_a}^*\hat c_{b,-m_b}^*\hat c_{d,-m_d}\hat c_{c m_c}|\Phi^{I_2;N_2}_{M_2,K_2;\kappa_2}\rangle.
$$

For its reduced form, use the scalar operator

$$
\hat X^{(0),\mathrm{pair}}_{\lambda;ab,cd}\equiv\sum_{\mu=-\lambda}^{\lambda}[\hat c_a^*\otimes\hat c_b^*]_{\lambda\mu}\left([\hat c_c^*\otimes\hat c_d^*]_{\lambda\mu}\right)^\dagger=\sum_{\mu,m_a,m_b,m_c,m_d}\langle j_a m_a;j_b,-m_b|\lambda\mu\rangle\langle j_c m_c;j_d,-m_d|\lambda\mu\rangle\hat c_{a m_a}^*\hat c_{b,-m_b}^*\hat c_{d,-m_d}\hat c_{c m_c}.
$$

Then

$$
\boxed{\langle\Phi^{I_1;N_1}_{K_1;\kappa_1}\Vert\hat H_\lambda^{\mathrm{pair}}\Vert\Phi^{I_2;N_2}_{K_2;\kappa_2}\rangle=-\frac{G_\lambda}{4(2\lambda+1)}\sum_{a,b,c,d}\langle a\Vert\hat q_\lambda\Vert b\rangle\langle c\Vert\hat q_\lambda\Vert d\rangle^*\langle\Phi^{I_1;N_1}_{K_1;\kappa_1}\Vert\hat X^{(0),\mathrm{pair}}_{\lambda;ab,cd}\Vert\Phi^{I_2;N_2}_{K_2;\kappa_2}\rangle}.
$$
