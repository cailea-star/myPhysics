# Axial Projection

Use the projection operators and transition kernels derived in [HFB Pfaffian Algorithm](../../notes/methods/hfb_pfaffian.md). All quantities refer to a single fermion species, without an isospin index.

Let $|\Phi_1\rangle$ and $|\Phi_2\rangle$ be normalized HFB vacua in the same single-particle basis, with Bogoliubov matrices $(U_1,V_1)$ and $(U_2,V_2)$. Matrix elements use $\langle\Phi_1|$ on the left and $|\Phi_2\rangle$ on the right.

### Particle-Number Projection

Let $\hat N$ denote the particle-number operator and $N$ the target particle number. For the gauge angle $\varphi\in[0,2\pi)$,

$$
\hat P^N=\frac{1}{2\pi}\int_0^{2\pi}e^{-iN\varphi}e^{i\varphi\hat N}\,\mathrm d\varphi,\qquad |\Phi_i^N\rangle=\hat P^N|\Phi_i\rangle,\qquad i\in\{1,2\}.
$$

The projected overlap kernel is

$$
\boxed{\mathcal N^N_{12}=\langle\Phi_1|\hat P^N|\Phi_2\rangle=\frac{1}{2\pi}\int_0^{2\pi}e^{-iN\varphi}\langle\Phi_1|e^{i\varphi\hat N}|\Phi_2\rangle\,\mathrm d\varphi}.
$$

For a number-conserving operator $\hat O$ and $\mathcal N^N_{12}\ne0$,

$$
\frac{\langle\Phi_1^N|\hat O|\Phi_2^N\rangle}{\langle\Phi_1^N|\Phi_2^N\rangle}=\frac{1}{2\pi\mathcal N^N_{12}}\int_0^{2\pi}e^{-iN\varphi}\langle\Phi_1|\hat Oe^{i\varphi\hat N}|\Phi_2\rangle\,\mathrm d\varphi.
$$

For identical vacua, this ratio reduces to the projected expectation value.

### Angular-Momentum Projection

Use the conventions in [SO(3) and SU(2)](../../notes/group_theory/so3_su2.md#wigner-d-functions), with $\gamma_{\max}=2\pi$ for integer $J$ and $4\pi$ for half-odd-integer $J$.

$$
\hat P^J_{MK}=\frac{2J+1}{4\pi\gamma_{\max}}\int_0^{2\pi}\mathrm d\alpha\int_0^\pi\sin\beta\,\mathrm d\beta\int_0^{\gamma_{\max}}\mathrm d\gamma\,D^{J*}_{MK}(\alpha,\beta,\gamma)\hat R(\alpha,\beta,\gamma).
$$

For $i\in\{1,2\}$ and $M,K\in\{-J,\ldots,J\}$, define the unnormalized projected component

$$
|\Phi^J_{i,MK}\rangle=\hat P^J_{MK}|\Phi_i\rangle.
$$

It satisfies

$$
\hat J^2|\Phi^J_{i,MK}\rangle=\hbar^2J(J+1)|\Phi^J_{i,MK}\rangle,\qquad \hat J_z|\Phi^J_{i,MK}\rangle=\hbar M|\Phi^J_{i,MK}\rangle.
$$

At fixed $J,M$, the overlap kernel is

$$
\mathcal N^J_{1K_1,2K_2}=\langle\Phi^J_{1,MK_1}|\Phi^J_{2,MK_2}\rangle=\langle\Phi_1|\hat P^J_{K_1K_2}|\Phi_2\rangle=\frac{2J+1}{4\pi\gamma_{\max}}\int\mathrm d\Omega\,D^{J*}_{K_1K_2}(\Omega)\langle\Phi_1|\hat R(\Omega)|\Phi_2\rangle,
$$

where $\Omega=(\alpha,\beta,\gamma)$ and $\mathrm d\Omega=\mathrm d\alpha\,\sin\beta\,\mathrm d\beta\,\mathrm d\gamma$ over the preceding ranges.

$$
\hat R(\alpha,\beta,\gamma)=e^{-i\alpha\hat J_z/\hbar}e^{-i\beta\hat J_y/\hbar}e^{-i\gamma\hat J_z/\hbar},\qquad D^J_{MK}(\alpha,\beta,\gamma)=e^{-iM\alpha}d^J_{MK}(\beta)e^{-iK\gamma}.
$$

### Transformation Matrices

Let $\{|\phi_a\rangle\}_{a=1}^{N_{\mathrm{sp}}}$ be an orthonormal axial basis containing both signs of $\Omega_a$, with $\hat J_z|\phi_a\rangle=\hbar\Omega_a|\phi_a\rangle$. For $g=(\varphi,\Omega)$, $\Omega=(\alpha,\beta,\gamma)$, define

$$
\hat G(g)=e^{i\varphi\hat N}\hat R(\Omega).
$$

Using [Axial-Basis Rotation Matrix](axial_rotation.md#axial-basis-rotation-matrix), the restricted single-particle matrix is

$$
D_{ab}(g)=e^{i\varphi}e^{-i\alpha\Omega_a}d_{ab}(\beta)e^{-i\gamma\Omega_b},\qquad d_{ab}(\beta)=\langle\phi_a|e^{-i\beta\hat J_y/\hbar}|\phi_b\rangle.
$$

If the retained space is closed under rotation, $D(g)$ is unitary and

$$
\hat G(g)\hat c_b^*\hat G(g)^\dagger=\sum_{a=1}^{N_{\mathrm{sp}}}D_{ab}(g)\hat c_a^*.
$$

For invertible $U_1,U_2$, define

$$
Z_1 = (V_1U_1^{-1})^*,\qquad Z_2 = (V_2U_2^{-1})^*.
$$

For invertible $D(g)$, use

$$
\boxed{U_2(g)=D(g)^{-\dagger}U_2,\qquad V_2(g)=D(g)^*V_2,\qquad Z_2(g)=D(g)Z_2D(g)^T}.
$$

For unitary $D(g)$, $D(g)^{-\dagger}=D(g)$. Otherwise, $U_2(g),V_2(g)$ are auxiliary matrices, and $Z_2(g)$ describes the transformed component with no particles outside the retained space.

The overlap auxiliary matrix is

$$
\boxed{\mathcal A(g)=U_1^TU_2(g)^*+V_1^TV_2(g)^*}.
$$

### Overlap and Transition Densities

Let $\nu_a\equiv\langle0|\Phi_a\rangle$ for $a=1,2$, $I_{N_{\mathrm{sp}}}$ denote the identity matrix, and $s_{N_{\mathrm{sp}}}\equiv(-1)^{N_{\mathrm{sp}}(N_{\mathrm{sp}}+1)/2}$. Define the overlap kernel by

$$
n_{12}(g)\equiv\langle\Phi_1|\hat G(g)|\Phi_2\rangle.
$$

Substituting the transformed Thouless matrix into [Quasiparticle Vacuum Overlap](../../notes/methods/hfb_pfaffian.md#quasiparticle-vacuum-overlap) gives

$$
\boxed{n_{12}(g)=\nu_1^*\nu_2s_{N_{\mathrm{sp}}}\operatorname{pf}\begin{pmatrix}Z_2(g)&-I_{N_{\mathrm{sp}}}\\I_{N_{\mathrm{sp}}}&-Z_1^*\end{pmatrix}}.
$$

Normalization fixes $|\nu_a|^2=|\det U_a|$, but not the relative vacuum phase. The Pfaffian fixes the overlap for a specified single-particle ordering and vacuum-phase convention.

For $n_{12}(g)\ne0$, the formulas in [Transition Densities](../../notes/methods/hfb_pfaffian.md#transition-densities) apply with

$$
\mathcal A(g)=U_1^TU_2(g)^*+V_1^TV_2(g)^*,
$$

$$
\boxed{\rho(g)=V_2(g)^*\mathcal A(g)^{-1}V_1^T},
$$

$$
\boxed{\kappa(g)=V_2(g)^*\mathcal A(g)^{-1}U_1^T,\qquad \bar\kappa(g)=U_2(g)^*\mathcal A(g)^{-1}V_1^T}.
$$

The overlap formula remains valid at $n_{12}(g)=0$, where normalized transition densities are undefined.

### Time-Reversal Symmetry

For the two axial, unblocked reference vacua, impose

$$
\hat{\mathcal T}|\Phi_a\rangle=|\Phi_a\rangle,\qquad a=1,2.
$$

Use the basis phases and paired quasiparticle columns defined in [Kramers-Representative Formulation](../../notes/methods/hfb_axial.md#kramers-representative-formulation). With particle rows ordered as $(\alpha,\bar\alpha)$ and quasiparticle columns as $(\mu,\bar\mu)$,

$$
\boxed{U_a=\begin{pmatrix}U_a^+&0\\0&\operatorname{diag}(\eta_\alpha)(U_a^+)^*\end{pmatrix},\qquad V_a=\begin{pmatrix}0&-\operatorname{diag}(\eta_\alpha)(V_a^-)^*\\V_a^-&0\end{pmatrix}}.
$$

Here $\Omega_\alpha>0$, $\Omega_{\bar\alpha}=-\Omega_\alpha$, and $\eta_\alpha=\pm1$ are the time-reversal phases. The independent amplitudes $U_a^+,V_a^-$ are generally complex.

For invertible $U_a$, axial symmetry restricts the Thouless pairs to opposite angular-momentum projections:

$$
(Z_a)_{\alpha\beta}\ne0\quad\Longrightarrow\quad\Omega_\alpha+\Omega_\beta=0.
$$

Each pair therefore carries zero $J_z$, giving

$$
\boxed{\hat J_z|\Phi_a\rangle=0,\qquad K_a=0}.
$$

This does not imply $\hat{\boldsymbol J}^{\,2}|\Phi_a\rangle=0$. The overlap under rotation reduces to

$$
\boxed{\langle\Phi_1|\hat R(\alpha,\beta,\gamma)|\Phi_2\rangle=\langle\Phi_1|e^{-i\beta\hat J_y/\hbar}|\Phi_2\rangle}.
$$

For the unrotated overlap, the full quasiparticle overlap matrix becomes

$$
\mathcal A=U_1^TU_2^*+V_1^TV_2^*=\begin{pmatrix}(U_1^+)^T(U_2^+)^*+(V_1^-)^T(V_2^-)^*&0\\0&(U_1^+)^\dagger U_2^++(V_1^-)^\dagger V_2^-\end{pmatrix}.
$$

Define the Kramers-representative matrix by

$$
\mathcal A^+\equiv(U_1^+)^T(U_2^+)^*+(V_1^-)^T(V_2^-)^*.
$$

Then

$$
\mathcal A=\begin{pmatrix}\mathcal A^+&0\\0&(\mathcal A^+)^*\end{pmatrix}.
$$

With $\nu_a\equiv\langle0|\Phi_a\rangle$, the [Quasiparticle Vacuum Overlap](../../notes/methods/hfb_pfaffian.md#quasiparticle-vacuum-overlap) reduces to

$$
\boxed{\langle\Phi_1|\Phi_2\rangle=\frac{\nu_1^*\nu_2}{\det U_1^+\,(\det U_2^+)^*}\det\mathcal A^+}.
$$

For positive vacuum phases,

$$
\nu_a=\sqrt{|\det U_a|}=|\det U_a^+|.
$$

If the representative amplitudes are real, define $s_a\equiv\operatorname{sgn}\det U_a^+$. Then

$$
\boxed{\langle\Phi_1|\Phi_2\rangle=s_1s_2\det\mathcal A^+}.
$$
