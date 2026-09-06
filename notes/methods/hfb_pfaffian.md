# HFB Pfaffian Algorithm

### Quasiparticle Vacuum Representation

Let $N_{\mathrm{sp}}$ denote the number of single-particle states. For invertible $U$, use the normalized Thouless vacuum derived in [Quasiparticle Vacuum](hfb.md#quasiparticle-vacuum):

$$
\boxed{|\Phi\rangle=\nu\exp\left[\frac12\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}Z_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta^\dagger\right]|0\rangle}.
$$

Here,

$$
\boxed{\nu\equiv\langle0|\Phi\rangle,\qquad Z\equiv(VU^{-1})^*=V^*U^{-*},\qquad Z^T=-Z}.
$$

Define the pair-creation operator by

$$
\hat A\equiv\frac12\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}Z_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta^\dagger.
$$

The Thouless vacuum expands as

$$
|\Phi\rangle=\nu\exp(\hat A)|0\rangle=\nu\sum_{p=0}^{\lfloor N_{\mathrm{sp}}/2\rfloor}\frac{\hat A^p}{p!}|0\rangle.
$$

Since $\hat A$ creates two particles,

$$
|\Phi^{(2p+1)}\rangle=0.
$$

The zero-particle component is

$$
|\Phi^{(0)}\rangle=\nu|0\rangle.
$$

The two-particle component is

$$
|\Phi^{(2)}\rangle=\nu\sum_{\alpha<\beta}Z_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta^\dagger|0\rangle.
$$

The four-particle component is

$$
|\Phi^{(4)}\rangle=\nu\sum_{\alpha<\beta<\gamma<\delta}\left(Z_{\alpha\beta}Z_{\gamma\delta}-Z_{\alpha\gamma}Z_{\beta\delta}+Z_{\alpha\delta}Z_{\beta\gamma}\right)\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\gamma^\dagger\hat c_\delta^\dagger|0\rangle.
$$

The three coefficients are the complete pairings of four orbitals, with signs fixed by fermionic permutations.

For a positive integer $p$, let $X\in\mathbb C^{2p\times2p}$ be skew-symmetric. Let $S_{2p}$ denote the permutation group of $2p$ elements. Define the Pfaffian by

$$
\boxed{\operatorname{pf}(X)\equiv\frac{1}{2^p p!}\sum_{\sigma\in S_{2p}}\operatorname{sgn}(\sigma)\prod_{r=1}^{p}X_{\sigma(2r-1),\sigma(2r)}}.
$$

The lowest-order Pfaffians are

$$
\operatorname{pf}(X_\varnothing)=1,
$$

$$
\operatorname{pf}\begin{pmatrix}0&X_{12}\\-X_{12}&0\end{pmatrix}=X_{12},
$$

$$
\operatorname{pf}\begin{pmatrix}0&X_{12}&X_{13}&X_{14}\\-X_{12}&0&X_{23}&X_{24}\\-X_{13}&-X_{23}&0&X_{34}\\-X_{14}&-X_{24}&-X_{34}&0\end{pmatrix}=X_{12}X_{34}-X_{13}X_{24}+X_{14}X_{23}.
$$

The two- and four-particle coefficients are therefore Pfaffians of the corresponding principal submatrices of $Z$.

The general $2p$-particle component is

$$
|\Phi^{(2p)}\rangle=\nu\sum_{1\leq\alpha_1<\cdots<\alpha_{2p}\leq N_{\mathrm{sp}}}\operatorname{pf}\left[\left(Z_{\alpha_r\alpha_s}\right)_{r,s=1}^{2p}\right]\hat c_{\alpha_1}^\dagger\cdots\hat c_{\alpha_{2p}}^\dagger|0\rangle.
$$

Therefore,

$$
\boxed{|\Phi\rangle=\nu\sum_{p=0}^{\lfloor N_{\mathrm{sp}}/2\rfloor}\sum_{1\leq\alpha_1<\cdots<\alpha_{2p}\leq N_{\mathrm{sp}}}\operatorname{pf}\left[\left(Z_{\alpha_r\alpha_s}\right)_{r,s=1}^{2p}\right]\hat c_{\alpha_1}^\dagger\cdots\hat c_{\alpha_{2p}}^\dagger|0\rangle}.
$$

### Quantum-Number Restoration

##### Particle-Number Projection

Let $\hat N$ denote the particle-number operator, $N,N'\in\mathbb Z_{\geq0}$ denote particle numbers, and $\varphi\in[0,2\pi)$ denote the gauge angle. The projector is

$$
\boxed{\hat P^N=\frac{1}{2\pi}\int_0^{2\pi}e^{i\varphi(\hat N-N)}\,\mathrm d\varphi}.
$$

The unnormalized projected state is

$$
|\Phi_N\rangle=\hat P^N|\Phi\rangle.
$$

With the Kronecker delta $\delta_{NN'}$,

$$
(\hat P^N)^\dagger=\hat P^N,\qquad \hat P^N\hat P^{N'}=\delta_{NN'}\hat P^N,\qquad \hat N\hat P^N=N\hat P^N.
$$

##### Angular-Momentum Projection

Let $\hat{\boldsymbol J}$ denote total angular momentum in units $\hbar=1$, with $J\geq0$ and $M,K\in\{-J,-J+1,\ldots,J\}$. Use the conventions in [SO(3) and SU(2)](../group_theory/so3_su2.md#wigner-d-functions):

- $SO(3)$ with $\gamma_{\max}=2\pi$ for integer $J$;
- $SU(2)$ with $\gamma_{\max}=4\pi$ for half-odd-integer $J$.

For $\Omega=(\alpha,\beta,\gamma)$,

$$
\hat R(\Omega)=e^{-i\alpha\hat J_z}e^{-i\beta\hat J_y}e^{-i\gamma\hat J_z}.
$$

The angular-momentum projection operator is

$$
\boxed{\hat P^J_{MK}=\frac{2J+1}{4\pi\gamma_{\max}}\int_0^{2\pi}\mathrm d\alpha\int_0^\pi\sin\beta\,\mathrm d\beta\int_0^{\gamma_{\max}}\mathrm d\gamma\,D^{J*}_{MK}(\Omega)\hat R(\Omega)}.
$$

The unnormalized projected state is

$$
|\Phi^J_{MK}\rangle=\hat P^J_{MK}|\Phi\rangle.
$$

The projector satisfies

$$
(\hat P^J_{MK_1})^\dagger=\hat P^J_{K_1M},\qquad \hat P^J_{K_1M}\hat P^J_{MK_2}=\hat P^J_{K_1K_2}.
$$

Therefore,

$$
\boxed{\mathcal N^J_{K_1K_2}=\langle\Phi^J_{MK_1}|\Phi^J_{MK_2}\rangle=\langle\Phi|\hat P^J_{K_1K_2}|\Phi\rangle}.
$$

##### Transformed Quasiparticle Vacuum

Use $|\Phi_0\rangle\equiv|\Phi\rangle$, $\nu_0\equiv\nu$, $Z_0\equiv Z$, and $(U_0,V_0)\equiv(U,V)$. Let $g$ denote a gauge transformation, a rotation, or their product:

$$
\hat G(g)\in\left\{e^{i\varphi\hat N},\hat R(\Omega),e^{i\varphi\hat N}\hat R(\Omega)\right\},\qquad |\Phi(g)\rangle\equiv\hat G(g)|\Phi_0\rangle.
$$

Vacuum invariance, $\hat G(g)|0\rangle=|0\rangle$, implies

$$
\nu(g)\equiv\langle0|\Phi(g)\rangle=\langle0|\hat G(g)|\Phi_0\rangle=\nu_0.
$$

Let $D(g)$ denote the single-particle representation:

$$
\hat G(g)\hat c_\alpha^\dagger\hat G^\dagger(g)=\sum_{\beta=1}^{N_{\mathrm{sp}}}\hat c_\beta^\dagger D_{\beta\alpha}(g).
$$

The Thouless and Bogoliubov matrices transform as

$$
\boxed{Z(g)=D(g)Z_0D^T(g),\qquad U(g)=D(g)U_0,\qquad V(g)=D^*(g)V_0}.
$$

Thus,

$$
|\Phi(g)\rangle=\nu_0\exp\left[\frac12\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}Z_{\alpha\beta}(g)\hat c_\alpha^\dagger\hat c_\beta^\dagger\right]|0\rangle.
$$

For $D(\varphi,\Omega)=e^{i\varphi}D(\Omega)$,

$$
\boxed{Z(\varphi,\Omega)=e^{2i\varphi}D(\Omega)Z_0D^T(\Omega)}.
$$

### Quasiparticle Vacuum Overlap

Orthogonality between distinct occupation patterns gives

$$
\langle\Phi_0|\Phi(g)\rangle=\nu_0^*\nu(g)\sum_{p=0}^{\lfloor N_{\mathrm{sp}}/2\rfloor}\sum_{1\leq\alpha_1<\cdots<\alpha_{2p}\leq N_{\mathrm{sp}}}\operatorname{pf}\left[\left((Z_0)_{\alpha_r\alpha_s}^*\right)_{r,s=1}^{2p}\right]\operatorname{pf}\left[\left([Z(g)]_{\alpha_r\alpha_s}\right)_{r,s=1}^{2p}\right].
$$

Let $I_{N_{\mathrm{sp}}}$ denote the identity matrix. Define

$$
\boxed{\mathbb X(g)\equiv\begin{pmatrix}Z(g)&-I_{N_{\mathrm{sp}}}\\I_{N_{\mathrm{sp}}}&-Z_0^*\end{pmatrix},\qquad s_{N_{\mathrm{sp}}}\equiv(-1)^{N_{\mathrm{sp}}(N_{\mathrm{sp}}+1)/2}}.
$$

The Pfaffian minor-summation identity yields

$$
\boxed{\langle\Phi_0|\Phi(g)\rangle=\nu_0^*\nu(g)s_{N_{\mathrm{sp}}}\operatorname{pf}[\mathbb X(g)]}.
$$

Since $\operatorname{pf}(X)^2=\det X$,

$$
\operatorname{pf}[\mathbb X(g)]^2=\det\mathbb X(g)=\det[I_{N_{\mathrm{sp}}}-Z_0^*Z(g)].
$$

Define the quasiparticle overlap matrix by

$$
\boxed{\mathcal A(g)\equiv U_0^TU(g)^*+V_0^TV(g)^*=U_0^T[I_{N_{\mathrm{sp}}}-Z_0^*Z(g)]U(g)^*}.
$$

Normalization gives $|\nu_0|^2=|\det U_0|$. Since $\nu(g)=\nu_0$ and $U(g)=D(g)U_0$,

$$
\det\mathcal A(g)=|\det U_0|^2\det D(g)^*\det[I_{N_{\mathrm{sp}}}-Z_0^*Z(g)]=\det D(g)^*\langle\Phi_0|\Phi(g)\rangle^2.
$$

Unitarity of $D(g)$ gives the Onishi relation in this convention:

$$
\boxed{\langle\Phi_0|\Phi(g)\rangle^2=\det D(g)\det\mathcal A(g),\qquad \langle\Phi_0|\Phi(g)\rangle=\pm\sqrt{\det D(g)\det\mathcal A(g)}}.
$$

The square root leaves a sign ambiguity; the Pfaffian fixes the branch for a fixed single-particle ordering and vacuum-phase convention.

### One-Body Transition Densities

For $\langle\Phi_0|\Phi(g)\rangle\ne0$, define the normal transition density:

$$
\rho_{\alpha\beta}(g)\equiv\frac{\langle\Phi_0|\hat c_\beta^\dagger\hat c_\alpha|\Phi(g)\rangle}{\langle\Phi_0|\Phi(g)\rangle}.
$$

Use the Thouless relations derived in [Quasiparticle Vacuum](hfb.md#quasiparticle-vacuum):

$$
\hat c_\alpha|\Phi(g)\rangle=\sum_{\gamma=1}^{N_{\mathrm{sp}}}Z_{\alpha\gamma}(g)\hat c_\gamma^\dagger|\Phi(g)\rangle,\qquad \langle\Phi_0|\hat c_\beta^\dagger=\sum_{\delta=1}^{N_{\mathrm{sp}}}(Z_0^*)_{\beta\delta}\langle\Phi_0|\hat c_\delta.
$$

Applying the bra relation, anticommuting the annihilation operators, and applying the ket relation gives

$$
\langle\Phi_0|\hat c_\beta^\dagger\hat c_\alpha|\Phi(g)\rangle=-\sum_{\delta=1}^{N_{\mathrm{sp}}}(Z_0^*)_{\beta\delta}\langle\Phi_0|\hat c_\alpha\hat c_\delta|\Phi(g)\rangle=-\sum_{\gamma,\delta=1}^{N_{\mathrm{sp}}}(Z_0^*)_{\beta\delta}Z_{\delta\gamma}(g)\langle\Phi_0|\hat c_\alpha\hat c_\gamma^\dagger|\Phi(g)\rangle.
$$

Anticommuting $\hat c_\gamma^\dagger$ to the left and dividing by the overlap yields

$$
\rho_{\alpha\beta}(g)=-\sum_{\gamma,\delta=1}^{N_{\mathrm{sp}}}[\delta_{\alpha\gamma}-\rho_{\alpha\gamma}(g)]Z_{\delta\gamma}(g)(Z_0^*)_{\beta\delta}.
$$

In matrix form,

$$
\rho(g)=[\rho(g)-I_{N_{\mathrm{sp}}}]Z^T(g)Z_0^\dagger=[I_{N_{\mathrm{sp}}}-\rho(g)]Z(g)Z_0^\dagger.
$$

Solving for $\rho(g)$ gives

$$
\boxed{\rho(g)=Z(g)Z_0^\dagger[I_{N_{\mathrm{sp}}}+Z(g)Z_0^\dagger]^{-1}=V(g)^*\mathcal A(g)^{-1}V_0^T}.
$$

For a Hermitian single-particle matrix $O$,

$$
\hat O\equiv\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}O_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta,\qquad \frac{\langle\Phi_0|\hat O|\Phi(g)\rangle}{\langle\Phi_0|\Phi(g)\rangle}=\operatorname{Tr}[O\rho(g)].
$$

### Pairing Transition Densities

For $\langle\Phi_0|\Phi(g)\rangle\ne0$, define

$$
\kappa_{\alpha\beta}(g)\equiv\frac{\langle\Phi_0|\hat c_\beta\hat c_\alpha|\Phi(g)\rangle}{\langle\Phi_0|\Phi(g)\rangle},\qquad \bar\kappa_{\alpha\beta}(g)\equiv\frac{\langle\Phi_0|\hat c_\beta^\dagger\hat c_\alpha^\dagger|\Phi(g)\rangle}{\langle\Phi_0|\Phi(g)\rangle}.
$$

Using the Thouless relations from [One-Body Transition Densities](#one-body-transition-densities), the ket relation gives

$$
\kappa_{\alpha\beta}(g)=\sum_{\gamma=1}^{N_{\mathrm{sp}}}Z_{\alpha\gamma}(g)\frac{\langle\Phi_0|\hat c_\beta\hat c_\gamma^\dagger|\Phi(g)\rangle}{\langle\Phi_0|\Phi(g)\rangle}=\sum_{\gamma=1}^{N_{\mathrm{sp}}}Z_{\alpha\gamma}(g)[\delta_{\beta\gamma}-\rho_{\beta\gamma}(g)].
$$

Similarly, the bra relation gives

$$
\bar\kappa_{\alpha\beta}(g)=\sum_{\gamma=1}^{N_{\mathrm{sp}}}(Z_0^*)_{\beta\gamma}\frac{\langle\Phi_0|\hat c_\gamma\hat c_\alpha^\dagger|\Phi(g)\rangle}{\langle\Phi_0|\Phi(g)\rangle}=\sum_{\gamma=1}^{N_{\mathrm{sp}}}(Z_0^*)_{\beta\gamma}[\delta_{\gamma\alpha}-\rho_{\gamma\alpha}(g)].
$$

In matrix form,

$$
\kappa(g)=Z(g)[I_{N_{\mathrm{sp}}}-\rho(g)]^T,\qquad \bar\kappa(g)=[I_{N_{\mathrm{sp}}}-\rho(g)]^TZ_0^\dagger.
$$

The normal transition density satisfies

$$
[I_{N_{\mathrm{sp}}}-\rho(g)]^T=[I_{N_{\mathrm{sp}}}+Z_0^\dagger Z(g)]^{-1}.
$$

Therefore,

$$
\boxed{\kappa(g)=Z(g)[I_{N_{\mathrm{sp}}}+Z_0^\dagger Z(g)]^{-1}=V(g)^*\mathcal A(g)^{-1}U_0^T},
$$

$$
\boxed{\bar\kappa(g)=[I_{N_{\mathrm{sp}}}+Z_0^\dagger Z(g)]^{-1}Z_0^\dagger=U(g)^*\mathcal A(g)^{-1}V_0^T}.
$$

Both pairing densities are skew-symmetric:

$$
\kappa(g)^T=-\kappa(g),\qquad \bar\kappa(g)^T=-\bar\kappa(g).
$$

At the group identity $e$, $\bar\kappa(e)=-\kappa(e)^*$; this conjugation relation does not generally hold for distinct vacua.

### Two-Body Transition Densities

For $\langle\Phi_0|\Phi(g)\rangle\ne0$, define the two-body transition density:

$$
\rho^{(2)}_{\alpha\beta;\gamma\delta}(g)\equiv\frac{\langle\Phi_0|\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma|\Phi(g)\rangle}{\langle\Phi_0|\Phi(g)\rangle}.
$$

Apply the bra Thouless relation from [One-Body Transition Densities](#one-body-transition-densities):

$$
\langle\Phi_0|\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma|\Phi(g)\rangle=\sum_{\mu=1}^{N_{\mathrm{sp}}}(Z_0^*)_{\alpha\mu}\langle\Phi_0|\hat c_\mu\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma|\Phi(g)\rangle.
$$

Anticommuting $\hat c_\mu$ to the right gives

$$
\langle\Phi_0|\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma|\Phi(g)\rangle=(Z_0^*)_{\alpha\beta}\langle\Phi_0|\hat c_\delta\hat c_\gamma|\Phi(g)\rangle-\sum_{\mu=1}^{N_{\mathrm{sp}}}(Z_0^*)_{\alpha\mu}\langle\Phi_0|\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma\hat c_\mu|\Phi(g)\rangle.
$$

Apply the ket Thouless relation:

$$
\langle\Phi_0|\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma|\Phi(g)\rangle=(Z_0^*)_{\alpha\beta}\langle\Phi_0|\hat c_\delta\hat c_\gamma|\Phi(g)\rangle-\sum_{\mu,\nu=1}^{N_{\mathrm{sp}}}(Z_0^*)_{\alpha\mu}Z_{\mu\nu}(g)\langle\Phi_0|\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma\hat c_\nu^\dagger|\Phi(g)\rangle.
$$

Anticommuting $\hat c_\nu^\dagger$ back to the left gives

$$
\langle\Phi_0|\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma|\Phi(g)\rangle=(Z_0^*)_{\alpha\beta}\langle\Phi_0|\hat c_\delta\hat c_\gamma|\Phi(g)\rangle-\sum_{\mu,\nu=1}^{N_{\mathrm{sp}}}(Z_0^*)_{\alpha\mu}Z_{\mu\nu}(g)\left[\delta_{\gamma\nu}\langle\Phi_0|\hat c_\beta^\dagger\hat c_\delta|\Phi(g)\rangle-\delta_{\delta\nu}\langle\Phi_0|\hat c_\beta^\dagger\hat c_\gamma|\Phi(g)\rangle-\langle\Phi_0|\hat c_\nu^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma|\Phi(g)\rangle\right].
$$

Divide by the overlap and collect the two-body terms, using $Z_0^*=-Z_0^\dagger$:

$$
\sum_{\nu=1}^{N_{\mathrm{sp}}}[I_{N_{\mathrm{sp}}}+Z_0^\dagger Z(g)]_{\alpha\nu}\rho^{(2)}_{\nu\beta;\gamma\delta}(g)=-(Z_0^\dagger)_{\alpha\beta}\kappa_{\gamma\delta}(g)+[Z_0^\dagger Z(g)]_{\alpha\gamma}\rho_{\delta\beta}(g)-[Z_0^\dagger Z(g)]_{\alpha\delta}\rho_{\gamma\beta}(g).
$$

The preceding density formulas imply

$$
[I_{N_{\mathrm{sp}}}+Z_0^\dagger Z(g)]^{-1}Z_0^\dagger=\bar\kappa(g),\qquad [I_{N_{\mathrm{sp}}}+Z_0^\dagger Z(g)]^{-1}Z_0^\dagger Z(g)=\rho(g)^T.
$$

For fixed $\beta,\gamma,\delta$, multiplying by the inverse matrix therefore yields

$$
\boxed{\rho^{(2)}_{\alpha\beta;\gamma\delta}(g)=\rho_{\gamma\alpha}(g)\rho_{\delta\beta}(g)-\rho_{\delta\alpha}(g)\rho_{\gamma\beta}(g)-\bar\kappa_{\alpha\beta}(g)\kappa_{\gamma\delta}(g)}.
$$

These three contractions form the fourth-order Pfaffian:

$$
\rho^{(2)}_{\alpha\beta;\gamma\delta}(g)=\operatorname{pf}\begin{pmatrix}0&-\bar\kappa_{\alpha\beta}(g)&\rho_{\delta\alpha}(g)&\rho_{\gamma\alpha}(g)\\\bar\kappa_{\alpha\beta}(g)&0&\rho_{\delta\beta}(g)&\rho_{\gamma\beta}(g)\\-\rho_{\delta\alpha}(g)&-\rho_{\delta\beta}(g)&0&\kappa_{\gamma\delta}(g)\\-\rho_{\gamma\alpha}(g)&-\rho_{\gamma\beta}(g)&-\kappa_{\gamma\delta}(g)&0\end{pmatrix}.
$$

Using the unsymmetrized matrix elements defined in [Second Quantization](../quantum_mechanics/second_quantization.md), the two-body operator and its matrix element are

$$
\hat V\equiv\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}v_{\alpha\beta;\gamma\delta}\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma,\qquad \frac{\langle\Phi_0|\hat V|\Phi(g)\rangle}{\langle\Phi_0|\Phi(g)\rangle}=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}v_{\alpha\beta;\gamma\delta}\rho^{(2)}_{\alpha\beta;\gamma\delta}(g).
$$
