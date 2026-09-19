# Nuclear Particle-Number and Angular-Momentum Projection

### Neutron–Proton Product States

Assume no neutron–proton mixing in the Bogoliubov transformations. The left and right reference vacua are

$$
\boxed{|\Phi_1\rangle=|\Phi_{1,\nu}\rangle\otimes|\Phi_{1,\pi}\rangle,\qquad |\Phi_2\rangle=|\Phi_{2,\nu}\rangle\otimes|\Phi_{2,\pi}\rangle.}
$$

Here $\nu$ and $\pi$ denote neutrons and protons. For $\tau\in\{\nu,\pi\}$, each species vacuum is normalized:

$$
\langle\Phi_{1,\tau}|\Phi_{1,\tau}\rangle=\langle\Phi_{2,\tau}|\Phi_{2,\tau}\rangle=1.
$$

The single-particle indices $\alpha,\beta$ refer to a common basis for the left and right states. The particle operators satisfy

$$
\boxed{\{c_{\tau,\alpha},c_{\tau',\beta}^\dagger\}=\delta_{\tau\tau'}\delta_{\alpha\beta},\qquad \{c_{\tau,\alpha},c_{\tau',\beta}\}=0.}
$$

Use neutron-before-proton ordering for creation operators. The tensor-product notation is understood with this fermionic ordering convention.

For either reference vacuum, suppress its left or right label and denote its Bogoliubov matrices by $U_\tau,V_\tau$. With operator vectors arranged as columns,

$$
\boxed{\hat{\boldsymbol\beta}_\tau=U_\tau^\dagger\hat{\boldsymbol c}_\tau+V_\tau^\dagger(\hat{\boldsymbol c}_\tau)^*,\qquad \hat{\boldsymbol c}_\tau=U_\tau\hat{\boldsymbol\beta}_\tau+V_\tau^*(\hat{\boldsymbol\beta}_\tau)^*.}
$$

For operator columns, $*$ denotes componentwise Hermitian conjugation; for matrices, it denotes complex conjugation. Restoring the reference-state labels,

$$
\hat{\boldsymbol\beta}_{1,\tau}|\Phi_{1,\tau}\rangle=0,\qquad \hat{\boldsymbol\beta}_{2,\tau}|\Phi_{2,\tau}\rangle=0.
$$

The particle operators are common to both reference states; the quasiparticle operators and Bogoliubov matrices depend on the reference vacuum.

Let $r,s$ denote the neutron and proton quasiparticle counts. The labels $\nu_i,\pi_j$ specify both species and quasiparticle orbitals, ordered increasingly within each species. Define the nuclear configuration $\kappa=(\nu_1,\ldots,\nu_r;\pi_1,\ldots,\pi_s)$ by

$$
\boxed{|\Phi_\kappa\rangle=\beta_{\nu_1}^\dagger\cdots\beta_{\nu_r}^\dagger\beta_{\pi_1}^\dagger\cdots\beta_{\pi_s}^\dagger|\Phi\rangle.}
$$

An empty species list contributes no operator; both lists empty give $|\Phi_\varnothing\rangle=|\Phi\rangle$. The adjoint reverses the complete operator order:

$$
\langle\Phi_\kappa|=\langle\Phi|\beta_{\pi_s}\cdots\beta_{\pi_1}\beta_{\nu_r}\cdots\beta_{\nu_1}.
$$

The left and right configurations are denoted by $|\Phi_{1;\kappa_1}\rangle$ and $|\Phi_{2;\kappa_2}\rangle$, constructed from their respective vacua using $\beta_{1,\tau}$ and $\beta_{2,\tau}$. Their quasiparticle counts may differ.

### Particle-Number and Angular-Momentum Projectors

For $\tau\in\{\nu,\pi\}$, define the species particle-number operator by

$$
\hat N_\tau=\sum_\alpha c_{\tau,\alpha}^\dagger c_{\tau,\alpha}.
$$

Let $N,Z$ denote the target neutron and proton numbers. Their independent projectors are

$$
P_\nu^N=\frac{1}{2\pi}\int_0^{2\pi}d\varphi_\nu\,e^{-i\varphi_\nu(\hat N_\nu-N)},\qquad P_\pi^Z=\frac{1}{2\pi}\int_0^{2\pi}d\varphi_\pi\,e^{-i\varphi_\pi(\hat N_\pi-Z)}.
$$

The total angular momentum is $\hat{\boldsymbol I}=\hat{\boldsymbol I}_\nu+\hat{\boldsymbol I}_\pi$. For Euler angles $\Omega=(\alpha,\beta,\gamma)$, define

$$
R(\Omega)=e^{-i\alpha\hat I_z/\hbar}e^{-i\beta\hat I_y/\hbar}e^{-i\gamma\hat I_z/\hbar}=R_\nu(\Omega)\otimes R_\pi(\Omega).
$$

Both species share the same Euler angles. Here $\alpha,\beta$ denote angles rather than single-particle indices.

Let $I$ denote the total angular-momentum quantum number and $M,K=-I,\ldots,I$. Define $D^I_{MK}(\Omega)=\langle IM|R(\Omega)|IK\rangle$. Using the $SU(2)$ integration domain, valid for integer and half-integer $I$,

$$
\boxed{P^I_{MK}=\frac{2I+1}{16\pi^2}\int_0^{2\pi}d\alpha\int_0^\pi\sin\beta\,d\beta\int_0^{4\pi}d\gamma\,D^{I*}_{MK}(\Omega)R(\Omega).}
$$

Particle-number projection commutes with spatial rotation:

$$
[P_\nu^N,P_\pi^Z]=0,\qquad [P_\nu^N,R(\Omega)]=[P_\pi^Z,R(\Omega)]=0.
$$

For a nuclear quasiparticle configuration $|\Phi_\kappa\rangle$, define the unnormalized projected state by

$$
\boxed{|NZIM;\kappa K\rangle=P_\nu^N P_\pi^Z P^I_{MK}|\Phi_\kappa\rangle.}
$$

This construction uses two independent gauge angles and one shared set of Euler angles. The angular-momentum projector acts on the total nuclear state.

For an operator $O$ with definite neutron- and proton-number changes $\Delta N,\Delta Z$,

$$
[\hat N_\nu,O]=\Delta N\,O,\qquad [\hat N_\pi,O]=\Delta Z\,O.
$$

Let $(N_1,Z_1)$ and $(N_2,Z_2)$ denote the left and right target particle numbers. Moving the projectors through $O$ gives

$$
P_\nu^{N_1}P_\pi^{Z_1}O=O\,P_\nu^{N_1-\Delta N}P_\pi^{Z_1-\Delta Z}.
$$

Consequently,

$$
\boxed{P_\nu^{N_1}P_\pi^{Z_1}O\,P_\nu^{N_2}P_\pi^{Z_2}=\delta_{N_1,N_2+\Delta N}\delta_{Z_1,Z_2+\Delta Z}\,O\,P_\nu^{N_2}P_\pi^{Z_2}.}
$$

For allowed transitions, only the two right-state gauge integrations are required. A species-conserving Hamiltonian has $\Delta N=\Delta Z=0$; the charge-exchange operator $c_{\pi,\alpha}^\dagger c_{\nu,\beta}$ has $\Delta N=-1,\Delta Z=+1$.

### Overlap Kernels

Let $g=(\varphi_\nu,\varphi_\pi,\Omega)$ collect the two gauge angles and the shared Euler angles. Define the species transformations and their nuclear product by

$$
G_\nu(g)=e^{-i\varphi_\nu\hat N_\nu}R_\nu(\Omega),\qquad G_\pi(g)=e^{-i\varphi_\pi\hat N_\pi}R_\pi(\Omega),\qquad G(g)=G_\nu(g)\otimes G_\pi(g).
$$

The reference-vacuum overlaps are

$$
n_\nu(g)=\langle\Phi_{1,\nu}|G_\nu(g)|\Phi_{2,\nu}\rangle,\qquad n_\pi(g)=\langle\Phi_{1,\pi}|G_\pi(g)|\Phi_{2,\pi}\rangle.
$$

Since both transformations preserve species particle-number parity,

$$
\boxed{n(g)\equiv\langle\Phi_1|G(g)|\Phi_2\rangle=n_\nu(g)n_\pi(g).}
$$

For the following configuration factorization, assume even particle-number parity for each reference species vacuum. Write the left and right quasiparticle lists explicitly as

$$
\kappa_1=\bigl(\nu_{i_1},\ldots,\nu_{i_{r_1}};\pi_{j_1},\ldots,\pi_{j_{s_1}}\bigr),\qquad \kappa_2=\bigl(\nu_{k_1},\ldots,\nu_{k_{r_2}};\pi_{l_1},\ldots,\pi_{l_{s_2}}\bigr).
$$

Here $i,j$ index the selected left neutron and proton orbitals, while $k,l$ index the corresponding right orbitals. Indices increase strictly within each list; left and right lists may share orbitals.

Their neutron and proton sublists are denoted by $\kappa_{1,\nu},\kappa_{1,\pi}$ and $\kappa_{2,\nu},\kappa_{2,\pi}$. Define the unnormalized species configuration kernels by

$$
\mathcal N^\nu_{\kappa_{1,\nu}\kappa_{2,\nu}}(g)=\langle\Phi_{1,\nu}|\beta_{1,\nu_{i_{r_1}}}\cdots\beta_{1,\nu_{i_1}}G_\nu(g)\beta_{2,\nu_{k_1}}^\dagger\cdots\beta_{2,\nu_{k_{r_2}}}^\dagger|\Phi_{2,\nu}\rangle,
$$

$$
\mathcal N^\pi_{\kappa_{1,\pi}\kappa_{2,\pi}}(g)=\langle\Phi_{1,\pi}|\beta_{1,\pi_{j_{s_1}}}\cdots\beta_{1,\pi_{j_1}}G_\pi(g)\beta_{2,\pi_{l_1}}^\dagger\cdots\beta_{2,\pi_{l_{s_2}}}^\dagger|\Phi_{2,\pi}\rangle.
$$

With the neutron-before-proton ordering fixed above, the nuclear configuration kernel factorizes as

$$
\boxed{\mathcal N_{\kappa_1\kappa_2}(g)\equiv\langle\Phi_{1;\kappa_1}|G(g)|\Phi_{2;\kappa_2}\rangle=\mathcal N^\nu_{\kappa_{1,\nu}\kappa_{2,\nu}}(g)\mathcal N^\pi_{\kappa_{1,\pi}\kappa_{2,\pi}}(g).}
$$

Parity conservation gives

$$
\mathcal N^\nu_{\kappa_{1,\nu}\kappa_{2,\nu}}(g)=0\quad\text{if }r_1+r_2\text{ is odd},\qquad \mathcal N^\pi_{\kappa_{1,\pi}\kappa_{2,\pi}}(g)=0\quad\text{if }s_1+s_2\text{ is odd}.
$$

The quasiparticle counts need not be equal: configurations containing zero and two quasiparticles may have a nonzero rotated overlap.

Factorization holds for any selected neutron–proton configurations. Evaluate both species kernels at the same Euler angles and multiply before angular-momentum integration.

### Contractions

Apply [Transition Densities](hfb_pfaffian.md#transition-densities) separately to $\tau=\nu,\pi$, using the transformed right vacua:

$$
|\Phi_{2,\tau}(g)\rangle=G_\tau(g)|\Phi_{2,\tau}\rangle,\qquad |\Phi_2(g)\rangle=G(g)|\Phi_2\rangle.
$$

The right quasiparticle operators transform consistently:

$$
\beta_{2,\tau,\mu}(g)=G_\tau(g)\beta_{2,\tau,\mu}G_\tau^\dagger(g),\qquad \beta_{2,\tau,\mu}(g)|\Phi_{2,\tau}(g)\rangle=0.
$$

Use $U_{1,\tau},V_{1,\tau}$ and $U_{2,\tau}(g),V_{2,\tau}(g)$ in the referenced formulas. For nuclear and species operators $X,X_\tau$, assume $n_\nu(g),n_\pi(g)\ne0$ and define

$$
\langle X\rangle_g=\frac{\langle\Phi_1|X|\Phi_2(g)\rangle}{n_\nu(g)n_\pi(g)},\qquad \langle X_\tau\rangle_{\tau,g}=\frac{\langle\Phi_{1,\tau}|X_\tau|\Phi_{2,\tau}(g)\rangle}{n_\tau(g)}.
$$

Let $a_{\tau,i}$ be a species-$\tau$ linear fermionic operator at position $i$ in the operator string. It includes particle, left quasiparticle, and transformed right quasiparticle operators and their adjoints. For even-parity species reference vacua,

$$
\boxed{\langle a_{\tau,i}a_{\tau',j}\rangle_g=\delta_{\tau\tau'}\,\langle a_{\tau,i}a_{\tau,j}\rangle_{\tau,g}.}
$$

Each species block uses the eleven referenced contractions. Cross-species contractions vanish; charge-exchange kernels between excited configurations can remain nonzero through same-species contractions.

### One-Body Operator Kernels

Following [One-Body Matrix Elements between Multiquasiparticle States](hfb_pfaffian.md#one-body-matrix-elements-between-multiquasiparticle-states), use hats for operators and $*$ for scalar creation operators. For single-particle matrix elements $\langle\tau_1,\alpha_1|\hat O|\tau_2,\alpha_2\rangle$, define

$$
\hat O^{(1)}=\sum_{\tau_1,\tau_2=\nu,\pi}\sum_{\alpha_1\alpha_2}\langle\tau_1,\alpha_1|\hat O|\tau_2,\alpha_2\rangle\,\hat c_{\tau_1,\alpha_1}^*\hat c_{\tau_2,\alpha_2}.
$$

For fixed configurations $\kappa_1,\kappa_2$, let $L=r_1+s_1+r_2+s_2+2$. Define the ordered operator column

$$
\hat{\boldsymbol a}(\tau_1,\alpha_1;\tau_2,\alpha_2;g)=\bigl(\hat\beta_{1,\pi_{j_{s_1}}},\ldots,\hat\beta_{1,\pi_{j_1}},\hat\beta_{1,\nu_{i_{r_1}}},\ldots,\hat\beta_{1,\nu_{i_1}},\hat c_{\tau_1,\alpha_1}^*,\hat c_{\tau_2,\alpha_2},\hat\beta_{2,\nu_{k_1}}^*(g),\ldots,\hat\beta_{2,\nu_{k_{r_2}}}^*(g),\hat\beta_{2,\pi_{l_1}}^*(g),\ldots,\hat\beta_{2,\pi_{l_{s_2}}}^*(g)\bigr)^T.
$$

Empty quasiparticle strings are omitted. The right operators are transformed by $G(g)$. For even $L$, construct $S(\tau_1,\alpha_1;\tau_2,\alpha_2;g)\in\mathbb C^{L\times L}$ using [Contractions](#contractions):

$$
S_{ij}(\tau_1,\alpha_1;\tau_2,\alpha_2;g)=\langle\hat a_i(\tau_1,\alpha_1;\tau_2,\alpha_2;g)\hat a_j(\tau_1,\alpha_1;\tau_2,\alpha_2;g)\rangle_g\quad(i<j),\qquad S=-S^T.
$$

The dependence on $\kappa_1,\kappa_2$ is implicit. The [generalized Wick theorem](hfb_pfaffian.md#generalized-wick-theorem) gives

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat O^{(1)}G(g)|\Phi_{2;\kappa_2}\rangle=n_\nu(g)n_\pi(g)\sum_{\tau_1,\tau_2=\nu,\pi}\sum_{\alpha_1\alpha_2}\langle\tau_1,\alpha_1|\hat O|\tau_2,\alpha_2\rangle\,\operatorname{pf}\!\left[S(\tau_1,\alpha_1;\tau_2,\alpha_2;g)\right].}
$$

For odd $L$, the matrix element vanishes. Cross-species entries of $S$ vanish; the operator ordering fixes all fermionic signs. Species-conserving channels require $r_1+r_2$ and $s_1+s_2$ both even; charge-exchange channels require both odd.

### Two-Body Operator Kernels

Use the configurations and operator conventions of [One-Body Operator Kernels](#one-body-operator-kernels). For unsymmetrized two-body matrix elements, define

$$
\hat O^{(2)}=\frac12\sum_{\tau_1,\tau_2,\tau_3,\tau_4=\nu,\pi}\sum_{\alpha_1\alpha_2\alpha_3\alpha_4}\langle\tau_1,\alpha_1;\tau_2,\alpha_2|\hat O|\tau_3,\alpha_3;\tau_4,\alpha_4\rangle\,\hat c_{\tau_1,\alpha_1}^*\hat c_{\tau_2,\alpha_2}^*\hat c_{\tau_4,\alpha_4}\hat c_{\tau_3,\alpha_3}.
$$

Here $|\tau_1,\alpha_1;\tau_2,\alpha_2\rangle=|\tau_1,\alpha_1\rangle\otimes|\tau_2,\alpha_2\rangle$. The annihilation order is reversed relative to the ket labels.

For fixed configurations $\kappa_1,\kappa_2$, let $L=r_1+s_1+r_2+s_2+4$. Define the ordered operator column

$$
\hat{\boldsymbol a}(\tau_1,\alpha_1;\tau_2,\alpha_2;\tau_3,\alpha_3;\tau_4,\alpha_4;g)=\bigl(\hat\beta_{1,\pi_{j_{s_1}}},\ldots,\hat\beta_{1,\pi_{j_1}},\hat\beta_{1,\nu_{i_{r_1}}},\ldots,\hat\beta_{1,\nu_{i_1}},\hat c_{\tau_1,\alpha_1}^*,\hat c_{\tau_2,\alpha_2}^*,\hat c_{\tau_4,\alpha_4},\hat c_{\tau_3,\alpha_3},\hat\beta_{2,\nu_{k_1}}^*(g),\ldots,\hat\beta_{2,\nu_{k_{r_2}}}^*(g),\hat\beta_{2,\pi_{l_1}}^*(g),\ldots,\hat\beta_{2,\pi_{l_{s_2}}}^*(g)\bigr)^T.
$$

Empty quasiparticle strings are omitted. For even $L$, construct $S(\tau_1,\alpha_1;\tau_2,\alpha_2;\tau_3,\alpha_3;\tau_4,\alpha_4;g)\in\mathbb C^{L\times L}$ using [Contractions](#contractions):

$$
S_{ij}=\langle\hat a_i\hat a_j\rangle_g\quad(i<j),\qquad S=-S^T.
$$

The dependence on $\kappa_1,\kappa_2$ is implicit. The [generalized Wick theorem](hfb_pfaffian.md#generalized-wick-theorem) gives

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat O^{(2)}G(g)|\Phi_{2;\kappa_2}\rangle=\frac12 n_\nu(g)n_\pi(g)\sum_{\tau_1,\tau_2,\tau_3,\tau_4=\nu,\pi}\sum_{\alpha_1\alpha_2\alpha_3\alpha_4}\langle\tau_1,\alpha_1;\tau_2,\alpha_2|\hat O|\tau_3,\alpha_3;\tau_4,\alpha_4\rangle\,\operatorname{pf}\!\left[S(\tau_1,\alpha_1;\tau_2,\alpha_2;\tau_3,\alpha_3;\tau_4,\alpha_4;g)\right].}
$$

For odd $L$, the matrix element vanishes. Cross-species contractions vanish; the specified operator order fixes all fermionic signs.

Each term changes the species particle numbers by

$$
\Delta N=\delta_{\tau_1\nu}+\delta_{\tau_2\nu}-\delta_{\tau_3\nu}-\delta_{\tau_4\nu},\qquad \Delta Z=-\Delta N.
$$

A species-conserving interaction retains only terms with $\Delta N=\Delta Z=0$. The factor $1/2$ assumes the full ordered species and orbital sums with unsymmetrized matrix elements; antisymmetrized matrix elements instead require $1/4$.

### Projected Matrix Elements

For fixed configurations $\kappa_1,\kappa_2$, define the kernel integral

$$
\mathcal I^{NZI}_{K_1K_2}[\hat O]\equiv\frac{2I+1}{16\pi^2(2\pi)^2}\int d\Omega\int_0^{2\pi}d\varphi_\nu\int_0^{2\pi}d\varphi_\pi\,D^{I*}_{K_1K_2}(\Omega)e^{i(N\varphi_\nu+Z\varphi_\pi)}\langle\Phi_{1;\kappa_1}|\hat O\,G(g)|\Phi_{2;\kappa_2}\rangle.
$$

Here $d\Omega=d\alpha\,\sin\beta\,d\beta\,d\gamma$, with the $SU(2)$ domain defined above. Equivalently,

$$
\mathcal I^{NZI}_{K_1K_2}[\hat O]=\langle\Phi_{1;\kappa_1}|\hat O\,P_\nu^N P_\pi^Z P^I_{K_1K_2}|\Phi_{2;\kappa_2}\rangle.
$$

Let $\mathbf 1$ denote the identity operator. The overlap between projected states is

$$
\boxed{\langle N_1Z_1I_1M_1;\kappa_1K_1|N_2Z_2I_2M_2;\kappa_2K_2\rangle=\delta_{N_1N_2}\delta_{Z_1Z_2}\delta_{I_1I_2}\delta_{M_1M_2}\,\mathcal I^{N_2Z_2I_2}_{K_1K_2}[\mathbf 1].}
$$

For a rotationally invariant, species-conserving Hamiltonian $\hat H$,

$$
[\hat H,\hat{\boldsymbol I}]=0,\qquad [\hat H,\hat N_\nu]=[\hat H,\hat N_\pi]=0,
$$

and

$$
\boxed{\langle N_1Z_1I_1M_1;\kappa_1K_1|\hat H|N_2Z_2I_2M_2;\kappa_2K_2\rangle=\delta_{N_1N_2}\delta_{Z_1Z_2}\delta_{I_1I_2}\delta_{M_1M_2}\,\mathcal I^{N_2Z_2I_2}_{K_1K_2}[\hat H].}
$$

Evaluate the integrands using the overlap, one-body, and two-body kernels defined above.

For an irreducible tensor $\hat T_{\lambda\mu}$ of angular rank $\lambda$, use the Wigner–Eckart convention

$$
\langle N_1Z_1I_1M_1;\kappa_1K_1|\hat T_{\lambda\mu}|N_2Z_2I_2M_2;\kappa_2K_2\rangle=\frac{\langle I_2M_2;\lambda\mu|I_1M_1\rangle}{\sqrt{2I_1+1}}\langle N_1Z_1I_1;\kappa_1K_1\Vert\hat T_\lambda\Vert N_2Z_2I_2;\kappa_2K_2\rangle.
$$

Assume each component has definite species-number changes:

$$
[\hat N_\nu,\hat T_{\lambda\mu}]=\Delta N\,\hat T_{\lambda\mu},\qquad [\hat N_\pi,\hat T_{\lambda\mu}]=\Delta Z\,\hat T_{\lambda\mu}.
$$

Using the [projected tensor identity](hfb_projection.md#irreducible-tensor-operators-and-projected-matrix-elements),

$$
\boxed{\langle N_1Z_1I_1;\kappa_1K_1\Vert\hat T_\lambda\Vert N_2Z_2I_2;\kappa_2K_2\rangle=\delta_{N_1,N_2+\Delta N}\delta_{Z_1,Z_2+\Delta Z}\sqrt{2I_1+1}\sum_{\mu'=-\lambda}^{\lambda}\langle I_2,K_1-\mu';\lambda\mu'|I_1K_1\rangle\,\mathcal I^{N_2Z_2I_2}_{K_1-\mu',K_2}[\hat T_{\lambda\mu'}].}
$$

Only terms with $|K_1-\mu'|\leq I_2$ contribute. These matrix elements refer to unnormalized projected configurations; normalization and configuration mixing require the projected overlap matrix.
