# Beta Decay

### Relativistic Kinematics

Let $P$, $D$, $e^-$, and $\bar\nu_e$ denote the parent nucleus, daughter nucleus, electron, and electron antineutrino, with masses $M_P$, $M_D$, $m_e$, and $m_\nu$.

The decay process is

$$
P\rightarrow D+e^-+\bar\nu_e.
$$

For $m_\nu=0$, the decay energy is

$$
Q_{\beta^-}=(M_P-M_D-m_e)c^2=T_D+T_e+T_\nu.
$$

Four-momentum conservation and the on-shell condition are

$$
p_P^\mu=p_D^\mu+p_e^\mu+p_\nu^\mu,\qquad p^\mu=(E,\mathbf p c),\qquad p^\mu p_\mu=E^2-\mathbf p^2c^2=m^2c^4.
$$

In the parent-nucleus rest frame, $p_P^\mu=(M_Pc^2,\mathbf 0)$. Define the squared invariant energy of the daughter nucleus and neutrino by

$$
s_{D\nu}\equiv(p_D+p_\nu)^2=(p_P-p_e)^2=M_P^2c^4+m_e^2c^4-2M_Pc^2E_e.
$$

Therefore,

$$
E_e=\frac{M_P^2c^4+m_e^2c^4-s_{D\nu}}{2M_Pc^2}.
$$

The maximum electron energy occurs at the minimum invariant mass

$$
s_{D\nu}^{\min}=(M_D+m_\nu)^2c^4.
$$

For $m_\nu=0$,

$$
E_e^{\max}=\frac{M_P^2+m_e^2-M_D^2}{2M_P}c^2.
$$

The corresponding maximum kinetic energy is

$$
T_e^{\max}=E_e^{\max}-m_ec^2=\frac{(M_P-m_e)^2-M_D^2}{2M_P}c^2.
$$

Using $M_P=M_D+m_e+Q_{\beta^-}/c^2$,

$$
\boxed{T_e^{\max}=\frac{Q_{\beta^-}(Q_{\beta^-}+2M_Dc^2)}{2M_Pc^2}\simeq Q_{\beta^-}}.
$$

The electron kinetic-energy range is

$$
0\leq T_e\leq T_e^{\max}.
$$

### Beta-Decay Modes

Let $A=Z+N$ be the mass number. Denote nuclear and neutral-atomic masses by $M_{\mathrm{nuc}}$ and $M_{\mathrm{atom}}$. Electronic binding-energy differences are neglected unless stated otherwise.

**Free-neutron decay.**

$$
n\rightarrow p+e^-+\bar\nu_e.
$$

For negligible neutrino mass,

$$
Q_n=(m_n-m_p-m_e)c^2\simeq0.782\,\mathrm{MeV}.
$$

**Beta-minus decay.**

$$
{}^{A}_{Z}P_N\rightarrow{}^{A}_{Z+1}D_{N-1}+e^-+\bar\nu_e.
$$

$$
\boxed{Q_{\beta^-}=\left[M_{\mathrm{nuc}}(P)-M_{\mathrm{nuc}}(D)-m_e\right]c^2\simeq\left[M_{\mathrm{atom}}(P)-M_{\mathrm{atom}}(D)\right]c^2}.
$$

**Beta-plus decay.**

$$
{}^{A}_{Z}P_N\rightarrow{}^{A}_{Z-1}D_{N+1}+e^++\nu_e.
$$

$$
\boxed{Q_{\beta^+}=\left[M_{\mathrm{nuc}}(P)-M_{\mathrm{nuc}}(D)-m_e\right]c^2\simeq\left[M_{\mathrm{atom}}(P)-M_{\mathrm{atom}}(D)-2m_e\right]c^2}.
$$

**Electron capture.** Let $e_x^-$ denote an electron captured from atomic shell $x$, with binding energy $B_x$.

$$
{}^{A}_{Z}P_N+e_x^-\rightarrow{}^{A}_{Z-1}D_{N+1}+\nu_e.
$$

The atomic mass difference is

$$
Q_{\mathrm{EC}}=\left[M_{\mathrm{atom}}(P)-M_{\mathrm{atom}}(D)\right]c^2.
$$

For a daughter-nucleus excitation energy $E_{\mathrm{exc}}$, the neutrino energy is approximately

$$
\boxed{E_\nu\simeq Q_{\mathrm{EC}}-B_x-E_{\mathrm{exc}}}.
$$

### Beta-Decay Phase Space

Fermi’s golden rule gives

$$
d\lambda_\beta=\frac{2\pi}{\hbar}|M_{fi}|^2dn(e,\nu),
$$

where $M_{fi}$ is the transition matrix element and $dn(e,\nu)$ is the energy-constrained final-state element.

In the parent-nucleus rest frame,

$$
\mathbf p_P=\mathbf p_e+\mathbf p_\nu+\mathbf p_D=\mathbf 0,\qquad E_P=E_e+E_\nu+E_D.
$$

Neglecting the daughter recoil, $E_D\simeq M_Dc^2$. Define

$$
E_0=E_e+E_\nu=E_P-E_D\simeq(M_P-M_D)c^2.
$$

For $\beta^-$ decay, let $\nu\equiv\bar\nu_e$. The same phase-space result applies to $\beta^+$ decay under $e^-\rightarrow e^+$ and $\bar\nu_e\rightarrow\nu_e$.

Let $V$ be the box-normalization volume. Let $\mathbf r_x$ and $\mathbf p_x$ denote the position and momentum of particle $x=e,\nu$, with $p_x=|\mathbf p_x|$.

The daughter momentum is fixed by momentum conservation, $\mathbf p_D=-(\mathbf p_e+\mathbf p_\nu)$, and is not an independent integration variable. The electron–neutrino final-state element is

$$
dn(e,\nu)\equiv\frac{d^3\mathbf r_e\,d^3\mathbf p_e}{(2\pi\hbar)^3}\frac{d^3\mathbf r_\nu\,d^3\mathbf p_\nu}{(2\pi\hbar)^3}\delta(E_0-E_e-E_\nu).
$$

Keeping $E_e$ as the spectral variable,

$$
dn(e,\nu)=\frac{p_e^2dp_e}{(2\pi\hbar)^6}\int_Vd^3\mathbf r_e\int_Vd^3\mathbf r_\nu\int_0^\infty p_\nu^2dp_\nu\int d\Omega_{\mathbf p_e}\int d\Omega_{\mathbf p_\nu}\,\delta(E_0-E_e-E_\nu).
$$

For a particle with energy $E$, momentum magnitude $p$, and mass $m$,

$$
E^2=p^2c^2+m^2c^4,\qquad E\,dE=pc^2\,dp.
$$

Therefore,

$$
\int_0^\infty p_\nu^2dp_\nu\,\delta(E_0-E_e-E_\nu)=\frac{p_\nu E_\nu}{c^2},\qquad p_e^2dp_e=\frac{p_eE_e}{c^2}\,dE_e.
$$

The spatial and angular integrations are

$$
\int_Vd^3\mathbf r_e\int_Vd^3\mathbf r_\nu=V^2,\qquad \int d\Omega_{\mathbf p_e}\int d\Omega_{\mathbf p_\nu}=(4\pi)^2.
$$

The electron momentum is

$$
p_e=\frac{1}{c}\sqrt{E_e^2-m_e^2c^4}.
$$

Thus,

$$
dn(e,\nu)=\frac{V^2}{4\pi^4\hbar^6c^4}p_\nu E_\nu\,p_eE_e\,dE_e=\frac{V^2}{4\pi^4\hbar^6c^5}(E_0-E_e)^2\sqrt{1-\frac{m_\nu^2c^4}{(E_0-E_e)^2}}\,p_eE_e\,dE_e.
$$

Let $\varphi_e$ be the outgoing electron wave function, $Z_D$ the daughter charge number, $\alpha=e^2/(4\pi\epsilon_0\hbar c)$ the fine-structure constant, and $v_e=p_ec^2/E_e$. In the point-Coulomb approximation, the electron distortion is described by the Fermi function

$$
|\varphi_e(0)|^2\simeq\frac{1}{V}F(Z_D,p_e),\qquad F(Z_D,p_e)=\frac{2\pi\eta}{1-e^{-2\pi\eta}},\qquad \eta=\pm\frac{Z_D\alpha c}{v_e}.
$$

The upper and lower signs correspond to $\beta^-$ and $\beta^+$ decay.

Define the reduced matrix element $M_{fi}'$ by

$$
|M_{fi}|^2=\frac{1}{V^2}F(Z_D,p_e)|M_{fi}'|^2.
$$

The electron energy spectrum is

$$
\frac{d\lambda_\beta}{dE_e}=\frac{|M_{fi}'|^2}{2\pi^3\hbar^7c^5}F(Z_D,p_e)p_eE_e(E_0-E_e)^2\sqrt{1-\frac{m_\nu^2c^4}{(E_0-E_e)^2}}.
$$

Within the recoil-free approximation,

$$
m_ec^2\leq E_e\leq E_0-m_\nu c^2.
$$

The total decay constant and half-life are

$$
\lambda_\beta=\int_{m_ec^2}^{E_0-m_\nu c^2}\frac{d\lambda_\beta}{dE_e}\,dE_e,\qquad T_{1/2}=\frac{\ln2}{\lambda_\beta}.
$$

### Transition Matrix Elements

Let $\Psi_P$ and $\Psi_D$ be the parent and daughter many-body nuclear wave functions. The transition matrix element is

$$
M_{fi}=\langle\Psi_D,e^-,\bar\nu_e|\hat H_{\mathrm{int}}|\Psi_P\rangle.
$$

For nucleon coordinates $\mathbf r_a$ and lepton coordinates $\mathbf r_e,\mathbf r_\nu$, let $\varphi_{\bar\nu_e}$ be the outgoing antineutrino wave function.

$$
M_{fi}=\int d^3\mathbf r_e\,d^3\mathbf r_\nu\prod_{a=1}^{A}d^3\mathbf r_a\,\Psi_D^*(\{\mathbf r_a\})\varphi_e^*(\mathbf r_e)\varphi_{\bar\nu_e}^*(\mathbf r_\nu)\hat H_{\mathrm{int}}\Psi_P(\{\mathbf r_a\}).
$$

In the spinless non-relativistic model, use the contact interaction with schematic coupling $g$,

$$
\hat H_{\mathrm{int}}=g\,\delta^3(\mathbf r_n-\mathbf r_p)\delta^3(\mathbf r_n-\mathbf r_e)\delta^3(\mathbf r_n-\mathbf r_\nu)\hat O_\beta,
$$

where $\hat O_\beta$ transforms the active neutron into a proton.

Separate the participating nucleon from the $A-1$ spectator nucleons:

$$
\Psi_P=\phi_P(\mathbf r_n)_n\Phi_s(\{\mathbf r_s\}),\qquad \Psi_D=\phi_D(\mathbf r_p)_p\Phi_s(\{\mathbf r_s\}),
$$

where $\phi_P$ and $\phi_D$ describe the participating neutron and proton, while $\Phi_s$ describes the unchanged spectators. Their overlap is

$$
\int\prod_s d^3\mathbf r_s\,\Phi_s^*\Phi_s=1.
$$

Evaluating the delta functions gives

$$
M_{fi}=g\int d^3\mathbf r\,\phi_D^*(\mathbf r)_p\varphi_e^*(\mathbf r)\varphi_{\bar\nu_e}^*(\mathbf r)\hat O_\beta\phi_P(\mathbf r)_n.
$$

After extracting the external-state normalization and Coulomb factor used above, let $\mathbf k_e$, $\mathbf k_{\bar\nu_e}$, and $\mathbf k_D$ be the outgoing wave vectors, with $\mathbf p_x=\hbar\mathbf k_x$, and define

$$
\mathbf k\equiv\mathbf k_e+\mathbf k_{\bar\nu_e}=-\mathbf k_D,\qquad k=|\mathbf k|.
$$

Using the resulting plane-wave phase,

$$
M_{fi}'=g\int d^3\mathbf r\,e^{i\mathbf k\cdot\mathbf r}\phi_D^*(\mathbf r)_p\hat O_\beta\phi_P(\mathbf r)_n.
$$

Let $L_\beta$ and $M_\beta$ be the orbital angular momentum carried by the lepton pair and its projection. Let $j_{L_\beta}$ and $Y_{L_\beta M_\beta}$ be the spherical Bessel function and spherical harmonic. Hats denote unit vectors. The plane-wave expansion is

$$
e^{i\mathbf k\cdot\mathbf r}=4\pi\sum_{L_\beta,M_\beta}i^{L_\beta}j_{L_\beta}(kr)Y_{L_\beta M_\beta}(\hat{\mathbf k})Y_{L_\beta M_\beta}^*(\hat{\mathbf r}).
$$

Therefore,

$$
M_{fi}'=4\pi g\sum_{L_\beta,M_\beta}i^{L_\beta}Y_{L_\beta M_\beta}(\hat{\mathbf k})\int d^3\mathbf r\,\phi_D^*(\mathbf r)_p\hat O_\beta j_{L_\beta}(kr)Y_{L_\beta M_\beta}^*(\hat{\mathbf r})\phi_P(\mathbf r)_n.
$$

In the spectator approximation, let $L_P$ and $L_D$ be the orbital angular momenta of the active-particle states $\phi_P$ and $\phi_D$, and let $\pi_P$ and $\pi_D$ be the parent and daughter parities. For a spatially scalar $\hat O_\beta$, the Wigner–Eckart theorem gives

$$
\boxed{|L_P-L_D|\leq L_\beta\leq L_P+L_D,\qquad \pi_P=\pi_D(-1)^{L_\beta}}.
$$

For $kr\ll1$,

$$
j_{L_\beta}(kr)\simeq\frac{(kr)^{L_\beta}}{(2L_\beta+1)!!}.
$$

Higher $L_\beta$ contributions are therefore rapidly suppressed.

### Fermi and Gamow–Teller Transitions

Let $J_P$ and $J_D$ be the parent and daughter nuclear spins. The electron and antineutrino are spin-$1/2$ particles. Let $\mathbf L_\beta$, $\mathbf S_\beta$, and $\mathbf J_\beta$ denote the orbital, intrinsic-spin, and total angular momenta carried by the lepton pair, with $S_\beta=0,1$.

$$
\boxed{\mathbf J_\beta=\mathbf L_\beta+\mathbf S_\beta,\qquad \mathbf J_P=\mathbf J_D+\mathbf J_\beta,\qquad \pi_P=\pi_D(-1)^{L_\beta}}.
$$

For allowed beta decay, $L_\beta=0$ and the parity is unchanged.

Let $\tau_\beta(a)$ be the isospin operator that converts the $a$th nucleon according to the decay mode, and let $\boldsymbol\sigma(a)$ be its Pauli spin operator.

**Fermi transitions.** The lepton spins couple to $S_\beta=0$. The nuclear operator is

$$
\hat O_{\mathrm F}=\sum_{a=1}^{A}\tau_\beta(a).
$$

The reduced nuclear matrix element is

$$
M_{\mathrm F}=\langle\Psi_D\|\hat O_{\mathrm F}\|\Psi_P\rangle.
$$

Since $\hat O_{\mathrm F}$ is a scalar in spin space,

$$
\boxed{J_P=J_D,\qquad \pi_P=\pi_D}.
$$

An allowed $0^+\rightarrow0^+$ transition is therefore purely Fermi.

**Gamow–Teller transitions.** The lepton spins couple to $S_\beta=1$. The nuclear operator is

$$
\hat{\mathbf O}_{\mathrm{GT}}=\sum_{a=1}^{A}\boldsymbol\sigma(a)\tau_\beta(a).
$$

The reduced nuclear matrix element is

$$
M_{\mathrm{GT}}=\langle\Psi_D\|\hat{\mathbf O}_{\mathrm{GT}}\|\Psi_P\rangle.
$$

Since $\hat{\mathbf O}_{\mathrm{GT}}$ has spin rank one,

$$
\boxed{|J_P-J_D|\leq1,\qquad \pi_P=\pi_D,\qquad 0\not\rightarrow0}.
$$

Let $g_{\mathrm F}$ and $g_{\mathrm{GT}}$ be the full low-energy Fermi and Gamow–Teller coupling strengths that replace the schematic $g$ when intrinsic spin is included. For the unpolarized angle-integrated allowed rate, average over the initial nuclear substates and sum over the final nuclear substates and lepton spins:

$$
\boxed{\overline{|M_{fi}'|^2}=\frac{g_{\mathrm F}^2|M_{\mathrm F}|^2+g_{\mathrm{GT}}^2|M_{\mathrm{GT}}|^2}{2J_P+1}}.
$$

### Allowed and Forbidden Transitions

Define the nuclear spin change by

$$
\Delta J\equiv|J_P-J_D|.
$$

Angular-momentum conservation requires

$$
\boxed{\Delta J\leq J_\beta\leq J_P+J_D}.
$$

**Allowed transitions.**

$$
L_\beta=0,\qquad \pi_P=\pi_D.
$$

$$
S_\beta=0:\qquad J_\beta=0,\qquad \Delta J=0\quad\text{for Fermi}.
$$

$$
S_\beta=1:\qquad J_\beta=1,\qquad \Delta J=0,1,\qquad 0\not\rightarrow0\quad\text{for Gamow--Teller}.
$$

**First-forbidden transitions.**

$$
L_\beta=1,\qquad \pi_P=-\pi_D.
$$

$$
S_\beta=0:\qquad J_\beta=1,\qquad \Delta J=0,1,\qquad 0\not\rightarrow0\quad\text{for Fermi}.
$$

$$
S_\beta=1:\qquad J_\beta=0,1,2,\qquad \Delta J=0,1,2\quad\text{for Gamow--Teller}.
$$

**Higher partial waves.** For $L_\beta\geq2$,

$$
\pi_P=\pi_D(-1)^{L_\beta}.
$$

$$
S_\beta=0:\qquad J_\beta=L_\beta,\qquad \Delta J=0,1,\ldots,L_\beta\quad\text{for Fermi}.
$$

$$
S_\beta=1:\qquad J_\beta=L_\beta-1,L_\beta,L_\beta+1,\qquad \Delta J=0,1,\ldots,L_\beta+1\quad\text{for Gamow--Teller}.
$$

**Leading forbidden contribution.** Let $R$ be the nuclear radius. For $kR\ll1$,

$$
j_{L_\beta}(kR)\simeq\frac{(kR)^{L_\beta}}{(2L_\beta+1)!!}.
$$

The smallest $L_\beta$ allowed by the angular-momentum and parity selection rules gives the forbidden order. Define

$$
n\equiv\min_{\mathrm{allowed}}L_\beta\!\left(J_P^{\pi_P}\rightarrow J_D^{\pi_D}\right).
$$

Thus,

$$
n=0\quad\text{allowed},\qquad n=1\quad\text{first-forbidden},\qquad n=2\quad\text{second-forbidden}.
$$

For first-forbidden transitions,

$$
n=1:\qquad \Delta J=0,1\quad\text{non-unique},\qquad \Delta J=2\quad\text{unique}.
$$

For higher forbidden transitions,

$$
n\geq2:\qquad \Delta J=n\quad\text{non-unique},\qquad \Delta J=n+1\quad\text{unique}.
$$

A forbidden transition is called unique when the selection rules permit only one leading-order nuclear tensor operator. For $\Delta J=n+1$, only the maximal Gamow--Teller coupling contributes:

$$
L_\beta=n,\qquad S_\beta=1,\qquad J_\beta=n+1.
$$

Otherwise, several leading-order nuclear operators may contribute, and the transition is non-unique. Here, “unique” refers to the leading operator, not to the number of final states.

### Dirac Vertex Structures

Let $d$, $u$, $n$, $p$, $e^-$, and $\bar\nu_e$ denote the down quark, up quark, neutron, proton, electron, and electron antineutrino. Let $W^{-*}$ denote an off-shell charged weak boson. Beta decay proceeds through

$$
d\rightarrow u+W^{-*}\rightarrow u+e^-+\bar\nu_e,\qquad n\rightarrow p+W^{-*}\rightarrow p+e^-+\bar\nu_e.
$$

Let $q^\mu=(\Delta E,\mathbf q c)$ be the transferred four-momentum, $M_W$ the $W$-boson mass, $g_2$ the $SU(2)_{\mathrm L}$ gauge coupling, and $G_{\mathrm F}$ the Fermi constant. Since $|q^2|\ll(M_Wc^2)^2$,

$$
\frac{1}{q^2-(M_Wc^2)^2}\simeq-\frac{1}{(M_Wc^2)^2},\qquad \frac{G_{\mathrm F}}{\sqrt2}=\frac{g_2^2(\hbar c)^3}{8(M_Wc^2)^2}.
$$

Low-energy $W$ exchange therefore reduces to a local four-fermion interaction.

Let $\psi_n$, $\psi_p$, $\psi_e$, and $\psi_{\nu_e}$ denote the neutron, proton, electron, and electron-neutrino fields, with $\bar\psi\equiv\psi^\dagger\gamma^0$. Define

$$
\gamma^5\equiv i\gamma^0\gamma^1\gamma^2\gamma^3,\qquad \sigma^{\mu\nu}\equiv\frac{i}{2}[\gamma^\mu,\gamma^\nu].
$$

The independent Dirac structures are

$$
\Gamma_{\mathrm S}=1,\qquad \Gamma_{\mathrm P}=\gamma^5,\qquad \Gamma_{\mathrm V}=\gamma^\mu,\qquad \Gamma_{\mathrm A}=\gamma^\mu\gamma^5,\qquad \Gamma_{\mathrm T}=\sigma^{\mu\nu}.
$$

Let $V_{ud}$ be the CKM matrix element and $C_i,C_i'$ dimensionless coupling coefficients. The general local beta-decay Hamiltonian density is

$$
\boxed{\mathcal H_\beta^{\mathrm{general}}=\frac{G_{\mathrm F}V_{ud}}{\sqrt2}\sum_{i=\mathrm S,\mathrm P,\mathrm V,\mathrm A,\mathrm T}\left(\bar\psi_p\Gamma_i\psi_n\right)\left[\bar\psi_e\Gamma_i(C_i+C_i'\gamma^5)\psi_{\nu_e}\right]+\mathrm{h.c.}}
$$

In this convention, the tree-level Standard Model imposes

$$
C_{\mathrm S}=C_{\mathrm S}'=C_{\mathrm P}=C_{\mathrm P}'=C_{\mathrm T}=C_{\mathrm T}'=0,\qquad C_{\mathrm V}'=-C_{\mathrm V},\qquad C_{\mathrm A}'=-C_{\mathrm A}.
$$

Thus, the Standard Model retains only the left-chiral vector and axial-vector interactions.

### Standard-Model Beta-Decay Matrix Element

Let $\psi_u$ and $\psi_d$ denote the up- and down-quark fields. The Standard Model quark and lepton charged-current operators are

$$
J_q^\mu=\bar\psi_u\gamma^\mu(1-\gamma^5)\psi_d,\qquad J_\ell^\mu=\bar\psi_e\gamma^\mu(1-\gamma^5)\psi_{\nu_e}.
$$

The quark current must be matched onto composite nucleons. Let $M_N$ be the nucleon mass and let $g_{\mathrm V}$, $g_{\mathrm M}$, $g_{\mathrm A}$, and $g_{\mathrm P}$ denote the vector, weak-magnetism, axial-vector, and induced-pseudoscalar form factors. The nucleon current is

$$
J_{\mathrm h}^\mu=\bar\psi_p\Gamma_{\mathrm h}^\mu(q)\psi_n,
$$

with

$$
\Gamma_{\mathrm h}^\mu(q)=g_{\mathrm V}(q^2)\gamma^\mu+\frac{ig_{\mathrm M}(q^2)}{2M_Nc^2}\sigma^{\mu\nu}q_\nu-g_{\mathrm A}(q^2)\gamma^\mu\gamma^5-\frac{g_{\mathrm P}(q^2)}{2M_Nc^2}q^\mu\gamma^5.
$$

The weak-magnetism and induced-pseudoscalar terms describe the momentum dependence of the composite-nucleon current. They do not represent fundamental tensor or pseudoscalar four-fermion interactions.

The Standard Model effective Hamiltonian density is

$$
\boxed{\mathcal H_\beta^{\mathrm{SM}}=\frac{G_{\mathrm F}V_{ud}}{\sqrt2}J_{\mathrm h}^\mu J_{\ell\mu}+\mathrm{h.c.}}.
$$

Let $\varphi_e(\mathbf r)$ and $\varphi_{\bar\nu_e}(\mathbf r)$ be the outgoing electron and antineutrino spinor wave functions. The corresponding lepton-current matrix element is

$$
j_\ell^\mu(\mathbf r)\equiv\langle e^-,\bar\nu_e|J_\ell^\mu(\mathbf r)|0\rangle=\bar\varphi_e(\mathbf r)\gamma^\mu(1-\gamma^5)\varphi_{\bar\nu_e}(\mathbf r).
$$

Thus, $J_\ell^\mu$ is a field operator, while $j_\ell^\mu$ is a complex four-vector function for specified lepton momenta and spin states. The nuclear beta-decay matrix element is

$$
\boxed{\mathcal M_{fi}=\frac{G_{\mathrm F}V_{ud}}{\sqrt2}\int d^3\mathbf r\,j_{\ell\mu}(\mathbf r)\langle\Psi_D|J_{\mathrm h}^\mu(\mathbf r)|\Psi_P\rangle}.
$$

Let $g_{\mathrm V}\equiv g_{\mathrm V}(0)$ and $g_{\mathrm A}\equiv g_{\mathrm A}(0)$. To leading order in the nucleon velocity and $q/(M_Nc^2)$,

$$
J_{\mathrm h}^0(\mathbf r)\simeq g_{\mathrm V}\sum_{a=1}^A\tau_\beta(a)\delta^3(\mathbf r-\mathbf r_a),
$$

$$
\mathbf J_{\mathrm h}(\mathbf r)\simeq-g_{\mathrm A}\sum_{a=1}^A\boldsymbol\sigma(a)\tau_\beta(a)\delta^3(\mathbf r-\mathbf r_a).
$$

Using $g_{\mu\nu}=\operatorname{diag}(1,-1,-1,-1)$ gives

$$
\mathcal M_{fi}\simeq\frac{G_{\mathrm F}V_{ud}}{\sqrt2}\sum_{a=1}^A\left\langle\Psi_D\left|g_{\mathrm V}j_\ell^0(\mathbf r_a)\tau_\beta(a)+g_{\mathrm A}\mathbf j_\ell(\mathbf r_a)\cdot\boldsymbol\sigma(a)\tau_\beta(a)\right|\Psi_P\right\rangle.
$$

For $qR\ll1$, the lepton current is approximately constant across the nucleus. The matrix element reduces to the Fermi and Gamow--Teller operators defined above:

$$
\mathcal M_{fi}\simeq\frac{G_{\mathrm F}V_{ud}}{\sqrt2}\left\langle\Psi_D\left|g_{\mathrm V}j_\ell^0\hat O_{\mathrm F}+g_{\mathrm A}\mathbf j_\ell\cdot\hat{\mathbf O}_{\mathrm{GT}}\right|\Psi_P\right\rangle.
$$

For finite $qR$, the lepton wave functions and nuclear current are expanded into angular-momentum channels. Let

$$
X\in\{\mathrm V,\mathrm A,\mathrm{WM},\mathrm{IP}\}
$$

label the vector, axial-vector, weak-magnetism, and induced-pseudoscalar contributions. For each allowed $(L_\beta,S_\beta,J_\beta)$ channel, the nuclear input is

$$
M_{XJ_\beta}(q)\equiv\langle\Psi_D,J_D\|\hat O_{XJ_\beta}(q)\|\Psi_P,J_P\rangle.
$$

The orbital and spin labels $L_\beta,S_\beta$ enter the lepton angular functions and their coupling to the nuclear tensor of rank $J_\beta$. For fixed nuclear substates and lepton spin states, write

$$
\mathcal M_{fi}(M_P,M_D;s_e,s_{\bar\nu})=\sum_{L_\beta,S_\beta,J_\beta,X}\mathcal M_{L_\beta S_\beta J_\beta X}(M_P,M_D;s_e,s_{\bar\nu}).
$$

For an unpolarized parent and unobserved final spins,

$$
\boxed{\overline{|\mathcal M_{fi}|^2}=\frac{1}{2J_P+1}\sum_{M_P,M_D}\sum_{s_e,s_{\bar\nu}}\left|\mathcal M_{fi}(M_P,M_D;s_e,s_{\bar\nu})\right|^2}.
$$

In the massless-neutrino limit, the antineutrino sum contains only the right-helicity state.
