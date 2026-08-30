# Beta Decay

### Relativistic Kinematics

Use natural units with $c=1$. Let $P$, $D$, $e^-$, and $\bar\nu_e$ denote the parent nucleus, daughter nucleus, electron, and electron antineutrino, with masses $M_P$, $M_D$, $m_e$, and $m_\nu$.

The decay process is

$$
P\rightarrow D+e^-+\bar\nu_e.
$$

For $m_\nu=0$, the decay energy is

$$
Q_{\beta^-}=M_P-M_D-m_e=T_D+T_e+T_\nu.
$$

Four-momentum conservation and the on-shell condition are

$$
p_P^\mu=p_D^\mu+p_e^\mu+p_\nu^\mu,\qquad p^\mu p_\mu=E^2-\mathbf p^2=m^2.
$$

In the parent-nucleus rest frame, $p_P^\mu=(M_P,\mathbf 0)$. Define the invariant mass of the daughter nucleus and neutrino by

$$
s_{D\nu}=(p_D+p_\nu)^2=(p_P-p_e)^2=M_P^2+m_e^2-2M_PE_e.
$$

Therefore,

$$
E_e=\frac{M_P^2+m_e^2-s_{D\nu}}{2M_P}.
$$

The maximum electron energy occurs at the minimum invariant mass

$$
s_{D\nu}^{\min}=(M_D+m_\nu)^2.
$$

For $m_\nu=0$,

$$
E_e^{\max}=\frac{M_P^2+m_e^2-M_D^2}{2M_P}.
$$

The corresponding maximum kinetic energy is

$$
T_e^{\max}=E_e^{\max}-m_e=\frac{(M_P-m_e)^2-M_D^2}{2M_P}.
$$

Using $M_P=Q_{\beta^-}+M_D+m_e$,

$$
\boxed{T_e^{\max}=\frac{Q_{\beta^-}(Q_{\beta^-}+2M_D)}{2M_P}\simeq Q_{\beta^-}}.
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
Q_n=m_n-m_p-m_e\simeq0.782\,\mathrm{MeV}.
$$

**Beta-minus decay.**

$$
{}^{A}_{Z}P_N\rightarrow{}^{A}_{Z+1}D_{N-1}+e^-+\bar\nu_e.
$$

$$
\boxed{Q_{\beta^-}=M_{\mathrm{nuc}}(P)-M_{\mathrm{nuc}}(D)-m_e\simeq M_{\mathrm{atom}}(P)-M_{\mathrm{atom}}(D)}.
$$

**Beta-plus decay.**

$$
{}^{A}_{Z}P_N\rightarrow{}^{A}_{Z-1}D_{N+1}+e^++\nu_e.
$$

$$
\boxed{Q_{\beta^+}=M_{\mathrm{nuc}}(P)-M_{\mathrm{nuc}}(D)-m_e\simeq M_{\mathrm{atom}}(P)-M_{\mathrm{atom}}(D)-2m_e}.
$$

**Electron capture.** Let $e_x^-$ denote an electron captured from atomic shell $x$, with binding energy $B_x$.

$$
{}^{A}_{Z}P_N+e_x^-\rightarrow{}^{A}_{Z-1}D_{N+1}+\nu_e.
$$

The atomic mass difference is

$$
Q_{\mathrm{EC}}=M_{\mathrm{atom}}(P)-M_{\mathrm{atom}}(D).
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

Neglecting the daughter recoil, $E_D\simeq M_D$. Define

$$
E_0=E_e+E_\nu=E_P-E_D\simeq M_P-M_D.
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
E^2=p^2+m^2,\qquad E\,dE=p\,dp.
$$

Therefore,

$$
\int_0^\infty p_\nu^2dp_\nu\,\delta(E_0-E_e-E_\nu)=p_\nu E_\nu,\qquad p_e^2dp_e=p_eE_e\,dE_e.
$$

The spatial and angular integrations are

$$
\int_Vd^3\mathbf r_e\int_Vd^3\mathbf r_\nu=V^2,\qquad \int d\Omega_{\mathbf p_e}\int d\Omega_{\mathbf p_\nu}=(4\pi)^2.
$$

The electron momentum is

$$
p_e=\sqrt{E_e^2-m_e^2}.
$$

Thus,

$$
dn(e,\nu)=\frac{V^2}{4\pi^4\hbar^6}p_\nu E_\nu\,p_eE_e\,dE_e=\frac{V^2}{4\pi^4\hbar^6}(E_0-E_e)^2\sqrt{1-\frac{m_\nu^2}{(E_0-E_e)^2}}\,p_eE_e\,dE_e.
$$

Let $\varphi_e$ be the outgoing electron wave function, $Z_D$ the daughter charge number, $\alpha$ the fine-structure constant, and $v_e=p_e/E_e$. In the point-Coulomb approximation, the electron distortion is described by the Fermi function

$$
|\varphi_e(0)|^2\simeq\frac{1}{V}F(Z_D,p_e),\qquad F(Z_D,p_e)=\frac{2\pi\eta}{1-e^{-2\pi\eta}},\qquad \eta=\pm\frac{Z_D\alpha}{v_e}.
$$

The upper and lower signs correspond to $\beta^-$ and $\beta^+$ decay.

Define the reduced matrix element $M_{fi}'$ by

$$
|M_{fi}|^2=\frac{1}{V^2}F(Z_D,p_e)|M_{fi}'|^2.
$$

The electron energy spectrum is

$$
\frac{d\lambda_\beta}{dE_e}=\frac{|M_{fi}'|^2}{2\pi^3\hbar^7}F(Z_D,p_e)p_eE_e(E_0-E_e)^2\sqrt{1-\frac{m_\nu^2}{(E_0-E_e)^2}}.
$$

Within the recoil-free approximation,

$$
m_e\leq E_e\leq E_0-m_\nu.
$$

The total decay constant and half-life are

$$
\lambda_\beta=\int_{m_e}^{E_0-m_\nu}\frac{d\lambda_\beta}{dE_e}\,dE_e,\qquad T_{1/2}=\frac{\ln2}{\lambda_\beta}.
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

### Beta Decay in the Standard Model

Use $\hbar=c=1$ in this section. Beta decay is mediated by the charged weak current through a virtual $W$ boson:

$$
d\rightarrow u+W^{-*}\rightarrow u+e^-+\bar\nu_e,\qquad n\rightarrow p+W^{-*}\rightarrow p+e^-+\bar\nu_e.
$$

Let $q^\mu$ be the transferred four-momentum, $M_W$ the $W$-boson mass, $g_2$ the $SU(2)_L$ gauge coupling, and $G_{\mathrm F}$ the Fermi constant. For nuclear beta decay, $|q^2|\ll M_W^2$, so

$$
\frac{1}{q^2-M_W^2}\simeq-\frac{1}{M_W^2},\qquad \frac{G_{\mathrm F}}{\sqrt2}=\frac{g_2^2}{8M_W^2}.
$$

The interaction range $R_W\sim M_W^{-1}$ is much smaller than the nuclear radius. Low-energy $W$ exchange therefore reduces to a local four-fermion interaction.

Let $\gamma^\mu$ be the Dirac matrices,

$$
\gamma^5\equiv i\gamma^0\gamma^1\gamma^2\gamma^3,\qquad \sigma^{\mu\nu}\equiv\frac{i}{2}[\gamma^\mu,\gamma^\nu].
$$

The independent Lorentz structures are

$$
\Gamma_{\mathrm S}=1,\qquad \Gamma_{\mathrm P}=\gamma^5,\qquad \Gamma_{\mathrm V}=\gamma^\mu,\qquad \Gamma_{\mathrm A}=\gamma^\mu\gamma^5,\qquad \Gamma_{\mathrm T}=\sigma^{\mu\nu}.
$$

Let $\psi_n$, $\psi_p$, $\psi_e$, and $\psi_{\nu_e}$ denote the neutron, proton, electron, and electron-neutrino fields, with $\bar\psi\equiv\psi^\dagger\gamma^0$. Let $C_i$ and $C_i'$ be general coupling coefficients, and let $\mathrm{h.c.}$ denote the Hermitian conjugate. A general local Lorentz-invariant beta interaction is

$$
\mathcal H_\beta^{\mathrm{general}}=\frac{G_{\mathrm F}}{\sqrt2}\sum_{i=\mathrm S,\mathrm P,\mathrm V,\mathrm A,\mathrm T}\left(\bar\psi_p\Gamma_i\psi_n\right)\left[\bar\psi_e\Gamma_i(C_i+C_i'\gamma^5)\psi_{\nu_e}\right]+\mathrm{h.c.}
$$

This expression is a model-independent Lorentz decomposition, not the Standard Model result. In the convention used above, the tree-level Standard Model imposes

$$
C_{\mathrm S}=C_{\mathrm S}'=C_{\mathrm P}=C_{\mathrm P}'=C_{\mathrm T}=C_{\mathrm T}'=0,\qquad C_{\mathrm V}'=-C_{\mathrm V},\qquad C_{\mathrm A}'=-C_{\mathrm A}.
$$

Let $\psi_u$ and $\psi_d$ denote the up- and down-quark fields, and let $J_q^\mu$ and $J_\ell^\mu$ denote the quark and lepton charged-current operators. The Standard Model fixes their $V-A$ structure:

$$
J_q^\mu=\bar\psi_u\gamma^\mu(1-\gamma^5)\psi_d,\qquad J_\ell^\mu=\bar\psi_e\gamma^\mu(1-\gamma^5)\psi_{\nu_e}.
$$

The coefficients of the vector and axial-vector terms are fixed at the fundamental-field level. The quark current must then be matched onto composite nucleon states. Let $g_{\mathrm V}$ and $g_{\mathrm A}$ be the vector and axial-vector nucleon coupling constants. At small momentum transfer,

$$
J_q^\mu\quad\xrightarrow{\text{nucleon matching}}\quad J_{\mathrm h}^\mu\simeq\bar\psi_p\gamma^\mu(g_{\mathrm V}-g_{\mathrm A}\gamma^5)\psi_n.
$$

Thus, $g_{\mathrm V}$ and $g_{\mathrm A}$ occur in the nucleon current but not in the fundamental quark or lepton currents. Let $V_{ud}$ be the CKM matrix element connecting the $u$ and $d$ quarks. In the low-energy normalization used above,

$$
g_{\mathrm F}\equiv\frac{G_{\mathrm F}V_{ud}}{\sqrt2}g_{\mathrm V},\qquad g_{\mathrm{GT}}\equiv\frac{G_{\mathrm F}V_{ud}}{\sqrt2}g_{\mathrm A}.
$$

Weak-magnetism and induced-pseudoscalar terms are omitted in this leading approximation. The Standard Model effective interaction is

$$
\boxed{\mathcal H_\beta^{\mathrm{SM}}=\frac{G_{\mathrm F}V_{ud}}{\sqrt2}J_{\mathrm h}^\mu J_{\ell\mu}+\mathrm{h.c.}}.
$$

Let $u_e$ and $v_{\bar\nu_e}$ denote the outgoing electron and antineutrino spinors. In relativistic normalization, the $\beta^-$-decay amplitude is

$$
\boxed{\mathcal M_{fi}=\frac{G_{\mathrm F}V_{ud}}{\sqrt2}\left[\bar u_e\gamma_\mu(1-\gamma^5)v_{\bar\nu_e}\right]\langle\Psi_D|J_{\mathrm h}^\mu|\Psi_P\rangle}.
$$

In the non-relativistic nucleon limit, the leading components of the hadronic current reduce to

$$
\boxed{\int d^3\mathbf r\,J_{\mathrm h}^0\longrightarrow g_{\mathrm V}\hat O_{\mathrm F},\qquad \int d^3\mathbf r\,\mathbf J_{\mathrm h}\longrightarrow-g_{\mathrm A}\hat{\mathbf O}_{\mathrm{GT}}}.
$$

Fermi transitions arise from the time component of the vector current, while Gamow--Teller transitions arise from the spatial component of the axial-vector current.

The factor $1-\gamma^5$ selects the left-chiral charged current. In the ultrarelativistic limit, beta decay emits predominantly left-helicity electrons and right-helicity antineutrinos. Since parity interchanges left- and right-chiral components, the $V-A$ interaction violates parity.
