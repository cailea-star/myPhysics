### General Scattering Theory

##### Asymptotic States and Scattering Matrix

Take $t_0=0$ as the reference time, where all pictures coincide. Let $H=H_0+V$ be time independent. In the [interaction picture](quantum_dynamics.md#interaction-picture),

$$
|\psi_I(t)\rangle=U_I(t,0)|\psi_0\rangle,\qquad |\psi_0\rangle=|\psi_S(0)\rangle=|\psi_I(0)\rangle.
$$

The inverse evolution brings the interaction-picture state back to the reference time:

$$
\boxed{\Omega(t)=U_I(0,t),\qquad \Omega(t)|\psi_I(t)\rangle=|\psi_0\rangle}.
$$

Define the incoming asymptotic state and the corresponding Møller operator by

$$
\Omega_+=\lim_{t\to-\infty}\Omega(t),\qquad |\psi_{\mathrm{in}}\rangle=\lim_{t\to-\infty}|\psi_I(t)\rangle,\qquad \Omega_+|\psi_{\mathrm{in}}\rangle=|\psi_0\rangle.
$$

Define the outgoing asymptotic state and the corresponding Møller operator by

$$
\Omega_-=\lim_{t\to+\infty}\Omega(t),\qquad |\psi_{\mathrm{out}}\rangle=\lim_{t\to+\infty}|\psi_I(t)\rangle,\qquad \Omega_-|\psi_{\mathrm{out}}\rangle=|\psi_0\rangle.
$$

For scattering wave packets, the interaction is asymptotically negligible:

$$
V_I(t)=e^{iH_0t/\hbar}Ve^{-iH_0t/\hbar},\qquad \lim_{t\to-\infty}\left\|V_I(t)|\psi_{\mathrm{in}}\rangle\right\|=0,\qquad \lim_{t\to+\infty}\left\|V_I(t)|\psi_{\mathrm{out}}\rangle\right\|=0.
$$

For a free eigenstate $H_0|k\rangle=E_k|k\rangle$, define the incoming and outgoing scattering states at the reference time by

$$
\boxed{|\psi_k^+\rangle=\Omega_+|k\rangle,\qquad |\psi_k^-\rangle=\Omega_-|k\rangle}.
$$

If the Hamiltonian is time independent, the asymptotic free states and the interacting scattering state have the same energy. The intertwining relation gives

$$
H\Omega_\pm=\Omega_\pm H_0,\qquad H|\psi_k^\pm\rangle=E_k|\psi_k^\pm\rangle.
$$

The Møller operators preserve the continuum normalization:

$$
\langle\psi_{k'}^\pm|\psi_k^\pm\rangle=\langle k'|k\rangle.
$$

Define the scattering matrix by

$$
S=\Omega_-^\dagger\Omega_+,\qquad |\psi_{\mathrm{out}}\rangle=S|\psi_{\mathrm{in}}\rangle,\qquad \langle k'|S|k\rangle=\langle\psi_{k'}^-|\psi_k^+\rangle.
$$

Using the [interaction-picture Dyson equation](quantum_dynamics.md#interaction-picture),

$$
U_I(t,t_0)=I-\frac{i}{\hbar}\int_{t_0}^{t}dt'\,V_I(t')U_I(t',t_0),
$$

$$
S=U_I(+\infty,-\infty)=I-\frac{i}{\hbar}\int_{-\infty}^{+\infty}dt\,V_I(t)U_I(t,-\infty).
$$

##### Transition Operator and S-Matrix Elements

The Møller operator defines the transition operator on each energy shell by

$$
\boxed{T^\pm(E_k)|k\rangle\equiv V\Omega_\pm|k\rangle=V|\psi_k^\pm\rangle}.
$$

For an incoming free state, write $T(E)\equiv T^+(E)$.

Expanding the interaction-picture operators gives

$$
\langle k'|V_I(t)U_I(t,-\infty)|k\rangle=\langle k'|\Big[U_0^\dagger(t,0)VU_0(t,0)\Big]\Big[U_I(t,0)\Omega_+\Big]|k\rangle.
$$

$$
=\langle k'|\Big[U_0^\dagger(t,0)VU_0(t,0)\Big]\Big[U_0^\dagger(t,0)U(t,0)\Big]\Omega_+|k\rangle=\langle k'|U_0^\dagger(t,0)VU(t,0)|\psi_k^+\rangle.
$$

Using $H_0|k'\rangle=E_{k'}|k'\rangle$ and $H|\psi_k^+\rangle=E_k|\psi_k^+\rangle$,

$$
=e^{i(E_{k'}-E_k)t/\hbar}\langle k'|V|\psi_k^+\rangle=e^{i(E_{k'}-E_k)t/\hbar}\langle k'|T(E_k)|k\rangle.
$$

Therefore,

$$
\langle k'|S|k\rangle=\langle k'|k\rangle-\frac{i}{\hbar}\int_{-\infty}^{+\infty}dt\,e^{i(E_{k'}-E_k)t/\hbar}\langle k'|T(E_k)|k\rangle.
$$

Using

$$
\int_{-\infty}^{+\infty}dt\,e^{i(E_{k'}-E_k)t/\hbar}=2\pi\hbar\,\delta(E_{k'}-E_k),
$$

gives

$$
\boxed{\langle k'|S|k\rangle=\langle k'|k\rangle-2\pi i\,\delta(E_{k'}-E_k)\langle k'|T(E_k)|k\rangle}.
$$

##### Lippmann–Schwinger Equation

Subtracting the free equation from the interacting equation gives

$$
(E_k-H_0)|\psi_k^\pm\rangle=V|\psi_k^\pm\rangle,\qquad (E_k-H_0)|k\rangle=0.
$$

Therefore,

$$
(E_k-H_0)\left(|\psi_k^\pm\rangle-|k\rangle\right)=V|\psi_k^\pm\rangle.
$$

Define the free [quantum Green function](../mathematical_physics/green_functions.md#quantum-green-function) by

$$
G_0^\pm(E)=\left(E-H_0\pm i0\right)^{-1}.
$$

Here $+$ selects the outgoing scattered wave, while $-$ selects the incoming scattered wave. The scattered correction can be written as

$$
|\psi_k^\pm\rangle=|k\rangle+G_0^\pm(E_k)V|\psi_k^\pm\rangle=|k\rangle+G_0^\pm(E_k)T^\pm(E_k)|k\rangle.
$$

Substitution gives

$$
T^\pm(E_k)|k\rangle=V|k\rangle+VG_0^\pm(E_k)T^\pm(E_k)|k\rangle,
$$

and therefore

$$
\boxed{T^\pm(E)=V+VG_0^\pm(E)T^\pm(E)}.
$$

### Potential Scattering

##### Scattering Amplitude and Cross Section

For elastic scattering by a localized potential,

$$
H_0=\frac{\mathbf p^2}{2M},\qquad H=H_0+V(\mathbf r),\qquad E_k=\frac{\hbar^2k^2}{2M},
$$

where $M$ is the reduced mass. The probability current is

$$
\mathbf j[\psi]=\frac{\hbar}{M}\operatorname{Im}\left(\psi^*\nabla\psi\right).
$$

Let $j_{\mathrm{in}}$ be the incident flux and $j_{\mathrm{sc}}=\hat{\mathbf r}\cdot\mathbf j[\psi_{\mathrm{sc}}]$ the outward scattered current. Define $d\sigma$ as the effective incident area corresponding to the scattered flux through $dA=r^2d\Omega$:

$$
j_{\mathrm{in}}\,d\sigma=j_{\mathrm{sc}}\,dA,\qquad dA=r^2d\Omega.
$$

Therefore,

$$
\boxed{\frac{d\sigma}{d\Omega}\equiv\lim_{r\to\infty}\frac{r^2j_{\mathrm{sc}}(\mathbf r)}{j_{\mathrm{in}}},\qquad \sigma_{\mathrm{tot}}=\int d\Omega\,\frac{d\sigma}{d\Omega}}.
$$

Use the plane-wave normalization

$$
\langle\mathbf r|\mathbf k\rangle=\frac{e^{i\mathbf k\cdot\mathbf r}}{(2\pi)^{3/2}},\qquad \langle\mathbf k'|\mathbf k\rangle=\delta^{(3)}(\mathbf k'-\mathbf k).
$$

Insert the momentum completeness relation into the Lippmann–Schwinger equation:

$$
\boxed{\psi_{\mathbf k}^+(\mathbf r)=\langle\mathbf r|\psi_{\mathbf k}^+\rangle=\langle\mathbf r|\mathbf k\rangle+\int d^3q\,\langle\mathbf r|\mathbf q\rangle\langle\mathbf q|G_0^+(E_k)T(E_k)|\mathbf k\rangle}.
$$

Using the momentum-space representation of the [quantum Green function](../mathematical_physics/green_functions.md#quantum-green-function),

$$
\langle\mathbf q|G_0^+(E)|\mathbf q'\rangle=\frac{\delta^{(3)}(\mathbf q-\mathbf q')}{E-E_q+i0},\qquad E_q=\frac{\hbar^2q^2}{2M},
$$

gives

$$
\psi_{\mathbf k}^+(\mathbf r)=\frac{e^{i\mathbf k\cdot\mathbf r}}{(2\pi)^{3/2}}+\int\frac{d^3q}{(2\pi)^{3/2}}\,\frac{e^{i\mathbf q\cdot\mathbf r}}{E_k-E_q+i0}\langle\mathbf q|T(E_k)|\mathbf k\rangle.
$$

In the far-field region, the transition matrix element is approximated by its on-shell value:

$$
r\to\infty:\qquad \langle\mathbf q|T(E_k)|\mathbf k\rangle\simeq\langle\mathbf k'|T(E_k)|\mathbf k\rangle,\qquad \mathbf k'=k\hat{\mathbf r}.
$$

Using

$$
E_k-E_q=\frac{\hbar^2}{2M}(k^2-q^2),\qquad\int\frac{d^3q}{(2\pi)^3}\frac{e^{i\mathbf q\cdot\mathbf r}}{k^2-q^2+i0}=-\frac{e^{ikr}}{4\pi r},
$$

gives

$$
\left.\psi_{\mathbf k}^+(\mathbf r)\right|_{r\to\infty}=\frac{1}{(2\pi)^{3/2}}\left[e^{i\mathbf k\cdot\mathbf r}-\frac{(2\pi)^2M}{\hbar^2}\langle\mathbf k'|T(E_k)|\mathbf k\rangle\frac{e^{ikr}}{r}\right].
$$

Define the scattering amplitude as the coefficient of the outgoing spherical wave:

$$
\boxed{\left.\psi_{\mathbf k}^+(\mathbf r)\right|_{r\to\infty}=\frac{1}{(2\pi)^{3/2}}\left[e^{i\mathbf k\cdot\mathbf r}+f(\mathbf k',\mathbf k)\frac{e^{ikr}}{r}\right],\qquad f(\mathbf k',\mathbf k)=-\frac{(2\pi)^2M}{\hbar^2}\langle\mathbf k'|T(E_k)|\mathbf k\rangle}.
$$

The incident and scattered currents are

$$
j_{\mathrm{in}}=\frac{\hbar k}{M(2\pi)^3},\qquad r^2j_{\mathrm{sc}}=\frac{\hbar k}{M(2\pi)^3}|f(\mathbf k',\mathbf k)|^2.
$$

Therefore, the differential cross section and optical theorem are

$$
\boxed{\frac{d\sigma}{d\Omega}=|f(\mathbf k',\mathbf k)|^2,\qquad \sigma_{\mathrm{tot}}=\frac{4\pi}{k}\operatorname{Im}f(\mathbf k,\mathbf k)}.
$$

Using $T(E_k)|\mathbf k\rangle=V|\psi_{\mathbf k}^+\rangle$,

$$
\boxed{f(\mathbf k',\mathbf k)=-\frac{(2\pi)^2M}{\hbar^2}\langle\mathbf k'|T(E_k)|\mathbf k\rangle=-\frac{\sqrt{2\pi}M}{\hbar^2}\int d^3r\,e^{-i\mathbf k'\cdot\mathbf r}V(\mathbf r)\psi_{\mathbf k}^+(\mathbf r)}.
$$

##### Born Approximation

Iterating the Lippmann–Schwinger equation gives the Born series:

$$
|\psi_{\mathbf k}^+\rangle=\sum_{n=0}^{\infty}\left[G_0^+(E_k)V\right]^n|\mathbf k\rangle,\qquad T(E)=V\sum_{n=0}^{\infty}\left[G_0^+(E)V\right]^n=V+VG_0^+(E)V+\cdots.
$$

The first Born approximation retains the leading term:

$$
|\psi_{\mathbf k}^+\rangle\simeq|\mathbf k\rangle,\qquad T_{\mathrm B}(E)=V.
$$

Therefore,

$$
f_{\mathrm B}(\mathbf k',\mathbf k)=-\frac{(2\pi)^2M}{\hbar^2}\langle\mathbf k'|V|\mathbf k\rangle=-\frac{M}{2\pi\hbar^2}\int d^3r\,e^{-i\mathbf Q\cdot\mathbf r}V(\mathbf r),
$$

where the momentum transfer is

$$
\mathbf Q=\mathbf k'-\mathbf k,\qquad Q=2k\sin\frac{\theta}{2},\qquad \cos\theta=\hat{\mathbf k}'\cdot\hat{\mathbf k}.
$$

For a central potential $V(\mathbf r)=V(r)$,

$$
\int d^3r\,e^{-i\mathbf Q\cdot\mathbf r}V(r)=4\pi\int_0^\infty dr\,r^2V(r)\frac{\sin(Qr)}{Qr}.
$$

Hence,

$$
\boxed{f_{\mathrm B}(Q)=-\frac{2M}{\hbar^2}\int_0^\infty dr\,r^2V(r)\frac{\sin(Qr)}{Qr},\qquad \frac{d\sigma_{\mathrm B}}{d\Omega}=|f_{\mathrm B}(Q)|^2}.
$$

The Born series requires

$$
\left\|G_0^+(E)V\right\|\ll1,
$$

and is typically reliable for weak short-range potentials or high incident energies.

### Partial-Wave Scattering

##### Phase Shifts and Partial Cross Sections

For elastic scattering by a real central potential, energy and angular momentum are conserved. Let

$$
E=\frac{\hbar^2k^2}{2M}.
$$

The scattering operator is diagonal in the angular-momentum basis. Elastic unitarity gives

$$
S|E l m\rangle=S_l(E)|E l m\rangle,\qquad |S_l(E)|=1.
$$

Define the phase shift $\delta_l(E)$ by

$$
\boxed{S|E l m\rangle=S_l(E)|E l m\rangle=e^{2i\delta_l(E)}|E l m\rangle}.
$$

Define the energy-direction states by

$$
|E,\hat{\mathbf k}\rangle=\frac{\sqrt{Mk}}{\hbar}|\mathbf k\rangle,\qquad \langle E',\hat{\mathbf k}'|E,\hat{\mathbf k}\rangle=\delta(E'-E)\delta(\Omega'-\Omega).
$$

For elastic scattering,

$$
|\mathbf k'|=|\mathbf k|=k,\qquad \cos\theta=\hat{\mathbf k}'\cdot\hat{\mathbf k}.
$$

The fixed-energy angular kernels are

$$
\langle E',\hat{\mathbf k}'|S|E,\hat{\mathbf k}\rangle=\delta(E'-E)S_E(\hat{\mathbf k}',\hat{\mathbf k}),\qquad T_E(\hat{\mathbf k}',\hat{\mathbf k})=\langle E,\hat{\mathbf k}'|T(E)|E,\hat{\mathbf k}\rangle,
$$

and satisfy

$$
S_E(\hat{\mathbf k}',\hat{\mathbf k})=\delta(\Omega'-\Omega)-2\pi i\,T_E(\hat{\mathbf k}',\hat{\mathbf k}).
$$

For a central potential, the fixed-energy kernels are rotationally invariant. Using the [angular expansions](../mathematical_physics/spherical_coordinates.md#angular-expansions),

$$
S_E(\hat{\mathbf k}',\hat{\mathbf k})=\sum_{l=0}^{\infty}\frac{2l+1}{4\pi}S_l(E)P_l(\cos\theta),\qquad T_E(\hat{\mathbf k}',\hat{\mathbf k})=\sum_{l=0}^{\infty}\frac{2l+1}{4\pi}T_l(E)P_l(\cos\theta),
$$

where $T_l(E)=\langle E l m|T(E)|E l m\rangle$. Since the identity operator contributes $1$ to every partial-wave channel,

$$
S_l(E)=1-2\pi i\,T_l(E).
$$

The elastic scattering amplitude is

$$
\boxed{f(\theta)=-\frac{(2\pi)^2}{k}T_E(\hat{\mathbf k}',\hat{\mathbf k})=\left.f(\mathbf k',\mathbf k)\right|_{|\mathbf k'|=|\mathbf k|=k,\ \hat{\mathbf k}'\cdot\hat{\mathbf k}=\cos\theta}}.
$$

Its partial-wave expansion is

$$
f(\theta)=\sum_{l=0}^{\infty}(2l+1)f_l(k)P_l(\cos\theta),\qquad f_l(k)=-\frac{\pi}{k}T_l(E)=\frac{S_l(E)-1}{2ik}=\frac{e^{i\delta_l(E)}\sin\delta_l(E)}{k}.
$$

Therefore,

$$
\boxed{f(\theta)=\frac{1}{2ik}\sum_{l=0}^{\infty}(2l+1)\left[e^{2i\delta_l(E)}-1\right]P_l(\cos\theta)}.
$$

Using $d\sigma/d\Omega=|f(\theta)|^2$ and the [Legendre orthogonality relation](../mathematical_physics/spherical_coordinates.md#angular-expansions), the total cross section decomposes into independent partial waves:

$$
\sigma_{\mathrm{tot}}=\int d\Omega\,|f(\theta)|^2=4\pi\sum_{l=0}^{\infty}(2l+1)|f_l(k)|^2\equiv\sum_{l=0}^{\infty}\sigma_l.
$$

Therefore,

$$
\boxed{\sigma_l=4\pi(2l+1)|f_l(k)|^2=\frac{4\pi}{k^2}(2l+1)\sin^2\delta_l(E)}.
$$

##### Phase-Shift Calculation

For a short-range central potential, let $u_l(r)$ be the reduced radial wavefunction. Each partial wave satisfies

$$
\left[-\frac{\hbar^2}{2M}\frac{d^2}{dr^2}+\frac{\hbar^2l(l+1)}{2Mr^2}+V(r)\right]u_l(r)=E\,u_l(r),\qquad u_l(r)\underset{r\to0}{\propto}r^{l+1}.
$$

Choose a matching radius $r_m$ such that $V(r)\simeq0$ for $r\geq r_m$. Define the Riccati–Bessel functions by

$$
\hat j_l(x)=xj_l(x),\qquad \hat n_l(x)=xn_l(x).
$$

The exterior solution is

$$
u_l(r)=A_l\left[\cos\delta_l(E)\,\hat j_l(kr)-\sin\delta_l(E)\,\hat n_l(kr)\right]\underset{r\to\infty}{\longrightarrow}A_l\sin\left(kr-\frac{l\pi}{2}+\delta_l(E)\right).
$$

Integrate the regular solution outward and define its logarithmic derivative at $r_m$ by

$$
L_l(E)=\frac{u_l'(r_m)}{u_l(r_m)}.
$$

Matching the logarithmic derivatives gives

$$
\boxed{\tan\delta_l(E)=\frac{k\hat j_l'(kr_m)-L_l(E)\hat j_l(kr_m)}{k\hat n_l'(kr_m)-L_l(E)\hat n_l(kr_m)}}.
$$

Here the primes on $\hat j_l$ and $\hat n_l$ denote derivatives with respect to their arguments. Integration to infinity is unnecessary; the extracted phase shift should be stable under changes of $r_m$ within the exterior region.

$$
V(r)\longrightarrow u_l(r_m),u_l'(r_m)\longrightarrow\delta_l(E)\longrightarrow S_l(E)\longrightarrow f_l(k)\longrightarrow\sigma_l.
$$

##### Resonant Scattering

An isolated resonance appears when the phase shift changes rapidly through $\pi/2$:

$$
\delta_l(E_R)=\frac{\pi}{2},\qquad \cot\delta_l(E_R)=0,
$$

where $E_R$ is the resonance energy. Neglecting a slowly varying background, parameterize the phase shift near $E_R$ by

$$
\boxed{\cot\delta_l(E)=\frac{E_R-E}{\Gamma/2}},
$$

where $\Gamma$ is the resonance width. Using the phase-shift parameterization,

$$
S_l(E)=e^{2i\delta_l(E)}=\frac{\cot\delta_l(E)+i}{\cot\delta_l(E)-i}=\frac{E-E_R-i\Gamma/2}{E-E_R+i\Gamma/2}.
$$

The pole of $S_l(E)$ is located at

$$
E=E_R-\frac{i\Gamma}{2},\qquad \tau=\frac{\hbar}{\Gamma},
$$

where $\tau$ is the resonance lifetime. The phase-shift parameterization gives

$$
\sin^2\delta_l(E)=\frac{(\Gamma/2)^2}{(E-E_R)^2+(\Gamma/2)^2},
$$

and therefore

$$
\boxed{\sigma_l(E)=\frac{4\pi}{k^2}(2l+1)\frac{(\Gamma/2)^2}{(E-E_R)^2+(\Gamma/2)^2}}.
$$

At the resonance energy,

$$
\delta_l(E_R)=\frac{\pi}{2},\qquad \sigma_l(E_R)=\frac{4\pi}{k^2}(2l+1).
$$

Thus,

$$
\delta_l(E)\text{ crosses }\frac{\pi}{2}\longrightarrow S_l(E)\text{ has a pole}\longrightarrow\sigma_l(E)\text{ has a resonance peak}.
$$
