# Axial Gaussian Gogny Matrix Elements

Let $a\in\{1,2,3,4\}$ label a one-particle state with spatial wave function $\phi_a(\boldsymbol r)$, spin projection $\sigma_a$, and isospin $q_a$. Let $\delta_{ab}$ denote the Kronecker delta.

For $g=1,\ldots,N_g$, define

$$
f_g(r_{12})=\exp\left(-\frac{r_{12}^2}{\mu_g^2}\right),\qquad r_{12}=\lvert\boldsymbol r_1-\boldsymbol r_2\rvert.
$$

The finite-range Gogny interaction is

$$
V_{\mathrm G}=\sum_{g=1}^{N_g}\left(W_g+B_gP_\sigma-H_gP_\tau-M_gP_\sigma P_\tau\right)f_g(r_{12}),
$$

where $P_\sigma$ and $P_\tau$ exchange spin and isospin, respectively. The spatial integrals follow [Axial Gaussian Matrix Elements](axial_gaussian_kernel.md).

Each channel is antisymmetrized by

$$
\bar v_{12;34}=\langle12|V_{\mathrm G}|34\rangle-\langle12|V_{\mathrm G}|43\rangle.
$$

The same-isospin configuration satisfies $q_1=q_2=q_3=q_4$. The direct cross-isospin configuration satisfies $q_1=q_3\ne q_2=q_4$. The exchanged cross-isospin configuration follows from $3\leftrightarrow4$.

### $W$ Channel

The $W$ channel contains no exchange operator:

$$
V_W=\sum_{g=1}^{N_g}W_gf_g(r_{12}).
$$

Define its spatial matrix element by

$$
V^W_{12,34}=\sum_{g=1}^{N_g}W_g\int\mathrm d^3r_1\,\mathrm d^3r_2\,\phi_1^*(\boldsymbol r_1)\phi_2^*(\boldsymbol r_2)f_g(r_{12})\phi_3(\boldsymbol r_1)\phi_4(\boldsymbol r_2).
$$

Antisymmetrization gives

$$
\boxed{\bar v^W_{12;34}=\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}\delta_{q_1q_3}\delta_{q_2q_4}V^W_{12,34}-\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}\delta_{q_1q_4}\delta_{q_2q_3}V^W_{12,43}}.
$$

For $q_1=q_2=q_3=q_4$,

$$
\boxed{\bar v^W_{\mathrm{same}}=\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}V^W_{12,34}-\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}V^W_{12,43}}.
$$

For $q_1=q_3\ne q_2=q_4$,

$$
\boxed{\bar v^W_{\mathrm{cross}}=\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}V^W_{12,34}}.
$$

### $B$ Channel

The $B$ channel exchanges spin:

$$
V_B=\sum_{g=1}^{N_g}B_gP_\sigma f_g(r_{12}).
$$

Define its spatial matrix element by

$$
V^B_{12,34}=\sum_{g=1}^{N_g}B_g\int\mathrm d^3r_1\,\mathrm d^3r_2\,\phi_1^*(\boldsymbol r_1)\phi_2^*(\boldsymbol r_2)f_g(r_{12})\phi_3(\boldsymbol r_1)\phi_4(\boldsymbol r_2).
$$

Antisymmetrization gives

$$
\boxed{\bar v^B_{12;34}=\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}\delta_{q_1q_3}\delta_{q_2q_4}V^B_{12,34}-\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}\delta_{q_1q_4}\delta_{q_2q_3}V^B_{12,43}}.
$$

For $q_1=q_2=q_3=q_4$,

$$
\boxed{\bar v^B_{\mathrm{same}}=\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}V^B_{12,34}-\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}V^B_{12,43}}.
$$

For $q_1=q_3\ne q_2=q_4$,

$$
\boxed{\bar v^B_{\mathrm{cross}}=\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}V^B_{12,34}}.
$$

### $H$ Channel

The $H$ channel exchanges isospin and carries a negative coefficient:

$$
V_H=-\sum_{g=1}^{N_g}H_gP_\tau f_g(r_{12}).
$$

Define its unsigned spatial matrix element by

$$
V^H_{12,34}=\sum_{g=1}^{N_g}H_g\int\mathrm d^3r_1\,\mathrm d^3r_2\,\phi_1^*(\boldsymbol r_1)\phi_2^*(\boldsymbol r_2)f_g(r_{12})\phi_3(\boldsymbol r_1)\phi_4(\boldsymbol r_2).
$$

Antisymmetrization gives

$$
\boxed{\bar v^H_{12;34}=-\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}\delta_{q_1q_4}\delta_{q_2q_3}V^H_{12,34}+\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}\delta_{q_1q_3}\delta_{q_2q_4}V^H_{12,43}}.
$$

For $q_1=q_2=q_3=q_4$,

$$
\boxed{\bar v^H_{\mathrm{same}}=-\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}V^H_{12,34}+\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}V^H_{12,43}}.
$$

For $q_1=q_3\ne q_2=q_4$,

$$
\boxed{\bar v^H_{\mathrm{cross}}=\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}V^H_{12,43}}.
$$

### $M$ Channel

The $M$ channel exchanges both spin and isospin:

$$
V_M=-\sum_{g=1}^{N_g}M_gP_\sigma P_\tau f_g(r_{12}).
$$

Define its unsigned spatial matrix element by

$$
V^M_{12,34}=\sum_{g=1}^{N_g}M_g\int\mathrm d^3r_1\,\mathrm d^3r_2\,\phi_1^*(\boldsymbol r_1)\phi_2^*(\boldsymbol r_2)f_g(r_{12})\phi_3(\boldsymbol r_1)\phi_4(\boldsymbol r_2).
$$

Antisymmetrization gives

$$
\boxed{\bar v^M_{12;34}=-\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}\delta_{q_1q_4}\delta_{q_2q_3}V^M_{12,34}+\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}\delta_{q_1q_3}\delta_{q_2q_4}V^M_{12,43}}.
$$

For $q_1=q_2=q_3=q_4$,

$$
\boxed{\bar v^M_{\mathrm{same}}=-\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}V^M_{12,34}+\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}V^M_{12,43}}.
$$

For $q_1=q_3\ne q_2=q_4$,

$$
\boxed{\bar v^M_{\mathrm{cross}}=\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}V^M_{12,43}}.
$$

Summing the four channels gives

$$
\boxed{\bar v_{\mathrm{same}}=\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}\left(V^W_{12,34}-V^B_{12,43}-V^H_{12,34}+V^M_{12,43}\right)+\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}\left(-V^W_{12,43}+V^B_{12,34}+V^H_{12,43}-V^M_{12,34}\right)}.
$$

For the direct cross-isospin configuration,

$$
\boxed{\bar v_{\mathrm{cross}}=\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}\left(V^W_{12,34}+V^M_{12,43}\right)+\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}\left(V^B_{12,34}+V^H_{12,43}\right)}.
$$

These matrix elements supply the representative $\Gamma$ and $\Delta$ fields defined in [Axial HFB](../../notes/methods/hfb_axial.md#representative-gamma-and-delta-fields).
