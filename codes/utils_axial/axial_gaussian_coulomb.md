# Axial Gaussian Coulomb Matrix Elements

### Gaussian Coulomb Expansion

Let $\boldsymbol{r}_1$ and $\boldsymbol{r}_2$ denote the particle coordinates, and define $r_{12}=\lvert\boldsymbol{r}_1-\boldsymbol{r}_2\rvert>0$. The Coulomb interaction with coupling $e^2$ is

$$
v_{\mathrm{C}}(r_{12})=\frac{e^2}{r_{12}}.
$$

Let $\pi$ denote the circle constant and $\alpha\geq0$ an inverse-length variable. The Coulomb kernel satisfies

$$
\frac{1}{r_{12}}=\frac{2}{\sqrt{\pi}}\int_0^\infty\exp(-\alpha^2r_{12}^2)\,\mathrm{d}\alpha.
$$

Let $b_z$ and $b_r$ denote the axial and transverse oscillator lengths, and define $b_{\max}=\max(b_z,b_r)$.

For $t\in(0,1)$, apply

$$
\alpha=\frac{t}{b_{\max}\sqrt{1-t^2}},\qquad \mathrm{d}\alpha=\frac{\mathrm{d}t}{b_{\max}(1-t^2)^{3/2}}.
$$

The Coulomb interaction becomes

$$
v_{\mathrm{C}}(r_{12})=\frac{2e^2}{\sqrt{\pi}b_{\max}}\int_0^1(1-t^2)^{-3/2}\exp\left[-\frac{t^2}{1-t^2}\left(\frac{r_{12}}{b_{\max}}\right)^2\right]\,\mathrm{d}t.
$$

Let $N_g$ be the number of Gaussian terms. Let $x_g$ and $w_g$, for $g=1,\ldots,N_g$, denote the Gauss–Legendre nodes and weights on $[-1,1]$; see [Gaussian Quadrature](../utils/integration_gauss.md). Define

$$
t_g=\frac{x_g+1}{2},\qquad \omega_g=\frac{w_g}{2}.
$$

The discrete Coulomb weights are

$$
W_g=\frac{2e^2\omega_g}{\sqrt{\pi}b_{\max}(1-t_g^2)^{3/2}}.
$$

The implementation uses $N_g=9$. The finite Gaussian expansion is

$$
\boxed{v_{\mathrm{C}}(r_{12})\approx\sum_{g=1}^{N_g}W_g\exp\left[-\frac{t_g^2}{1-t_g^2}\left(\frac{r_{12}}{b_{\max}}\right)^2\right]}.
$$

Let $\phi_a(\boldsymbol{r})$, with $a\in\{1,2,3,4\}$, denote normalized axial harmonic-oscillator states. A superscript $^*$ denotes complex conjugation. The spatial Coulomb matrix element is

$$
V_{12,34}=\int\mathrm{d}^3r_1\,\mathrm{d}^3r_2\,\phi_1^*(\boldsymbol{r}_1)\phi_2^*(\boldsymbol{r}_2)v_{\mathrm{C}}(r_{12})\phi_3(\boldsymbol{r}_1)\phi_4(\boldsymbol{r}_2).
$$

For each $g$, let $G^z_{12,34}(g)$ and $G^r_{12,34}(g)$ denote the axial and transverse matrix elements of the corresponding Gaussian term; see [Axial Gaussian Matrix Elements](axial_gaussian_kernel.md). The Coulomb matrix element is

$$
\boxed{V_{12,34}\approx\sum_{g=1}^{N_g}W_gG^z_{12,34}(g)G^r_{12,34}(g)}.
$$

Parity and angular-projection selection rules follow from the canonical Gaussian matrix elements.

### Kramers-Sector Coulomb Elements

Coulomb contributes only to the proton field $\Gamma$; Coulomb pairing is omitted.

Let $\alpha,\beta,\gamma,\delta$ denote positive-$\Omega$ representatives and barred indices their [Kramers partners](../../notes/methods/hfb_axial.md#kramers-representative-formulation). Define $\sigma_a=2\Sigma_a\in\{-1,+1\}$, $\sigma_{\bar a}=-\sigma_a$, and the Kronecker delta $\delta_{xy}$.

$$
\bar v_{12;34}=\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}V_{12,34}-\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}V_{12,43}.
$$

The first sign labels the external $\Gamma$ sector; the second labels the contracted $\rho$ sector.

$$
\boxed{\bar v^{++}=\bar v_{\alpha\gamma;\beta\delta},\qquad \bar v^{+-}=\bar v_{\alpha\bar\gamma;\beta\bar\delta},\qquad \bar v^{-+}=\bar v_{\bar\alpha\gamma;\bar\beta\delta},\qquad \bar v^{--}=\bar v_{\bar\alpha\bar\gamma;\bar\beta\bar\delta}}.
$$

For the real Coulomb interaction, simultaneous time reversal gives

$$
\boxed{\bar v^{-+}=\bar v^{+-},\qquad \bar v^{--}=\bar v^{++}}.
$$

Density contractions and time-reversal phases follow [Axial HFB](../../notes/methods/hfb_axial.md#representative-gamma-and-delta-fields).
