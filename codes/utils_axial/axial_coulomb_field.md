# Axial Direct Coulomb Field

The direct Coulomb field is evaluated on the coordinate grid defined by the [Axial Harmonic-Oscillator Basis](axial_basis.md).

### Axial Coulomb Kernel

Let $\boldsymbol r_{\mathrm{s}}$ and $\boldsymbol r_{\mathrm{t}}$ denote the source and target positions. The direct Coulomb field is

$$
V_{\mathrm C}^{\mathrm{dir}}(\boldsymbol r_{\mathrm{t}})=e^2\int\frac{\rho_p(\boldsymbol r_{\mathrm{s}})}{|\boldsymbol r_{\mathrm{t}}-\boldsymbol r_{\mathrm{s}}|}\,\mathrm d^3r_{\mathrm{s}}.
$$

The Coulomb denominator has the Gaussian representation

$$
\frac{1}{|\boldsymbol r_{\mathrm{t}}-\boldsymbol r_{\mathrm{s}}|}=\frac{2}{\sqrt{\pi}}\int_0^\infty e^{-u^2|\boldsymbol r_{\mathrm{t}}-\boldsymbol r_{\mathrm{s}}|^2}\,\mathrm du.
$$

Under axial symmetry, take $\varphi_{\mathrm{t}}=0$ and write $s=(z_{\mathrm{s}},r_{\mathrm{s}})$ and $t=(z_{\mathrm{t}},r_{\mathrm{t}})$. Then

$$
V_{\mathrm C}^{\mathrm{dir}}(t)=\frac{2e^2}{\sqrt{\pi}}\int_{-\infty}^{\infty}\mathrm dz_{\mathrm{s}}\int_0^\infty r_{\mathrm{s}}\,\mathrm dr_{\mathrm{s}}\,\rho_p(s)\int_0^\infty\mathrm du\int_0^{2\pi}\mathrm d\varphi_{\mathrm{s}}\,e^{-u^2[(z_{\mathrm{s}}-z_{\mathrm{t}})^2+r_{\mathrm{s}}^2+r_{\mathrm{t}}^2-2r_{\mathrm{s}}r_{\mathrm{t}}\cos\varphi_{\mathrm{s}}]}.
$$

The azimuthally averaged part defines the bare axial kernel,

$$
\mathcal K_{\mathrm C}(s,t)=\frac{2e^2}{\sqrt{\pi}}\int_0^\infty\mathrm du\,\frac{1}{2\pi}\int_0^{2\pi}\mathrm d\varphi_{\mathrm{s}}\,e^{-u^2[(z_{\mathrm{s}}-z_{\mathrm{t}})^2+r_{\mathrm{s}}^2+r_{\mathrm{t}}^2-2r_{\mathrm{s}}r_{\mathrm{t}}\cos\varphi_{\mathrm{s}}]}.
$$

Define

$$
I_0^{\mathrm e}(y)=\frac{1}{2\pi}\int_0^{2\pi}e^{-y(1-\cos\varphi)}\,\mathrm d\varphi=e^{-y}I_0(y).
$$

After integrating over $\varphi_{\mathrm{s}}$,

$$
\boxed{\mathcal K_{\mathrm C}(s,t)=\frac{2e^2}{\sqrt{\pi}}\int_0^\infty e^{-u^2[(z_{\mathrm{s}}-z_{\mathrm{t}})^2+(r_{\mathrm{s}}-r_{\mathrm{t}})^2]}I_0^{\mathrm e}(2u^2r_{\mathrm{s}}r_{\mathrm{t}})\,\mathrm du.}
$$

### Semi-Infinite Mapping

The semi-infinite interval is mapped to $x\in[0,1)$ by

$$
u(x)=\frac{x}{b\sqrt{1-x^2}},\qquad \frac{\mathrm du}{\mathrm dx}=\frac{1}{b(1-x^2)^{3/2}}.
$$

Therefore, the bare kernel becomes

$$
\boxed{\mathcal K_{\mathrm C}(s,t)=\frac{2e^2}{\sqrt{\pi}}\int_0^1\frac{\mathrm dx}{b(1-x^2)^{3/2}}e^{-u(x)^2[(z_{\mathrm{s}}-z_{\mathrm{t}})^2+(r_{\mathrm{s}}-r_{\mathrm{t}})^2]}I_0^{\mathrm e}(2u(x)^2r_{\mathrm{s}}r_{\mathrm{t}}).}
$$

The scale $b=50\,\mathrm{fm}$ controls the mapping without changing the integral. The finite integral is evaluated by Gauss–Legendre quadrature.

### Reflection Symmetry

For a reflection-reduced grid,

$$
\rho_p(-z,r)=\rho_p(z,r).
$$

The axial exponential is replaced by

$$
e^{-u^2(z_{\mathrm{s}}-z_{\mathrm{t}})^2}\longrightarrow\frac12\left[e^{-u^2(z_{\mathrm{s}}-z_{\mathrm{t}})^2}+e^{-u^2(z_{\mathrm{s}}+z_{\mathrm{t}})^2}\right].
$$

The factor $1/2$ compensates for the doubled positive-$z$ quadrature weights.
