# Kramers-Representative Formulation

This note applies time-reversal symmetry to the [Axial HFB Equation](hfb_axial.md#axial-hfb-equation).

### Time-Reversal Transformations

Time reversal leaves position unchanged and reverses momentum and angular momentum:

$$
\hat{\mathcal T}\hat{\mathbf r}\hat{\mathcal T}^{-1}=\hat{\mathbf r},\qquad \hat{\mathcal T}\hat{\mathbf p}\hat{\mathcal T}^{-1}=-\hat{\mathbf p},\qquad \hat{\mathcal T}\hat{\mathbf J}\hat{\mathcal T}^{-1}=-\hat{\mathbf J}.
$$

For a spin-$s$ fermion, the time-reversal operator is antiunitary and satisfies:

$$
\hat{\mathcal T}\hat S_z\hat{\mathcal T}^{-1}=-\hat S_z,\qquad \hat{\mathcal T}^2=(-1)^{2s}=-1.
$$

Choose real basis phases such that

$$
\hat{\mathcal T}|\phi_\alpha\rangle=\eta_\alpha|\phi_{\bar\alpha}\rangle,\qquad \hat{\mathcal T}|\phi_{\bar\alpha}\rangle=-\eta_\alpha|\phi_\alpha\rangle,\qquad \eta_\alpha=\pm1.
$$

The corresponding creation operators satisfy

$$
\hat{\mathcal T}\hat c_\alpha^*\hat{\mathcal T}^{-1}=\eta_\alpha\hat c_{\bar\alpha}^*,\qquad \hat{\mathcal T}\hat c_{\bar\alpha}^*\hat{\mathcal T}^{-1}=-\eta_\alpha\hat c_\alpha^*,
$$

with the same relations for the annihilation operators.

For a time-even one-body operator $\hat A$, define

$$
A_{\alpha\beta}\equiv\langle\phi_\alpha|\hat A|\phi_\beta\rangle.
$$

Antiunitarity gives

$$
\eta_\alpha\eta_\beta A_{\bar\alpha\bar\beta}=\langle\hat{\mathcal T}\phi_\alpha|\hat A|\hat{\mathcal T}\phi_\beta\rangle=\langle\phi_\alpha|\hat{\mathcal T}^{-1}\hat A\hat{\mathcal T}|\phi_\beta\rangle^*=A_{\alpha\beta}^*.
$$

For spin-$1/2$, a concrete representation is

$$
\hat{\mathcal T}=-i\sigma_yK,\qquad [\hat{\mathcal T}\phi_\alpha](\mathbf r)=-i\sigma_y\phi_\alpha^*(\mathbf r),
$$

and therefore

$$
\hat{\mathcal T}|\uparrow\rangle=|\downarrow\rangle,\qquad \hat{\mathcal T}|\downarrow\rangle=-|\uparrow\rangle.
$$

### Time-Reversal Reduction of the HFB Equation

For an unblocked even-even HFB vacuum, impose

$$
\hat{\mathcal T}|\Phi\rangle=|\Phi\rangle.
$$

The normal and pairing blocks satisfy

$$
h_{\bar\alpha\bar\beta}=\eta_\alpha\eta_\beta h_{\alpha\beta}^*,\qquad \Delta_{\alpha\bar\beta}=-\eta_\alpha\eta_\beta\Delta_{\bar\alpha\beta}^*.
$$

The blocks defined in [Solutions of the HFB Block Equations](hfb_axial.md#solutions-of-the-hfb-block-equations) become

$$
\mathcal H^+=\begin{pmatrix}h_{\alpha\beta}-\lambda\delta_{\alpha\beta}&\Delta_{\alpha\bar\beta}\\\eta_\alpha\eta_\beta\Delta_{\alpha\bar\beta}&-\eta_\alpha\eta_\beta h_{\alpha\beta}+\lambda\delta_{\alpha\beta}\end{pmatrix},\qquad \mathcal H^-=\begin{pmatrix}\eta_\alpha\eta_\beta h_{\alpha\beta}^*-\lambda\delta_{\alpha\beta}&-\eta_\alpha\eta_\beta\Delta_{\alpha\bar\beta}^*\\-\Delta_{\alpha\bar\beta}^*&-h_{\alpha\beta}^*+\lambda\delta_{\alpha\beta}\end{pmatrix}.
$$

Define the phase transformation

$$
\mathcal S_\eta\equiv\begin{pmatrix}\operatorname{diag}(\eta_\alpha)&0\\0&-\operatorname{diag}(\eta_\alpha)\end{pmatrix},\qquad \mathcal S_\eta^2=I.
$$

The two blocks satisfy

$$
\boxed{\mathcal H^-=\mathcal S_\eta\mathcal H^{+*}\mathcal S_\eta}.
$$

### Solutions of the Time-Reversal-Reduced HFB Equation

Using the eigenvectors defined in [Solutions of the HFB Block Equations](hfb_axial.md#solutions-of-the-hfb-block-equations), complex conjugation gives

$$
\mathcal H^{+*}X_\mu^{+*}=E_\mu^+X_\mu^{+*}.
$$

Multiplying by $\mathcal S_\eta$ yields

$$
\mathcal S_\eta\mathcal H^{+*}X_\mu^{+*}=E_\mu^+\mathcal S_\eta X_\mu^{+*}.
$$

Using $\mathcal S_\eta^2=I$,

$$
\mathcal S_\eta\mathcal H^{+*}\mathcal S_\eta\mathcal S_\eta X_\mu^{+*}=E_\mu^+\mathcal S_\eta X_\mu^{+*}.
$$

Since $\mathcal S_\eta\mathcal H^{+*}\mathcal S_\eta=\mathcal H^-$,

$$
\mathcal H^-\mathcal S_\eta X_\mu^{+*}=E_\mu^+\mathcal S_\eta X_\mu^{+*}.
$$

Choose the paired eigenvectors and their phases such that

$$
\boxed{X_{\bar\mu}^-=\mathcal S_\eta X_\mu^{+*},\qquad E_{\bar\mu}^-=E_\mu^+}.
$$

Collecting the paired positive-energy solutions gives

$$
\boxed{\begin{pmatrix}U^-\\V^+\end{pmatrix}=\begin{pmatrix}\operatorname{diag}(\eta_\alpha)&0\\0&-\operatorname{diag}(\eta_\alpha)\end{pmatrix}\begin{pmatrix}(U^+)^*\\(V^-)^*\end{pmatrix}}.
$$

Only $\mathcal H^+$ needs to be solved. With particle rows ordered as $(\alpha,\bar\alpha)$ and quasiparticle columns paired as $(\mu,\bar\mu)$, the full Bogoliubov matrices become

$$
\boxed{U=\begin{pmatrix}U^+&0\\0&\operatorname{diag}(\eta_\alpha)(U^+)^*\end{pmatrix},\qquad V=\begin{pmatrix}0&-\operatorname{diag}(\eta_\alpha)(V^-)^*\\V^-&0\end{pmatrix}}.
$$

The amplitudes are generally complex.

Using the density definitions in [Density Blocks](hfb_axial.md#density-blocks) and the Kramers amplitudes,

$$
\boxed{\rho=V^*V^T=\begin{pmatrix}\rho^+&0\\0&\operatorname{diag}(\eta_\alpha)(\rho^+)^*\operatorname{diag}(\eta_\alpha)\end{pmatrix}},
$$

$$
\boxed{\kappa=V^*U^T=\begin{pmatrix}0&\kappa^{+-}\\-\operatorname{diag}(\eta_\alpha)(\kappa^{+-})^*\operatorname{diag}(\eta_\alpha)&0\end{pmatrix}}.
$$

Combining time-reversal invariance with $\kappa^{-+}=-(\kappa^{+-})^T$ gives

$$
\kappa^{+-}_{\alpha\beta}=\eta_\alpha\eta_\beta\bigl[(\kappa^{+-})^\dagger\bigr]_{\alpha\beta}.
$$

### Representative $\Gamma$ and $\Delta$ Fields

Let $N_{\mathrm{sp}}$ denote the full single-particle dimension. Using the field definitions in [Hartree–Fock–Bogoliubov Method](hfb.md),

$$
\Gamma_{\alpha\beta}=\sum_{\gamma,\delta=1}^{N_{\mathrm{sp}}}\bar v_{\alpha\gamma;\beta\delta}\rho_{\delta\gamma},\qquad \alpha,\beta=1,\ldots,N_{\mathrm{sp}}.
$$

Substituting the density blocks from [Solutions of the Time-Reversal-Reduced HFB Equation](#solutions-of-the-time-reversal-reduced-hfb-equation) gives

$$
\boxed{\Gamma_{\alpha\beta}=\sum_{\gamma,\delta\in\{\alpha\mid\Omega_\alpha>0\}}\left[\bar v_{\alpha\gamma;\beta\delta}\rho_{\delta\gamma}^++\eta_\gamma\eta_\delta\bar v_{\alpha\bar\gamma;\beta\bar\delta}(\rho_{\delta\gamma}^+)^*\right]},\qquad \alpha,\beta\in\{\alpha\mid\Omega_\alpha>0\}.
$$

The remaining blocks satisfy

$$
\Gamma_{\alpha\bar\beta}=\Gamma_{\bar\alpha\beta}=0,\qquad \Gamma_{\bar\alpha\bar\beta}=\eta_\alpha\eta_\beta\Gamma_{\alpha\beta}^*,\qquad \alpha,\beta\in\{\alpha\mid\Omega_\alpha>0\}.
$$

The pairing field in the full single-particle space is

$$
\Delta_{\alpha\beta}=\frac12\sum_{\gamma,\delta=1}^{N_{\mathrm{sp}}}\bar v_{\alpha\beta;\gamma\delta}\kappa_{\gamma\delta},\qquad \alpha,\beta=1,\ldots,N_{\mathrm{sp}}.
$$

Its two cross-sector contributions give

$$
\Delta_{\alpha\bar\beta}=\frac12\sum_{\gamma,\delta\in\{\alpha\mid\Omega_\alpha>0\}}\left[\bar v_{\alpha\bar\beta;\gamma\bar\delta}\kappa_{\gamma\delta}^{+-}+\bar v_{\alpha\bar\beta;\bar\gamma\delta}\kappa_{\gamma\delta}^{-+}\right],\qquad \alpha,\beta\in\{\alpha\mid\Omega_\alpha>0\}.
$$

Using $\kappa^{-+}=-(\kappa^{+-})^T$, exchanging $\gamma,\delta$ in the second term, and applying the antisymmetry of $\bar v$ gives

$$
\boxed{\Delta_{\alpha\bar\beta}=\sum_{\gamma,\delta\in\{\alpha\mid\Omega_\alpha>0\}}\bar v_{\alpha\bar\beta;\gamma\bar\delta}\kappa_{\gamma\delta}^{+-}},\qquad \alpha,\beta\in\{\alpha\mid\Omega_\alpha>0\}.
$$

The remaining blocks satisfy

$$
\Delta_{\alpha\beta}=\Delta_{\bar\alpha\bar\beta}=0,\qquad \Delta_{\bar\alpha\beta}=-(\Delta^T)_{\bar\alpha\beta},\qquad \alpha,\beta\in\{\alpha\mid\Omega_\alpha>0\}.
$$
