# BCS Reduction

This note specializes the [Axial HFB Equation](hfb_axial.md#axial-hfb-equation) through simultaneous diagonalization and subsequent time-reversal reduction. All matrices are real.

### HFB Equation

For each $(\Omega,\pi)$ block, let $\alpha,\beta$ label positive-$\Omega$ states and $\bar\alpha,\bar\beta$ their opposite-$\Omega$ partners. This labeling does not assume time-reversal symmetry. The supplied single-particle field is $h=h_0+\Gamma$, and $\lambda$ is the chemical potential.

$$
\sum_\beta\begin{pmatrix}h_{\alpha\beta}-\lambda\delta_{\alpha\beta}&0&0&\Delta_{\alpha\bar\beta}\\0&h_{\bar\alpha\bar\beta}-\lambda\delta_{\alpha\beta}&\Delta_{\bar\alpha\beta}&0\\0&-\Delta_{\alpha\bar\beta}&-h_{\alpha\beta}+\lambda\delta_{\alpha\beta}&0\\-\Delta_{\bar\alpha\beta}&0&0&-h_{\bar\alpha\bar\beta}+\lambda\delta_{\alpha\beta}\end{pmatrix}\begin{pmatrix}U_{\beta\mu}^+\\U_{\beta\mu}^-\\V_{\beta\mu}^+\\V_{\beta\mu}^-\end{pmatrix}=E_\mu\begin{pmatrix}U_{\alpha\mu}^+\\U_{\alpha\mu}^-\\V_{\alpha\mu}^+\\V_{\alpha\mu}^-\end{pmatrix}.
$$

The $(U^+,V^-)$ and $(U^-,V^+)$ sectors decouple.

### Simultaneous Diagonalization

Assume that real orthonormal transformations $f_{\alpha i}^+$ and $f_{\alpha i}^-$ simultaneously diagonalize the normal fields and the pairing field. Define the single-particle energies $\epsilon_i^+,\epsilon_i^-$ and pairing matrix elements $\Delta_{i\bar i}$ by

$$
\sum_{\alpha\beta}f_{\alpha i}^+h_{\alpha\beta}f_{\beta j}^+=\epsilon_i^+\delta_{ij},\qquad \sum_{\alpha\beta}f_{\alpha i}^-h_{\bar\alpha\bar\beta}f_{\beta j}^-=\epsilon_i^-\delta_{ij},\qquad \sum_{\alpha\beta}f_{\alpha i}^+\Delta_{\alpha\bar\beta}f_{\beta j}^-=\Delta_{i\bar i}\delta_{ij}.
$$

The transformations satisfy

$$
\sum_\alpha f_{\alpha i}^+f_{\alpha j}^+=\delta_{ij},\qquad \sum_\alpha f_{\alpha i}^-f_{\alpha j}^-=\delta_{ij}.
$$

Pairing antisymmetry gives the remaining block:

$$
\sum_{\alpha\beta}f_{\alpha i}^-\Delta_{\bar\alpha\beta}f_{\beta j}^+=\Delta_{\bar i i}\delta_{ij}=-\Delta_{i\bar i}\delta_{ij}.
$$

This simultaneous diagonalization is the BCS assumption; diagonalizing $h$ alone does not ensure diagonal pairing.

For an eigenvector supported on the pair $(i,\bar i)$, define the amplitudes by

$$
U_{\alpha i}^+=f_{\alpha i}^+u_i^+,\qquad U_{\alpha i}^-=f_{\alpha i}^-u_i^-,
$$

$$
V_{\alpha i}^+=f_{\alpha i}^+v_i^+,\qquad V_{\alpha i}^-=f_{\alpha i}^-v_i^-.
$$

Each pair $(i,\bar i)$ then defines an independent BCS equation:

$$
\boxed{\begin{pmatrix}\epsilon_i^+-\lambda&0&0&\Delta_{i\bar i}\\0&\epsilon_i^--\lambda&\Delta_{\bar i i}&0\\0&-\Delta_{i\bar i}&-(\epsilon_i^+-\lambda)&0\\-\Delta_{\bar i i}&0&0&-(\epsilon_i^--\lambda)\end{pmatrix}\begin{pmatrix}u_i^+\\u_i^-\\v_i^+\\v_i^-\end{pmatrix}=E_i\begin{pmatrix}u_i^+\\u_i^-\\v_i^+\\v_i^-\end{pmatrix}.}
$$

### Time-Reversal Reduction

Choose the single-particle eigenstates as time-reversal partners with real phases $\eta_i=\pm1$:

$$
\hat{\mathcal T}|i\rangle=\eta_i|\bar i\rangle,\qquad \hat{\mathcal T}|\bar i\rangle=-\eta_i|i\rangle.
$$

Time-reversal invariance and pairing antisymmetry give, respectively,

$$
\epsilon_i^-=\epsilon_i^+,\qquad \Delta_{\bar i i}=-\Delta_{i\bar i}.
$$

The two BCS sectors are equivalent. Choose the time-reversed partner amplitudes as

$$
u_i^-=\eta_i u_i^+,\qquad v_i^+=-\eta_i v_i^-.
$$

Retain the $(u_i^+,v_i^-)$ representative equation:

$$
\boxed{\begin{pmatrix}\epsilon_i^+-\lambda&\Delta_{i\bar i}\\\Delta_{i\bar i}&-(\epsilon_i^+-\lambda)\end{pmatrix}\begin{pmatrix}u_i^+\\v_i^-\end{pmatrix}=E_i\begin{pmatrix}u_i^+\\v_i^-\end{pmatrix}.}
$$

For $E_i>0$, the normalized positive-energy solution satisfies

$$
E_i=\sqrt{(\epsilon_i^+-\lambda)^2+(\Delta_{i\bar i})^2}.
$$

$$
(u_i^+)^2=\frac12\left(1+\frac{\epsilon_i^+-\lambda}{E_i}\right),\qquad (v_i^-)^2=\frac12\left(1-\frac{\epsilon_i^+-\lambda}{E_i}\right).
$$

Consequently,

$$
u_i^+v_i^-=\frac{\Delta_{i\bar i}}{2E_i},\qquad\epsilon_i^+=\lambda+E_i\Big[(u_i^+)^2-(v_i^-)^2\Big].
$$

At zero temperature, $\rho=V^*V^T$ and $\kappa=V^*U^T$ give the following matrices for each pair $(i,\bar i)$:

$$
\rho_i=\begin{pmatrix}0&v_i^+\\v_i^-&0\end{pmatrix}\begin{pmatrix}0&v_i^-\\v_i^+&0\end{pmatrix}=\begin{pmatrix}(v_i^+)^2&0\\0&(v_i^-)^2\end{pmatrix}=\begin{pmatrix}\rho_i^+&0\\0&\rho_i^-\end{pmatrix}.
$$

The pairing tensor is

$$
\kappa_i=\begin{pmatrix}0&v_i^+\\v_i^-&0\end{pmatrix}\begin{pmatrix}u_i^+&0\\0&u_i^-\end{pmatrix}=\begin{pmatrix}0&v_i^+u_i^-\\v_i^-u_i^+&0\end{pmatrix}=\begin{pmatrix}0&\kappa_i^{+-}\\\kappa_i^{-+}&0\end{pmatrix}.
$$

Time-reversal symmetry implies $\rho_i^+=\rho_i^-$ and $\kappa_i^{+-}=-\kappa_i^{-+}=-u_i^+v_i^-$. The total occupation of each pair is

$$
\rho_i^++\rho_i^-=2(v_i^-)^2.
$$

### Constant Pairing Strength

Let $v$ denote antisymmetrized two-body matrix elements. Pairing antisymmetry gives

$$
\Delta_{i\bar i}=\frac12\sum_j\left(v_{i\bar i j\bar j}\kappa_j^{+-}+v_{i\bar i\bar j j}\kappa_j^{-+}\right)=\sum_jv_{i\bar i j\bar j}\kappa_j^{+-}.
$$

Within a pairing space $\mathcal W$ near the Fermi surface, approximate the pair-scattering interaction by a constant attractive strength $G>0$:

$$
v_{i\bar i j\bar j}\approx-G\eta_i\eta_j,\qquad \Delta_{i\bar i}=\eta_i\Delta,\qquad \Delta=-G\sum_{j\in\mathcal W}\eta_j\kappa_j^{+-}.
$$

Using $\kappa_j^{+-}=-u_j^+v_j^-$, the self-consistency equations are

$$
\boxed{\Delta=G\sum_{j\in\mathcal W}\eta_j u_j^+v_j^-,\qquad N=\sum_i(\rho_i^++\rho_i^-)=2\sum_i(v_i^-)^2.}
$$

For $\Delta\ne0$, the gap equation is equivalently

$$
1=G\sum_{j\in\mathcal W}\frac{1}{2E_j},\qquad E_j=\sqrt{(\epsilon_j^+-\lambda)^2+\Delta^2}.
$$

Each sum counts one representative per pair. The number sum covers the full model space; pairing vanishes outside $\mathcal W$. Determine $\lambda$ and $\Delta$ from the two self-consistency equations, allowing $\Delta=0$ when pairing collapses.
