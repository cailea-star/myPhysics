# BCS Reduction

Starting from the real [Kramers-representative formulation](hfb_kramers.md), use the subscript $\mathrm{nat}$ for matrices in the natural basis.

### Natural-Basis Transformation

The real symmetric density matrix admits an orthogonal diagonalization:

$$
\boxed{\rho^{++}=S\rho^{++}_{\mathrm{nat}}S^T,\qquad \rho^{++}_{\mathrm{nat}}=S^T\rho^{++}S,\qquad S^TS=I.}
$$

To use the same transformation for both Kramers sectors while preserving their phases, assume

$$
S\operatorname{diag}(\eta_\alpha)=\operatorname{diag}(\eta_\alpha)S.
$$

### BCS Diagonalization

The BCS approximation assumes that the pairing interaction acts predominantly between time-reversed partners. In the natural basis, retain only pairing-field matrix elements connecting each orbital to its Kramers partner, and additionally assume that the normal field is diagonal. The amplitudes then admit a diagonal representation with a compatible quasiparticle transformation.

For an orthogonal quasiparticle transformation $Q$,

$$
\boxed{U^+_{\mathrm{nat}}=S^TU^+Q=\operatorname{diag}(u_\alpha),\qquad V^-_{\mathrm{nat}}=S^TV^-Q=\operatorname{diag}(v_\alpha),\qquad Q^TQ=I.}
$$

The normal and pairing fields transform as

$$
h^{++}_{\mathrm{nat}}=S^Th^{++}S,\qquad h^{--}_{\mathrm{nat}}=S^Th^{--}S,
$$

$$
\Delta^{+-}_{\mathrm{nat}}=S^T\Delta^{+-}S,\qquad \Delta^{-+}_{\mathrm{nat}}=S^T\Delta^{-+}S.
$$

The transformed fields are assumed diagonal:

$$
(h^{++}_{\mathrm{nat}})_{\alpha\beta}=(h^{++}_{\mathrm{nat}})_{\alpha\alpha}\delta_{\alpha\beta},\qquad (h^{--}_{\mathrm{nat}})_{\alpha\beta}=(h^{--}_{\mathrm{nat}})_{\alpha\alpha}\delta_{\alpha\beta},
$$

$$
(\Delta^{+-}_{\mathrm{nat}})_{\alpha\beta}=(\Delta^{+-}_{\mathrm{nat}})_{\alpha\alpha}\delta_{\alpha\beta},\qquad (\Delta^{-+}_{\mathrm{nat}})_{\alpha\beta}=(\Delta^{-+}_{\mathrm{nat}})_{\alpha\alpha}\delta_{\alpha\beta}.
$$

The representative HFB matrix therefore transforms as

$$
\boxed{\begin{pmatrix}S^T&0\\0&S^T\end{pmatrix}\mathcal H^+\begin{pmatrix}S&0\\0&S\end{pmatrix}=\begin{pmatrix}h^{++}_{\mathrm{nat}}-\lambda I&\Delta^{+-}_{\mathrm{nat}}\\-\Delta^{-+}_{\mathrm{nat}}&-h^{--}_{\mathrm{nat}}+\lambda I\end{pmatrix}.}
$$

All four blocks on the right are diagonal.

Each $\alpha$ defines an independent HFB block:

$$
\boxed{\mathcal H^+_\alpha(\mathrm{nat})=\begin{pmatrix}(h^{++}_{\mathrm{nat}})_{\alpha\alpha}-\lambda&(\Delta^{+-}_{\mathrm{nat}})_{\alpha\alpha}\\-(\Delta^{-+}_{\mathrm{nat}})_{\alpha\alpha}&-(h^{--}_{\mathrm{nat}})_{\alpha\alpha}+\lambda\end{pmatrix}.}
$$
