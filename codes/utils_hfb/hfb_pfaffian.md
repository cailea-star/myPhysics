# HFB Pfaffian Algorithm

### Quasiparticle Vacuum Representation

Let $N_{\mathrm{sp}}$ denote the number of single-particle states and $|0\rangle$ the particle vacuum. Use the Bogoliubov convention in [Quasiparticle Vacuum](hfb_full.md#quasiparticle-vacuum), with operator vectors arranged as columns:

$$
\hat{\boldsymbol\beta}=U^\dagger\hat{\boldsymbol c}+V^\dagger\hat{\boldsymbol c}^*,\qquad \hat{\boldsymbol c}=U\hat{\boldsymbol\beta}+V^*\hat{\boldsymbol\beta}^*.
$$

The normalized quasiparticle vacuum satisfies

$$
\hat{\boldsymbol\beta}|\Phi\rangle=0,\qquad \langle\Phi|\Phi\rangle=1.
$$

For invertible $U$, its Thouless representation is

$$
\boxed{|\Phi\rangle=\nu\exp\left[\frac12\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}Z_{\alpha\beta}\hat c_\alpha^*\hat c_\beta^*\right]|0\rangle}.
$$

Here,

$$
\boxed{\nu\equiv\langle0|\Phi\rangle,\qquad Z\equiv(VU^{-1})^*=V^*U^{-*},\qquad Z^T=-Z}.
$$

Normalization fixes the magnitude of $\nu$, while its phase specifies the overall vacuum phase:

$$
|\nu|^2=|\det U|.
$$

Define the pair-creation operator by

$$
\hat A\equiv\frac12\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}Z_{\alpha\beta}\hat c_\alpha^*\hat c_\beta^*.
$$

The Thouless vacuum expands as

$$
|\Phi\rangle=\nu\exp(\hat A)|0\rangle=\nu\sum_{p=0}^{\lfloor N_{\mathrm{sp}}/2\rfloor}\frac{\hat A^p}{p!}|0\rangle.
$$

Let $|\Phi^{(n)}\rangle$ denote its $n$-particle component. Since $\hat A$ creates two particles,

$$
|\Phi^{(2p+1)}\rangle=0.
$$

The zero-particle component is

$$
|\Phi^{(0)}\rangle=\nu|0\rangle.
$$

The two-particle component is

$$
|\Phi^{(2)}\rangle=\nu\sum_{\alpha<\beta}Z_{\alpha\beta}\hat c_\alpha^*\hat c_\beta^*|0\rangle.
$$

The four-particle component is

$$
|\Phi^{(4)}\rangle=\nu\sum_{\alpha<\beta<\gamma<\delta}\left(Z_{\alpha\beta}Z_{\gamma\delta}-Z_{\alpha\gamma}Z_{\beta\delta}+Z_{\alpha\delta}Z_{\beta\gamma}\right)\hat c_\alpha^*\hat c_\beta^*\hat c_\gamma^*\hat c_\delta^*|0\rangle.
$$

The three terms represent all pairings of four orbitals, with signs fixed by fermionic permutations.

For a positive integer $p$, let $X\in\mathbb C^{2p\times2p}$ be skew-symmetric and $S_{2p}$ the permutation group of $2p$ elements. Define the Pfaffian by

$$
\boxed{\operatorname{pf}(X)\equiv\frac{1}{2^p p!}\sum_{\sigma\in S_{2p}}\operatorname{sgn}(\sigma)\prod_{r=1}^{p}X_{\sigma(2r-1),\sigma(2r)}}.
$$

With $X_\varnothing$ denoting the empty matrix, the lowest-order Pfaffians are

$$
\operatorname{pf}(X_\varnothing)=1,
$$

$$
\operatorname{pf}\begin{pmatrix}0&X_{12}\\-X_{12}&0\end{pmatrix}=X_{12},
$$

$$
\operatorname{pf}\begin{pmatrix}0&X_{12}&X_{13}&X_{14}\\-X_{12}&0&X_{23}&X_{24}\\-X_{13}&-X_{23}&0&X_{34}\\-X_{14}&-X_{24}&-X_{34}&0\end{pmatrix}=X_{12}X_{34}-X_{13}X_{24}+X_{14}X_{23}.
$$

The two- and four-particle coefficients are therefore Pfaffians of the corresponding principal submatrices of $Z$.

The general $2p$-particle component is

$$
|\Phi^{(2p)}\rangle=\nu\sum_{1\leq\alpha_1<\cdots<\alpha_{2p}\leq N_{\mathrm{sp}}}\operatorname{pf}\left[\left(Z_{\alpha_r\alpha_s}\right)_{r,s=1}^{2p}\right]\hat c_{\alpha_1}^*\cdots\hat c_{\alpha_{2p}}^*|0\rangle.
$$

Therefore,

$$
\boxed{|\Phi\rangle=\nu\sum_{p=0}^{\lfloor N_{\mathrm{sp}}/2\rfloor}\sum_{1\leq\alpha_1<\cdots<\alpha_{2p}\leq N_{\mathrm{sp}}}\operatorname{pf}\left[\left(Z_{\alpha_r\alpha_s}\right)_{r,s=1}^{2p}\right]\hat c_{\alpha_1}^*\cdots\hat c_{\alpha_{2p}}^*|0\rangle}.
$$

### Quasiparticle Vacuum Overlap

Let $|\Phi_1\rangle$ and $|\Phi_2\rangle$ be normalized Thouless vacua in the same ordered single-particle basis, with parameters $(U_a,V_a,Z_a,\nu_a)$ for $a=1,2$ and invertible $U_1,U_2$.

Orthogonality between distinct occupation patterns gives

$$
\langle\Phi_1|\Phi_2\rangle=\nu_1^*\nu_2\sum_{p=0}^{\lfloor N_{\mathrm{sp}}/2\rfloor}\sum_{1\leq\alpha_1<\cdots<\alpha_{2p}\leq N_{\mathrm{sp}}}\operatorname{pf}\left[\left((Z_1)^*_{\alpha_r\alpha_s}\right)_{r,s=1}^{2p}\right]\operatorname{pf}\left[\left((Z_2)_{\alpha_r\alpha_s}\right)_{r,s=1}^{2p}\right].
$$

Let $I_{N_{\mathrm{sp}}}$ denote the identity matrix. Define

$$
\mathbb X\equiv\begin{pmatrix}Z_2&-I_{N_{\mathrm{sp}}}\\I_{N_{\mathrm{sp}}}&-Z_1^*\end{pmatrix},\qquad s_{N_{\mathrm{sp}}}\equiv(-1)^{N_{\mathrm{sp}}(N_{\mathrm{sp}}+1)/2}.
$$

The Pfaffian minor-summation identity yields

$$
\boxed{\langle\Phi_1|\Phi_2\rangle=\nu_1^*\nu_2s_{N_{\mathrm{sp}}}\operatorname{pf}(\mathbb X)}.
$$

Normalization gives

$$
\boxed{|\nu_1|^2=|\det U_1|,\qquad |\nu_2|^2=|\det U_2|}.
$$

Since $\operatorname{pf}(X)^2=\det X$,

$$
\operatorname{pf}(\mathbb X)^2=\det\mathbb X=\det(I_{N_{\mathrm{sp}}}-Z_1^*Z_2).
$$

Define the quasiparticle overlap matrix by

$$
\boxed{\mathcal A\equiv U_1^TU_2^*+V_1^TV_2^*=U_1^T(I_{N_{\mathrm{sp}}}-Z_1^*Z_2)U_2^*}.
$$

Taking determinants gives

$$
\det\mathcal A=\det U_1\,\det U_2^*\,\det(I_{N_{\mathrm{sp}}}-Z_1^*Z_2).
$$

Combining this with the squared Pfaffian overlap yields

$$
\boxed{\langle\Phi_1|\Phi_2\rangle^2=\frac{(\nu_1^*\nu_2)^2}{\det U_1\,\det U_2^*}\det\mathcal A}.
$$

Normalization fixes the magnitudes of $\nu_1,\nu_2$, but not the relative vacuum phase. The square root leaves a sign ambiguity; the Pfaffian fixes the overlap for a specified single-particle ordering and vacuum-phase convention.

These overlap formulas remain valid when $\langle\Phi_1|\Phi_2\rangle=0$. Nonzero overlap is required only for normalized transition contractions and the inverse $\mathcal A^{-1}$.

### Transition Densities

##### One-Body Transition Densities

Assume $\langle\Phi_1|\Phi_2\rangle\ne0$, with all operator vectors arranged as columns. Define the normal transition density by

$$
\rho^T\equiv\frac{\langle\Phi_1|\hat{\boldsymbol c}^*\hat{\boldsymbol c}^T|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}.
$$

The Bogoliubov transformations are

$$
\hat{\boldsymbol c}=U_a\hat{\boldsymbol\beta}_a+V_a^*\hat{\boldsymbol\beta}_a^*,\qquad \hat{\boldsymbol c}^*=V_a\hat{\boldsymbol\beta}_a+U_a^*\hat{\boldsymbol\beta}_a^*,\qquad a=1,2.
$$

Using the left transformation for $\hat{\boldsymbol c}^*$, the right transformation for $\hat{\boldsymbol c}$, and the vacuum conditions

$$
\langle\Phi_1|\hat{\boldsymbol\beta}_1^*=0,\qquad \hat{\boldsymbol\beta}_2|\Phi_2\rangle=0,
$$

gives

$$
\rho^T=V_1\frac{\langle\Phi_1|\hat{\boldsymbol\beta}_1\hat{\boldsymbol\beta}_2^\dagger|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}V_2^\dagger.
$$

The two quasiparticle bases satisfy

$$
\hat{\boldsymbol\beta}_1^*=\mathcal B\hat{\boldsymbol\beta}_2+\mathcal A\hat{\boldsymbol\beta}_2^*,
$$

where

$$
\mathcal A=U_1^TU_2^*+V_1^TV_2^*,\qquad \mathcal B\equiv V_1^TU_2+U_1^TV_2.
$$

The canonical anticommutation relations, quasiparticle transformation, and vacuum conditions give

$$
\langle\Phi_1|\Phi_2\rangle I_{N_{\mathrm{sp}}}=\langle\Phi_1|\hat{\boldsymbol\beta}_1\hat{\boldsymbol\beta}_1^\dagger|\Phi_2\rangle=\langle\Phi_1|\hat{\boldsymbol\beta}_1\hat{\boldsymbol\beta}_2^\dagger|\Phi_2\rangle\mathcal A^T.
$$

Therefore,

$$
\boxed{\rho=V_2^*\mathcal A^{-1}V_1^T}.
$$

##### Pairing Transition Densities

Assume $\langle\Phi_1|\Phi_2\rangle\ne0$. Define the pairing transition densities by

$$
\kappa^T\equiv\frac{\langle\Phi_1|\hat{\boldsymbol c}\hat{\boldsymbol c}^T|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle},\qquad \bar\kappa^T\equiv\frac{\langle\Phi_1|\hat{\boldsymbol c}^*\hat{\boldsymbol c}^\dagger|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}.
$$

Using the left Bogoliubov transformation for the first operator vector and the right transformation for the second, the vacuum conditions give

$$
\kappa^T=U_1\frac{\langle\Phi_1|\hat{\boldsymbol\beta}_1\hat{\boldsymbol\beta}_2^\dagger|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}V_2^\dagger,\qquad \bar\kappa^T=V_1\frac{\langle\Phi_1|\hat{\boldsymbol\beta}_1\hat{\boldsymbol\beta}_2^\dagger|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}U_2^\dagger.
$$

The mixed quasiparticle matrix element satisfies the relation derived in [One-Body Transition Densities](#one-body-transition-densities):

$$
\langle\Phi_1|\Phi_2\rangle I_{N_{\mathrm{sp}}}=\langle\Phi_1|\hat{\boldsymbol\beta}_1\hat{\boldsymbol\beta}_2^\dagger|\Phi_2\rangle\mathcal A^T.
$$

Therefore,

$$
\boxed{\kappa=V_2^*\mathcal A^{-1}U_1^T,\qquad \bar\kappa=U_2^*\mathcal A^{-1}V_1^T}.
$$

Fermionic anticommutation implies

$$
\kappa^T=-\kappa,\qquad \bar\kappa^T=-\bar\kappa.
$$

For identical vacua with identical Bogoliubov matrices, $\mathcal A=I_{N_{\mathrm{sp}}}$, recovering

$$
\kappa=V^*U^T,\qquad \bar\kappa=U^*V^T=-\kappa^*.
$$

The relation $\bar\kappa=-\kappa^*$ does not generally hold for distinct vacua.

##### Other Contractions

Assume $\langle\Phi_1|\Phi_2\rangle\ne0$ and define

$$
\langle\hat X\rangle_{12}\equiv\frac{\langle\Phi_1|\hat X|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}.
$$

- Quasiparticle operators. Creation operators form the column $\hat{\boldsymbol\beta}_a^*$ and the row $\hat{\boldsymbol\beta}_a^\dagger=(\hat{\boldsymbol\beta}_a^*)^T$.

For $|\Phi_1\rangle$,

$$
\hat{\boldsymbol\beta}_1=U_1^\dagger\hat{\boldsymbol c}+V_1^\dagger\hat{\boldsymbol c}^*,\qquad \hat{\boldsymbol\beta}_1^\dagger=\hat{\boldsymbol c}^\dagger U_1+\hat{\boldsymbol c}^TV_1.
$$

For $|\Phi_2\rangle$,

$$
\hat{\boldsymbol\beta}_2=U_2^\dagger\hat{\boldsymbol c}+V_2^\dagger\hat{\boldsymbol c}^*,\qquad \hat{\boldsymbol\beta}_2^\dagger=\hat{\boldsymbol c}^\dagger U_2+\hat{\boldsymbol c}^TV_2.
$$

- All contractions follow by substitution using $\mathcal A^{-T}\equiv(\mathcal A^{-1})^T$:

$$
\langle\hat{\boldsymbol c}^*\hat{\boldsymbol c}^T\rangle_{12}=\rho^T=V_1\mathcal A^{-T}V_2^\dagger,\qquad \langle\hat{\boldsymbol c}\hat{\boldsymbol c}^\dagger\rangle_{12}=I_{N_{\mathrm{sp}}}-\rho=U_1\mathcal A^{-T}U_2^\dagger,
$$

$$
\langle\hat{\boldsymbol c}\hat{\boldsymbol c}^T\rangle_{12}=\kappa^T=U_1\mathcal A^{-T}V_2^\dagger,\qquad \langle\hat{\boldsymbol c}^*\hat{\boldsymbol c}^\dagger\rangle_{12}=\bar\kappa^T=V_1\mathcal A^{-T}U_2^\dagger.
$$

- Quasiparticle–quasiparticle contractions.

For $|\Phi_1\rangle$,

$$
\langle\hat{\boldsymbol\beta}_1^*\hat{\boldsymbol\beta}_1^T\rangle_{12}=0,\qquad \langle\hat{\boldsymbol\beta}_1\hat{\boldsymbol\beta}_1^\dagger\rangle_{12}=I_{N_{\mathrm{sp}}},\qquad \langle\hat{\boldsymbol\beta}_1^*\hat{\boldsymbol\beta}_1^\dagger\rangle_{12}=0.
$$

$$
\langle\hat{\boldsymbol\beta}_1\hat{\boldsymbol\beta}_1^T\rangle_{12}=U_1^\dagger\kappa^TU_1^*+U_1^\dagger(I_{N_{\mathrm{sp}}}-\rho)V_1^*+V_1^\dagger\rho^TU_1^*+V_1^\dagger\bar\kappa^TV_1^*=\mathcal A^{-T}(V_2^\dagger U_1^*+U_2^\dagger V_1^*).
$$

For $|\Phi_2\rangle$,

$$
\langle\hat{\boldsymbol\beta}_2^*\hat{\boldsymbol\beta}_2^T\rangle_{12}=0,\qquad \langle\hat{\boldsymbol\beta}_2\hat{\boldsymbol\beta}_2^\dagger\rangle_{12}=I_{N_{\mathrm{sp}}},\qquad \langle\hat{\boldsymbol\beta}_2\hat{\boldsymbol\beta}_2^T\rangle_{12}=0.
$$

$$
\langle\hat{\boldsymbol\beta}_2^*\hat{\boldsymbol\beta}_2^\dagger\rangle_{12}=V_2^T\kappa^TV_2+V_2^T(I_{N_{\mathrm{sp}}}-\rho)U_2+U_2^T\rho^TV_2+U_2^T\bar\kappa^TU_2=(V_2^TU_1+U_2^TV_1)\mathcal A^{-T}.
$$

Between $|\Phi_1\rangle$ and $|\Phi_2\rangle$,

$$
\langle\hat{\boldsymbol\beta}_1\hat{\boldsymbol\beta}_2^\dagger\rangle_{12}=U_1^\dagger\kappa^TV_2+U_1^\dagger(I_{N_{\mathrm{sp}}}-\rho)U_2+V_1^\dagger\rho^TV_2+V_1^\dagger\bar\kappa^TU_2=\mathcal A^{-T}.
$$

$$
\langle\hat{\boldsymbol\beta}_1^*\hat{\boldsymbol\beta}_2^T\rangle_{12}=\langle\hat{\boldsymbol\beta}_1\hat{\boldsymbol\beta}_2^T\rangle_{12}=\langle\hat{\boldsymbol\beta}_1^*\hat{\boldsymbol\beta}_2^\dagger\rangle_{12}=0.
$$

- Quasiparticle–particle contractions.

For $|\Phi_1\rangle$, quasiparticle operators are placed on the left:

$$
\langle\hat{\boldsymbol\beta}_1^*\hat{\boldsymbol c}^T\rangle_{12}=0,\qquad \langle\hat{\boldsymbol\beta}_1^*\hat{\boldsymbol c}^\dagger\rangle_{12}=0,
$$

$$
\langle\hat{\boldsymbol\beta}_1\hat{\boldsymbol c}^T\rangle_{12}=U_1^\dagger\kappa^T+V_1^\dagger\rho^T=\mathcal A^{-T}V_2^\dagger,\qquad \langle\hat{\boldsymbol\beta}_1\hat{\boldsymbol c}^\dagger\rangle_{12}=U_1^\dagger(I_{N_{\mathrm{sp}}}-\rho)+V_1^\dagger\bar\kappa^T=\mathcal A^{-T}U_2^\dagger.
$$

For $|\Phi_2\rangle$, quasiparticle operators are placed on the right:

$$
\langle\hat{\boldsymbol c}\hat{\boldsymbol\beta}_2^\dagger\rangle_{12}=\kappa^TV_2+(I_{N_{\mathrm{sp}}}-\rho)U_2=U_1\mathcal A^{-T},\qquad \langle\hat{\boldsymbol c}^*\hat{\boldsymbol\beta}_2^\dagger\rangle_{12}=\rho^TV_2+\bar\kappa^TU_2=V_1\mathcal A^{-T},
$$

$$
\langle\hat{\boldsymbol c}\hat{\boldsymbol\beta}_2^T\rangle_{12}=0,\qquad \langle\hat{\boldsymbol c}^*\hat{\boldsymbol\beta}_2^T\rangle_{12}=0.
$$

### Generalized Wick Theorem

Assume $\langle\Phi_1|\Phi_2\rangle\ne0$. Let $\hat s_1,\ldots,\hat s_{2p}$ be linear combinations of particle creation and annihilation operators. Define

$$
\langle\hat s_1\cdots\hat s_{2p}\rangle_{12}\equiv\frac{\langle\Phi_1|\hat s_1\cdots\hat s_{2p}|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}.
$$

For the specified operator ordering, define the contraction matrix $S\in\mathbb C^{2p\times2p}$ by

$$
S_{ij}\equiv\langle\hat s_i\hat s_j\rangle_{12}\quad(i<j),\qquad S=-S^T.
$$

The generalized Wick theorem gives

$$
\boxed{\langle\hat s_1\cdots\hat s_{2p}\rangle_{12}=\operatorname{pf}(S)=\operatorname{pf}\begin{pmatrix}0&S_{12}&\cdots&S_{1,2p}\\-S_{12}&0&\cdots&S_{2,2p}\\\vdots&\vdots&\ddots&\vdots\\-S_{1,2p}&-S_{2,2p}&\cdots&0\end{pmatrix}}.
$$

Both Thouless vacua have even particle-number parity, so

$$
\langle\hat s_1\cdots\hat s_{2p+1}\rangle_{12}=0.
$$

Let the ordered insertion group be

$$
\hat X=(\hat x_1,\ldots,\hat x_{N_X}).
$$

For configurations

$$
\kappa_1=(\mu_1,\ldots,\mu_{r_1}),\qquad \kappa_2=(\nu_1,\ldots,\nu_{r_2}),
$$

the excited states are

$$
|\Phi_{1;\kappa_1}\rangle=\hat\beta_{1,\mu_1}^*\cdots\hat\beta_{1,\mu_{r_1}}^*|\Phi_1\rangle,\qquad |\Phi_{2;\kappa_2}\rangle=\hat\beta_{2,\nu_1}^*\cdots\hat\beta_{2,\nu_{r_2}}^*|\Phi_2\rangle.
$$

The full ordered operator column is

$$
\hat{\boldsymbol s}=\begin{pmatrix}\hat\beta_{1,\mu_{r_1}}&\cdots&\hat\beta_{1,\mu_1}&\hat x_1&\cdots&\hat x_{N_X}&\hat\beta_{2,\nu_1}^*&\cdots&\hat\beta_{2,\nu_{r_2}}^*\end{pmatrix}^T,\qquad L=r_1+N_X+r_2.
$$

In matrix elements, $\hat X$ denotes the ordered product of its entries. For even $L$,

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat X|\Phi_{2;\kappa_2}\rangle\equiv\langle\Phi_{1;\kappa_1}|\hat x_1\cdots\hat x_{N_X}|\Phi_{2;\kappa_2}\rangle=\langle\Phi_1|\Phi_2\rangle\operatorname{pf}(S)}.
$$

Here $S_{ij}=\langle\hat s_i\hat s_j\rangle_{12}$ for $i<j$, with $S=-S^T$. For odd $L$, the matrix element vanishes.

##### One-Body Matrix Elements between Multiquasiparticle States

For single-particle matrix elements $O_{\alpha\beta}=\langle\alpha|\hat O|\beta\rangle$,

$$
\hat O^{(1)}=\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}O_{\alpha\beta}\hat c_\alpha^*\hat c_\beta.
$$

Choose the ordered insertion group

$$
\hat X(\alpha,\beta)=(\hat c_\alpha^*,\hat c_\beta).
$$

Let $S(\alpha,\beta)$ be the corresponding contraction matrix. For even $r_1+r_2$,

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat O^{(1)}|\Phi_{2;\kappa_2}\rangle=\langle\Phi_1|\Phi_2\rangle\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}O_{\alpha\beta}\operatorname{pf}[S(\alpha,\beta)]}.
$$

For odd $r_1+r_2$, the matrix element vanishes. For $r_1=r_2=0$,

$$
\langle\Phi_1|\hat O^{(1)}|\Phi_2\rangle=\langle\Phi_1|\Phi_2\rangle\operatorname{Tr}(O\rho).
$$

##### Two-Body Matrix Elements between Multiquasiparticle States

For unsymmetrized matrix elements $V_{\alpha\beta\gamma\delta}=\langle\alpha\otimes\beta|\hat O|\gamma\otimes\delta\rangle$,

$$
\hat O^{(2)}=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}V_{\alpha\beta\gamma\delta}\hat c_\alpha^*\hat c_\beta^*\hat c_\delta\hat c_\gamma.
$$

Choose the ordered insertion group

$$
\hat X(\alpha,\beta,\gamma,\delta)=(\hat c_\alpha^*,\hat c_\beta^*,\hat c_\delta,\hat c_\gamma).
$$

Let $S(\alpha,\beta,\gamma,\delta)$ be the corresponding contraction matrix. For even $r_1+r_2$,

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat O^{(2)}|\Phi_{2;\kappa_2}\rangle=\frac12\langle\Phi_1|\Phi_2\rangle\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}V_{\alpha\beta\gamma\delta}\operatorname{pf}[S(\alpha,\beta,\gamma,\delta)]}.
$$

For odd $r_1+r_2$, the matrix element vanishes. For $r_1=r_2=0$,

$$
\operatorname{pf}[S(\alpha,\beta,\gamma,\delta)]=\rho_{\gamma\alpha}\rho_{\delta\beta}-\rho_{\delta\alpha}\rho_{\gamma\beta}-\bar\kappa_{\alpha\beta}\kappa_{\gamma\delta}.
$$

See [HFB Quantum-Number Projection](hfb_projection.md) for quantum-number restoration.
