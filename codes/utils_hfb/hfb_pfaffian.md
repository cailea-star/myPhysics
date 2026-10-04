# HFB Pfaffian Algorithm

## Quasiparticle Vacuum Representation

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

## Quasiparticle Vacuum Overlap

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

## Transition Densities

### One-Body Transition Densities

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

### Pairing Transition Densities

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

### Other Contractions

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

## Generalized Wick Theorem

Assume $\langle\Phi_1|\Phi_2\rangle\ne0$. Let $\hat s_1,\ldots,\hat s_{2p}$ be linear combinations of particle creation and annihilation operators. Define

$$
\langle\hat s_1\cdots\hat s_{2p}\rangle_{12}\equiv\frac{\langle\Phi_1|\hat s_1\cdots\hat s_{2p}|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}.
$$

For the specified operator ordering, define the contraction matrix $S\in\mathbb C^{2p\times2p}$ by

$$
S_{ij}=\langle\hat s_i\hat s_j\rangle_{12}\quad(i<j),\qquad S=\begin{pmatrix}0&S_{12}&\cdots&S_{1,2p}\\-S_{12}&0&\cdots&S_{2,2p}\\\vdots&\vdots&\ddots&\vdots\\-S_{1,2p}&-S_{2,2p}&\cdots&0\end{pmatrix}.
$$

With $\mathfrak S_{2p}$ the permutation group of $2p$ indices, the generalized Wick theorem gives

$$
\boxed{\langle\hat s_1\cdots\hat s_{2p}\rangle_{12}=\operatorname{pf}(S)=\frac{1}{2^p p!}\sum_{\sigma\in\mathfrak S_{2p}}\operatorname{sgn}(\sigma)\prod_{k=1}^{p}S_{\sigma(2k-1),\,\sigma(2k)}}.
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

### Creation and Annihilation Matrix Elements between Multiquasiparticle States

For single-particle coefficients $O_\alpha^{(+)}$ and $O_\alpha^{(-)}$,

$$
\hat O^{(+)}=\sum_{\alpha=1}^{N_{\mathrm{sp}}}O_\alpha^{(+)}\hat c_\alpha^*,\qquad \hat O^{(-)}=\sum_{\alpha=1}^{N_{\mathrm{sp}}}O_\alpha^{(-)}\hat c_\alpha.
$$

Choose the respective ordered insertion groups

$$
\hat X^{(+)}(\alpha)=(\hat c_\alpha^*),\qquad \hat X^{(-)}(\alpha)=(\hat c_\alpha).
$$

Let $S^{(\pm)}(\alpha)$ be the corresponding contraction matrices. For odd $r_1+r_2$,

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat O^{(\pm)}|\Phi_{2;\kappa_2}\rangle=\langle\Phi_1|\Phi_2\rangle\sum_{\alpha=1}^{N_{\mathrm{sp}}}O_\alpha^{(\pm)}\operatorname{pf}[S^{(\pm)}(\alpha)]}.
$$

For even $r_1+r_2$, the matrix element vanishes. For $r_1=r_2=0$,

$$
\langle\Phi_1|\hat O^{(\pm)}|\Phi_2\rangle=0.
$$

### One-Body Matrix Elements between Multiquasiparticle States

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

### Two-Body Matrix Elements between Multiquasiparticle States

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

## Pfaffian Minors and Kernel Expansion

### Pfaffian Minor Expansion

For an even-order antisymmetric matrix $S$, expansion along index $p$ gives

$$
\boxed{\operatorname{pf}(S)=\sum_{i_1\ne p}(-1)^{i_1+p+1}\operatorname{sgn}(p-i_1)\,S_{i_1p}\operatorname{pf}(S_{\hat i_1\hat p})}.
$$

Hatted indices denote deleted rows and columns, with the remaining order preserved.

The left and right configurations are specified by

$$
\kappa_1=(\mu_1,\ldots,\mu_{r_1}),\qquad \kappa_2=(\nu_1,\ldots,\nu_{r_2}).
$$

The pure quasiparticle column and its contraction matrix are

$$
\hat{\mathbf f}=(\hat\beta_{1;\mu_{r_1}},\ldots,\hat\beta_{1;\mu_1},\hat\beta_{2;\nu_1}^{*},\ldots,\hat\beta_{2;\nu_{r_2}}^{*})^T,
$$

$$
F_{i_1i_2}=\langle\hat f_{i_1}\hat f_{i_2}\rangle_{12}\quad(i_1<i_2),\qquad F_{i_2i_1}=-F_{i_1i_2},\qquad F_{i_1i_1}=0.
$$

Define the signed Pfaffian minors

$$
F^0=\operatorname{pf}(F),
$$

$$
F^1_{i_1}=(-1)^{i_1}\operatorname{pf}(F_{\hat i_1}),
$$

$$
F^2_{i_1i_2}=(-1)^{i_1+i_2}\operatorname{pf}(F_{\hat i_1\hat i_2}),\qquad i_1<i_2,
$$

$$
F^4_{i_1i_2i_3i_4}=(-1)^{i_1+i_2+i_3+i_4}\operatorname{pf}(F_{\hat i_1\hat i_2\hat i_3\hat i_4}),\qquad i_1<i_2<i_3<i_4.
$$

We use $\operatorname{pf}(\varnothing)=1$. Superscripts count removed quasiparticle positions, not matrix powers. Even $r_1+r_2$ uses $F^0,F^2,F^4$; odd $r_1+r_2$ uses $F^1$.

### Overlap Matrix Elements

With no particle-operator insertions,

$$
\hat{\boldsymbol s}=\hat{\mathbf f}=(\hat\beta_{1;\mu_{r_1}},\ldots,\hat\beta_{1;\mu_1},\hat\beta_{2;\nu_1}^{*},\ldots,\hat\beta_{2;\nu_{r_2}}^{*})^T,\qquad S=F.
$$

For even $r_1+r_2$,

$$
\boxed{\langle\Phi_{1;\kappa_1}|\Phi_{2;\kappa_2}\rangle=\langle\Phi_1|\Phi_2\rangle\,\operatorname{pf}(S)=\langle\Phi_1|\Phi_2\rangle F^0}.
$$

For odd $r_1+r_2$, the overlap vanishes. For the vacuum configurations $r_1=r_2=0$, $F^0=\operatorname{pf}(\varnothing)=1$.

### Creation and Annihilation Matrix Elements

For a single fermionic insertion,

$$
\hat{\boldsymbol s}=(\hat\beta_{1;\mu_{r_1}},\ldots,\hat\beta_{1;\mu_1},\hat x,\hat\beta_{2;\nu_1}^{*},\ldots,\hat\beta_{2;\nu_{r_2}}^{*})^T,\qquad p=r_1+1.
$$

This subsection uses positions in the complete operator column. Define

$$
F^1_{i_1}=(-1)^{i_1}\operatorname{pf}(S_{\hat i_1\hat p}),\qquad i_1\ne p.
$$

Expansion along the insertion position $p$ gives

$$
\operatorname{pf}(S)=(-1)^{p+1}\sum_{i_1\ne p}\operatorname{sgn}(p-i_1)\,F^1_{i_1}\operatorname{pf}\begin{pmatrix}0&S_{i_1p}\\-S_{i_1p}&0\end{pmatrix}.
$$

Therefore,

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat x|\Phi_{2;\kappa_2}\rangle=\langle\Phi_1|\Phi_2\rangle\,\operatorname{pf}(S)}.
$$

Creation uses $\hat x=\hat c_\alpha^*$; annihilation uses $\hat x=\hat c_\alpha$. Both share these minors. This expansion applies to odd $r_1+r_2$; for even $r_1+r_2$, the matrix element vanishes.

### One-Body Matrix Elements

For two fermionic insertions,

$$
\hat{\boldsymbol s}=(\hat\beta_{1;\mu_{r_1}},\ldots,\hat\beta_{1;\mu_1},\hat x_1,\hat x_2,\hat\beta_{2;\nu_1}^{*},\ldots,\hat\beta_{2;\nu_{r_2}}^{*})^T,\qquad p_1=r_1+1,\quad p_2=r_1+2.
$$

This subsection uses positions in the complete operator column, with $i_1,i_2\notin\{p_1,p_2\}$. The corresponding minors are

$$
F^0=\operatorname{pf}(S_{\hat p_1\hat p_2}),
$$

$$
F^2_{i_1i_2}=(-1)^{i_1+i_2}\operatorname{pf}(S_{\hat i_1\hat i_2\hat p_1\hat p_2}),\qquad i_1<i_2.
$$

Right quasiparticle positions increase by $2$, leaving the sign unchanged. Expansion gives

$$
\operatorname{pf}(S)=F^0\operatorname{pf}\begin{pmatrix}0&S_{p_1p_2}\\-S_{p_1p_2}&0\end{pmatrix}-\sum_{\substack{i_1<i_2\\i_1,i_2\notin\{p_1,p_2\}}}F^2_{i_1i_2}\operatorname{pf}\begin{pmatrix}0&0&S_{i_1p_1}&S_{i_1p_2}\\0&0&S_{i_2p_1}&S_{i_2p_2}\\-S_{i_1p_1}&-S_{i_2p_1}&0&0\\-S_{i_1p_2}&-S_{i_2p_2}&0&0\end{pmatrix}.
$$

Therefore,

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat x_1\hat x_2|\Phi_{2;\kappa_2}\rangle=\langle\Phi_1|\Phi_2\rangle\,\operatorname{pf}(S)}.
$$

One-body matrix elements use $\hat x_1=\hat c_\alpha^*$ and $\hat x_2=\hat c_\beta$. Different $\alpha,\beta$ share these minors. This expansion applies to even $r_1+r_2$; for odd $r_1+r_2$, the matrix element vanishes.

### Two-Body Matrix Elements

For four fermionic insertions,

$$
\hat{\boldsymbol s}=(\hat\beta_{1;\mu_{r_1}},\ldots,\hat\beta_{1;\mu_1},\hat x_1,\hat x_2,\hat x_3,\hat x_4,\hat\beta_{2;\nu_1}^{*},\ldots,\hat\beta_{2;\nu_{r_2}}^{*})^T,\qquad p_j=r_1+j,\quad j=1,2,3,4.
$$

This subsection uses positions in the complete operator column. All quasiparticle indices exclude $p_1,p_2,p_3,p_4$. Define

$$
F^0=\operatorname{pf}(S_{\hat p_1\hat p_2\hat p_3\hat p_4}),
$$

$$
F^2_{i_1i_2}=(-1)^{i_1+i_2}\operatorname{pf}(S_{\hat i_1\hat i_2\hat p_1\hat p_2\hat p_3\hat p_4}),\qquad i_1<i_2,
$$

$$
F^4_{i_1i_2i_3i_4}=(-1)^{i_1+i_2+i_3+i_4}\operatorname{pf}(S_{\hat i_1\hat i_2\hat i_3\hat i_4\hat p_1\hat p_2\hat p_3\hat p_4}),\qquad i_1<i_2<i_3<i_4.
$$

Right quasiparticle positions increase by $4$, leaving the signs unchanged. Expansion gives

$$
\begin{aligned}\operatorname{pf}(S)={}&F^0\operatorname{pf}\begin{pmatrix}0&S_{p_1p_2}&S_{p_1p_3}&S_{p_1p_4}\\-S_{p_1p_2}&0&S_{p_2p_3}&S_{p_2p_4}\\-S_{p_1p_3}&-S_{p_2p_3}&0&S_{p_3p_4}\\-S_{p_1p_4}&-S_{p_2p_4}&-S_{p_3p_4}&0\end{pmatrix}\\[4pt]&+\sum_{i_1<i_2}F^2_{i_1i_2}\sum_{1\le j_1<j_2\le4}(-1)^{j_1+j_2}\operatorname{pf}\begin{pmatrix}0&S_{p_{j_3}p_{j_4}}\\-S_{p_{j_3}p_{j_4}}&0\end{pmatrix}\operatorname{pf}\begin{pmatrix}0&0&S_{i_1p_{j_1}}&S_{i_1p_{j_2}}\\0&0&S_{i_2p_{j_1}}&S_{i_2p_{j_2}}\\-S_{i_1p_{j_1}}&-S_{i_2p_{j_1}}&0&0\\-S_{i_1p_{j_2}}&-S_{i_2p_{j_2}}&0&0\end{pmatrix}\\[4pt]&+\sum_{i_1<i_2<i_3<i_4}F^4_{i_1i_2i_3i_4}\operatorname{pf}\begin{pmatrix}0&0&0&0&S_{i_1p_1}&S_{i_1p_2}&S_{i_1p_3}&S_{i_1p_4}\\0&0&0&0&S_{i_2p_1}&S_{i_2p_2}&S_{i_2p_3}&S_{i_2p_4}\\0&0&0&0&S_{i_3p_1}&S_{i_3p_2}&S_{i_3p_3}&S_{i_3p_4}\\0&0&0&0&S_{i_4p_1}&S_{i_4p_2}&S_{i_4p_3}&S_{i_4p_4}\\-S_{i_1p_1}&-S_{i_2p_1}&-S_{i_3p_1}&-S_{i_4p_1}&0&0&0&0\\-S_{i_1p_2}&-S_{i_2p_2}&-S_{i_3p_2}&-S_{i_4p_2}&0&0&0&0\\-S_{i_1p_3}&-S_{i_2p_3}&-S_{i_3p_3}&-S_{i_4p_3}&0&0&0&0\\-S_{i_1p_4}&-S_{i_2p_4}&-S_{i_3p_4}&-S_{i_4p_4}&0&0&0&0\end{pmatrix}.\end{aligned}
$$

Here $j_3<j_4$ label the two insertions complementary to $j_1,j_2$. Therefore,

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat x_1\hat x_2\hat x_3\hat x_4|\Phi_{2;\kappa_2}\rangle=\langle\Phi_1|\Phi_2\rangle\,\operatorname{pf}(S)}.
$$

Two-body matrix elements use

$$
(\hat x_1,\hat x_2,\hat x_3,\hat x_4)=(\hat c_\alpha^*,\hat c_\beta^*,\hat c_\delta,\hat c_\gamma).
$$

Different $\alpha,\beta,\gamma,\delta$ share these minors. This expansion applies to even $r_1+r_2$; for odd $r_1+r_2$, the matrix element vanishes.

## Separable Two-Body Interactions

Contractions and signed minors follow the preceding sections. Indices $i_1,i_2,\ldots$ label positions in the pure quasiparticle column $\hat{\mathbf f}$.

### Q-Type Interaction

For complex $Q$ and real $g_Q$,

$$
\hat Q=\sum_{\alpha\beta}Q_{\alpha\beta}\hat c_\alpha^*\hat c_\beta,\qquad \hat Q^\dagger=\sum_{\alpha\beta}Q_{\beta\alpha}^*\hat c_\alpha^*\hat c_\beta.
$$

$$
\boxed{\hat H_Q=g_Q\hat Q^\dagger\hat Q=\hat H_Q^{(1)}+\hat H_Q^{(2)}.}
$$

$$
\hat H_Q^{(1)}=\sum_{\alpha\delta}h^Q_{\alpha\delta}\hat c_\alpha^*\hat c_\delta,\qquad h^Q=g_QQ^\dagger Q.
$$

$$
\hat H_Q^{(2)}=g_Q\sum_{\alpha\beta\gamma\delta}Q_{\beta\alpha}^*Q_{\gamma\delta}\hat c_\alpha^*\hat c_\gamma^*\hat c_\delta\hat c_\beta.
$$

Define the insertions

$$
\hat x_1=\sum_\alpha Q_{\beta\alpha}^*\hat c_\alpha^*,\qquad \hat x_2=\sum_\gamma Q_{\gamma\delta}\hat c_\gamma^*,\qquad \hat x_3=\hat c_\delta,\qquad \hat x_4=\hat c_\beta.
$$

Here $\hat x_1,\hat x_4$ carry $\beta$; $\hat x_2,\hat x_3$ carry $\delta$. Orbital arguments are implicit below.

$$
\hat H_Q^{(2)}=g_Q\sum_{\beta\delta}\hat x_1\hat x_2\hat x_3\hat x_4.
$$

**$T^0$**

$$
\boxed{T_Q^0=\sum_{\beta\delta}\left[\langle\hat x_1\hat x_2\rangle_{12}\langle\hat x_3\hat x_4\rangle_{12}-\langle\hat x_1\hat x_3\rangle_{12}\langle\hat x_2\hat x_4\rangle_{12}+\langle\hat x_1\hat x_4\rangle_{12}\langle\hat x_2\hat x_3\rangle_{12}\right].}
$$

**$T^2$**

Arrange the insertion contractions in configuration-chain order:

$$
\eta^j=\begin{pmatrix}\bigl(\langle\hat\beta_{1;\mu_{r_1+1-i}}\hat x_j\rangle_{12}\bigr)_{i=1}^{r_1}\\-\bigl(\langle\hat x_j\hat\beta_{2;\nu_i}^{*}\rangle_{12}\bigr)_{i=1}^{r_2}\end{pmatrix},\qquad j=1,2,3,4.
$$

$$
W^{j_1j_2}_{i_1i_2}=\det\begin{pmatrix}\eta^{j_1}_{i_1}&\eta^{j_2}_{i_1}\\\eta^{j_1}_{i_2}&\eta^{j_2}_{i_2}\end{pmatrix},\qquad W^{j_1j_2}_{i_2i_1}=-W^{j_1j_2}_{i_1i_2}.
$$

$$
\boxed{\begin{aligned}T^2_{Q;i_1i_2}=\sum_{\beta\delta}\Big[&-W^{12}_{i_1i_2}\langle\hat x_3\hat x_4\rangle_{12}+W^{13}_{i_1i_2}\langle\hat x_2\hat x_4\rangle_{12}-W^{14}_{i_1i_2}\langle\hat x_2\hat x_3\rangle_{12}\\&-W^{23}_{i_1i_2}\langle\hat x_1\hat x_4\rangle_{12}+W^{24}_{i_1i_2}\langle\hat x_1\hat x_3\rangle_{12}-W^{34}_{i_1i_2}\langle\hat x_1\hat x_2\rangle_{12}\Big].\end{aligned}}
$$

This convention enters the kernel as $-F^2T_Q^2$; the signed minors remain unchanged.

**$T^4$**

$$
\Pi^{14}_{i_1i_2}=\sum_\beta W^{14}_{i_1i_2},\qquad \Pi^{23}_{i_1i_2}=\sum_\delta W^{23}_{i_1i_2}.
$$

$$
(\Pi^{14})^T=-\Pi^{14},\qquad (\Pi^{23})^T=-\Pi^{23}.
$$

$$
\boxed{\begin{aligned}T^4_{Q;i_1i_2i_3i_4}={}&\Pi^{14}_{i_1i_2}\Pi^{23}_{i_3i_4}-\Pi^{14}_{i_1i_3}\Pi^{23}_{i_2i_4}+\Pi^{14}_{i_1i_4}\Pi^{23}_{i_2i_3}\\&+\Pi^{14}_{i_2i_3}\Pi^{23}_{i_1i_4}-\Pi^{14}_{i_2i_4}\Pi^{23}_{i_1i_3}+\Pi^{14}_{i_3i_4}\Pi^{23}_{i_1i_2}.\end{aligned}}
$$

The coefficients $T_Q^0,T_Q^2,T_Q^4$ exclude $g_Q$, the signed minors, and the reference overlap.

### P-Type Interaction

For complex antisymmetric $P$ and real $g_P$,

$$
P^T=-P.
$$

$$
\hat P^\dagger=\frac12\sum_{\alpha\beta}P_{\alpha\beta}\hat c_\alpha^*\hat c_\beta^*,\qquad \hat P=\frac12\sum_{\gamma\delta}P_{\gamma\delta}^*\hat c_\delta\hat c_\gamma.
$$

$$
\boxed{\hat H_P=g_P\hat P^\dagger\hat P=\frac{g_P}{4}\sum_{\alpha\beta\gamma\delta}P_{\alpha\beta}P_{\gamma\delta}^*\hat c_\alpha^*\hat c_\beta^*\hat c_\delta\hat c_\gamma.}
$$

Define the insertions

$$
\hat x_1=\sum_\alpha P_{\alpha\beta}\hat c_\alpha^*,\qquad \hat x_2=\hat c_\beta^*,\qquad \hat x_3=\hat c_\delta,\qquad \hat x_4=\sum_\gamma P_{\gamma\delta}^*\hat c_\gamma.
$$

Here $\hat x_1,\hat x_2$ carry $\beta$; $\hat x_3,\hat x_4$ carry $\delta$.

$$
\hat H_P=\frac{g_P}{4}\sum_{\beta\delta}\hat x_1\hat x_2\hat x_3\hat x_4.
$$

**$T^0$**

$$
\boxed{T_P^0=\sum_{\beta\delta}\left[\langle\hat x_1\hat x_2\rangle_{12}\langle\hat x_3\hat x_4\rangle_{12}-\langle\hat x_1\hat x_3\rangle_{12}\langle\hat x_2\hat x_4\rangle_{12}+\langle\hat x_1\hat x_4\rangle_{12}\langle\hat x_2\hat x_3\rangle_{12}\right].}
$$

**$T^2$**

Use the same definitions of $\eta^j$ and $W^{j_1j_2}$, with the P-type insertions.

$$
\boxed{\begin{aligned}T^2_{P;i_1i_2}=\sum_{\beta\delta}\Big[&-W^{12}_{i_1i_2}\langle\hat x_3\hat x_4\rangle_{12}+W^{13}_{i_1i_2}\langle\hat x_2\hat x_4\rangle_{12}-W^{14}_{i_1i_2}\langle\hat x_2\hat x_3\rangle_{12}\\&-W^{23}_{i_1i_2}\langle\hat x_1\hat x_4\rangle_{12}+W^{24}_{i_1i_2}\langle\hat x_1\hat x_3\rangle_{12}-W^{34}_{i_1i_2}\langle\hat x_1\hat x_2\rangle_{12}\Big].\end{aligned}}
$$

**$T^4$**

$$
\Pi^{12}_{i_1i_2}=\sum_\beta W^{12}_{i_1i_2},\qquad \Pi^{34}_{i_1i_2}=\sum_\delta W^{34}_{i_1i_2}.
$$

$$
(\Pi^{12})^T=-\Pi^{12},\qquad (\Pi^{34})^T=-\Pi^{34}.
$$

$$
\boxed{\begin{aligned}T^4_{P;i_1i_2i_3i_4}={}&\Pi^{12}_{i_1i_2}\Pi^{34}_{i_3i_4}-\Pi^{12}_{i_1i_3}\Pi^{34}_{i_2i_4}+\Pi^{12}_{i_1i_4}\Pi^{34}_{i_2i_3}\\&+\Pi^{12}_{i_2i_3}\Pi^{34}_{i_1i_4}-\Pi^{12}_{i_2i_4}\Pi^{34}_{i_1i_3}+\Pi^{12}_{i_3i_4}\Pi^{34}_{i_1i_2}.\end{aligned}}
$$

The coefficients $T_P^0,T_P^2,T_P^4$ exclude $g_P/4$, the signed minors, and the reference overlap.

### Configuration Kernel Assembly

Extend $F^2,F^4$ antisymmetrically; repeated indices vanish. All indices below run over the configuration chain.

$$
\sum_{i_1<i_2}F^2_{i_1i_2}T^2_{i_1i_2}=\frac12\sum_{i_1i_2}F^2_{i_1i_2}T^2_{i_1i_2}.
$$

$$
\sum_{i_1<i_2<i_3<i_4}F^4_{i_1i_2i_3i_4}T^4_{Q;i_1i_2i_3i_4}=\frac14\sum_{i_1i_2i_3i_4}F^4_{i_1i_2i_3i_4}\Pi^{14}_{i_1i_2}\Pi^{23}_{i_3i_4}.
$$

$$
\sum_{i_1<i_2<i_3<i_4}F^4_{i_1i_2i_3i_4}T^4_{P;i_1i_2i_3i_4}=\frac14\sum_{i_1i_2i_3i_4}F^4_{i_1i_2i_3i_4}\Pi^{12}_{i_1i_2}\Pi^{34}_{i_3i_4}.
$$

For the Q-type two-body contribution,

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat H_Q^{(2)}|\Phi_{2;\kappa_2}\rangle=g_Q\langle\Phi_1|\Phi_2\rangle\left[F^0T_Q^0-\frac12\sum_{i_1i_2}F^2_{i_1i_2}T^2_{Q;i_1i_2}+\frac14\sum_{i_1i_2i_3i_4}F^4_{i_1i_2i_3i_4}\Pi^{14}_{i_1i_2}\Pi^{23}_{i_3i_4}\right].}
$$

The full $g_Q\hat Q^\dagger\hat Q$ kernel additionally includes

$$
\langle\Phi_{1;\kappa_1}|\hat H_Q^{(1)}|\Phi_{2;\kappa_2}\rangle=\sum_{\alpha\delta}h^Q_{\alpha\delta}\langle\Phi_{1;\kappa_1}|\hat c_\alpha^*\hat c_\delta|\Phi_{2;\kappa_2}\rangle.
$$

For the P-type interaction,

$$
\boxed{\langle\Phi_{1;\kappa_1}|\hat H_P|\Phi_{2;\kappa_2}\rangle=\frac{g_P}{4}\langle\Phi_1|\Phi_2\rangle\left[F^0T_P^0-\frac12\sum_{i_1i_2}F^2_{i_1i_2}T^2_{P;i_1i_2}+\frac14\sum_{i_1i_2i_3i_4}F^4_{i_1i_2i_3i_4}\Pi^{12}_{i_1i_2}\Pi^{34}_{i_3i_4}\right].}
$$

The outer $1/4$ comes from the pair definitions; the inner $1/4$ comes from antisymmetric contraction.

For fixed reference states and configurations, the same $F^0,F^2,F^4$ serve all Q-type and P-type terms. For odd $r_1+r_2$, these kernels vanish.
