# HFB Pfaffian Algorithm

### Quasiparticle Vacuum Representation

Let $N_{\mathrm{sp}}$ denote the number of single-particle states and $|0\rangle$ the particle vacuum. Use the Bogoliubov convention in [Quasiparticle Vacuum](hfb.md#quasiparticle-vacuum), with operator vectors arranged as columns:

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

Normalization fixes $|\nu_a|^2=|\det U_a|$, but not the relative vacuum phase. The square root leaves a sign ambiguity; the Pfaffian fixes the overlap for a specified single-particle ordering and vacuum-phase convention.

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

### Generalized Wick Theorem

For $\langle\Phi_1|\Phi_2\rangle\ne0$, define the normalized transition matrix element by

$$
\langle\hat X\rangle_{12}\equiv\frac{\langle\Phi_1|\hat X|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}.
$$

Let $\hat a_1,\ldots,\hat a_{2p}$ be linear combinations of the particle creation and annihilation operators. For their specified ordering, define the skew-symmetric contraction matrix $S\in\mathbb C^{2p\times2p}$ by

$$
S_{ij}\equiv\langle\hat a_i\hat a_j\rangle_{12}\quad(i<j),\qquad S_{ji}\equiv-S_{ij},\qquad S_{ii}\equiv0.
$$

The lower-triangular entries are defined by antisymmetry, not by reversed operator products.

The generalized Wick theorem expresses the transition matrix element as a sum over all complete pairings, with fermionic permutation signs:

$$
\boxed{\langle\hat a_1\cdots\hat a_{2p}\rangle_{12}=\operatorname{pf}(S)}.
$$

For four operators,

$$
\langle\hat a_1\hat a_2\hat a_3\hat a_4\rangle_{12}=S_{12}S_{34}-S_{13}S_{24}+S_{14}S_{23}.
$$

In particular,

$$
\rho_{\alpha\beta\gamma\delta}\equiv\langle\hat c_\alpha^*\hat c_\beta^*\hat c_\delta\hat c_\gamma\rangle_{12}=\rho_{\gamma\alpha}\rho_{\delta\beta}-\rho_{\delta\alpha}\rho_{\gamma\beta}-\bar\kappa_{\alpha\beta}\kappa_{\gamma\delta}.
$$

Both Thouless vacua have even particle-number parity, so products of an odd number of these operators have vanishing transition matrix elements:

$$
\langle\hat a_1\cdots\hat a_{2p+1}\rangle_{12}=0.
$$

The elementary contractions are evaluated in [Transition Densities](#transition-densities).

##### One-Body Matrix Elements between Multiquasiparticle States

Let $N_{\mathrm{qp}}^{(a)}$ denote the number of excited quasiparticles built on $|\Phi_a\rangle$, with $a=1,2$. Define

$$
|\Phi_{1;\boldsymbol\mu^{(1)}}\rangle\equiv\hat\beta_{1,\mu_1^{(1)}}^*\cdots\hat\beta_{1,\mu_{N_{\mathrm{qp}}^{(1)}}^{(1)}}^*|\Phi_1\rangle,\qquad |\Phi_{2;\boldsymbol\mu^{(2)}}\rangle\equiv\hat\beta_{2,\mu_1^{(2)}}^*\cdots\hat\beta_{2,\mu_{N_{\mathrm{qp}}^{(2)}}^{(2)}}^*|\Phi_2\rangle.
$$

The mode labels satisfy

$$
1\leq\mu_1^{(1)}<\cdots<\mu_{N_{\mathrm{qp}}^{(1)}}^{(1)}\leq N_{\mathrm{sp}},\qquad 1\leq\mu_1^{(2)}<\cdots<\mu_{N_{\mathrm{qp}}^{(2)}}^{(2)}\leq N_{\mathrm{sp}}.
$$

Taking the adjoint reverses the operator order:

$$
\langle\Phi_{1;\boldsymbol\mu^{(1)}}|=\langle\Phi_1|\hat\beta_{1,\mu_{N_{\mathrm{qp}}^{(1)}}^{(1)}}\cdots\hat\beta_{1,\mu_1^{(1)}},\qquad \langle\Phi_{2;\boldsymbol\mu^{(2)}}|=\langle\Phi_2|\hat\beta_{2,\mu_{N_{\mathrm{qp}}^{(2)}}^{(2)}}\cdots\hat\beta_{2,\mu_1^{(2)}}.
$$

For zero excited quasiparticles, the corresponding operator chain is omitted.

For a single-particle matrix $O$ with entries $O(\alpha,\beta)$, define the one-body operator by

$$
\hat O^{(1)}\equiv\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}O(\alpha,\beta)\hat c_\alpha^*\hat c_\beta.
$$

Its matrix element between the excited states is

$$
\langle\Phi_{1;\boldsymbol\mu^{(1)}}|\hat O^{(1)}|\Phi_{2;\boldsymbol\mu^{(2)}}\rangle=\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}O(\alpha,\beta)\langle\Phi_1|\hat\beta_{1,\mu_{N_{\mathrm{qp}}^{(1)}}^{(1)}}\cdots\hat\beta_{1,\mu_1^{(1)}}\hat c_\alpha^*\hat c_\beta\hat\beta_{2,\mu_1^{(2)}}^*\cdots\hat\beta_{2,\mu_{N_{\mathrm{qp}}^{(2)}}^{(2)}}^*|\Phi_2\rangle.
$$

Assume $\langle\Phi_1|\Phi_2\rangle\ne0$. For each $(\alpha,\beta)$, define the ordered operator column of length $L$ by

$$
\hat{\boldsymbol a}(\alpha,\beta)\equiv\begin{pmatrix}\hat\beta_{1,\mu_{N_{\mathrm{qp}}^{(1)}}^{(1)}}&\cdots&\hat\beta_{1,\mu_1^{(1)}}&\hat c_\alpha^*&\hat c_\beta&\hat\beta_{2,\mu_1^{(2)}}^*&\cdots&\hat\beta_{2,\mu_{N_{\mathrm{qp}}^{(2)}}^{(2)}}^*\end{pmatrix}^T,\qquad L\equiv N_{\mathrm{qp}}^{(1)}+N_{\mathrm{qp}}^{(2)}+2.
$$

Define the skew-symmetric contraction matrix $S(\alpha,\beta)\in\mathbb C^{L\times L}$ by

$$
S_{ij}(\alpha,\beta)\equiv\frac{\langle\Phi_1|\hat a_i(\alpha,\beta)\hat a_j(\alpha,\beta)|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}\quad(i<j),\qquad S_{ji}(\alpha,\beta)=-S_{ij}(\alpha,\beta),\qquad S_{ii}(\alpha,\beta)=0.
$$

All entries reduce to the elementary contractions in [Transition Densities](#transition-densities) using

$$
\hat{\boldsymbol\beta}_1=U_1^\dagger\hat{\boldsymbol c}+V_1^\dagger\hat{\boldsymbol c}^*,\qquad \hat{\boldsymbol\beta}_2^*=V_2^T\hat{\boldsymbol c}+U_2^T\hat{\boldsymbol c}^*,
$$

and the density matrices

$$
\frac{\langle\Phi_1|\hat{\boldsymbol c}^*\hat{\boldsymbol c}^T|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}=\rho^T,\qquad \frac{\langle\Phi_1|\hat{\boldsymbol c}\hat{\boldsymbol c}^\dagger|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}=I_{N_{\mathrm{sp}}}-\rho,
$$

and the pairing matrices

$$
\frac{\langle\Phi_1|\hat{\boldsymbol c}\hat{\boldsymbol c}^T|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}=\kappa^T,\qquad \frac{\langle\Phi_1|\hat{\boldsymbol c}^*\hat{\boldsymbol c}^\dagger|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}=\bar\kappa^T.
$$

For even $N_{\mathrm{qp}}^{(1)}+N_{\mathrm{qp}}^{(2)}$, the generalized Wick theorem gives

$$
\boxed{\langle\Phi_{1;\boldsymbol\mu^{(1)}}|\hat O^{(1)}|\Phi_{2;\boldsymbol\mu^{(2)}}\rangle=\langle\Phi_1|\Phi_2\rangle\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}O(\alpha,\beta)\operatorname{pf}[S(\alpha,\beta)]}.
$$

For odd $N_{\mathrm{qp}}^{(1)}+N_{\mathrm{qp}}^{(2)}$, the matrix element vanishes by particle-number parity.

When $N_{\mathrm{qp}}^{(1)}=N_{\mathrm{qp}}^{(2)}=0$, the contraction matrix has only one independent entry, $S_{12}(\alpha,\beta)=\rho_{\beta\alpha}$, recovering

$$
\langle\Phi_1|\hat O^{(1)}|\Phi_2\rangle=\langle\Phi_1|\Phi_2\rangle\operatorname{Tr}(O\rho).
$$

##### Two-Body Matrix Elements between Multiquasiparticle States

Use the excited states and operator ordering defined in [One-Body Matrix Elements between Multiquasiparticle States](#one-body-matrix-elements-between-multiquasiparticle-states). For unsymmetrized two-body matrix elements $O(\alpha,\beta,\gamma,\delta)$, define

$$
\hat O^{(2)}\equiv\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}O(\alpha,\beta,\gamma,\delta)\hat c_\alpha^*\hat c_\beta^*\hat c_\delta\hat c_\gamma.
$$

Its matrix element between the excited states is

$$
\langle\Phi_{1;\boldsymbol\mu^{(1)}}|\hat O^{(2)}|\Phi_{2;\boldsymbol\mu^{(2)}}\rangle=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}O(\alpha,\beta,\gamma,\delta)\langle\Phi_1|\hat\beta_{1,\mu_{N_{\mathrm{qp}}^{(1)}}^{(1)}}\cdots\hat\beta_{1,\mu_1^{(1)}}\hat c_\alpha^*\hat c_\beta^*\hat c_\delta\hat c_\gamma\hat\beta_{2,\mu_1^{(2)}}^*\cdots\hat\beta_{2,\mu_{N_{\mathrm{qp}}^{(2)}}^{(2)}}^*|\Phi_2\rangle.
$$

Assume $\langle\Phi_1|\Phi_2\rangle\ne0$. For each $(\alpha,\beta,\gamma,\delta)$, define the ordered operator column by

$$
\hat{\boldsymbol a}(\alpha,\beta,\gamma,\delta)\equiv\begin{pmatrix}\hat\beta_{1,\mu_{N_{\mathrm{qp}}^{(1)}}^{(1)}}&\cdots&\hat\beta_{1,\mu_1^{(1)}}&\hat c_\alpha^*&\hat c_\beta^*&\hat c_\delta&\hat c_\gamma&\hat\beta_{2,\mu_1^{(2)}}^*&\cdots&\hat\beta_{2,\mu_{N_{\mathrm{qp}}^{(2)}}^{(2)}}^*\end{pmatrix}^T.
$$

Its length is

$$
L\equiv N_{\mathrm{qp}}^{(1)}+N_{\mathrm{qp}}^{(2)}+4.
$$

Define the skew-symmetric contraction matrix $S(\alpha,\beta,\gamma,\delta)\in\mathbb C^{L\times L}$ by

$$
S_{ij}(\alpha,\beta,\gamma,\delta)\equiv\frac{\langle\Phi_1|\hat a_i(\alpha,\beta,\gamma,\delta)\hat a_j(\alpha,\beta,\gamma,\delta)|\Phi_2\rangle}{\langle\Phi_1|\Phi_2\rangle}\quad(i<j),
$$

$$
S_{ji}(\alpha,\beta,\gamma,\delta)=-S_{ij}(\alpha,\beta,\gamma,\delta),\qquad S_{ii}(\alpha,\beta,\gamma,\delta)=0.
$$

Its entries are evaluated from the same Bogoliubov transformations and elementary contractions used for the one-body matrix elements.

For even $N_{\mathrm{qp}}^{(1)}+N_{\mathrm{qp}}^{(2)}$, the generalized Wick theorem gives

$$
\boxed{\langle\Phi_{1;\boldsymbol\mu^{(1)}}|\hat O^{(2)}|\Phi_{2;\boldsymbol\mu^{(2)}}\rangle=\frac12\langle\Phi_1|\Phi_2\rangle\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}O(\alpha,\beta,\gamma,\delta)\operatorname{pf}[S(\alpha,\beta,\gamma,\delta)]}.
$$

For odd $N_{\mathrm{qp}}^{(1)}+N_{\mathrm{qp}}^{(2)}$, the matrix element vanishes by particle-number parity.

When $N_{\mathrm{qp}}^{(1)}=N_{\mathrm{qp}}^{(2)}=0$, the contraction matrix reduces to

$$
S(\alpha,\beta,\gamma,\delta)=\begin{pmatrix}0&-\bar\kappa_{\alpha\beta}&\rho_{\delta\alpha}&\rho_{\gamma\alpha}\\\bar\kappa_{\alpha\beta}&0&\rho_{\delta\beta}&\rho_{\gamma\beta}\\-\rho_{\delta\alpha}&-\rho_{\delta\beta}&0&\kappa_{\gamma\delta}\\-\rho_{\gamma\alpha}&-\rho_{\gamma\beta}&-\kappa_{\gamma\delta}&0\end{pmatrix}.
$$

Therefore,

$$
\operatorname{pf}[S(\alpha,\beta,\gamma,\delta)]=\rho_{\gamma\alpha}\rho_{\delta\beta}-\rho_{\delta\alpha}\rho_{\gamma\beta}-\bar\kappa_{\alpha\beta}\kappa_{\gamma\delta},
$$

recovering

$$
\langle\Phi_1|\hat O^{(2)}|\Phi_2\rangle=\frac12\langle\Phi_1|\Phi_2\rangle\sum_{\alpha,\beta,\gamma,\delta=1}^{N_{\mathrm{sp}}}O(\alpha,\beta,\gamma,\delta)\rho_{\alpha\beta\gamma\delta}.
$$

### Quantum-Number Restoration

Use the two normalized Thouless vacua defined in [Quasiparticle Vacuum Overlap](#quasiparticle-vacuum-overlap).

##### Particle-Number Projection

Let $\hat N$ denote the particle-number operator, $N,N'\in\mathbb Z_{\geq0}$ particle numbers, and $\varphi\in[0,2\pi)$ the gauge angle. The projector is

$$
\boxed{\hat P^N=\frac{1}{2\pi}\int_0^{2\pi}e^{i\varphi(\hat N-N)}\,\mathrm d\varphi}.
$$

The unnormalized projected state is

$$
|\Phi_a^N\rangle=\hat P^N|\Phi_a\rangle,\qquad a=1,2.
$$

The projector satisfies

$$
(\hat P^N)^\dagger=\hat P^N,\qquad \hat P^N\hat P^{N'}=\delta_{NN'}\hat P^N,\qquad \hat N\hat P^N=N\hat P^N.
$$

##### Angular-Momentum Projection

Let $\hat{\boldsymbol J}$ denote total angular momentum, with quantum numbers $J\geq0$ and $M,K\in\{-J,-J+1,\ldots,J\}$. Use the conventions in [SO(3) and SU(2)](../group_theory/so3_su2.md#wigner-d-functions):

- $SO(3)$ with $\gamma_{\max}=2\pi$ for integer $J$;
- $SU(2)$ with $\gamma_{\max}=4\pi$ for half-odd-integer $J$.

For Euler angles $\Omega=(\alpha,\beta,\gamma)$,

$$
\hat R(\Omega)=e^{-i\alpha\hat J_z/\hbar}e^{-i\beta\hat J_y/\hbar}e^{-i\gamma\hat J_z/\hbar}.
$$

The angular-momentum projection operator is

$$
\boxed{\hat P^J_{MK}=\frac{2J+1}{4\pi\gamma_{\max}}\int_0^{2\pi}\mathrm d\alpha\int_0^\pi\sin\beta\,\mathrm d\beta\int_0^{\gamma_{\max}}\mathrm d\gamma\,D^{J*}_{MK}(\Omega)\hat R(\Omega)}.
$$

The even-particle-number Thouless vacua considered here contain only integer-$J$ components.

The unnormalized projected state is

$$
|\Phi^J_{a,MK}\rangle=\hat P^J_{MK}|\Phi_a\rangle,\qquad a=1,2.
$$

The projection operators satisfy

$$
(\hat P^J_{MK_1})^\dagger=\hat P^J_{K_1M},\qquad \hat P^J_{K_1M}\hat P^J_{MK_2}=\hat P^J_{K_1K_2}.
$$

##### Transformed Quasiparticle Vacuum

Let $g$ denote a gauge transformation, a rotation, or their product:

$$
\hat G(g)\in\left\{e^{i\varphi\hat N},\hat R(\Omega),e^{i\varphi\hat N}\hat R(\Omega)\right\},\qquad |\Phi_2(g)\rangle\equiv\hat G(g)|\Phi_2\rangle.
$$

Vacuum invariance, $\hat G(g)|0\rangle=|0\rangle$, implies

$$
\nu_2(g)\equiv\langle0|\Phi_2(g)\rangle=\nu_2.
$$

Assume the single-particle space is closed under the transformation. Let $D(g)$ denote its unitary representation:

$$
\hat G(g)\hat c_\alpha^*\hat G^\dagger(g)=\sum_{\beta=1}^{N_{\mathrm{sp}}}\hat c_\beta^*D_{\beta\alpha}(g).
$$

The Bogoliubov and Thouless matrices transform as

$$
\boxed{U_2(g)=D(g)U_2,\qquad V_2(g)=D^*(g)V_2,\qquad Z_2(g)=D(g)Z_2D^T(g)}.
$$

Thus,

$$
|\Phi_2(g)\rangle=\nu_2\exp\left[\frac12\sum_{\alpha,\beta=1}^{N_{\mathrm{sp}}}[Z_2(g)]_{\alpha\beta}\hat c_\alpha^*\hat c_\beta^*\right]|0\rangle.
$$

For $D(\varphi,\Omega)=e^{i\varphi}D(\Omega)$,

$$
\boxed{Z_2(\varphi,\Omega)=e^{2i\varphi}D(\Omega)Z_2D^T(\Omega)}.
$$

##### Projected Overlaps and Matrix Elements

Define the overlap kernel by

$$
n_{12}(g)\equiv\langle\Phi_1|\hat G(g)|\Phi_2\rangle.
$$

Substituting the transformed Thouless matrix into [Quasiparticle Vacuum Overlap](#quasiparticle-vacuum-overlap) gives

$$
\boxed{n_{12}(g)=\nu_1^*\nu_2s_{N_{\mathrm{sp}}}\operatorname{pf}\begin{pmatrix}Z_2(g)&-I_{N_{\mathrm{sp}}}\\I_{N_{\mathrm{sp}}}&-Z_1^*\end{pmatrix}}.
$$

The particle-number projected overlap is

$$
\boxed{\langle\Phi_1^N|\Phi_2^N\rangle=\langle\Phi_1|\hat P^N|\Phi_2\rangle=\frac{1}{2\pi}\int_0^{2\pi}e^{-iN\varphi}n_{12}(\varphi)\,\mathrm d\varphi}.
$$

The angular-momentum projected overlap is

$$
\mathcal N^J_{1K_1,2K_2}\equiv\langle\Phi^J_{1,MK_1}|\Phi^J_{2,MK_2}\rangle=\langle\Phi_1|\hat P^J_{K_1K_2}|\Phi_2\rangle,
$$

$$
\boxed{\mathcal N^J_{1K_1,2K_2}=\frac{2J+1}{4\pi\gamma_{\max}}\int_0^{2\pi}\mathrm d\alpha\int_0^\pi\sin\beta\,\mathrm d\beta\int_0^{\gamma_{\max}}\mathrm d\gamma\,D^{J*}_{K_1K_2}(\Omega)n_{12}(\Omega)}.
$$

For $n_{12}(g)\ne0$, the formulas in [Transition Densities](#transition-densities) apply with

$$
\mathcal A(g)=U_1^TU_2(g)^*+V_1^TV_2(g)^*,
$$

$$
\rho(g)=V_2(g)^*\mathcal A(g)^{-1}V_1^T,\qquad \kappa(g)=V_2(g)^*\mathcal A(g)^{-1}U_1^T,\qquad \bar\kappa(g)=U_2(g)^*\mathcal A(g)^{-1}V_1^T.
$$

For an operator $\hat O$, define its unnormalized kernel by

$$
o_{12}(g)\equiv\langle\Phi_1|\hat O\hat G(g)|\Phi_2\rangle.
$$

For a one-body operator with single-particle matrix $O$,

$$
o_{12}(g)=n_{12}(g)\operatorname{Tr}[O\rho(g)].
$$

Higher-order operator kernels follow from [Generalized Wick Theorem](#generalized-wick-theorem). At zero overlap, normalized transition densities are undefined; the unnormalized kernels must be evaluated without dividing by $n_{12}(g)$.

For a particle-number-conserving operator, $[\hat O,\hat N]=0$,

$$
\boxed{\langle\Phi_1^N|\hat O|\Phi_2^N\rangle=\frac{1}{2\pi}\int_0^{2\pi}e^{-iN\varphi}o_{12}(\varphi)\,\mathrm d\varphi}.
$$

For a rotational scalar operator, $[\hat O,\hat J_i]=0$ for $i=x,y,z$,

$$
\boxed{\langle\Phi^J_{1,MK_1}|\hat O|\Phi^J_{2,MK_2}\rangle=\frac{2J+1}{4\pi\gamma_{\max}}\int_0^{2\pi}\mathrm d\alpha\int_0^\pi\sin\beta\,\mathrm d\beta\int_0^{\gamma_{\max}}\mathrm d\gamma\,D^{J*}_{K_1K_2}(\Omega)o_{12}(\Omega)}.
$$
