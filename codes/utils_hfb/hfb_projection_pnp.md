# Particle-Number Projection of Q-Type and P-Type Interactions

## Conventions and Basic Contractions

At fixed gauge angle $\phi$ and rotation angles $\Omega$,

$$
|\Phi_2(\phi,\Omega)\rangle=e^{-i\phi\hat N}\hat R(\Omega)|\Phi_2\rangle.
$$

$$
\boxed{\langle\hat A\hat B\rangle(\phi,\Omega)=\frac{\langle\Phi_1|\hat A\hat B|\Phi_2(\phi,\Omega)\rangle}{\langle\Phi_1|\Phi_2(\phi,\Omega)\rangle}.}
$$

The same convention applies to longer operator products.

### Indices and Operator Order

- Single-particle orbitals: $\alpha,\beta,\gamma,\delta$.
- Quasiparticle orbitals: $\mu,\nu$.
- Reference states: $1,2$.
- Insertion positions: $i=1,2,3,4$.

The right quasiparticle operators transform with the reference state:

$$
\hat\beta_{2;\nu}^\dagger(\phi,\Omega)=e^{-i\phi\hat N}\hat R(\Omega)\hat\beta_{2;\nu}^\dagger\hat R^\dagger(\Omega)e^{i\phi\hat N}.
$$

Below, $\hat\beta_{2;\nu}^\dagger$ denotes this transformed operator.

For configurations containing $r_1,r_2$ quasiparticles,

$$
\mu_1<\cdots<\mu_{r_1},\qquad \nu_1<\cdots<\nu_{r_2}.
$$

The insertion kernel is

$$
\boxed{\langle\Phi_1|\Phi_2(\phi,\Omega)\rangle\left\langle\hat\beta_{1;\mu_{r_1}}\cdots\hat\beta_{1;\mu_1}\,\hat x_1\hat x_2\hat x_3\hat x_4\,\hat\beta_{2;\nu_1}^\dagger(\phi,\Omega)\cdots\hat\beta_{2;\nu_{r_2}}^\dagger(\phi,\Omega)\right\rangle(\phi,\Omega).}
$$

Hermitian conjugation reverses the left configuration. Each $\hat x_{i,\alpha}$ is a complex linear combination of particle creation and annihilation operators; its definition and orbital argument depend on the interaction.

### Basic Contractions

$$
\rho_{\alpha\beta}=\langle\hat c_\alpha^\dagger\hat c_\beta\rangle(\phi,\Omega),\qquad \kappa_{\alpha\beta}=\langle\hat c_\alpha\hat c_\beta\rangle(\phi,\Omega),\qquad \bar\kappa_{\alpha\beta}=\langle\hat c_\alpha^\dagger\hat c_\beta^\dagger\rangle(\phi,\Omega).
$$

$$
\kappa^T=-\kappa,\qquad \bar\kappa^T=-\bar\kappa,\qquad \langle\hat c_\alpha\hat c_\beta^\dagger\rangle(\phi,\Omega)=\delta_{\alpha\beta}-\rho_{\beta\alpha}.
$$

Neither $\rho^\dagger=\rho$ nor $\bar\kappa=\kappa^*$ is assumed.

$$
\langle\hat\beta_{1;\mu}\hat\beta_{1;\nu}\rangle(\phi,\Omega),\qquad \langle\hat\beta_{1;\mu}\hat\beta_{2;\nu}^\dagger(\phi,\Omega)\rangle(\phi,\Omega),\qquad \langle\hat\beta_{2;\mu}^\dagger(\phi,\Omega)\hat\beta_{2;\nu}^\dagger(\phi,\Omega)\rangle(\phi,\Omega).
$$

Insertion contractions:

$$
\boxed{\eta^{Li}_{\mu\alpha}(\phi,\Omega)=\langle\hat\beta_{1;\mu}\hat x_{i,\alpha}\rangle(\phi,\Omega),\qquad \eta^{iR}_{\alpha\nu}(\phi,\Omega)=\langle\hat x_{i,\alpha}\hat\beta_{2;\nu}^\dagger(\phi,\Omega)\rangle(\phi,\Omega).}
$$

The superscripts identify operator order; the subscripts follow that order.

### Wick–Pfaffian Convention

Write the ordered chain as $\hat A_1,\ldots,\hat A_{2m}$. Contractions follow chain order; the lower triangle is fixed by antisymmetry:

$$
S_{ij}=\langle\hat A_i\hat A_j\rangle(\phi,\Omega),\qquad S=-S^T.
$$

Thus the blocks opposite $\eta^{Li}(\phi,\Omega)$ and $\eta^{iR}(\phi,\Omega)$ are their negative transposes.

$$
\boxed{\langle\hat A_1\cdots\hat A_{2m}\rangle(\phi,\Omega)=\operatorname{pf}(S)},\qquad \operatorname{pf}(\varnothing)=1.
$$

## Q-Type Interaction

### Operator

For complex $Q$,

$$
\hat Q=\sum_{\alpha\beta}Q_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta,\qquad \hat Q^\dagger=\sum_{\alpha\beta}Q^*_{\beta\alpha}\hat c_\alpha^\dagger\hat c_\beta.
$$

$$
\boxed{\hat H_Q=g_Q\hat Q^\dagger\hat Q=\hat H_Q^{(1)}+\hat H_Q^{(2)}.}
$$

Using the fermionic anticommutation relation,

$$
\hat H_Q^{(1)}=g_Q\sum_{\alpha\delta}(Q^\dagger Q)_{\alpha\delta}\hat c_\alpha^\dagger\hat c_\delta.
$$

$$
\hat H_Q^{(2)}=g_Q\sum_{\alpha\beta\gamma\delta}Q^*_{\beta\alpha}Q_{\gamma\delta}\hat c_\alpha^\dagger\hat c_\gamma^\dagger\hat c_\delta\hat c_\beta.
$$

The real coefficient $g_Q$ includes the interaction sign and strength.

### One-Body Contribution

The one-body matrix is

$$
h^Q_{\alpha\delta}=g_Q(Q^\dagger Q)_{\alpha\delta}.
$$

Its configuration kernel is

$$
\mathcal K_Q^{(1)}(\phi,\Omega)=\sum_{\alpha\delta}h^Q_{\alpha\delta}\langle\Phi_1|\hat\beta_{1;\mu_{r_1}}\cdots\hat\beta_{1;\mu_1}\hat c_\alpha^\dagger\hat c_\delta\hat\beta^\dagger_{2;\nu_1}(\phi,\Omega)\cdots\hat\beta^\dagger_{2;\nu_{r_2}}(\phi,\Omega)|\Phi_2(\phi,\Omega)\rangle.
$$

### Insertions

For the two-body contribution,

$$
\hat x_{1,\beta}=\sum_\alpha Q^*_{\beta\alpha}\hat c_\alpha^\dagger,\qquad \hat x_{2,\delta}=\sum_\gamma Q_{\gamma\delta}\hat c_\gamma^\dagger,\qquad \hat x_{3,\delta}=\hat c_\delta,\qquad \hat x_{4,\beta}=\hat c_\beta.
$$

$$
\hat H_Q^{(2)}=g_Q\sum_{\beta\delta}\hat x_{1,\beta}\hat x_{2,\delta}\hat x_{3,\delta}\hat x_{4,\beta}.
$$

### $T^0$

The internal contractions are

$$
\langle\hat x_{1,\beta}\hat x_{2,\delta}\rangle(\phi,\Omega)=(Q^*\bar\kappa Q)_{\beta\delta}.
$$

$$
\langle\hat x_{1,\beta}\hat x_{3,\delta}\rangle(\phi,\Omega)=(Q^*\rho)_{\beta\delta},\qquad \langle\hat x_{1,\beta}\hat x_{4,\beta}\rangle(\phi,\Omega)=(Q^*\rho)_{\beta\beta}.
$$

$$
\langle\hat x_{2,\delta}\hat x_{3,\delta}\rangle(\phi,\Omega)=(Q^T\rho)_{\delta\delta},\qquad \langle\hat x_{2,\delta}\hat x_{4,\beta}\rangle(\phi,\Omega)=(Q^T\rho)_{\delta\beta}.
$$

$$
\langle\hat x_{3,\delta}\hat x_{4,\beta}\rangle(\phi,\Omega)=\kappa_{\delta\beta}.
$$

$$
\boxed{T^0(\phi,\Omega)=\sum_{\beta\delta}\left[\langle\hat x_{1,\beta}\hat x_{2,\delta}\rangle(\phi,\Omega)\langle\hat x_{3,\delta}\hat x_{4,\beta}\rangle(\phi,\Omega)-\langle\hat x_{1,\beta}\hat x_{3,\delta}\rangle(\phi,\Omega)\langle\hat x_{2,\delta}\hat x_{4,\beta}\rangle(\phi,\Omega)+\langle\hat x_{1,\beta}\hat x_{4,\beta}\rangle(\phi,\Omega)\langle\hat x_{2,\delta}\hat x_{3,\delta}\rangle(\phi,\Omega)\right].}
$$

### $T^2$

Combine the left/right insertion contractions:

$$
\eta^i_{\nu\alpha}(\phi,\Omega)=\begin{pmatrix}\eta^{Li}_{\nu\alpha}(\phi,\Omega)\\-\eta^{iR}_{\alpha\nu}(\phi,\Omega)\end{pmatrix}=\begin{pmatrix}\langle\hat\beta_{1;\nu}\hat x_{i,\alpha}\rangle(\phi,\Omega)\\-\langle\hat x_{i,\alpha}\hat\beta^\dagger_{2;\nu}(\phi,\Omega)\rangle(\phi,\Omega)\end{pmatrix}.
$$

Blockwise notation: $\nu$ runs within the corresponding left/right block.

$$
W^{ij}_{\nu\nu';\alpha\beta}(\phi,\Omega)=\eta^i_{\nu\alpha}(\phi,\Omega)\eta^j_{\nu'\beta}(\phi,\Omega)-\eta^j_{\nu\beta}(\phi,\Omega)\eta^i_{\nu'\alpha}(\phi,\Omega).
$$

$$
\boxed{\begin{aligned}T^2_{\nu\nu'}(\phi,\Omega)=\sum_{\beta\delta}\Big[&-W^{12}_{\nu\nu';\beta\delta}(\phi,\Omega)\langle\hat x_{3,\delta}\hat x_{4,\beta}\rangle(\phi,\Omega)\\&+W^{13}_{\nu\nu';\beta\delta}(\phi,\Omega)\langle\hat x_{2,\delta}\hat x_{4,\beta}\rangle(\phi,\Omega)\\&-W^{14}_{\nu\nu';\beta\beta}(\phi,\Omega)\langle\hat x_{2,\delta}\hat x_{3,\delta}\rangle(\phi,\Omega)\\&-W^{23}_{\nu\nu';\delta\delta}(\phi,\Omega)\langle\hat x_{1,\beta}\hat x_{4,\beta}\rangle(\phi,\Omega)\\&+W^{24}_{\nu\nu';\delta\beta}(\phi,\Omega)\langle\hat x_{1,\beta}\hat x_{3,\delta}\rangle(\phi,\Omega)\\&-W^{34}_{\nu\nu';\delta\beta}(\phi,\Omega)\langle\hat x_{1,\beta}\hat x_{2,\delta}\rangle(\phi,\Omega)\Big].\end{aligned}}
$$

### $T^4$

$$
\Pi^{14}_{\nu\nu'}(\phi,\Omega)=\sum_\beta\left[\eta^1_{\nu\beta}(\phi,\Omega)\eta^4_{\nu'\beta}(\phi,\Omega)-\eta^1_{\nu'\beta}(\phi,\Omega)\eta^4_{\nu\beta}(\phi,\Omega)\right].
$$

$$
\Pi^{23}_{\nu\nu'}(\phi,\Omega)=\sum_\delta\left[\eta^2_{\nu\delta}(\phi,\Omega)\eta^3_{\nu'\delta}(\phi,\Omega)-\eta^2_{\nu'\delta}(\phi,\Omega)\eta^3_{\nu\delta}(\phi,\Omega)\right].
$$

$$
\bigl(\Pi^{14}(\phi,\Omega)\bigr)^T=-\Pi^{14}(\phi,\Omega),\qquad \bigl(\Pi^{23}(\phi,\Omega)\bigr)^T=-\Pi^{23}(\phi,\Omega).
$$

$$
\boxed{\begin{aligned}T^4_{\nu_1\nu_2\nu_3\nu_4}(\phi,\Omega)={}&\Pi^{14}_{\nu_1\nu_2}(\phi,\Omega)\Pi^{23}_{\nu_3\nu_4}(\phi,\Omega)\\&-\Pi^{14}_{\nu_1\nu_3}(\phi,\Omega)\Pi^{23}_{\nu_2\nu_4}(\phi,\Omega)\\&+\Pi^{14}_{\nu_1\nu_4}(\phi,\Omega)\Pi^{23}_{\nu_2\nu_3}(\phi,\Omega)\\&+\Pi^{14}_{\nu_2\nu_3}(\phi,\Omega)\Pi^{23}_{\nu_1\nu_4}(\phi,\Omega)\\&-\Pi^{14}_{\nu_2\nu_4}(\phi,\Omega)\Pi^{23}_{\nu_1\nu_3}(\phi,\Omega)\\&+\Pi^{14}_{\nu_3\nu_4}(\phi,\Omega)\Pi^{23}_{\nu_1\nu_2}(\phi,\Omega).\end{aligned}}
$$

$T^0,T^2,T^4$ describe the four-operator contribution; they exclude $g_Q$, configuration factors $F$, reference overlap, and projection weights.

The full kernel includes both contributions:

$$
\boxed{\mathcal K_Q(\phi,\Omega)=\mathcal K_Q^{(1)}(\phi,\Omega)+\mathcal K_Q^{(2)}(\phi,\Omega).}
$$

## P-Type Interaction

### Operator

For complex antisymmetric $P$,

$$
P^T=-P.
$$

$$
\hat P^\dagger=\frac12\sum_{\alpha\beta}P_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta^\dagger,\qquad \hat P=\frac12\sum_{\gamma\delta}P^*_{\gamma\delta}\hat c_\delta\hat c_\gamma.
$$

$$
\boxed{\hat H_P=g_P\hat P^\dagger\hat P=\frac{g_P}{4}\sum_{\alpha\beta\gamma\delta}P_{\alpha\beta}P^*_{\gamma\delta}\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma.}
$$

The real coefficient $g_P$ includes the interaction sign and strength. The product is already normally ordered; no one-body contribution occurs.

### Insertions

$$
\hat x_{1,\beta}=\sum_\alpha P_{\alpha\beta}\hat c_\alpha^\dagger,\qquad \hat x_{2,\beta}=\hat c_\beta^\dagger,\qquad \hat x_{3,\delta}=\hat c_\delta,\qquad \hat x_{4,\delta}=\sum_\gamma P^*_{\gamma\delta}\hat c_\gamma.
$$

$$
\hat H_P=\frac{g_P}{4}\sum_{\beta\delta}\hat x_{1,\beta}\hat x_{2,\beta}\hat x_{3,\delta}\hat x_{4,\delta}.
$$

### $T^0$

The internal contractions are

$$
\langle\hat x_{1,\beta}\hat x_{2,\beta}\rangle(\phi,\Omega)=(P^T\bar\kappa)_{\beta\beta}.
$$

$$
\langle\hat x_{1,\beta}\hat x_{3,\delta}\rangle(\phi,\Omega)=(P^T\rho)_{\beta\delta},\qquad \langle\hat x_{1,\beta}\hat x_{4,\delta}\rangle(\phi,\Omega)=(P^T\rho P^*)_{\beta\delta}.
$$

$$
\langle\hat x_{2,\beta}\hat x_{3,\delta}\rangle(\phi,\Omega)=\rho_{\beta\delta},\qquad \langle\hat x_{2,\beta}\hat x_{4,\delta}\rangle(\phi,\Omega)=(\rho P^*)_{\beta\delta}.
$$

$$
\langle\hat x_{3,\delta}\hat x_{4,\delta}\rangle(\phi,\Omega)=(\kappa P^*)_{\delta\delta}.
$$

$$
\boxed{T^0(\phi,\Omega)=\sum_{\beta\delta}\left[\langle\hat x_{1,\beta}\hat x_{2,\beta}\rangle(\phi,\Omega)\langle\hat x_{3,\delta}\hat x_{4,\delta}\rangle(\phi,\Omega)-\langle\hat x_{1,\beta}\hat x_{3,\delta}\rangle(\phi,\Omega)\langle\hat x_{2,\beta}\hat x_{4,\delta}\rangle(\phi,\Omega)+\langle\hat x_{1,\beta}\hat x_{4,\delta}\rangle(\phi,\Omega)\langle\hat x_{2,\beta}\hat x_{3,\delta}\rangle(\phi,\Omega)\right].}
$$

### $T^2$

Combine the left/right insertion contractions:

$$
\eta^i_{\nu\alpha}(\phi,\Omega)=\begin{pmatrix}\eta^{Li}_{\nu\alpha}(\phi,\Omega)\\-\eta^{iR}_{\alpha\nu}(\phi,\Omega)\end{pmatrix}=\begin{pmatrix}\langle\hat\beta_{1;\nu}\hat x_{i,\alpha}\rangle(\phi,\Omega)\\-\langle\hat x_{i,\alpha}\hat\beta^\dagger_{2;\nu}(\phi,\Omega)\rangle(\phi,\Omega)\end{pmatrix}.
$$

Blockwise notation: $\nu$ runs within the corresponding left/right block.

$$
W^{ij}_{\nu\nu';\alpha\beta}(\phi,\Omega)=\eta^i_{\nu\alpha}(\phi,\Omega)\eta^j_{\nu'\beta}(\phi,\Omega)-\eta^j_{\nu\beta}(\phi,\Omega)\eta^i_{\nu'\alpha}(\phi,\Omega).
$$

$$
\boxed{\begin{aligned}T^2_{\nu\nu'}(\phi,\Omega)=\sum_{\beta\delta}\Big[&-W^{12}_{\nu\nu';\beta\beta}(\phi,\Omega)\langle\hat x_{3,\delta}\hat x_{4,\delta}\rangle(\phi,\Omega)\\&+W^{13}_{\nu\nu';\beta\delta}(\phi,\Omega)\langle\hat x_{2,\beta}\hat x_{4,\delta}\rangle(\phi,\Omega)\\&-W^{14}_{\nu\nu';\beta\delta}(\phi,\Omega)\langle\hat x_{2,\beta}\hat x_{3,\delta}\rangle(\phi,\Omega)\\&-W^{23}_{\nu\nu';\beta\delta}(\phi,\Omega)\langle\hat x_{1,\beta}\hat x_{4,\delta}\rangle(\phi,\Omega)\\&+W^{24}_{\nu\nu';\beta\delta}(\phi,\Omega)\langle\hat x_{1,\beta}\hat x_{3,\delta}\rangle(\phi,\Omega)\\&-W^{34}_{\nu\nu';\delta\delta}(\phi,\Omega)\langle\hat x_{1,\beta}\hat x_{2,\beta}\rangle(\phi,\Omega)\Big].\end{aligned}}
$$

### $T^4$

$$
\Pi^{12}_{\nu\nu'}(\phi,\Omega)=\sum_\beta\left[\eta^1_{\nu\beta}(\phi,\Omega)\eta^2_{\nu'\beta}(\phi,\Omega)-\eta^1_{\nu'\beta}(\phi,\Omega)\eta^2_{\nu\beta}(\phi,\Omega)\right].
$$

$$
\Pi^{34}_{\nu\nu'}(\phi,\Omega)=\sum_\delta\left[\eta^3_{\nu\delta}(\phi,\Omega)\eta^4_{\nu'\delta}(\phi,\Omega)-\eta^3_{\nu'\delta}(\phi,\Omega)\eta^4_{\nu\delta}(\phi,\Omega)\right].
$$

$$
\bigl(\Pi^{12}(\phi,\Omega)\bigr)^T=-\Pi^{12}(\phi,\Omega),\qquad \bigl(\Pi^{34}(\phi,\Omega)\bigr)^T=-\Pi^{34}(\phi,\Omega).
$$

$$
\boxed{\begin{aligned}T^4_{\nu_1\nu_2\nu_3\nu_4}(\phi,\Omega)={}&\Pi^{12}_{\nu_1\nu_2}(\phi,\Omega)\Pi^{34}_{\nu_3\nu_4}(\phi,\Omega)\\&-\Pi^{12}_{\nu_1\nu_3}(\phi,\Omega)\Pi^{34}_{\nu_2\nu_4}(\phi,\Omega)\\&+\Pi^{12}_{\nu_1\nu_4}(\phi,\Omega)\Pi^{34}_{\nu_2\nu_3}(\phi,\Omega)\\&+\Pi^{12}_{\nu_2\nu_3}(\phi,\Omega)\Pi^{34}_{\nu_1\nu_4}(\phi,\Omega)\\&-\Pi^{12}_{\nu_2\nu_4}(\phi,\Omega)\Pi^{34}_{\nu_1\nu_3}(\phi,\Omega)\\&+\Pi^{12}_{\nu_3\nu_4}(\phi,\Omega)\Pi^{34}_{\nu_1\nu_2}(\phi,\Omega).\end{aligned}}
$$

$T^0,T^2,T^4$ describe the four-operator contribution; they exclude $g_P/4$, configuration factors $F$, reference overlap, and projection weights.

## Configuration Kernel Assembly

### External Quasiparticle Chain

For a configuration pair, define the ordered external chain:

$$
(\hat A_1,\ldots,\hat A_{r_1+r_2})=\left(\hat\beta_{1;\mu_{r_1}},\ldots,\hat\beta_{1;\mu_1},\hat\beta^\dagger_{2;\nu_1}(\phi,\Omega),\ldots,\hat\beta^\dagger_{2;\nu_{r_2}}(\phi,\Omega)\right).
$$

The external contraction matrix follows the Wick–Pfaffian convention:

$$
S_{ab}(\phi,\Omega)=\langle\hat A_a\hat A_b\rangle(\phi,\Omega),\qquad S=-S^T.
$$

Here $a,b,c,d$ label positions in the configuration chain.

### Configuration Factors

Let $S_{\widehat{a\cdots}}$ denote $S$ with the indicated rows and columns removed, preserving the remaining order.

$$
F^0(\phi,\Omega)=\operatorname{pf}S(\phi,\Omega).
$$

$$
F^2_{ab}(\phi,\Omega)=(-1)^{a+b-3}\operatorname{pf}S_{\widehat{ab}}(\phi,\Omega),\qquad a<b.
$$

$$
F^4_{abcd}(\phi,\Omega)=(-1)^{a+b+c+d-10}\operatorname{pf}S_{\widehat{abcd}}(\phi,\Omega),\qquad a<b<c<d.
$$

Chain positions start at $1$; $\operatorname{pf}(\varnothing)=1$. Extend $F^2,F^4$ antisymmetrically; repeated indices give zero.

### Antisymmetric Contractions

Restrict the quasiparticle indices of $T^2,T^4,\Pi$ to the configuration chain, retaining the left/right blocks.

$$
\sum_{a<b}F^2_{ab}(\phi,\Omega)T^2_{ab}(\phi,\Omega)=\frac12\sum_{ab}F^2_{ab}(\phi,\Omega)T^2_{ab}(\phi,\Omega).
$$

For the Q-type interaction,

$$
\sum_{a<b<c<d}F^4_{abcd}(\phi,\Omega)T^4_{Q;abcd}(\phi,\Omega)=\frac14\sum_{abcd}F^4_{abcd}(\phi,\Omega)\Pi^{14}_{ab}(\phi,\Omega)\Pi^{23}_{cd}(\phi,\Omega).
$$

For the P-type interaction,

$$
\sum_{a<b<c<d}F^4_{abcd}(\phi,\Omega)T^4_{P;abcd}(\phi,\Omega)=\frac14\sum_{abcd}F^4_{abcd}(\phi,\Omega)\Pi^{12}_{ab}(\phi,\Omega)\Pi^{34}_{cd}(\phi,\Omega).
$$

All unrestricted indices run from $1$ to $r_1+r_2$. The labels $Q,P$ distinguish the interaction-specific $T$ coefficients.

### Q-Type Kernel

The four-operator contribution is

$$
\boxed{\mathcal K_Q^{(2)}(\phi,\Omega)=g_Q\langle\Phi_1|\Phi_2(\phi,\Omega)\rangle\left[F^0(\phi,\Omega)T_Q^0(\phi,\Omega)+\frac12\sum_{ab}F^2_{ab}(\phi,\Omega)T^2_{Q;ab}(\phi,\Omega)+\frac14\sum_{abcd}F^4_{abcd}(\phi,\Omega)\Pi^{14}_{ab}(\phi,\Omega)\Pi^{23}_{cd}(\phi,\Omega)\right].}
$$

Add the one-body kernel defined above:

$$
\boxed{\mathcal K_Q(\phi,\Omega)=\mathcal K_Q^{(1)}(\phi,\Omega)+\mathcal K_Q^{(2)}(\phi,\Omega).}
$$

### P-Type Kernel

$$
\boxed{\mathcal K_P(\phi,\Omega)=\frac{g_P}{4}\langle\Phi_1|\Phi_2(\phi,\Omega)\rangle\left[F^0(\phi,\Omega)T_P^0(\phi,\Omega)+\frac12\sum_{ab}F^2_{ab}(\phi,\Omega)T^2_{P;ab}(\phi,\Omega)+\frac14\sum_{abcd}F^4_{abcd}(\phi,\Omega)\Pi^{12}_{ab}(\phi,\Omega)\Pi^{34}_{cd}(\phi,\Omega)\right].}
$$

The outer $1/4$ comes from the pair-operator definitions; the inner $1/4$ comes from antisymmetric contraction.

At fixed configuration pair and $(\phi,\Omega)$, the same $F^0,F^2,F^4$ serve every Q-type and P-type operator. For odd $r_1+r_2$, the kernels vanish.
