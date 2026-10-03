# Particle-Number Projection of Q-Type and P-Type Interactions

## Conventions and Basic Contractions

At gauge angle $\phi$ and Euler angles $\Omega$,

$$
|\Phi_2(\phi,\Omega)\rangle=e^{-i\phi\hat N}\hat R(\Omega)|\Phi_2\rangle.
$$

$$
\boxed{\langle\hat A\hat B\rangle(\phi,\Omega)=\frac{\langle\Phi_1|\hat A\hat B|\Phi_2(\phi,\Omega)\rangle}{\langle\Phi_1|\Phi_2(\phi,\Omega)\rangle}.}
$$

The same normalization applies to longer operator products.

### Indices and Operator Order

- Single-particle orbitals: $\alpha,\beta,\gamma,\delta$.
- Quasiparticle orbitals: $\mu,\nu$.
- Quasiparticle chain positions: $i_1,i_2,\ldots$.
- Insertion labels: $j=1,2,3,4$; full-chain positions: $p_j$.

The right quasiparticle operators transform as

$$
\hat\beta_{2;\nu}^{\dagger}(\phi,\Omega)=e^{-i\phi\hat N}\hat R(\Omega)\hat\beta_{2;\nu}^{\dagger}\hat R^{\dagger}(\Omega)e^{i\phi\hat N}.
$$

For $\mu_1<\cdots<\mu_{r_1}$ and $\nu_1<\cdots<\nu_{r_2}$, the ordered chain is

$$
\hat{\boldsymbol s}=(\hat\beta_{1;\mu_{r_1}},\ldots,\hat\beta_{1;\mu_1},\hat x_1,\hat x_2,\hat x_3,\hat x_4,\hat\beta_{2;\nu_1}^{\dagger}(\phi,\Omega),\ldots,\hat\beta_{2;\nu_{r_2}}^{\dagger}(\phi,\Omega))^T,\qquad p_j=r_1+j.
$$

Each $\hat x_{j,\alpha}$ is a linear combination of particle operators, defined separately for Q-type and P-type interactions.

### Basic Contractions

Particle contractions:

$$
\langle\hat c_\alpha^\dagger\hat c_\beta\rangle(\phi,\Omega),\qquad \langle\hat c_\alpha\hat c_\beta\rangle(\phi,\Omega),\qquad \langle\hat c_\alpha^\dagger\hat c_\beta^\dagger\rangle(\phi,\Omega).
$$

$$
\langle\hat c_\alpha\hat c_\beta^\dagger\rangle(\phi,\Omega)=\delta_{\alpha\beta}-\langle\hat c_\beta^\dagger\hat c_\alpha\rangle(\phi,\Omega).
$$

Quasiparticle contractions:

$$
\langle\hat\beta_{1;\mu}\hat\beta_{1;\nu}\rangle(\phi,\Omega),\qquad \langle\hat\beta_{1;\mu}\hat\beta_{2;\nu}^\dagger(\phi,\Omega)\rangle(\phi,\Omega),\qquad \langle\hat\beta_{2;\mu}^\dagger(\phi,\Omega)\hat\beta_{2;\nu}^\dagger(\phi,\Omega)\rangle(\phi,\Omega).
$$

Insertion contractions:

$$
\boxed{\eta^{Lj}_{\mu\alpha}(\phi,\Omega)=\langle\hat\beta_{1;\mu}\hat x_{j,\alpha}\rangle(\phi,\Omega),\qquad \eta^{jR}_{\alpha\nu}(\phi,\Omega)=\langle\hat x_{j,\alpha}\hat\beta_{2;\nu}^{\dagger}(\phi,\Omega)\rangle(\phi,\Omega).}
$$

Subscripts follow operator order; opposite matrix blocks are negative transposes.

### Wick–Pfaffian Convention

Upper-triangle entries follow the ordered chain:

$$
S_{ab}(\phi,\Omega)=\langle\hat s_a\hat s_b\rangle(\phi,\Omega)\quad(a<b),\qquad S=-S^T.
$$

$$
\boxed{\langle\hat s_1\cdots\hat s_{r_1+r_2+4}\rangle(\phi,\Omega)=\operatorname{pf}S(\phi,\Omega)},\qquad r_1+r_2\ \text{even},\qquad \operatorname{pf}(\varnothing)=1.
$$

The unnormalized insertion kernel is

$$
\boxed{\langle\Phi_1|\Phi_2(\phi,\Omega)\rangle\,\operatorname{pf}S(\phi,\Omega).}
$$

## Q-Type Interaction

### Operator

For complex $Q$ and real coupling $g_Q$,

$$
\hat Q=\sum_{\alpha\beta}Q_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta,\qquad \hat Q^\dagger=\sum_{\alpha\beta}Q_{\beta\alpha}^*\hat c_\alpha^\dagger\hat c_\beta.
$$

$$
\boxed{\hat H_Q=g_Q\hat Q^\dagger\hat Q=\hat H_Q^{(1)}+\hat H_Q^{(2)}.}
$$

$$
\hat H_Q^{(1)}=\sum_{\alpha\delta}h^Q_{\alpha\delta}\hat c_\alpha^\dagger\hat c_\delta,\qquad h^Q=g_QQ^\dagger Q.
$$

$$
\hat H_Q^{(2)}=g_Q\sum_{\alpha\beta\gamma\delta}Q_{\beta\alpha}^*Q_{\gamma\delta}\hat c_\alpha^\dagger\hat c_\gamma^\dagger\hat c_\delta\hat c_\beta.
$$

### Insertions

$$
\hat x_1=\sum_\alpha Q_{\beta\alpha}^*\hat c_\alpha^\dagger,\qquad \hat x_2=\sum_\gamma Q_{\gamma\delta}\hat c_\gamma^\dagger,\qquad \hat x_3=\hat c_\delta,\qquad \hat x_4=\hat c_\beta.
$$

The orbital arguments are implicit: $\hat x_1,\hat x_4$ carry $\beta$; $\hat x_2,\hat x_3$ carry $\delta$.

$$
\hat H_Q^{(2)}=g_Q\sum_{\beta\delta}\hat x_1\hat x_2\hat x_3\hat x_4.
$$

### $T^0$

$$
\boxed{T^0(\phi,\Omega)=\sum_{\beta\delta}\left[\langle\hat x_1\hat x_2\rangle(\phi,\Omega)\langle\hat x_3\hat x_4\rangle(\phi,\Omega)-\langle\hat x_1\hat x_3\rangle(\phi,\Omega)\langle\hat x_2\hat x_4\rangle(\phi,\Omega)+\langle\hat x_1\hat x_4\rangle(\phi,\Omega)\langle\hat x_2\hat x_3\rangle(\phi,\Omega)\right].}
$$

### $T^2$

Combine the left/right contraction columns:

$$
\eta^j(\phi,\Omega)=\begin{pmatrix}\bigl(\langle\hat\beta_{1;\mu}\hat x_j\rangle(\phi,\Omega)\bigr)_\mu\\-\bigl(\langle\hat x_j\hat\beta_{2;\nu}^\dagger(\phi,\Omega)\rangle(\phi,\Omega)\bigr)_\nu\end{pmatrix},\qquad j=1,2,3,4.
$$

Below, $\nu,\nu'$ label rows of the combined quasiparticle space. Orbital arguments remain implicit.

$$
W^{j_1j_2}_{\nu\nu'}(\phi,\Omega)=\det\begin{pmatrix}\eta^{j_1}_{\nu}(\phi,\Omega)&\eta^{j_2}_{\nu}(\phi,\Omega)\\\eta^{j_1}_{\nu'}(\phi,\Omega)&\eta^{j_2}_{\nu'}(\phi,\Omega)\end{pmatrix}.
$$

$$
\boxed{\begin{aligned}T^2_{\nu\nu'}(\phi,\Omega)=\sum_{\beta\delta}\Big[&-W^{12}_{\nu\nu'}(\phi,\Omega)\langle\hat x_3\hat x_4\rangle(\phi,\Omega)\\&+W^{13}_{\nu\nu'}(\phi,\Omega)\langle\hat x_2\hat x_4\rangle(\phi,\Omega)\\&-W^{14}_{\nu\nu'}(\phi,\Omega)\langle\hat x_2\hat x_3\rangle(\phi,\Omega)\\&-W^{23}_{\nu\nu'}(\phi,\Omega)\langle\hat x_1\hat x_4\rangle(\phi,\Omega)\\&+W^{24}_{\nu\nu'}(\phi,\Omega)\langle\hat x_1\hat x_3\rangle(\phi,\Omega)\\&-W^{34}_{\nu\nu'}(\phi,\Omega)\langle\hat x_1\hat x_2\rangle(\phi,\Omega)\Big].\end{aligned}}
$$

Here $T^2$ is the negative of the full two-index coefficient in the Pfaffian expansion; the signed minors $F^2$ are unchanged:

$$
-\sum_{i_1<i_2}F^2_{i_1i_2}T^2_{i_1i_2}=-\frac12\sum_{i_1i_2}F^2_{i_1i_2}T^2_{i_1i_2}.
$$

### $T^4$

$$
\Pi^{14}_{\nu\nu'}(\phi,\Omega)=\sum_\beta W^{14}_{\nu\nu'}(\phi,\Omega),\qquad \Pi^{23}_{\nu\nu'}(\phi,\Omega)=\sum_\delta W^{23}_{\nu\nu'}(\phi,\Omega).
$$

$$
[\Pi^{14}(\phi,\Omega)]^T=-\Pi^{14}(\phi,\Omega),\qquad [\Pi^{23}(\phi,\Omega)]^T=-\Pi^{23}(\phi,\Omega).
$$

$$
\boxed{\begin{aligned}T^4_{\nu_1\nu_2\nu_3\nu_4}(\phi,\Omega)={}&\Pi^{14}_{\nu_1\nu_2}(\phi,\Omega)\Pi^{23}_{\nu_3\nu_4}(\phi,\Omega)\\&-\Pi^{14}_{\nu_1\nu_3}(\phi,\Omega)\Pi^{23}_{\nu_2\nu_4}(\phi,\Omega)\\&+\Pi^{14}_{\nu_1\nu_4}(\phi,\Omega)\Pi^{23}_{\nu_2\nu_3}(\phi,\Omega)\\&+\Pi^{14}_{\nu_2\nu_3}(\phi,\Omega)\Pi^{23}_{\nu_1\nu_4}(\phi,\Omega)\\&-\Pi^{14}_{\nu_2\nu_4}(\phi,\Omega)\Pi^{23}_{\nu_1\nu_3}(\phi,\Omega)\\&+\Pi^{14}_{\nu_3\nu_4}(\phi,\Omega)\Pi^{23}_{\nu_1\nu_2}(\phi,\Omega).\end{aligned}}
$$

$T^0,T^2,T^4$ exclude $g_Q$, configuration minors, reference overlap, and projection weights.

$$
\boxed{\mathcal K_Q(\phi,\Omega)=\mathcal K_Q^{(1)}(\phi,\Omega)+\mathcal K_Q^{(2)}(\phi,\Omega).}
$$

## P-Type Interaction

### Operator

For complex antisymmetric $P$ and real coupling $g_P$,

$$
P^T=-P.
$$

$$
\hat P^\dagger=\frac12\sum_{\alpha\beta}P_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta^\dagger,\qquad \hat P=\frac12\sum_{\gamma\delta}P_{\gamma\delta}^*\hat c_\delta\hat c_\gamma.
$$

$$
\boxed{\hat H_P=g_P\hat P^\dagger\hat P=\frac{g_P}{4}\sum_{\alpha\beta\gamma\delta}P_{\alpha\beta}P_{\gamma\delta}^*\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma.}
$$

The product is normally ordered; no one-body contribution occurs.

### Insertions

$$
\hat x_1=\sum_\alpha P_{\alpha\beta}\hat c_\alpha^\dagger,\qquad \hat x_2=\hat c_\beta^\dagger,\qquad \hat x_3=\hat c_\delta,\qquad \hat x_4=\sum_\gamma P_{\gamma\delta}^*\hat c_\gamma.
$$

The orbital arguments are implicit: $\hat x_1,\hat x_2$ carry $\beta$; $\hat x_3,\hat x_4$ carry $\delta$.

$$
\hat H_P=\frac{g_P}{4}\sum_{\beta\delta}\hat x_1\hat x_2\hat x_3\hat x_4.
$$

### $T^0$

$$
\boxed{T^0(\phi,\Omega)=\sum_{\beta\delta}\left[\langle\hat x_1\hat x_2\rangle(\phi,\Omega)\langle\hat x_3\hat x_4\rangle(\phi,\Omega)-\langle\hat x_1\hat x_3\rangle(\phi,\Omega)\langle\hat x_2\hat x_4\rangle(\phi,\Omega)+\langle\hat x_1\hat x_4\rangle(\phi,\Omega)\langle\hat x_2\hat x_3\rangle(\phi,\Omega)\right].}
$$

### $T^2$

Combine the left/right contraction columns:

$$
\eta^j(\phi,\Omega)=\begin{pmatrix}\bigl(\langle\hat\beta_{1;\mu}\hat x_j\rangle(\phi,\Omega)\bigr)_\mu\\-\bigl(\langle\hat x_j\hat\beta_{2;\nu}^\dagger(\phi,\Omega)\rangle(\phi,\Omega)\bigr)_\nu\end{pmatrix},\qquad j=1,2,3,4.
$$

Below, $\nu,\nu'$ label rows of the combined quasiparticle space. Orbital arguments remain implicit.

$$
W^{j_1j_2}_{\nu\nu'}(\phi,\Omega)=\det\begin{pmatrix}\eta^{j_1}_{\nu}(\phi,\Omega)&\eta^{j_2}_{\nu}(\phi,\Omega)\\\eta^{j_1}_{\nu'}(\phi,\Omega)&\eta^{j_2}_{\nu'}(\phi,\Omega)\end{pmatrix}.
$$

$$
\boxed{\begin{aligned}T^2_{\nu\nu'}(\phi,\Omega)=\sum_{\beta\delta}\Big[&-W^{12}_{\nu\nu'}(\phi,\Omega)\langle\hat x_3\hat x_4\rangle(\phi,\Omega)\\&+W^{13}_{\nu\nu'}(\phi,\Omega)\langle\hat x_2\hat x_4\rangle(\phi,\Omega)\\&-W^{14}_{\nu\nu'}(\phi,\Omega)\langle\hat x_2\hat x_3\rangle(\phi,\Omega)\\&-W^{23}_{\nu\nu'}(\phi,\Omega)\langle\hat x_1\hat x_4\rangle(\phi,\Omega)\\&+W^{24}_{\nu\nu'}(\phi,\Omega)\langle\hat x_1\hat x_3\rangle(\phi,\Omega)\\&-W^{34}_{\nu\nu'}(\phi,\Omega)\langle\hat x_1\hat x_2\rangle(\phi,\Omega)\Big].\end{aligned}}
$$

As for Q-type interactions, $T^2$ is the negative of the full two-index coefficient in the Pfaffian expansion; $F^2$ is unchanged:

$$
-\sum_{i_1<i_2}F^2_{i_1i_2}(\phi,\Omega)T^2_{i_1i_2}(\phi,\Omega)=-\frac12\sum_{i_1i_2}F^2_{i_1i_2}(\phi,\Omega)T^2_{i_1i_2}(\phi,\Omega).
$$

### $T^4$

$$
\Pi^{12}_{\nu\nu'}(\phi,\Omega)=\sum_\beta W^{12}_{\nu\nu'}(\phi,\Omega),\qquad \Pi^{34}_{\nu\nu'}(\phi,\Omega)=\sum_\delta W^{34}_{\nu\nu'}(\phi,\Omega).
$$

$$
[\Pi^{12}(\phi,\Omega)]^T=-\Pi^{12}(\phi,\Omega),\qquad [\Pi^{34}(\phi,\Omega)]^T=-\Pi^{34}(\phi,\Omega).
$$

$$
\boxed{\begin{aligned}T^4_{\nu_1\nu_2\nu_3\nu_4}(\phi,\Omega)={}&\Pi^{12}_{\nu_1\nu_2}(\phi,\Omega)\Pi^{34}_{\nu_3\nu_4}(\phi,\Omega)\\&-\Pi^{12}_{\nu_1\nu_3}(\phi,\Omega)\Pi^{34}_{\nu_2\nu_4}(\phi,\Omega)\\&+\Pi^{12}_{\nu_1\nu_4}(\phi,\Omega)\Pi^{34}_{\nu_2\nu_3}(\phi,\Omega)\\&+\Pi^{12}_{\nu_2\nu_3}(\phi,\Omega)\Pi^{34}_{\nu_1\nu_4}(\phi,\Omega)\\&-\Pi^{12}_{\nu_2\nu_4}(\phi,\Omega)\Pi^{34}_{\nu_1\nu_3}(\phi,\Omega)\\&+\Pi^{12}_{\nu_3\nu_4}(\phi,\Omega)\Pi^{34}_{\nu_1\nu_2}(\phi,\Omega).\end{aligned}}
$$

$T^0,T^2,T^4$ exclude $g_P/4$, configuration minors, reference overlap, and projection weights.

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
F^2_{ab}(\phi,\Omega)=(-1)^{a+b}\operatorname{pf}S_{\widehat{ab}}(\phi,\Omega),\qquad a<b.
$$

$$
F^4_{abcd}(\phi,\Omega)=(-1)^{a+b+c+d}\operatorname{pf}S_{\widehat{abcd}}(\phi,\Omega),\qquad a<b<c<d.
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
\boxed{\mathcal K_Q^{(2)}(\phi,\Omega)=g_Q\langle\Phi_1|\Phi_2(\phi,\Omega)\rangle\left[F^0(\phi,\Omega)T_Q^0(\phi,\Omega)-\frac12\sum_{ab}F^2_{ab}(\phi,\Omega)T^2_{Q;ab}(\phi,\Omega)+\frac14\sum_{abcd}F^4_{abcd}(\phi,\Omega)\Pi^{14}_{ab}(\phi,\Omega)\Pi^{23}_{cd}(\phi,\Omega)\right].}
$$

Add the one-body kernel defined above:

$$
\boxed{\mathcal K_Q(\phi,\Omega)=\mathcal K_Q^{(1)}(\phi,\Omega)+\mathcal K_Q^{(2)}(\phi,\Omega).}
$$

### P-Type Kernel

$$
\boxed{\mathcal K_P(\phi,\Omega)=\frac{g_P}{4}\langle\Phi_1|\Phi_2(\phi,\Omega)\rangle\left[F^0(\phi,\Omega)T_P^0(\phi,\Omega)-\frac12\sum_{ab}F^2_{ab}(\phi,\Omega)T^2_{P;ab}(\phi,\Omega)+\frac14\sum_{abcd}F^4_{abcd}(\phi,\Omega)\Pi^{12}_{ab}(\phi,\Omega)\Pi^{34}_{cd}(\phi,\Omega)\right].}
$$

The outer $1/4$ comes from the pair-operator definitions; the inner $1/4$ comes from antisymmetric contraction.

At fixed configuration pair and $(\phi,\Omega)$, the same $F^0,F^2,F^4$ serve every Q-type and P-type operator. For odd $r_1+r_2$, the kernels vanish.
