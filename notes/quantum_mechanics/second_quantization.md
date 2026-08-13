# Second Quantization

### Bosonic Second Quantization

##### Single-Particle Modes and Field Operators

Notation:

- $i,j=1,\ldots,N$: particle labels.
- $\alpha,\beta=1,\ldots,M$: single-particle state labels.

Let $\hat h_0$ be a reference single-particle Hamiltonian with a complete orthonormal eigenbasis:

$$
\hat h_0|\phi_\alpha\rangle=\varepsilon_\alpha|\phi_\alpha\rangle,\qquad \langle\phi_\alpha|\phi_\beta\rangle=\delta_{\alpha\beta},\qquad \sum_{\alpha=1}^{M}|\phi_\alpha\rangle\langle\phi_\alpha|=\hat I_1.
$$

The coordinate representation of $|\phi_\alpha\rangle$ is $\phi_\alpha(x)=\langle x|\phi_\alpha\rangle$. In the one-particle sector, the mode and coordinate states are created from the vacuum by

$$
\hat a_\alpha^\dagger|0\rangle=|\phi_\alpha\rangle,\qquad \hat\psi^\dagger(x)|0\rangle=\sum_{\beta=1}^{M}\phi_\beta^*(x)\hat a_\beta^\dagger|0\rangle=\sum_{\beta=1}^{M}|\phi_\beta\rangle\langle\phi_\beta|x\rangle=|x\rangle.
$$

Hence,

$$
\langle0|\hat\psi(x)\hat a_\alpha^\dagger|0\rangle=\phi_\alpha(x).
$$

The bosonic field operators satisfy the equal-time commutation relations

$$
\boxed{[\hat\psi(x),\hat\psi^\dagger(x')]=\delta(x-x'),\qquad [\hat\psi(x),\hat\psi(x')]=[\hat\psi^\dagger(x),\hat\psi^\dagger(x')]=0}.
$$

Expand the field operators in the single-particle basis:

$$
\boxed{\hat\psi(x)=\sum_{\alpha=1}^{M}\phi_\alpha(x)\hat a_\alpha,\qquad \hat\psi^\dagger(x)=\sum_{\alpha=1}^{M}\phi_\alpha^*(x)\hat a_\alpha^\dagger}.
$$

Orthonormality gives

$$
\hat a_\alpha=\int dx\,\phi_\alpha^*(x)\hat\psi(x),\qquad \hat a_\alpha^\dagger=\int dx\,\phi_\alpha(x)\hat\psi^\dagger(x).
$$

Projecting the field commutator onto modes $\alpha$ and $\beta$ gives

$$
[\hat a_\alpha,\hat a_\beta^\dagger]=\int dx\,dx'\,\phi_\alpha^*(x)\phi_\beta(x')[\hat\psi(x),\hat\psi^\dagger(x')]=\delta_{\alpha\beta}.
$$

Therefore,

$$
\boxed{[\hat a_\alpha,\hat a_\beta^\dagger]=\delta_{\alpha\beta},\qquad [\hat a_\alpha,\hat a_\beta]=[\hat a_\alpha^\dagger,\hat a_\beta^\dagger]=0}.
$$

With $h_0(x,x')=\langle x|\hat h_0|x'\rangle$, the independent-particle Hamiltonian is

$$
\hat H_0=\sum_{i=1}^{N}\hat h_0(i)=\int dx\,dx'\,\hat\psi^\dagger(x)h_0(x,x')\hat\psi(x')=\sum_{\alpha=1}^{M}\varepsilon_\alpha\hat a_\alpha^\dagger\hat a_\alpha.
$$

##### Exchange Symmetry and Fock Space

Let $\mathcal H_1$ denote the single-particle Hilbert space. The state space of $N$ labeled particles is $\mathcal H_1^{\otimes N}$. For mode labels $\alpha_i$ assigned to particles $i$, define the ordered product state

$$
|\Phi(\alpha_1,\ldots,\alpha_N)\rangle=|\phi_{\alpha_1}(1)\rangle\otimes\cdots\otimes|\phi_{\alpha_N}(N)\rangle.
$$

These product states span $\mathcal H_1^{\otimes N}$. For a permutation $\sigma\in S_N$,

$$
\hat P_\sigma|\Phi(\alpha_1,\ldots,\alpha_N)\rangle=|\Phi(\alpha_{\sigma^{-1}(1)},\ldots,\alpha_{\sigma^{-1}(N)})\rangle.
$$

Bosonic states satisfy

$$
\boxed{\hat P_\sigma|\Psi_B\rangle=|\Psi_B\rangle}.
$$

Define the symmetrization projector by

$$
\hat{\mathcal S}_N=\frac{1}{N!}\sum_{\sigma\in S_N}\hat P_\sigma,\qquad \hat{\mathcal S}_N^\dagger=\hat{\mathcal S}_N,\qquad \hat{\mathcal S}_N^2=\hat{\mathcal S}_N.
$$

The normalized symmetric product state is

$$
|\Phi_B(\alpha_1,\ldots,\alpha_N)\rangle=\frac{\hat{\mathcal S}_N|\Phi(\alpha_1,\ldots,\alpha_N)\rangle}{\sqrt{\langle\Phi(\alpha_1,\ldots,\alpha_N)|\hat{\mathcal S}_N|\Phi(\alpha_1,\ldots,\alpha_N)\rangle}}.
$$

Since the creation operators commute, their order is irrelevant:

$$
[\hat a_{\alpha_i}^\dagger,\hat a_{\alpha_j}^\dagger]=0,\qquad |\Phi_B(\alpha_1,\ldots,\alpha_N)\rangle\propto\hat a_{\alpha_1}^\dagger\cdots\hat a_{\alpha_N}^\dagger|0\rangle.
$$

These states span the fixed-$N$ bosonic sector. Hence,

$$
|\Psi_B\rangle=\sum_{\alpha_1\leq\cdots\leq\alpha_N}C_{\alpha_1\cdots\alpha_N}|\Phi_B(\alpha_1,\ldots,\alpha_N)\rangle.
$$

For $N$ distinct orthonormal modes,

$$
\boxed{\langle x_1,\ldots,x_N|\Phi_B(\alpha_1,\ldots,\alpha_N)\rangle=\frac{1}{\sqrt{N!}}\operatorname{perm}\left[\phi_{\alpha_j}(x_i)\right]_{i,j=1}^{N}=\frac{1}{\sqrt{N!}}\sum_{\sigma\in S_N}\prod_{i=1}^{N}\phi_{\alpha_{\sigma(i)}}(x_i)}.
$$

The fixed-$N$ bosonic sector and the bosonic Fock space are

$$
\boxed{\mathcal H_B^{(N)}=\hat{\mathcal S}_N\mathcal H_1^{\otimes N},\qquad \mathcal F_B(\mathcal H_1)=\bigoplus_{N=0}^{\infty}\mathcal H_B^{(N)}}.
$$

where $\mathcal H_B^{(0)}=\mathbb C|0\rangle$.

##### Number Representation

Define the mode-number operator. The bosonic commutation relations give

$$
\hat n_\alpha=\hat a_\alpha^\dagger\hat a_\alpha,\qquad [\hat n_\alpha,\hat a_\beta^\dagger]=\delta_{\alpha\beta}\hat a_\beta^\dagger,\qquad [\hat n_\alpha,\hat a_\beta]=-\delta_{\alpha\beta}\hat a_\beta.
$$

For the single-particle state $|\phi_\beta\rangle=\hat a_\beta^\dagger|0\rangle$,

$$
\hat n_\alpha|\phi_\beta\rangle=\delta_{\alpha\beta}|\phi_\beta\rangle.
$$

Repeatedly applying the number-operator commutator gives

$$
\hat n_\alpha|\Phi_B(\alpha_1,\ldots,\alpha_N)\rangle=\left(\sum_{i=1}^{N}\delta_{\alpha\alpha_i}\right)|\Phi_B(\alpha_1,\ldots,\alpha_N)\rangle.
$$

Define the occupation numbers by

$$
n_\alpha=\sum_{i=1}^{N}\delta_{\alpha\alpha_i},\qquad \boldsymbol n=(n_1,\ldots,n_M),\qquad \sum_{\alpha=1}^{M}n_\alpha=N.
$$

$$
|\boldsymbol n\rangle_B\equiv|\Phi_B(\alpha_1,\ldots,\alpha_N)\rangle,\qquad \hat n_\alpha|\boldsymbol n\rangle_B=n_\alpha|\boldsymbol n\rangle_B.
$$

Their normalization follows from

$$
\left\|\hat a_\alpha^\dagger|\boldsymbol n\rangle_B\right\|^2=\langle\boldsymbol n|\hat a_\alpha\hat a_\alpha^\dagger|\boldsymbol n\rangle_B=n_\alpha+1,\qquad \left\|\hat a_\alpha|\boldsymbol n\rangle_B\right\|^2=\langle\boldsymbol n|\hat a_\alpha^\dagger\hat a_\alpha|\boldsymbol n\rangle_B=n_\alpha.
$$

Let $\boldsymbol e_\alpha$ increase $n_\alpha$ by one. Choosing the phases of the occupation-number states gives

$$
\boxed{\hat a_\alpha^\dagger|\boldsymbol n\rangle_B=\sqrt{n_\alpha+1}\,|\boldsymbol n+\boldsymbol e_\alpha\rangle_B,\qquad \hat a_\alpha|\boldsymbol n\rangle_B=\sqrt{n_\alpha}\,|\boldsymbol n-\boldsymbol e_\alpha\rangle_B}.
$$

Starting from the vacuum, the occupation-number states form a basis of the bosonic Fock space:

$$
\hat a_\alpha|0\rangle=0,\qquad |\boldsymbol n\rangle_B=\prod_{\alpha=1}^{M}\frac{(\hat a_\alpha^\dagger)^{n_\alpha}}{\sqrt{n_\alpha!}}|0\rangle,\qquad |\Psi_B\rangle=\sum_{\boldsymbol n}C_{\boldsymbol n}|\boldsymbol n\rangle_B.
$$

##### One- and Two-Body Operators

On the single-particle sector, $|\phi_\gamma\rangle=\hat a_\gamma^\dagger|0\rangle$, and

$$
\hat a_\alpha^\dagger\hat a_\beta|\phi_\gamma\rangle=\delta_{\beta\gamma}|\phi_\alpha\rangle=\left(|\phi_\alpha\rangle\langle\phi_\beta|\right)|\phi_\gamma\rangle.
$$

Therefore, on the single-particle and fixed-$N$ sectors,

$$
\left.\hat a_\alpha^\dagger\hat a_\beta\right|_{\mathcal H_1}=|\phi_\alpha\rangle\langle\phi_\beta|,\qquad \left.\hat a_\alpha^\dagger\hat a_\beta\right|_{\mathcal H_B^{(N)}}=\sum_{i=1}^{N}|\phi_\alpha(i)\rangle\langle\phi_\beta(i)|.
$$

For $o_{\alpha\beta}=\langle\phi_\alpha|\hat o|\phi_\beta\rangle$, inserting the single-particle completeness relation gives

$$
\hat O^{(1)}=\sum_{i=1}^{N}\hat o(i)=\sum_{i=1}^{N}\left[\sum_{\alpha,\beta=1}^{M}|\phi_\alpha(i)\rangle\langle\phi_\alpha(i)|\hat o(i)|\phi_\beta(i)\rangle\langle\phi_\beta(i)|\right]=\sum_{\alpha,\beta=1}^{M}o_{\alpha\beta}\sum_{i=1}^{N}|\phi_\alpha(i)\rangle\langle\phi_\beta(i)|=\sum_{\alpha,\beta=1}^{M}o_{\alpha\beta}\hat a_\alpha^\dagger\hat a_\beta.
$$

In coordinate space,

$$
\boxed{\hat O^{(1)}=\sum_{\alpha,\beta=1}^{M}o_{\alpha\beta}\hat a_\alpha^\dagger\hat a_\beta=\int dx\,dx'\,\hat\psi^\dagger(x)o(x,x')\hat\psi(x')}.
$$

For a symmetric two-particle operator $\hat v$, let $v_{\alpha\beta;\gamma\delta}=\langle\phi_\alpha\phi_\beta|\hat v|\phi_\gamma\phi_\delta\rangle$. Inserting the two-particle completeness relation gives

$$
\hat V^{(2)}=\frac12\sum_{i\neq j}\hat v(i,j)=\frac12\sum_{i\neq j}\left[\sum_{\alpha,\beta,\gamma,\delta=1}^{M}|\phi_\alpha(i)\rangle\langle\phi_\alpha(i)|\,|\phi_\beta(j)\rangle\langle\phi_\beta(j)|\hat v(i,j)|\phi_\gamma(i)\rangle\langle\phi_\gamma(i)|\,|\phi_\delta(j)\rangle\langle\phi_\delta(j)|\right].
$$

Using the two-particle transition operator,

$$
\hat V^{(2)}=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{M}v_{\alpha\beta;\gamma\delta}\sum_{i\neq j}|\phi_\alpha(i)\rangle\langle\phi_\gamma(i)|\,|\phi_\beta(j)\rangle\langle\phi_\delta(j)|=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{M}v_{\alpha\beta;\gamma\delta}\hat a_\alpha^\dagger\hat a_\beta^\dagger\hat a_\delta\hat a_\gamma.
$$

For a coordinate-space interaction $v(x,x')$,

$$
\boxed{\hat V^{(2)}=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{M}v_{\alpha\beta;\gamma\delta}\hat a_\alpha^\dagger\hat a_\beta^\dagger\hat a_\delta\hat a_\gamma=\frac12\int dx\,dx'\,\hat\psi^\dagger(x)\hat\psi^\dagger(x')v(x,x')\hat\psi(x')\hat\psi(x)}.
$$

For $h_{\alpha\beta}=\langle\phi_\alpha|\hat h|\phi_\beta\rangle$, the one- and two-body Hamiltonian is

$$
\boxed{\hat H=\sum_{\alpha,\beta=1}^{M}h_{\alpha\beta}\hat a_\alpha^\dagger\hat a_\beta+\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{M}v_{\alpha\beta;\gamma\delta}\hat a_\alpha^\dagger\hat a_\beta^\dagger\hat a_\delta\hat a_\gamma}.
$$

### Fermionic Second Quantization

##### Single-Particle Modes and Field Operators

Notation:

- $i,j=1,\ldots,N$: particle labels.
- $\alpha,\beta=1,\ldots,M$: single-particle state labels.

Let $\hat h_0$ be a reference single-particle Hamiltonian with a complete orthonormal eigenbasis:

$$
\hat h_0|\phi_\alpha\rangle=\varepsilon_\alpha|\phi_\alpha\rangle,\qquad \langle\phi_\alpha|\phi_\beta\rangle=\delta_{\alpha\beta},\qquad \sum_{\alpha=1}^{M}|\phi_\alpha\rangle\langle\phi_\alpha|=\hat I_1.
$$

The coordinate representation of $|\phi_\alpha\rangle$ is $\phi_\alpha(x)=\langle x|\phi_\alpha\rangle$. In the one-particle sector, the mode and coordinate states are created from the vacuum by

$$
\hat c_\alpha^\dagger|0\rangle=|\phi_\alpha\rangle,\qquad \hat\psi_F^\dagger(x)|0\rangle=\sum_{\beta=1}^{M}\phi_\beta^*(x)\hat c_\beta^\dagger|0\rangle=\sum_{\beta=1}^{M}|\phi_\beta\rangle\langle\phi_\beta|x\rangle=|x\rangle.
$$

Hence,

$$
\langle0|\hat\psi_F(x)\hat c_\alpha^\dagger|0\rangle=\phi_\alpha(x).
$$

With $\{\hat A,\hat B\}=\hat A\hat B+\hat B\hat A$, the fermionic field operators satisfy the equal-time anticommutation relations

$$
\boxed{\{\hat\psi_F(x),\hat\psi_F^\dagger(x')\}=\delta(x-x'),\qquad \{\hat\psi_F(x),\hat\psi_F(x')\}=\{\hat\psi_F^\dagger(x),\hat\psi_F^\dagger(x')\}=0}.
$$

Expand the field operators in the single-particle basis:

$$
\boxed{\hat\psi_F(x)=\sum_{\alpha=1}^{M}\phi_\alpha(x)\hat c_\alpha,\qquad \hat\psi_F^\dagger(x)=\sum_{\alpha=1}^{M}\phi_\alpha^*(x)\hat c_\alpha^\dagger}.
$$

Orthonormality gives

$$
\hat c_\alpha=\int dx\,\phi_\alpha^*(x)\hat\psi_F(x),\qquad \hat c_\alpha^\dagger=\int dx\,\phi_\alpha(x)\hat\psi_F^\dagger(x).
$$

Projecting the field anticommutator onto modes $\alpha$ and $\beta$ gives

$$
\{\hat c_\alpha,\hat c_\beta^\dagger\}=\int dx\,dx'\,\phi_\alpha^*(x)\phi_\beta(x')\{\hat\psi_F(x),\hat\psi_F^\dagger(x')\}=\delta_{\alpha\beta}.
$$

Therefore,

$$
\boxed{\{\hat c_\alpha,\hat c_\beta^\dagger\}=\delta_{\alpha\beta},\qquad \{\hat c_\alpha,\hat c_\beta\}=\{\hat c_\alpha^\dagger,\hat c_\beta^\dagger\}=0}.
$$

With $h_0(x,x')=\langle x|\hat h_0|x'\rangle$, the independent-particle Hamiltonian is

$$
\hat H_0=\sum_{i=1}^{N}\hat h_0(i)=\int dx\,dx'\,\hat\psi_F^\dagger(x)h_0(x,x')\hat\psi_F(x')=\sum_{\alpha=1}^{M}\varepsilon_\alpha\hat c_\alpha^\dagger\hat c_\alpha.
$$

##### Exchange Symmetry and Fock Space

Let $\mathcal H_1$ denote the single-particle Hilbert space. For mode labels $\alpha_i$ assigned to particles $i$, define

$$
|\Phi(\alpha_1,\ldots,\alpha_N)\rangle=|\phi_{\alpha_1}(1)\rangle\otimes\cdots\otimes|\phi_{\alpha_N}(N)\rangle.
$$

For a permutation $\sigma\in S_N$,

$$
\hat P_\sigma|\Phi(\alpha_1,\ldots,\alpha_N)\rangle=|\Phi(\alpha_{\sigma^{-1}(1)},\ldots,\alpha_{\sigma^{-1}(N)})\rangle,\qquad \boxed{\hat P_\sigma|\Psi_F\rangle=\operatorname{sgn}(\sigma)|\Psi_F\rangle}.
$$

Define the antisymmetrization projector by

$$
\hat{\mathcal A}_N=\frac{1}{N!}\sum_{\sigma\in S_N}\operatorname{sgn}(\sigma)\hat P_\sigma,\qquad \hat{\mathcal A}_N^\dagger=\hat{\mathcal A}_N,\qquad \hat{\mathcal A}_N^2=\hat{\mathcal A}_N.
$$

The normalized antisymmetric product state is

$$
|\Phi_F(\alpha_1,\ldots,\alpha_N)\rangle=\frac{\hat{\mathcal A}_N|\Phi(\alpha_1,\ldots,\alpha_N)\rangle}{\sqrt{\langle\Phi(\alpha_1,\ldots,\alpha_N)|\hat{\mathcal A}_N|\Phi(\alpha_1,\ldots,\alpha_N)\rangle}}.
$$

The anticommuting creation operators encode the exchange symmetry:

$$
\{\hat c_{\alpha_i}^\dagger,\hat c_{\alpha_j}^\dagger\}=0,\qquad |\Phi_F(\alpha_1,\ldots,\alpha_N)\rangle\propto\hat c_{\alpha_1}^\dagger\cdots\hat c_{\alpha_N}^\dagger|0\rangle,\qquad \alpha_1<\cdots<\alpha_N.
$$

If two mode labels coincide,

$$
\hat{\mathcal A}_N|\Phi(\alpha_1,\ldots,\alpha_N)\rangle=0,\qquad (\hat c_\alpha^\dagger)^2=0.
$$

The nonzero states span the fixed-$N$ fermionic sector:

$$
|\Psi_F\rangle=\sum_{\alpha_1<\cdots<\alpha_N}C_{\alpha_1\cdots\alpha_N}|\Phi_F(\alpha_1,\ldots,\alpha_N)\rangle.
$$

For distinct orthonormal modes,

$$
\boxed{\langle x_1,\ldots,x_N|\Phi_F(\alpha_1,\ldots,\alpha_N)\rangle=\frac{1}{\sqrt{N!}}\det\left[\phi_{\alpha_j}(x_i)\right]_{i,j=1}^{N}=\frac{1}{\sqrt{N!}}\sum_{\sigma\in S_N}\operatorname{sgn}(\sigma)\prod_{i=1}^{N}\phi_{\alpha_{\sigma(i)}}(x_i)}.
$$

The fixed-$N$ sector and fermionic Fock space are

$$
\boxed{\mathcal H_F^{(N)}=\hat{\mathcal A}_N\mathcal H_1^{\otimes N},\qquad \mathcal F_F(\mathcal H_1)=\bigoplus_{N=0}^{\infty}\mathcal H_F^{(N)}},\qquad \mathcal H_F^{(0)}=\mathbb C|0\rangle.
$$

##### Number Representation

Define the mode-number operator. The fermionic anticommutation relations give

$$
\hat n_\alpha=\hat c_\alpha^\dagger\hat c_\alpha,\qquad [\hat n_\alpha,\hat c_\beta^\dagger]=\delta_{\alpha\beta}\hat c_\beta^\dagger,\qquad [\hat n_\alpha,\hat c_\beta]=-\delta_{\alpha\beta}\hat c_\beta.
$$

For the single-particle state $|\phi_\beta\rangle=\hat c_\beta^\dagger|0\rangle$,

$$
\hat n_\alpha|\phi_\beta\rangle=\delta_{\alpha\beta}|\phi_\beta\rangle.
$$

Repeatedly applying the number-operator commutator gives

$$
\hat n_\alpha|\Phi_F(\alpha_1,\ldots,\alpha_N)\rangle=\left(\sum_{i=1}^{N}\delta_{\alpha\alpha_i}\right)|\Phi_F(\alpha_1,\ldots,\alpha_N)\rangle.
$$

Define the occupation numbers by

$$
n_\alpha=\sum_{i=1}^{N}\delta_{\alpha\alpha_i}\in\{0,1\},\qquad \boldsymbol n=(n_1,\ldots,n_M),\qquad \sum_{\alpha=1}^{M}n_\alpha=N.
$$

$$
|\boldsymbol n\rangle_F\equiv|\Phi_F(\alpha_1,\ldots,\alpha_N)\rangle,\qquad \hat n_\alpha|\boldsymbol n\rangle_F=n_\alpha|\boldsymbol n\rangle_F.
$$

Their normalization follows from

$$
\left\|\hat c_\alpha^\dagger|\boldsymbol n\rangle_F\right\|^2=\langle\boldsymbol n|\hat c_\alpha\hat c_\alpha^\dagger|\boldsymbol n\rangle_F=1-n_\alpha,\qquad \left\|\hat c_\alpha|\boldsymbol n\rangle_F\right\|^2=\langle\boldsymbol n|\hat c_\alpha^\dagger\hat c_\alpha|\boldsymbol n\rangle_F=n_\alpha.
$$

Let $\boldsymbol e_\alpha$ increase $n_\alpha$ by one and define

$$
\eta_\alpha(\boldsymbol n)=\sum_{\beta<\alpha}n_\beta.
$$

The factor $(-1)^{\eta_\alpha(\boldsymbol n)}$ follows from fermionic exchange antisymmetry: $\hat c_\alpha^{(\dagger)}$ must pass through $\eta_\alpha(\boldsymbol n)$ occupied modes to preserve the canonical ordering.

With the canonical mode ordering,

$$
\boxed{\hat c_\alpha^\dagger|\boldsymbol n\rangle_F=(-1)^{\eta_\alpha(\boldsymbol n)}\sqrt{1-n_\alpha}\,|\boldsymbol n+\boldsymbol e_\alpha\rangle_F,\qquad \hat c_\alpha|\boldsymbol n\rangle_F=(-1)^{\eta_\alpha(\boldsymbol n)}\sqrt{n_\alpha}\,|\boldsymbol n-\boldsymbol e_\alpha\rangle_F}.
$$

Starting from the vacuum, the occupation-number states form a basis of the fermionic Fock space:

$$
\hat c_\alpha|0\rangle=0,\qquad |\boldsymbol n\rangle_F=\prod_{\alpha=1}^{M}(\hat c_\alpha^\dagger)^{n_\alpha}|0\rangle,\qquad |\Psi_F\rangle=\sum_{\boldsymbol n}C_{\boldsymbol n}|\boldsymbol n\rangle_F,
$$

where the product is ordered by increasing $\alpha$.

##### One- and Two-Body Operators

On the single-particle sector, $|\phi_\gamma\rangle=\hat c_\gamma^\dagger|0\rangle$, and

$$
\hat c_\alpha^\dagger\hat c_\beta|\phi_\gamma\rangle=\delta_{\beta\gamma}|\phi_\alpha\rangle=\left(|\phi_\alpha\rangle\langle\phi_\beta|\right)|\phi_\gamma\rangle.
$$

Therefore, on the single-particle and fixed-$N$ sectors,

$$
\left.\hat c_\alpha^\dagger\hat c_\beta\right|_{\mathcal H_1}=|\phi_\alpha\rangle\langle\phi_\beta|,\qquad \left.\hat c_\alpha^\dagger\hat c_\beta\right|_{\mathcal H_F^{(N)}}=\sum_{i=1}^{N}|\phi_\alpha(i)\rangle\langle\phi_\beta(i)|.
$$

For $o_{\alpha\beta}=\langle\phi_\alpha|\hat o|\phi_\beta\rangle$, inserting the single-particle completeness relation gives

$$
\hat O^{(1)}=\sum_{i=1}^{N}\hat o(i)=\sum_{i=1}^{N}\left[\sum_{\alpha,\beta=1}^{M}|\phi_\alpha(i)\rangle\langle\phi_\alpha(i)|\hat o(i)|\phi_\beta(i)\rangle\langle\phi_\beta(i)|\right]=\sum_{\alpha,\beta=1}^{M}o_{\alpha\beta}\sum_{i=1}^{N}|\phi_\alpha(i)\rangle\langle\phi_\beta(i)|=\sum_{\alpha,\beta=1}^{M}o_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta.
$$

In coordinate space,

$$
\boxed{\hat O^{(1)}=\sum_{\alpha,\beta=1}^{M}o_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta=\int dx\,dx'\,\hat\psi_F^\dagger(x)o(x,x')\hat\psi_F(x')}.
$$

For a symmetric two-particle operator $\hat v$, let $v_{\alpha\beta;\gamma\delta}=\langle\phi_\alpha\phi_\beta|\hat v|\phi_\gamma\phi_\delta\rangle$. Inserting the two-particle completeness relation gives

$$
\hat V^{(2)}=\frac12\sum_{i\neq j}\hat v(i,j)=\frac12\sum_{i\neq j}\left[\sum_{\alpha,\beta,\gamma,\delta=1}^{M}|\phi_\alpha(i)\rangle\langle\phi_\alpha(i)|\,|\phi_\beta(j)\rangle\langle\phi_\beta(j)|\hat v(i,j)|\phi_\gamma(i)\rangle\langle\phi_\gamma(i)|\,|\phi_\delta(j)\rangle\langle\phi_\delta(j)|\right].
$$

Using the two-particle transition operator,

$$
\hat V^{(2)}=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{M}v_{\alpha\beta;\gamma\delta}\sum_{i\neq j}|\phi_\alpha(i)\rangle\langle\phi_\gamma(i)|\,|\phi_\beta(j)\rangle\langle\phi_\delta(j)|=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{M}v_{\alpha\beta;\gamma\delta}\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma.
$$

For a coordinate-space interaction $v(x,x')$,

$$
\boxed{\hat V^{(2)}=\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{M}v_{\alpha\beta;\gamma\delta}\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma=\frac12\int dx\,dx'\,\hat\psi_F^\dagger(x)\hat\psi_F^\dagger(x')v(x,x')\hat\psi_F(x')\hat\psi_F(x)}.
$$

For $h_{\alpha\beta}=\langle\phi_\alpha|\hat h|\phi_\beta\rangle$, the one- and two-body Hamiltonian is

$$
\boxed{\hat H=\sum_{\alpha,\beta=1}^{M}h_{\alpha\beta}\hat c_\alpha^\dagger\hat c_\beta+\frac12\sum_{\alpha,\beta,\gamma,\delta=1}^{M}v_{\alpha\beta;\gamma\delta}\hat c_\alpha^\dagger\hat c_\beta^\dagger\hat c_\delta\hat c_\gamma}.
$$
