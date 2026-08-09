# Spherical Harmonics

The spherical-coordinate conventions, $\Omega=(\theta,\phi)$, and orbital angular-momentum operators used below are defined in [Spherical Coordinates](spherical_coordinates.md).

### Scalar Spherical Harmonics

Let $l=0,1,\ldots$, $m=-l,\ldots,l$, and let $P_l(x)$ and $P_l^m(x)$ denote the Legendre polynomial and associated Legendre function. In the Condon–Shortley convention, for $m\geq0$,

$$
\boxed{Y_{lm}(\theta,\phi)=\sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}\,P_l^m(\cos\theta)e^{im\phi}}.
$$

The spherical harmonics satisfy:

- Eigenvalue equations

  $$
  \hat{\mathbf L}^2Y_{lm}=\hbar^2l(l+1)Y_{lm},\qquad\hat L_zY_{lm}=\hbar mY_{lm}.
  $$

- Conjugation and parity

  $$
  Y_{l,-m}(\Omega)=(-1)^mY_{lm}^*(\Omega),\qquad Y_{lm}(-\hat{\mathbf r})=(-1)^lY_{lm}(\hat{\mathbf r}).
  $$

- Orthonormality and completeness

  $$
  \int d\Omega\,Y_{lm}^*(\Omega)Y_{l'm'}(\Omega)=\delta_{ll'}\delta_{mm'},\qquad\sum_{l=0}^{\infty}\sum_{m=-l}^{l}Y_{lm}(\Omega)Y_{lm}^*(\Omega')=\delta(\Omega-\Omega').
  $$

- Ladder relations, with $\hat L_\pm=\hat L_x\pm i\hat L_y$

  $$
  \hat L_\pm Y_{lm}=\hbar\sqrt{l(l+1)-m(m\pm1)}\,Y_{l,m\pm1}.
  $$

- Addition theorem, where $\cos\gamma=\hat{\mathbf r}\cdot\hat{\mathbf r}'$

  $$
  \sum_{m=-l}^{l}Y_{lm}(\Omega)Y_{lm}^*(\Omega')=\frac{2l+1}{4\pi}P_l(\cos\gamma).
  $$

- Gaunt integral

  $$
  \int d\Omega\,Y_{l_1m_1}Y_{l_2m_2}Y_{l_3m_3}=\sqrt{\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi}}\begin{pmatrix}l_1&l_2&l_3\\0&0&0\end{pmatrix}\begin{pmatrix}l_1&l_2&l_3\\m_1&m_2&m_3\end{pmatrix}.
  $$

  The matrices are Wigner $3j$ symbols. The integral vanishes unless $m_1+m_2+m_3=0$, $|l_1-l_2|\leq l_3\leq l_1+l_2$, and $l_1+l_2+l_3$ is even.

A scalar function can be expanded as

$$
f(r,\Omega)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\frac{u_{lm}(r)}{r}Y_{lm}(\Omega),\qquad u_{lm}(r)=r\int d\Omega\,Y_{lm}^*(\Omega)f(r,\Omega).
$$

### Vector Spherical Harmonics

Let $\mathbf e_0$ and $\mathbf e_\pm$ denote the spherical vector basis,

$$
\mathbf e_0=\hat{\mathbf z},\qquad \mathbf e_\pm=\mp\frac{\hat{\mathbf x}\pm i\hat{\mathbf y}}{\sqrt2}.
$$

Coupling $Y_{lm}(\Omega)$ to this spin-$1$ basis defines

$$
\boxed{\mathbf Y_{JlM}(\Omega)=\sum_{m=-l}^{l}\sum_{\mu=-1}^{1}\langle l,m;1,\mu|J,M\rangle Y_{lm}(\Omega)\mathbf e_\mu},
$$

where $\langle l,m;1,\mu|J,M\rangle$ is a Clebsch–Gordan coefficient, and $\mathbf e_\mu$ denotes $\mathbf e_0$, $\mathbf e_+$, or $\mathbf e_-$ according to $\mu=0,+1,-1$:

$$
J=|l-1|,\ldots,l+1,\qquad M=-J,\ldots,J.
$$

Let $\hat{\mathbf S}$ act on the vector index and define $\hat{\mathbf J}=\hat{\mathbf L}+\hat{\mathbf S}$. The vector spherical harmonics satisfy:

- Eigenvalue equations

  $$
  \hat{\mathbf J}^2\mathbf Y_{JlM}=\hbar^2J(J+1)\mathbf Y_{JlM},\qquad \hat J_z\mathbf Y_{JlM}=\hbar M\mathbf Y_{JlM}.
  $$

  $$
  \hat{\mathbf L}^2\mathbf Y_{JlM}=\hbar^2l(l+1)\mathbf Y_{JlM},\qquad \hat{\mathbf S}^2\mathbf Y_{JlM}=2\hbar^2\mathbf Y_{JlM}.
  $$

- Conjugation and parity

  $$
  \mathbf Y_{JlM}^*(\Omega)=(-1)^{J+l+M+1}\mathbf Y_{Jl,-M}(\Omega),\qquad \mathbf Y_{JlM}(-\hat{\mathbf r})=(-1)^l\mathbf Y_{JlM}(\hat{\mathbf r}).
  $$

  The second relation transforms the argument while keeping the Cartesian basis fixed. A polar vector acquires one additional minus sign under the parity operator.

- Orthonormality and completeness

  $$
  \int d\Omega\,\mathbf Y_{J'l'M'}^*(\Omega)\cdot\mathbf Y_{JlM}(\Omega)=\delta_{JJ'}\delta_{ll'}\delta_{MM'}.
  $$

  $$
  \sum_{l=0}^{\infty}\sum_{J=|l-1|}^{l+1}\sum_{M=-J}^{J}\left[\mathbf Y_{JlM}(\Omega)\right]_a\left[\mathbf Y_{JlM}(\Omega')\right]_b^*=\delta_{ab}\delta(\Omega-\Omega').
  $$

- Ladder relations, with $\hat J_\pm=\hat J_x\pm i\hat J_y$

  $$
  \hat J_\pm\mathbf Y_{JlM}=\hbar\sqrt{J(J+1)-M(M\pm1)}\,\mathbf Y_{Jl,M\pm1}.
  $$

- For $J\geq1$, an alternative basis adapted to angular differential operators is

  $$
  \mathbf Y_{JM}^{(r)}=\hat{\mathbf r}Y_{JM},\qquad \mathbf\Psi_{JM}=\frac{\nabla_\Omega Y_{JM}}{\sqrt{J(J+1)}},\qquad \mathbf\Phi_{JM}=\frac{\hat{\mathbf r}\times\nabla_\Omega Y_{JM}}{\sqrt{J(J+1)}}.
  $$

For fixed $(J,M)$, the coupled and differential basis expansions are

$$
\mathbf F_{JM}(r,\Omega)=\sum_{l=|J-1|}^{J+1}F_{JlM}(r)\mathbf Y_{JlM}(\Omega).
$$

$$
\mathbf F_{JM}(r,\Omega)=F_{JM}^{(r)}(r)\mathbf Y_{JM}^{(r)}(\Omega)+F_{JM}^{(\Psi)}(r)\mathbf\Psi_{JM}(\Omega)+F_{JM}^{(\Phi)}(r)\mathbf\Phi_{JM}(\Omega).
$$

For $J\geq1$, the two bases span the same three-dimensional space and are related by a unitary change of basis. For $J=0$, only $\mathbf Y_{00}^{(r)}$ remains.

### Spinor Spherical Harmonics

Let $\hat I_2$ be the two-dimensional identity operator. The Pauli matrices are

$$
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad \sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad \sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$

For $i,j,k\in\{x,y,z\}$,

$$
\sigma_i\sigma_j=\delta_{ij}\hat I_2+i\sum_k\epsilon_{ijk}\sigma_k,
$$

where $\delta_{ij}$ and $\epsilon_{ijk}$ are the Kronecker delta and Levi–Civita symbol. The spin-$\frac12$ operators are

$$
\boxed{\hat S_i=\frac{\hbar}{2}\sigma_i,\qquad \hat{\mathbf S}^2=\frac{3}{4}\hbar^2\hat I_2},
$$

and satisfy

$$
[\hat S_i,\hat S_j]=i\hbar\sum_k\epsilon_{ijk}\hat S_k.
$$

Let $\chi_\mu$ denote the $\hat S_z$ eigenstate with $\mu=\pm\frac12$:

$$
\chi_{1/2}=\begin{pmatrix}1\\0\end{pmatrix},\qquad \chi_{-1/2}=\begin{pmatrix}0\\1\end{pmatrix},\qquad \hat S_z\chi_\mu=\hbar\mu\chi_\mu.
$$

On two-component functions, define the total angular momentum by

$$
\hat{\mathbf J}=\hat{\mathbf L}\otimes\hat I_2+\hat{\mathbf S}.
$$

Coupling $Y_{lm}(\Omega)$ to this basis defines

$$
\boxed{\mathcal Y_{JlM}(\Omega)=\sum_{m=-l}^{l}\sum_{\mu=\pm1/2}\left\langle l,m;\frac12,\mu\middle|J,M\right\rangle Y_{lm}(\Omega)\chi_\mu}.
$$

Here $J$ and $M$ are the total angular-momentum quantum numbers:

$$
J=l\pm\frac12,\qquad M=-J,-J+1,\ldots,J,
$$

with only $J=\frac12$ allowed for $l=0$. Since $M=m+\mu$, the two components are explicitly

$$
\mathcal Y_{JlM}=\begin{pmatrix}\left\langle l,M-\frac12;\frac12,\frac12\middle|J,M\right\rangle Y_{l,M-1/2}\\\left\langle l,M+\frac12;\frac12,-\frac12\middle|J,M\right\rangle Y_{l,M+1/2}\end{pmatrix}.
$$

In the Condon–Shortley convention, for $J=l+\frac12$,

$$
\mathcal Y_{l+1/2,l,M}=\begin{pmatrix}\sqrt{\frac{l+M+1/2}{2l+1}}\,Y_{l,M-1/2}\\\sqrt{\frac{l-M+1/2}{2l+1}}\,Y_{l,M+1/2}\end{pmatrix}.
$$

For $J=l-\frac12$ with $l\geq1$,

$$
\mathcal Y_{l-1/2,l,M}=\begin{pmatrix}-\sqrt{\frac{l-M+1/2}{2l+1}}\,Y_{l,M-1/2}\\\sqrt{\frac{l+M+1/2}{2l+1}}\,Y_{l,M+1/2}\end{pmatrix}.
$$

The spinor spherical harmonics satisfy:

- Eigenvalue equations

  $$
  \hat{\mathbf J}^2\mathcal Y_{JlM}=\hbar^2J(J+1)\mathcal Y_{JlM},\qquad \hat J_z\mathcal Y_{JlM}=\hbar M\mathcal Y_{JlM}.
  $$

  $$
  \hat{\mathbf L}^2\mathcal Y_{JlM}=\hbar^2l(l+1)\mathcal Y_{JlM},\qquad \hat{\mathbf S}^2\mathcal Y_{JlM}=\frac34\hbar^2\mathcal Y_{JlM}.
  $$

- Complex conjugation and parity

  $$
  \mathcal Y_{JlM}^*(\Omega)=(-1)^{l+J-M}i\sigma_y\mathcal Y_{Jl,-M}(\Omega),\qquad \mathcal Y_{JlM}(-\hat{\mathbf r})=(-1)^l\mathcal Y_{JlM}(\hat{\mathbf r}).
  $$

- Orthonormality and completeness

  $$
  \int d\Omega\,\mathcal Y_{J'l'M'}^\dagger(\Omega)\mathcal Y_{JlM}(\Omega)=\delta_{JJ'}\delta_{ll'}\delta_{MM'}.
  $$

  $$
  \sum_{l=0}^{\infty}\sum_{J=|l-1/2|}^{l+1/2}\sum_{M=-J}^{J}\mathcal Y_{JlM}(\Omega)\mathcal Y_{JlM}^\dagger(\Omega')=\hat I_2\delta(\Omega-\Omega').
  $$

- Ladder relations

  $$
  \hat J_\pm\mathcal Y_{JlM}=\hbar\sqrt{J(J+1)-M(M\pm1)}\,\mathcal Y_{Jl,M\pm1}.
  $$

A two-component field can be expanded as

$$
\psi(r,\Omega)=\sum_{l=0}^{\infty}\sum_{J=|l-1/2|}^{l+1/2}\sum_{M=-J}^{J}\frac{u_{JlM}(r)}{r}\mathcal Y_{JlM}(\Omega),\qquad u_{JlM}(r)=r\int d\Omega\,\mathcal Y_{JlM}^\dagger(\Omega)\psi(r,\Omega).
$$
