# Angular Momentum

### Angular Momentum Algebra

Let $\hat{\mathbf J}=(\hat J_x,\hat J_y,\hat J_z)$ have Hermitian components. Angular momentum is characterized by

$$
\boxed{[\hat J_i,\hat J_j]=i\hbar\sum_k\epsilon_{ijk}\hat J_k},\qquad i,j,k\in\{x,y,z\},
$$

where $\epsilon_{ijk}$ is the Levi–Civita symbol. Define the Casimir operator

$$
\hat J^2=\hat J_x^2+\hat J_y^2+\hat J_z^2,\qquad [\hat J^2,\hat J_i]=0.
$$

The ladder operators are

$$
\hat J_\pm=\hat J_x\pm i\hat J_y,\qquad \hat J_\pm^\dagger=\hat J_\mp,
$$

with commutation relations

$$
[\hat J_z,\hat J_\pm]=\pm\hbar\hat J_\pm,\qquad [\hat J_+,\hat J_-]=2\hbar\hat J_z,\qquad [\hat J^2,\hat J_\pm]=0.
$$

Consequently,

$$
\hat J^2=\hat J_z^2+\frac{1}{2}\left(\hat J_+\hat J_-+\hat J_-\hat J_+\right),
$$

$$
\hat J_+\hat J_-=\hat J^2-\hat J_z^2+\hbar\hat J_z,\qquad \hat J_-\hat J_+=\hat J^2-\hat J_z^2-\hbar\hat J_z.
$$

### Angular Momentum Eigenstates

Since $[\hat J^2,\hat J_z]=0$, choose normalized common eigenstates $|\lambda,m\rangle$ with dimensionless eigenvalues $\lambda$ and $m$:

$$
\hat J^2|\lambda,m\rangle=\hbar^2\lambda|\lambda,m\rangle,\qquad \hat J_z|\lambda,m\rangle=\hbar m|\lambda,m\rangle.
$$

The ladder operators preserve $\lambda$:

$$
\hat J^2\hat J_\pm|\lambda,m\rangle=\hbar^2\lambda\hat J_\pm|\lambda,m\rangle.
$$

For the raising operator,

$$
\hat J_z\hat J_+|\lambda,m\rangle=\hbar(m+1)\hat J_+|\lambda,m\rangle.
$$

$$
\left\|\hat J_+|\lambda,m\rangle\right\|^2=\langle\lambda,m|\hat J_-\hat J_+|\lambda,m\rangle=\hbar^2\left[\lambda-m(m+1)\right].
$$

Similarly, for the lowering operator,

$$
\hat J_z\hat J_-|\lambda,m\rangle=\hbar(m-1)\hat J_-|\lambda,m\rangle.
$$

$$
\left\|\hat J_-|\lambda,m\rangle\right\|^2=\langle\lambda,m|\hat J_+\hat J_-|\lambda,m\rangle=\hbar^2\left[\lambda-m(m-1)\right].
$$

At the upper and lower endpoints,

$$
\hat J_+|\lambda,m_{\max}\rangle=0,\qquad \lambda=m_{\max}(m_{\max}+1),
$$

$$
\hat J_-|\lambda,m_{\min}\rangle=0,\qquad \lambda=m_{\min}(m_{\min}-1).
$$

Writing $m_{\max}=j$ gives $m_{\min}=-j$ and $2j\in\mathbb N_0$. Therefore,

$$
\boxed{\hat J^2|j,m\rangle=\hbar^2j(j+1)|j,m\rangle,\qquad \hat J_z|j,m\rangle=\hbar m|j,m\rangle},
$$

where

$$
j=0,\frac12,1,\frac32,\ldots,\qquad m=-j,-j+1,\ldots,j.
$$

With the conventional phase choice,

$$
\boxed{\hat J_\pm|j,m\rangle=\hbar\sqrt{j(j+1)-m(m\pm1)}\,|j,m\pm1\rangle}.
$$

Each $j$ multiplet has dimension $2j+1$.

### Orbital Angular Momentum in Coordinate Space

For $\psi(\mathbf r)=\langle\mathbf r|\psi\rangle$,

$$
\hat{\mathbf p}=-i\hbar\nabla,\qquad \hat{\mathbf L}=\mathbf r\times\hat{\mathbf p}=-i\hbar\,\mathbf r\times\nabla.
$$

Let $r=|\mathbf r|$, $\hat{\mathbf r}=\mathbf r/r$, and $\Omega=(\theta,\phi)$. In spherical coordinates,

$$
\boxed{\hat L_z=-i\hbar\frac{\partial}{\partial\phi}}.
$$

The gradient is

$$
\boxed{\nabla=\hat{\mathbf r}\frac{\partial}{\partial r}-\frac{i}{\hbar r}\hat{\mathbf r}\times\hat{\mathbf L}}.
$$

The Laplacian is

$$
\boxed{\nabla^2=\frac{1}{r^2}\left[\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)-\frac{\hat{\mathbf L}^2}{\hbar^2}\right]}.
$$

The angular coordinate representation of $|l,m\rangle$ is

$$
\langle\Omega|l,m\rangle=Y_{lm}(\Omega).
$$

Therefore,

$$
\boxed{\hat{\mathbf L}^2Y_{lm}=\hbar^2l(l+1)Y_{lm},\qquad \hat L_zY_{lm}=\hbar mY_{lm}}.
$$

A wavefunction with fixed $(l,m)$ has the form

$$
\psi_{lm}(\mathbf r)=R_l(r)Y_{lm}(\Omega).
$$

The differential identities and properties of $Y_{lm}$ are collected in [Spherical Coordinates](../mathematical_physics/coordinate_spherical.md#differential-operators) and [Spherical Harmonics](../mathematical_physics/specialfunction_spherical_harmonics.md#scalar-spherical-harmonics).

### Angular Momentum in Matrix Form

For fixed $j$, use the ordered basis

$$
\mathcal B_j=\left\{|j,j\rangle,|j,j-1\rangle,\ldots,|j,-j\rangle\right\},
$$

with dimension $2j+1$. Define $(\hat J_i)_{mm'}=\langle j,m|\hat J_i|j,m'\rangle$, where $\delta_{ab}$ denotes the Kronecker delta. Then

$$
(\hat J^2)_{mm'}=\hbar^2j(j+1)\delta_{mm'},\qquad (\hat J_z)_{mm'}=\hbar m'\delta_{mm'},
$$

$$
\boxed{(\hat J_\pm)_{mm'}=\hbar\sqrt{j(j+1)-m'(m'\pm1)}\,\delta_{m,m'\pm1}}.
$$

The remaining components follow from

$$
\hat J_x=\frac{\hat J_++\hat J_-}{2},\qquad \hat J_y=\frac{\hat J_+-\hat J_-}{2i}.
$$

### Spin-$1/2$ in Matrix Form

##### Pauli Representation

In the ordered basis $\mathcal B_{1/2}$, let $\hat I_2$ denote the two-dimensional identity operator. The Pauli matrices are

$$
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad \sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad \sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$

The $\sigma_z$ eigenvectors are

$$
\chi_+=\begin{pmatrix}1\\0\end{pmatrix},\qquad \chi_-=\begin{pmatrix}0\\1\end{pmatrix},\qquad \sigma_z\chi_\pm=\pm\chi_\pm.
$$

The spin-$\frac12$ operator is

$$
\boxed{\hat{\mathbf S}=\frac{\hbar}{2}\boldsymbol{\sigma}},\qquad \boldsymbol{\sigma}=(\sigma_x,\sigma_y,\sigma_z),\qquad \hat S^2=\frac34\hbar^2\hat I_2.
$$

The Pauli matrices satisfy

$$
\sigma_i\sigma_j=\delta_{ij}\hat I_2+i\sum_k\epsilon_{ijk}\sigma_k.
$$

For ordinary vectors $\mathbf a$ and $\mathbf b$,

$$
\boxed{(\mathbf a\cdot\boldsymbol{\sigma})(\mathbf b\cdot\boldsymbol{\sigma})=(\mathbf a\cdot\mathbf b)\hat I_2+i\boldsymbol{\sigma}\cdot(\mathbf a\times\mathbf b)}.
$$

##### Pauli Operators in Local Frames

For a local orthonormal basis $\{\hat{\mathbf e}_a\}$, define

$$
\boxed{\sigma_a=\boldsymbol{\sigma}\cdot\hat{\mathbf e}_a}.
$$

These are spatial projections of the Cartesian Pauli vector, not new spin-space matrices.

For the local basis of [Spherical Coordinates](../mathematical_physics/coordinate_spherical.md#coordinate-differentials),

$$
\sigma_r=\begin{pmatrix}\cos\theta&e^{-i\phi}\sin\theta\\e^{i\phi}\sin\theta&-\cos\theta\end{pmatrix},\qquad \sigma_\theta=\begin{pmatrix}-\sin\theta&e^{-i\phi}\cos\theta\\e^{i\phi}\cos\theta&\sin\theta\end{pmatrix},\qquad \sigma_\phi=\begin{pmatrix}0&-ie^{-i\phi}\\ie^{i\phi}&0\end{pmatrix}.
$$

They satisfy

$$
\sigma_r\sigma_\theta=i\sigma_\phi,\qquad \sigma_\theta\sigma_\phi=i\sigma_r,\qquad \sigma_\phi\sigma_r=i\sigma_\theta.
$$

For the local basis of [Cylindrical Coordinates](../mathematical_physics/coordinate_cylindrical.md#coordinate-differentials),

$$
\sigma_{r_\perp}=\begin{pmatrix}0&e^{-i\phi}\\e^{i\phi}&0\end{pmatrix},\qquad \sigma_\phi=\begin{pmatrix}0&-ie^{-i\phi}\\ie^{i\phi}&0\end{pmatrix},\qquad \sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$

They satisfy

$$
\sigma_{r_\perp}\sigma_\phi=i\sigma_z,\qquad \sigma_\phi\sigma_z=i\sigma_{r_\perp},\qquad \sigma_z\sigma_{r_\perp}=i\sigma_\phi.
$$

### Angular Momentum Coupling

Let $\hat{\mathbf J}_1$ and $\hat{\mathbf J}_2$ be independent angular momenta. Define the total angular momentum on $\mathcal H=\mathcal H_{j_1}\otimes\mathcal H_{j_2}$ by

$$
[\hat J_{1\alpha},\hat J_{2\beta}]=0,\qquad \hat{\mathbf J}_3=\hat{\mathbf J}_1+\hat{\mathbf J}_2,\qquad \alpha,\beta\in\{x,y,z\}.
$$

The uncoupled basis is

$$
|j_1m_1;j_2m_2\rangle\equiv|j_1m_1\rangle\otimes|j_2m_2\rangle.
$$

The three angular momenta satisfy

$$
\hat J_1^2|j_1m_1;j_2m_2\rangle=\hbar^2j_1(j_1+1)|j_1m_1;j_2m_2\rangle,\qquad \hat J_{1z}|j_1m_1;j_2m_2\rangle=\hbar m_1|j_1m_1;j_2m_2\rangle.
$$

$$
\hat J_2^2|j_1m_1;j_2m_2\rangle=\hbar^2j_2(j_2+1)|j_1m_1;j_2m_2\rangle,\qquad \hat J_{2z}|j_1m_1;j_2m_2\rangle=\hbar m_2|j_1m_1;j_2m_2\rangle.
$$

$$
\hat J_3^2|(j_1j_2)j_3m_3\rangle=\hbar^2j_3(j_3+1)|(j_1j_2)j_3m_3\rangle,\qquad \hat J_{3z}|(j_1j_2)j_3m_3\rangle=\hbar m_3|(j_1j_2)j_3m_3\rangle.
$$

The allowed coupled quantum numbers are

$$
j_3=|j_1-j_2|,|j_1-j_2|+1,\ldots,j_1+j_2,\qquad m_3=-j_3,-j_3+1,\ldots,j_3,
$$

and a nonzero coupling coefficient requires $m_3=m_1+m_2$. The dimensions satisfy

$$
(2j_1+1)(2j_2+1)=\sum_{j_3=|j_1-j_2|}^{j_1+j_2}(2j_3+1).
$$

The uncoupled and coupled bases are related by

$$
\boxed{|(j_1j_2)j_3m_3\rangle=\sum_{m_1,m_2}|j_1m_1;j_2m_2\rangle\langle j_1m_1;j_2m_2|j_3m_3\rangle}.
$$

The inverse transformation is

$$
|j_1m_1;j_2m_2\rangle=\sum_{j_3,m_3}|(j_1j_2)j_3m_3\rangle\langle j_3m_3|j_1m_1;j_2m_2\rangle.
$$

The Condon–Shortley convention and the orthogonality relations of the Clebsch–Gordan coefficients are collected in [Wigner Symbols](../mathematical_physics/wigner_symbols.md).

Since

$$
\hat J_3^2=\hat J_1^2+\hat J_2^2+2\hat{\mathbf J}_1\cdot\hat{\mathbf J}_2,
$$

the scalar coupling operator is

$$
\boxed{\hat{\mathbf J}_1\cdot\hat{\mathbf J}_2=\frac12\left(\hat J_3^2-\hat J_1^2-\hat J_2^2\right)}.
$$

Therefore,

$$
\hat{\mathbf J}_1\cdot\hat{\mathbf J}_2|(j_1j_2)j_3m_3\rangle=\frac{\hbar^2}{2}\left[j_3(j_3+1)-j_1(j_1+1)-j_2(j_2+1)\right]|(j_1j_2)j_3m_3\rangle.
$$

The specialization $\hat{\mathbf J}_1=\hat{\mathbf L}$, $\hat{\mathbf J}_2=\hat{\mathbf S}$, and $\hat{\mathbf J}_3=\hat{\mathbf J}$ gives the algebraic basis of spin–orbit coupling.
