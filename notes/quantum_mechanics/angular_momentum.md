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

### Matrix Representation

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

For $j=\frac12$, using the basis $\mathcal B_{1/2}=\{|1/2,1/2\rangle,|1/2,-1/2\rangle\}$,

$$
\hat J_x=\frac{\hbar}{2}\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad \hat J_y=\frac{\hbar}{2}\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad \hat J_z=\frac{\hbar}{2}\begin{pmatrix}1&0\\0&-1\end{pmatrix},
$$

$$
\hat J^2=\frac{3}{4}\hbar^2\hat I_2,
$$

where $\hat I_2$ is the two-dimensional identity operator. The component matrices are $\hbar/2$ times the Pauli matrices.
