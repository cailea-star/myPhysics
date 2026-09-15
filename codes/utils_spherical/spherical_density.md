# Spherical Densities

### Basis and Block Matrices

Following the [spherical basis](spherical_basis.md), label each positive-branch representative by

$$
\alpha=\{n_\alpha,l_\alpha,j_\alpha,m_\alpha\},\qquad m_\alpha>0.
$$

Let $\hat{\mathbf r}=\boldsymbol r/r$ denote the radial unit vector. The positive-branch basis spinor and its negative-branch partner are

$$
\phi_\alpha^+(\boldsymbol r)=\phi_{n_\alpha l_\alpha}(r)\mathcal Y^{l_\alpha}_{j_\alpha m_\alpha}(\hat{\mathbf r}),\qquad \phi_\alpha^-(\boldsymbol r)=\phi_{n_\alpha l_\alpha}(r)\mathcal Y^{l_\alpha}_{j_\alpha,-m_\alpha}(\hat{\mathbf r}).
$$

The spinor spherical harmonics are

$$
\mathcal Y^{l_\alpha}_{j_\alpha m_\alpha}(\hat{\mathbf r})=\sum_{m_s=\pm1/2}\left\langle l_\alpha,m_\alpha-m_s;\frac12,m_s\middle|j_\alpha,m_\alpha\right\rangle Y_{l_\alpha,m_\alpha-m_s}(\hat{\mathbf r})\chi_{m_s},
$$

$$
\mathcal Y^{l_\alpha}_{j_\alpha,-m_\alpha}(\hat{\mathbf r})=\sum_{m_s=\pm1/2}\left\langle l_\alpha,-m_\alpha-m_s;\frac12,m_s\middle|j_\alpha,-m_\alpha\right\rangle Y_{l_\alpha,-m_\alpha-m_s}(\hat{\mathbf r})\chi_{m_s}.
$$

Here $\chi_{m_s}$ is the spin state with projection $m_s=\pm\frac12$. Terms with orbital projection outside $[-l_\alpha,l_\alpha]$ are omitted.

Both branches use the same representative label and radial function. Their total angular-momentum projections are opposite.

Let $\hat c_\alpha^\pm$ annihilate the corresponding branch state. The fields are

$$
\hat\psi(\boldsymbol r)=\sum_\alpha\left[\hat c_\alpha^+\phi_\alpha^+(\boldsymbol r)+\hat c_\alpha^-\phi_\alpha^-(\boldsymbol r)\right],
$$

$$
\hat\psi^\dagger(\boldsymbol r)=\sum_\alpha\left[(\hat c_\alpha^+)^\dagger(\phi_\alpha^+)^\dagger(\boldsymbol r)+(\hat c_\alpha^-)^\dagger(\phi_\alpha^-)^\dagger(\boldsymbol r)\right].
$$

Order all positive-branch states before their partners. The normal density matrix is

$$
\rho_{\alpha\beta}=\begin{pmatrix}\rho_{\alpha\beta}^{++}&\rho_{\alpha\beta}^{+-}\\\rho_{\alpha\beta}^{-+}&\rho_{\alpha\beta}^{--}\end{pmatrix}=\begin{pmatrix}\langle(\hat c_\beta^+)^\dagger\hat c_\alpha^+\rangle&\langle(\hat c_\beta^-)^\dagger\hat c_\alpha^+\rangle\\\langle(\hat c_\beta^+)^\dagger\hat c_\alpha^-\rangle&\langle(\hat c_\beta^-)^\dagger\hat c_\alpha^-\rangle\end{pmatrix}.
$$

Hermiticity implies

$$
(\rho^{++})^\dagger=\rho^{++},\qquad (\rho^{--})^\dagger=\rho^{--},\qquad \rho^{-+}=(\rho^{+-})^\dagger.
$$

The pairing tensor is

$$
\kappa_{\alpha\beta}=\begin{pmatrix}\kappa_{\alpha\beta}^{++}&\kappa_{\alpha\beta}^{+-}\\\kappa_{\alpha\beta}^{-+}&\kappa_{\alpha\beta}^{--}\end{pmatrix}=\begin{pmatrix}\langle\hat c_\beta^+\hat c_\alpha^+\rangle&\langle\hat c_\beta^-\hat c_\alpha^+\rangle\\\langle\hat c_\beta^+\hat c_\alpha^-\rangle&\langle\hat c_\beta^-\hat c_\alpha^-\rangle\end{pmatrix}.
$$

Fermionic antisymmetry implies

$$
\kappa^T=-\kappa,
$$

$$
(\kappa^{++})^T=-\kappa^{++},\qquad (\kappa^{--})^T=-\kappa^{--},\qquad \kappa^{-+}=-(\kappa^{+-})^T.
$$

No spherical, axial, or time-reversal symmetry is assumed; all four branch blocks of each tensor are retained.

### Axial Symmetry

Under a rotation about the $z$ axis,

$$
\hat R_z(\theta)\phi_\alpha^+=e^{-im_\alpha\theta}\phi_\alpha^+,\qquad \hat R_z(\theta)\phi_\alpha^-=e^{+im_\alpha\theta}\phi_\alpha^-.
$$

Axial invariance requires

$$
R_z(\theta)\rho R_z^\dagger(\theta)=\rho,\qquad R_z(\theta)\kappa R_z^T(\theta)=\kappa.
$$

Since $m_\alpha,m_\beta>0$, the normal density connects states within the same branch, while the pairing tensor connects opposite branches:

$$
\boxed{\rho=\begin{pmatrix}\rho^{++}&0\\0&\rho^{--}\end{pmatrix},\qquad \kappa=\begin{pmatrix}0&\kappa^{+-}\\\kappa^{-+}&0\end{pmatrix},\qquad \kappa^{-+}=-(\kappa^{+-})^T.}
$$

The remaining matrix elements satisfy

$$
(m_\alpha-m_\beta)\rho_{\alpha\beta}^{++}=0,\qquad (m_\alpha-m_\beta)\rho_{\alpha\beta}^{--}=0,\qquad (m_\alpha-m_\beta)\kappa_{\alpha\beta}^{+-}=0.
$$

Each block is labeled by the positive projection $m$:

$$
\mathcal B_m=\{\alpha:m_\alpha=m\},\qquad m>0.
$$

Different $N,j$ states may mix within the same block. For $\alpha,\beta\in\mathcal B_m$, define the normal-density combinations

$$
\boxed{R_{\alpha\beta}^{(+)}=\rho_{\alpha\beta}^{++}+\rho_{\alpha\beta}^{--},\qquad R_{\alpha\beta}^{(-)}=\rho_{\alpha\beta}^{++}-\rho_{\alpha\beta}^{--}.}
$$

The superscripts on $\rho$ denote the $+m$ and $-m$ branches; the parenthesized superscripts on $R$ denote their sum and difference. Axial symmetry alone imposes neither parity conservation nor time-reversal invariance.

### Particle Density

The particle density is

$$
\rho(\boldsymbol r)=\left\langle\hat\psi^\dagger(\boldsymbol r)\hat\psi(\boldsymbol r)\right\rangle=\sum_{\text{blocks}}\sum_{\alpha,\beta}\left[\rho_{\alpha\beta}^{++}(\phi_\beta^+)^\dagger(\boldsymbol r)\phi_\alpha^+(\boldsymbol r)+\rho_{\alpha\beta}^{--}(\phi_\beta^-)^\dagger(\boldsymbol r)\phi_\alpha^-(\boldsymbol r)\right].
$$

Write $\phi_\alpha(r)=\phi_{n_\alpha l_\alpha}(r)$. Contracting the spin components gives

$$
(\phi_\beta^+)^\dagger(\boldsymbol r)\phi_\alpha^+(\boldsymbol r)=\phi_\beta(r)\phi_\alpha(r)\left[\mathcal Y^{l_\beta}_{j_\beta m_\beta}(\hat{\mathbf r})\right]^\dagger\mathcal Y^{l_\alpha}_{j_\alpha m_\alpha}(\hat{\mathbf r})=\phi_\beta(r)\phi_\alpha(r)A_{\rho,\alpha\beta}^{++}(\hat{\mathbf r}),
$$

$$
(\phi_\beta^-)^\dagger(\boldsymbol r)\phi_\alpha^-(\boldsymbol r)=\phi_\beta(r)\phi_\alpha(r)\left[\mathcal Y^{l_\beta}_{j_\beta,-m_\beta}(\hat{\mathbf r})\right]^\dagger\mathcal Y^{l_\alpha}_{j_\alpha,-m_\alpha}(\hat{\mathbf r})=\phi_\beta(r)\phi_\alpha(r)A_{\rho,\alpha\beta}^{--}(\hat{\mathbf r}).
$$

The sums over $\alpha,\beta$ run within each axial block. The particle-density angular factors are

$$
A_{\rho,\alpha\beta}^{++}(\hat{\mathbf r})=\left[\mathcal Y^{l_\beta}_{j_\beta m_\beta}(\hat{\mathbf r})\right]^\dagger\mathcal Y^{l_\alpha}_{j_\alpha m_\alpha}(\hat{\mathbf r}),\qquad A_{\rho,\alpha\beta}^{--}(\hat{\mathbf r})=\left[\mathcal Y^{l_\beta}_{j_\beta,-m_\beta}(\hat{\mathbf r})\right]^\dagger\mathcal Y^{l_\alpha}_{j_\alpha,-m_\alpha}(\hat{\mathbf r}).
$$

Define their spherical-harmonic expansion coefficients by

$$
A_{\rho,\alpha\beta}^{++,LM}=\int d\hat{\mathbf r}\,Y_{LM}^{*}(\hat{\mathbf r})A_{\rho,\alpha\beta}^{++}(\hat{\mathbf r})=(-1)^M\left\langle(l_\beta s)j_\beta m_\beta\middle|Y_{L,-M}\middle|(l_\alpha s)j_\alpha m_\alpha\right\rangle,
$$

$$
A_{\rho,\alpha\beta}^{--,LM}=\int d\hat{\mathbf r}\,Y_{LM}^{*}(\hat{\mathbf r})A_{\rho,\alpha\beta}^{--}(\hat{\mathbf r})=(-1)^M\left\langle(l_\beta s)j_\beta,-m_\beta\middle|Y_{L,-M}\middle|(l_\alpha s)j_\alpha,-m_\alpha\right\rangle,
$$

where $s=\frac12$ and $Y_{LM}^{*}=(-1)^M Y_{L,-M}$. Reconstructing the angular factors gives

$$
A_{\rho,\alpha\beta}^{++}(\hat{\mathbf r})=\sum_{LM}A_{\rho,\alpha\beta}^{++,LM}Y_{LM}(\hat{\mathbf r}),\qquad A_{\rho,\alpha\beta}^{--}(\hat{\mathbf r})=\sum_{LM}A_{\rho,\alpha\beta}^{--,LM}Y_{LM}(\hat{\mathbf r}).
$$

The particle density therefore has the multipole expansion

$$
\rho(\boldsymbol r)=\sum_{\text{blocks}}\sum_{\alpha,\beta}\phi_\beta(r)\phi_\alpha(r)\left[\rho_{\alpha\beta}^{++}A_{\rho,\alpha\beta}^{++}(\hat{\mathbf r})+\rho_{\alpha\beta}^{--}A_{\rho,\alpha\beta}^{--}(\hat{\mathbf r})\right]=\sum_{LM}\rho_{LM}(r)Y_{LM}(\hat{\mathbf r}),
$$

with radial multipole components

$$
\boxed{\rho_{LM}(r)=\sum_{\text{blocks}}\sum_{\alpha,\beta}\phi_\beta(r)\phi_\alpha(r)\left[\rho_{\alpha\beta}^{++}A_{\rho,\alpha\beta}^{++,LM}+\rho_{\alpha\beta}^{--}A_{\rho,\alpha\beta}^{--,LM}\right].}
$$

### Kinetic Density

The kinetic density is

$$
\tau(\boldsymbol r)=\left\langle\nabla\hat\psi^\dagger(\boldsymbol r)\cdot\nabla\hat\psi(\boldsymbol r)\right\rangle=\sum_{\text{blocks}}\sum_{\alpha,\beta}\left[\rho_{\alpha\beta}^{++}\nabla(\phi_\beta^+)^\dagger(\boldsymbol r)\cdot\nabla\phi_\alpha^+(\boldsymbol r)+\rho_{\alpha\beta}^{--}\nabla(\phi_\beta^-)^\dagger(\boldsymbol r)\cdot\nabla\phi_\alpha^-(\boldsymbol r)\right].
$$

Following the [spherical-coordinate differential operators](../../notes/mathematical_physics/coordinate_spherical.md#differential-operators), use

$$
\nabla=\hat{\mathbf r}\,\partial_r+\frac1r\nabla_{\hat{\mathbf r}},\qquad \Delta_{\hat{\mathbf r}}\equiv\nabla_{\hat{\mathbf r}}^2=-\frac{\hat{\mathbf L}^{\,2}}{\hbar^2}.
$$

Radial and angular directions are orthogonal. With $\phi_\alpha(r)=\phi_{n_\alpha l_\alpha}(r)$, the branch contractions become

$$
\nabla(\phi_\beta^+)^\dagger\cdot\nabla\phi_\alpha^+=\partial_r\phi_\beta\,\partial_r\phi_\alpha\,A_{\rho,\alpha\beta}^{++}(\hat{\mathbf r})+\frac{\phi_\beta\phi_\alpha}{r^2}A_{\tau,\alpha\beta}^{++}(\hat{\mathbf r}),
$$

$$
\nabla(\phi_\beta^-)^\dagger\cdot\nabla\phi_\alpha^-=\partial_r\phi_\beta\,\partial_r\phi_\alpha\,A_{\rho,\alpha\beta}^{--}(\hat{\mathbf r})+\frac{\phi_\beta\phi_\alpha}{r^2}A_{\tau,\alpha\beta}^{--}(\hat{\mathbf r}).
$$

For angular spinors $u,v$ in a fixed spin basis, use the spherical Laplacian product rule,

$$
2(\nabla_{\hat{\mathbf r}}u)^\dagger\cdot\nabla_{\hat{\mathbf r}}v=\Delta_{\hat{\mathbf r}}(u^\dagger v)-(\Delta_{\hat{\mathbf r}}u)^\dagger v-u^\dagger\Delta_{\hat{\mathbf r}}v,
$$

and $\Delta_{\hat{\mathbf r}}\mathcal Y^l_{jm}=-l(l+1)\mathcal Y^l_{jm}$ to define the kinetic-density angular factors as

$$
A_{\tau,\alpha\beta}^{++}(\hat{\mathbf r})=\left[\nabla_{\hat{\mathbf r}}\mathcal Y^{l_\beta}_{j_\beta m_\beta}(\hat{\mathbf r})\right]^\dagger\cdot\nabla_{\hat{\mathbf r}}\mathcal Y^{l_\alpha}_{j_\alpha m_\alpha}(\hat{\mathbf r})=\frac12\left[\Delta_{\hat{\mathbf r}}+l_\alpha(l_\alpha+1)+l_\beta(l_\beta+1)\right]A_{\rho,\alpha\beta}^{++}(\hat{\mathbf r}),
$$

$$
A_{\tau,\alpha\beta}^{--}(\hat{\mathbf r})=\left[\nabla_{\hat{\mathbf r}}\mathcal Y^{l_\beta}_{j_\beta,-m_\beta}(\hat{\mathbf r})\right]^\dagger\cdot\nabla_{\hat{\mathbf r}}\mathcal Y^{l_\alpha}_{j_\alpha,-m_\alpha}(\hat{\mathbf r})=\frac12\left[\Delta_{\hat{\mathbf r}}+l_\alpha(l_\alpha+1)+l_\beta(l_\beta+1)\right]A_{\rho,\alpha\beta}^{--}(\hat{\mathbf r}).
$$

Using $\Delta_{\hat{\mathbf r}}Y_{LM}=-L(L+1)Y_{LM}$, define the kinetic-density multipole coefficients as

$$
A_{\tau,\alpha\beta}^{++,LM}=\frac{l_\alpha(l_\alpha+1)+l_\beta(l_\beta+1)-L(L+1)}{2}\,A_{\rho,\alpha\beta}^{++,LM},
$$

$$
A_{\tau,\alpha\beta}^{--,LM}=\frac{l_\alpha(l_\alpha+1)+l_\beta(l_\beta+1)-L(L+1)}{2}\,A_{\rho,\alpha\beta}^{--,LM}.
$$

Reconstructing the angular factors gives

$$
A_{\tau,\alpha\beta}^{++}(\hat{\mathbf r})=\sum_{LM}A_{\tau,\alpha\beta}^{++,LM}Y_{LM}(\hat{\mathbf r}),\qquad A_{\tau,\alpha\beta}^{--}(\hat{\mathbf r})=\sum_{LM}A_{\tau,\alpha\beta}^{--,LM}Y_{LM}(\hat{\mathbf r}).
$$

The kinetic density therefore has the multipole expansion

$$
\tau(\boldsymbol r)=\sum_{LM}\tau_{LM}(r)Y_{LM}(\hat{\mathbf r}),
$$

with radial multipole components

$$
\boxed{\tau_{LM}(r)=\sum_{\text{blocks}}\sum_{\alpha,\beta}\left\{\partial_r\phi_\beta(r)\,\partial_r\phi_\alpha(r)\left[\rho_{\alpha\beta}^{++}A_{\rho,\alpha\beta}^{++,LM}+\rho_{\alpha\beta}^{--}A_{\rho,\alpha\beta}^{--,LM}\right]+\frac{\phi_\beta(r)\phi_\alpha(r)}{r^2}\left[\rho_{\alpha\beta}^{++}A_{\tau,\alpha\beta}^{++,LM}+\rho_{\alpha\beta}^{--}A_{\tau,\alpha\beta}^{--,LM}\right]\right\}.}
$$
