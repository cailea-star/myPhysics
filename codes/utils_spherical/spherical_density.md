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

### Laplacian of Particle Density

Following the [spherical-coordinate differential operators](../../notes/mathematical_physics/coordinate_spherical.md#differential-operators), use

$$
\nabla^2=\frac1{r^2}\partial_r\!\left(r^2\partial_r\right)+\frac1{r^2}\Delta_{\hat{\mathbf r}}=\partial_r^2+\frac2r\partial_r+\frac1{r^2}\Delta_{\hat{\mathbf r}}.
$$

Applying this operator to the particle-density multipole expansion and using $\Delta_{\hat{\mathbf r}}Y_{LM}=-L(L+1)Y_{LM}$ gives

$$
\nabla^2\rho(\boldsymbol r)=\sum_{LM}\left[\partial_r^2+\frac2r\partial_r-\frac{L(L+1)}{r^2}\right]\rho_{LM}(r)Y_{LM}(\hat{\mathbf r}).
$$

The density-matrix elements and angular coefficients are independent of $r$. Substituting $\rho_{LM}(r)$ therefore yields

$$
\nabla^2\rho(\boldsymbol r)=\sum_{LM}Y_{LM}(\hat{\mathbf r})\sum_{\text{blocks}}\sum_{\alpha,\beta}\left\{\left[\partial_r^2+\frac2r\partial_r-\frac{L(L+1)}{r^2}\right]\bigl[\phi_\beta(r)\phi_\alpha(r)\bigr]\right\}\left[\rho_{\alpha\beta}^{++}A_{\rho,\alpha\beta}^{++,LM}+\rho_{\alpha\beta}^{--}A_{\rho,\alpha\beta}^{--,LM}\right].
$$

Expanding the radial derivatives gives

$$
\boxed{\nabla^2\rho(\boldsymbol r)=\sum_{LM}Y_{LM}(\hat{\mathbf r})\sum_{\text{blocks}}\sum_{\alpha,\beta}\left[\partial_r^2\phi_\beta\,\phi_\alpha+2\partial_r\phi_\beta\,\partial_r\phi_\alpha+\phi_\beta\,\partial_r^2\phi_\alpha+\frac2r\left(\partial_r\phi_\beta\,\phi_\alpha+\phi_\beta\,\partial_r\phi_\alpha\right)-\frac{L(L+1)}{r^2}\phi_\beta\phi_\alpha\right]\left[\rho_{\alpha\beta}^{++}A_{\rho,\alpha\beta}^{++,LM}+\rho_{\alpha\beta}^{--}A_{\rho,\alpha\beta}^{--,LM}\right].}
$$

All radial functions and their derivatives are evaluated at $r$.

### Current Density

The current density is defined without the factor $\hbar/M$ as

$$
\boldsymbol j(\boldsymbol r)=\frac1{2i}\left\langle\hat\psi^\dagger(\boldsymbol r)\nabla\hat\psi(\boldsymbol r)-[\nabla\hat\psi^\dagger(\boldsymbol r)]\hat\psi(\boldsymbol r)\right\rangle.
$$

Following the [spherical-coordinate differential operators](../../notes/mathematical_physics/coordinate_spherical.md#differential-operators), use

$$
\nabla=\hat{\mathbf r}\,\partial_r+\frac1r\nabla_{\hat{\mathbf r}},\qquad \nabla_{\hat{\mathbf r}}=-\frac{i}{\hbar}\hat{\mathbf r}\times\hat{\mathbf L}.
$$

With the real radial functions $\phi_\alpha(r)=\phi_{n_\alpha l_\alpha}(r)$, the two branch contributions separate as

$$
\frac{(\phi_\beta^+)^\dagger\nabla\phi_\alpha^+-[\nabla(\phi_\beta^+)^\dagger]\phi_\alpha^+}{2i}=\frac{\hat{\mathbf r}}{2i}\left(\phi_\beta\partial_r\phi_\alpha-\partial_r\phi_\beta\,\phi_\alpha\right)A_{\rho,\alpha\beta}^{++}+\frac{\phi_\beta\phi_\alpha}{r}\boldsymbol A_{j,\alpha\beta}^{++},
$$

$$
\frac{(\phi_\beta^-)^\dagger\nabla\phi_\alpha^--[\nabla(\phi_\beta^-)^\dagger]\phi_\alpha^-}{2i}=\frac{\hat{\mathbf r}}{2i}\left(\phi_\beta\partial_r\phi_\alpha-\partial_r\phi_\beta\,\phi_\alpha\right)A_{\rho,\alpha\beta}^{--}+\frac{\phi_\beta\phi_\alpha}{r}\boldsymbol A_{j,\alpha\beta}^{--}.
$$

Define the tangential current-density angular factors as

$$
\boldsymbol A_{j,\alpha\beta}^{++}=\frac1{2i}\left\{\left[\mathcal Y^{l_\beta}_{j_\beta m_\beta}\right]^\dagger\nabla_{\hat{\mathbf r}}\mathcal Y^{l_\alpha}_{j_\alpha m_\alpha}-\left[\nabla_{\hat{\mathbf r}}\mathcal Y^{l_\beta}_{j_\beta m_\beta}\right]^\dagger\mathcal Y^{l_\alpha}_{j_\alpha m_\alpha}\right\},
$$

$$
\boldsymbol A_{j,\alpha\beta}^{--}=\frac1{2i}\left\{\left[\mathcal Y^{l_\beta}_{j_\beta,-m_\beta}\right]^\dagger\nabla_{\hat{\mathbf r}}\mathcal Y^{l_\alpha}_{j_\alpha,-m_\alpha}-\left[\nabla_{\hat{\mathbf r}}\mathcal Y^{l_\beta}_{j_\beta,-m_\beta}\right]^\dagger\mathcal Y^{l_\alpha}_{j_\alpha,-m_\alpha}\right\}.
$$

All angular functions are evaluated at $\hat{\mathbf r}$. Following the [tangential vector spherical harmonics](../../notes/mathematical_physics/specialfunction_spherical_harmonics.md#tangential-vector-spherical-harmonics), expand

$$
\boldsymbol A_{j,\alpha\beta}^{++}(\hat{\mathbf r})=\sum_{L\geq1,M}\left[A_{j,\alpha\beta}^{++,LM}(\Psi)\boldsymbol\Psi_{LM}(\hat{\mathbf r})+A_{j,\alpha\beta}^{++,LM}(\Phi)\boldsymbol\Phi_{LM}(\hat{\mathbf r})\right],
$$

$$
\boldsymbol A_{j,\alpha\beta}^{--}(\hat{\mathbf r})=\sum_{L\geq1,M}\left[A_{j,\alpha\beta}^{--,LM}(\Psi)\boldsymbol\Psi_{LM}(\hat{\mathbf r})+A_{j,\alpha\beta}^{--,LM}(\Phi)\boldsymbol\Phi_{LM}(\hat{\mathbf r})\right].
$$

For either branch, let $u$ and $v$ denote the $\beta$ and $\alpha$ spinors in the fixed spin basis. The divergence and surface curl satisfy

$$
\nabla_{\hat{\mathbf r}}\cdot\frac{u^\dagger\nabla_{\hat{\mathbf r}}v-(\nabla_{\hat{\mathbf r}}u)^\dagger v}{2i}=\frac{u^\dagger\Delta_{\hat{\mathbf r}}v-(\Delta_{\hat{\mathbf r}}u)^\dagger v}{2i},
$$

$$
\operatorname{curl}_{\hat{\mathbf r}}\frac{u^\dagger\nabla_{\hat{\mathbf r}}v-(\nabla_{\hat{\mathbf r}}u)^\dagger v}{2i}=-\frac1{2\hbar}\nabla_{\hat{\mathbf r}}\cdot\left[u^\dagger\hat{\mathbf L}v+(\hat{\mathbf L}u)^\dagger v\right]=-\frac1\hbar(\nabla_{\hat{\mathbf r}}u)^\dagger\cdot\hat{\mathbf L}v.
$$

Here $\operatorname{curl}_{\hat{\mathbf r}}\mathbf V=-\nabla_{\hat{\mathbf r}}\cdot(\hat{\mathbf r}\times\mathbf V)$. Using $\Delta_{\hat{\mathbf r}}\mathcal Y^l_{jm}=-l(l+1)\mathcal Y^l_{jm}$, spherical integration by parts gives the gradient-type coefficients

$$
A_{j,\alpha\beta}^{++,LM}(\Psi)=\frac{l_\alpha(l_\alpha+1)-l_\beta(l_\beta+1)}{2i\sqrt{L(L+1)}}A_{\rho,\alpha\beta}^{++,LM},
$$

$$
A_{j,\alpha\beta}^{--,LM}(\Psi)=\frac{l_\alpha(l_\alpha+1)-l_\beta(l_\beta+1)}{2i\sqrt{L(L+1)}}A_{\rho,\alpha\beta}^{--,LM}.
$$

For the rotational-type coefficients, the solid-harmonic gradient relation and $\hat{\mathbf r}\cdot\hat{\mathbf L}=0$ give

$$
(\nabla_{\hat{\mathbf r}}Y_{LM})\cdot\hat{\mathbf L}=\sqrt{L(2L+1)}[Y_{L-1}\otimes\hat{\mathbf L}]_{LM}.
$$

The coupled tensor is explicitly

$$
[Y_{L-1}\otimes\hat{\mathbf L}]_{LM}=\sum_{\mu=-1}^{1}\langle L-1,M-\mu;1,\mu|LM\rangle Y_{L-1,M-\mu}\hat L_\mu,
$$

where $\hat L_0=\hat L_z$ and $\hat L_{\pm1}=\mp(\hat L_x\pm i\hat L_y)/\sqrt2$. Spherical integration by parts then yields

$$
A_{j,\alpha\beta}^{++,LM}(\Phi)=-(-1)^M\sqrt{\frac{2L+1}{L+1}}\left\langle(l_\beta s)j_\beta m_\beta\middle|\left[Y_{L-1}\otimes\frac{\hat{\mathbf L}}{\hbar}\right]_{L,-M}\middle|(l_\alpha s)j_\alpha m_\alpha\right\rangle,
$$

$$
A_{j,\alpha\beta}^{--,LM}(\Phi)=-(-1)^M\sqrt{\frac{2L+1}{L+1}}\left\langle(l_\beta s)j_\beta,-m_\beta\middle|\left[Y_{L-1}\otimes\frac{\hat{\mathbf L}}{\hbar}\right]_{L,-M}\middle|(l_\alpha s)j_\alpha,-m_\alpha\right\rangle.
$$

Here $s=\frac12$, and $\hat{\mathbf L}$ acts on the ket. Combining both branches within each axial block gives

$$
\boxed{\boldsymbol j(\boldsymbol r)=\sum_{\text{blocks}}\sum_{\alpha,\beta}\left\{\frac{\hat{\mathbf r}}{2i}\left(\phi_\beta\partial_r\phi_\alpha-\partial_r\phi_\beta\,\phi_\alpha\right)\left[\rho_{\alpha\beta}^{++}A_{\rho,\alpha\beta}^{++}(\hat{\mathbf r})+\rho_{\alpha\beta}^{--}A_{\rho,\alpha\beta}^{--}(\hat{\mathbf r})\right]+\frac{\phi_\beta\phi_\alpha}{r}\left[\rho_{\alpha\beta}^{++}\boldsymbol A_{j,\alpha\beta}^{++}(\hat{\mathbf r})+\rho_{\alpha\beta}^{--}\boldsymbol A_{j,\alpha\beta}^{--}(\hat{\mathbf r})\right]\right\}.}
$$

All radial functions and their derivatives are evaluated at $r$.

### Spin Density

The spin density is defined without the factor $\hbar/2$ as

$$
\boldsymbol s(\boldsymbol r)=\left\langle\hat\psi^\dagger(\boldsymbol r)\boldsymbol\sigma\hat\psi(\boldsymbol r)\right\rangle,
$$

where $\boldsymbol\sigma$ denotes the Pauli matrices. Under axial symmetry,

$$
\boldsymbol s(\boldsymbol r)=\sum_{\text{blocks}}\sum_{\alpha,\beta}\left[\rho_{\alpha\beta}^{++}(\phi_\beta^+)^\dagger\boldsymbol\sigma\phi_\alpha^++\rho_{\alpha\beta}^{--}(\phi_\beta^-)^\dagger\boldsymbol\sigma\phi_\alpha^-\right].
$$

With the real radial functions $\phi_\alpha(r)=\phi_{n_\alpha l_\alpha}(r)$, define the spin-density angular factors as

$$
\boldsymbol A_{s,\alpha\beta}^{++}(\hat{\mathbf r})=\left[\mathcal Y^{l_\beta}_{j_\beta m_\beta}(\hat{\mathbf r})\right]^\dagger\boldsymbol\sigma\mathcal Y^{l_\alpha}_{j_\alpha m_\alpha}(\hat{\mathbf r}),
$$

$$
\boldsymbol A_{s,\alpha\beta}^{--}(\hat{\mathbf r})=\left[\mathcal Y^{l_\beta}_{j_\beta,-m_\beta}(\hat{\mathbf r})\right]^\dagger\boldsymbol\sigma\mathcal Y^{l_\alpha}_{j_\alpha,-m_\alpha}(\hat{\mathbf r}).
$$

The spin density can have both radial and tangential components. Following the [coupled vector spherical harmonics](../../notes/mathematical_physics/specialfunction_spherical_harmonics.md#vector-spherical-harmonics), use

$$
\mathbf Y_{LM}^{\lambda}(\hat{\mathbf r})=\sum_{\mu=-1}^{1}\langle\lambda,M-\mu;1,\mu|L,M\rangle Y_{\lambda,M-\mu}(\hat{\mathbf r})\mathbf e_\mu,
$$

where $\lambda$ is the orbital rank, $L=|\lambda-1|,\ldots,\lambda+1$ is the total rank, and $M=-L,\ldots,L$. Expand the angular factors as

$$
\boldsymbol A_{s,\alpha\beta}^{++}(\hat{\mathbf r})=\sum_{\lambda LM}A_{s,\alpha\beta}^{++,\lambda LM}\mathbf Y_{LM}^{\lambda}(\hat{\mathbf r}),
$$

$$
\boldsymbol A_{s,\alpha\beta}^{--}(\hat{\mathbf r})=\sum_{\lambda LM}A_{s,\alpha\beta}^{--,\lambda LM}\mathbf Y_{LM}^{\lambda}(\hat{\mathbf r}).
$$

Using $[\mathbf Y_{LM}^{\lambda}]^*=(-1)^{L+\lambda+M+1}\mathbf Y_{L,-M}^{\lambda}$, the expansion coefficients become

$$
A_{s,\alpha\beta}^{++,\lambda LM}=\int d\hat{\mathbf r}\,[\mathbf Y_{LM}^{\lambda}]^*\cdot\boldsymbol A_{s,\alpha\beta}^{++}=(-1)^{L+\lambda+M+1}\left\langle(l_\beta s)j_\beta m_\beta\middle|[Y_\lambda\otimes\boldsymbol\sigma]_{L,-M}\middle|(l_\alpha s)j_\alpha m_\alpha\right\rangle,
$$

$$
A_{s,\alpha\beta}^{--,\lambda LM}=\int d\hat{\mathbf r}\,[\mathbf Y_{LM}^{\lambda}]^*\cdot\boldsymbol A_{s,\alpha\beta}^{--}=(-1)^{L+\lambda+M+1}\left\langle(l_\beta s)j_\beta,-m_\beta\middle|[Y_\lambda\otimes\boldsymbol\sigma]_{L,-M}\middle|(l_\alpha s)j_\alpha,-m_\alpha\right\rangle.
$$

Here $s=\frac12$. The coupled operator acts on the orbital and spin spaces separately:

$$
[Y_\lambda\otimes\boldsymbol\sigma]_{LM}=\sum_{\mu=-1}^{1}\langle\lambda,M-\mu;1,\mu|L,M\rangle Y_{\lambda,M-\mu}\sigma_\mu,
$$

with $\sigma_0=\sigma_z$ and $\sigma_{\pm1}=\mp(\sigma_x\pm i\sigma_y)/\sqrt2$.

Combining both branches gives the vector multipole expansion

$$
\boldsymbol s(\boldsymbol r)=\sum_{\text{blocks}}\sum_{\alpha,\beta}\phi_\beta(r)\phi_\alpha(r)\left[\rho_{\alpha\beta}^{++}\boldsymbol A_{s,\alpha\beta}^{++}(\hat{\mathbf r})+\rho_{\alpha\beta}^{--}\boldsymbol A_{s,\alpha\beta}^{--}(\hat{\mathbf r})\right]=\sum_{\lambda LM}s_{\lambda LM}(r)\mathbf Y_{LM}^{\lambda}(\hat{\mathbf r}),
$$

with radial coefficients

$$
\boxed{s_{\lambda LM}(r)=\sum_{\text{blocks}}\sum_{\alpha,\beta}\phi_\beta(r)\phi_\alpha(r)\left[\rho_{\alpha\beta}^{++}A_{s,\alpha\beta}^{++,\lambda LM}+\rho_{\alpha\beta}^{--}A_{s,\alpha\beta}^{--,\lambda LM}\right].}
$$
