# Cylindrical Densities

### Basis and Block Matrices

Following the [cylindrical basis](cylindrical_basis.md), label each positive-branch representative by

$$
\alpha=\{n_z,n_r,\Lambda,\Sigma\},\qquad \Sigma=\pm\frac12.
$$

Define the real spatial amplitude as

$$
\phi_\alpha(z,r_\perp)=\phi_{n_{z,\alpha}}(z)\phi_{n_{r,\alpha}}^{\Lambda_\alpha}(r_\perp).
$$

The positive-branch basis spinor and its negative-branch partner are

$$
\phi_\alpha^+(\mathbf r)=\phi_\alpha(z,r_\perp)\frac{e^{i\Lambda_\alpha\varphi}}{\sqrt{2\pi}}\chi_{\Sigma_\alpha},\qquad \phi_\alpha^-(\mathbf r)=\phi_\alpha(z,r_\perp)\frac{e^{-i\Lambda_\alpha\varphi}}{\sqrt{2\pi}}\chi_{-\Sigma_\alpha},
$$

where

$$
\chi_\uparrow=\begin{pmatrix}1\\0\end{pmatrix},\qquad \chi_\downarrow=\begin{pmatrix}0\\1\end{pmatrix}.
$$

Both branches use the same representative label and spatial amplitude. Their physical spin and orbital projections are opposite.

Let $\hat c_\alpha^\pm$ annihilate the corresponding branch state; the superscripts denote branches. The fields are

$$
\hat\psi(\mathbf r)=\sum_\alpha\left[\hat c_\alpha^+\phi_\alpha^+(\mathbf r)+\hat c_\alpha^-\phi_\alpha^-(\mathbf r)\right],\qquad \hat\psi^\dagger(\mathbf r)=\sum_\alpha\left[(\hat c_\alpha^+)^\dagger(\phi_\alpha^+)^\dagger(\mathbf r)+(\hat c_\alpha^-)^\dagger(\phi_\alpha^-)^\dagger(\mathbf r)\right].
$$

The normal density blocks and pairing blocks are

$$
\rho_{\alpha\beta}^{++}=\left\langle(\hat c_\beta^+)^\dagger\hat c_\alpha^+\right\rangle,\qquad \rho_{\alpha\beta}^{--}=\left\langle(\hat c_\beta^-)^\dagger\hat c_\alpha^-\right\rangle,\qquad \kappa_{\alpha\beta}^{+-}=\left\langle\hat c_\beta^-\hat c_\alpha^+\right\rangle,\qquad \kappa_{\alpha\beta}^{-+}=\left\langle\hat c_\beta^+\hat c_\alpha^-\right\rangle.
$$

Following the [axial block structure](../utils/hfb_axial.md#density-blocks), order positive-branch states before their partners within each block:

$$
\rho=\begin{pmatrix}\rho^{++}&0\\0&\rho^{--}\end{pmatrix},\qquad \kappa=\begin{pmatrix}0&\kappa^{+-}\\\kappa^{-+}&0\end{pmatrix}.
$$

Hermiticity and fermionic antisymmetry imply

$$
\rho^\dagger=\rho,\qquad \kappa^T=-\kappa,\qquad \kappa^{-+}=-(\kappa^{+-})^T.
$$

### Spin Groups

Within each positive-branch block with angular-momentum projection $\Omega>0$, separate the representative labels by spin:

$$
\{\alpha\}=\{\alpha^\uparrow\}\cup\{\alpha^\downarrow\}.
$$

Since $\Lambda_\alpha+\Sigma_\alpha=\Omega$,

$$
\Sigma_{\alpha^\uparrow}=\frac12,\qquad \Lambda_{\alpha^\uparrow}=\Omega-\frac12,\qquad \Sigma_{\alpha^\downarrow}=-\frac12,\qquad \Lambda_{\alpha^\downarrow}=\Omega+\frac12.
$$

The arrows label representative spins in both branches. The corresponding basis spinors have spin factors

$$
\phi_{\alpha^\uparrow}^+\propto\chi_\uparrow,\qquad \phi_{\alpha^\uparrow}^-\propto\chi_\downarrow,\qquad \phi_{\alpha^\downarrow}^+\propto\chi_\downarrow,\qquad \phi_{\alpha^\downarrow}^-\propto\chi_\uparrow.
$$

For fixed representative indices, both branches share the same real spatial amplitudes. Define their combined matrix elements by

$$
R_{\alpha^\uparrow,\beta^\uparrow}^{(\pm)}=\rho_{\alpha^\uparrow\beta^\uparrow}^{++}\pm\rho_{\alpha^\uparrow\beta^\uparrow}^{--},\qquad R_{\alpha^\downarrow,\beta^\downarrow}^{(\pm)}=\rho_{\alpha^\downarrow\beta^\downarrow}^{++}\pm\rho_{\alpha^\downarrow\beta^\downarrow}^{--},
$$

$$
R_{\alpha^\uparrow,\beta^\downarrow}^{(\pm)}=\rho_{\alpha^\uparrow\beta^\downarrow}^{++}\pm\rho_{\alpha^\uparrow\beta^\downarrow}^{--},\qquad R_{\alpha^\downarrow,\beta^\uparrow}^{(\pm)}=\rho_{\alpha^\downarrow\beta^\uparrow}^{++}\pm\rho_{\alpha^\downarrow\beta^\uparrow}^{--}.
$$

### Particle Density

The particle density is

$$
\rho(\mathbf r)=\left\langle\hat\psi^\dagger(\mathbf r)\hat\psi(\mathbf r)\right\rangle=\sum_{\alpha,\beta}\left[\rho_{\alpha\beta}^{++}(\phi_\beta^+)^\dagger(\mathbf r)\phi_\alpha^+(\mathbf r)+\rho_{\alpha\beta}^{--}(\phi_\beta^-)^\dagger(\mathbf r)\phi_\alpha^-(\mathbf r)\right].
$$

The two branch contractions are

$$
(\phi_\beta^+)^\dagger(\mathbf r)\phi_\alpha^+(\mathbf r)=\frac{\phi_\beta(z,r_\perp)\phi_\alpha(z,r_\perp)}{2\pi}e^{i(\Lambda_\alpha-\Lambda_\beta)\varphi}\delta_{\Sigma_\alpha,\Sigma_\beta},
$$

$$
(\phi_\beta^-)^\dagger(\mathbf r)\phi_\alpha^-(\mathbf r)=\frac{\phi_\beta(z,r_\perp)\phi_\alpha(z,r_\perp)}{2\pi}e^{-i(\Lambda_\alpha-\Lambda_\beta)\varphi}\delta_{\Sigma_\alpha,\Sigma_\beta}.
$$

Cross-group terms vanish by spin orthogonality. Within each spin group of a block, $\Lambda_\alpha=\Lambda_\beta$, so both phases equal one. The two branches therefore share the same spatial kernel and combine through $R_{\alpha\beta}^{(+)}=\rho_{\alpha\beta}^{++}+\rho_{\alpha\beta}^{--}$:

$$
\boxed{\rho(z,r_\perp)=\frac{1}{2\pi}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(+)}\phi_{\alpha^\uparrow}\phi_{\beta^\uparrow}+\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(+)}\phi_{\alpha^\downarrow}\phi_{\beta^\downarrow}\right].}
$$

All spatial amplitudes are evaluated at $(z,r_\perp)$, and arrows denote representative spins. The result is the three-dimensional local density, with the factor $1/(2\pi)$ supplied by the angular normalization of the basis.

### Kinetic Density

The kinetic density is

$$
\tau(\mathbf r)=\left\langle\nabla\hat\psi^\dagger(\mathbf r)\cdot\nabla\hat\psi(\mathbf r)\right\rangle=\sum_{\alpha,\beta}\left[\rho_{\alpha\beta}^{++}\nabla(\phi_\beta^+)^\dagger(\mathbf r)\cdot\nabla\phi_\alpha^+(\mathbf r)+\rho_{\alpha\beta}^{--}\nabla(\phi_\beta^-)^\dagger(\mathbf r)\cdot\nabla\phi_\alpha^-(\mathbf r)\right].
$$

In cylindrical coordinates,

$$
\nabla_z=\partial_z,\qquad \nabla_r=\partial_{r_\perp},\qquad \nabla_\varphi=\frac{1}{r_\perp}\partial_\varphi.
$$

The two branch contractions are

$$
\nabla(\phi_\beta^+)^\dagger(\mathbf r)\cdot\nabla\phi_\alpha^+(\mathbf r)=\frac{e^{i(\Lambda_\alpha-\Lambda_\beta)\varphi}}{2\pi}\left[\partial_z\phi_\beta\,\partial_z\phi_\alpha+\partial_{r_\perp}\phi_\beta\,\partial_{r_\perp}\phi_\alpha+\frac{\Lambda_\alpha\Lambda_\beta}{r_\perp^2}\phi_\beta\phi_\alpha\right]\delta_{\Sigma_\alpha,\Sigma_\beta},
$$

$$
\nabla(\phi_\beta^-)^\dagger(\mathbf r)\cdot\nabla\phi_\alpha^-(\mathbf r)=\frac{e^{-i(\Lambda_\alpha-\Lambda_\beta)\varphi}}{2\pi}\left[\partial_z\phi_\beta\,\partial_z\phi_\alpha+\partial_{r_\perp}\phi_\beta\,\partial_{r_\perp}\phi_\alpha+\frac{\Lambda_\alpha\Lambda_\beta}{r_\perp^2}\phi_\beta\phi_\alpha\right]\delta_{\Sigma_\alpha,\Sigma_\beta}.
$$

Cross-group terms vanish by spin orthogonality. Within each spin group of a block, $\Lambda_\alpha=\Lambda_\beta$, so both phases equal one. Reversing the branch changes both angular-derivative signs, leaving their product unchanged. The two branches therefore share the same derivative kernel and combine through $R_{\alpha\beta}^{(+)}=\rho_{\alpha\beta}^{++}+\rho_{\alpha\beta}^{--}$:

$$
\boxed{\tau(z,r_\perp)=\frac{1}{2\pi}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(+)}\left(\partial_z\phi_{\alpha^\uparrow}\partial_z\phi_{\beta^\uparrow}+\partial_{r_\perp}\phi_{\alpha^\uparrow}\partial_{r_\perp}\phi_{\beta^\uparrow}+\frac{\Lambda_{\alpha^\uparrow}^2}{r_\perp^2}\phi_{\alpha^\uparrow}\phi_{\beta^\uparrow}\right)+\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(+)}\left(\partial_z\phi_{\alpha^\downarrow}\partial_z\phi_{\beta^\downarrow}+\partial_{r_\perp}\phi_{\alpha^\downarrow}\partial_{r_\perp}\phi_{\beta^\downarrow}+\frac{\Lambda_{\alpha^\downarrow}^2}{r_\perp^2}\phi_{\alpha^\downarrow}\phi_{\beta^\downarrow}\right)\right].}
$$

All spatial amplitudes are evaluated at $(z,r_\perp)$, and arrows denote representative spins. The kinetic density $\tau$ excludes the kinetic-energy coefficient $\hbar^2/(2m)$.

### Density Laplacian

The density Laplacian follows from the product rule:

$$
\nabla^2\rho(\mathbf r)=\sum_{\alpha,\beta}\left\{\rho_{\alpha\beta}^{++}\left[(\nabla^2(\phi_\beta^+)^\dagger)\phi_\alpha^++2\nabla(\phi_\beta^+)^\dagger\cdot\nabla\phi_\alpha^++(\phi_\beta^+)^\dagger\nabla^2\phi_\alpha^+\right]+\rho_{\alpha\beta}^{--}\left[(\nabla^2(\phi_\beta^-)^\dagger)\phi_\alpha^-+2\nabla(\phi_\beta^-)^\dagger\cdot\nabla\phi_\alpha^-+(\phi_\beta^-)^\dagger\nabla^2\phi_\alpha^-\right]\right\}.
$$

In cylindrical coordinates,

$$
\nabla^2=\partial_z^2+\partial_{r_\perp}^2+\frac{1}{r_\perp}\partial_{r_\perp}+\frac{1}{r_\perp^2}\partial_\varphi^2.
$$

Define the spatial Laplacian amplitude

$$
L_\alpha=\left(\partial_z^2+\partial_{r_\perp}^2+\frac{1}{r_\perp}\partial_{r_\perp}-\frac{\Lambda_\alpha^2}{r_\perp^2}\right)\phi_\alpha(z,r_\perp).
$$

The two branch contractions are

$$
\nabla^2\left[(\phi_\beta^+)^\dagger(\mathbf r)\phi_\alpha^+(\mathbf r)\right]=\frac{e^{i(\Lambda_\alpha-\Lambda_\beta)\varphi}}{2\pi}\left[L_\beta\phi_\alpha+\phi_\beta L_\alpha+2\partial_z\phi_\beta\,\partial_z\phi_\alpha+2\partial_{r_\perp}\phi_\beta\,\partial_{r_\perp}\phi_\alpha+\frac{2\Lambda_\alpha\Lambda_\beta}{r_\perp^2}\phi_\beta\phi_\alpha\right]\delta_{\Sigma_\alpha,\Sigma_\beta},
$$

$$
\nabla^2\left[(\phi_\beta^-)^\dagger(\mathbf r)\phi_\alpha^-(\mathbf r)\right]=\frac{e^{-i(\Lambda_\alpha-\Lambda_\beta)\varphi}}{2\pi}\left[L_\beta\phi_\alpha+\phi_\beta L_\alpha+2\partial_z\phi_\beta\,\partial_z\phi_\alpha+2\partial_{r_\perp}\phi_\beta\,\partial_{r_\perp}\phi_\alpha+\frac{2\Lambda_\alpha\Lambda_\beta}{r_\perp^2}\phi_\beta\phi_\alpha\right]\delta_{\Sigma_\alpha,\Sigma_\beta}.
$$

Here the spatial amplitudes are evaluated at $(z,r_\perp)$, and arrows denote representative spins. Cross-group terms vanish by spin orthogonality. Within each spin group of a block, $\Lambda_\alpha=\Lambda_\beta$, so both azimuthal phases equal one. The two branches share the same spatial kernel and combine through $R_{\alpha\beta}^{(+)}=\rho_{\alpha\beta}^{++}+\rho_{\alpha\beta}^{--}$ defined in [Spin Groups](#spin-groups):

$$
\boxed{\nabla^2\rho(z,r_\perp)=\frac{1}{2\pi}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(+)}\left(L_{\alpha^\uparrow}\phi_{\beta^\uparrow}+\phi_{\alpha^\uparrow}L_{\beta^\uparrow}+2\partial_z\phi_{\alpha^\uparrow}\partial_z\phi_{\beta^\uparrow}+2\partial_{r_\perp}\phi_{\alpha^\uparrow}\partial_{r_\perp}\phi_{\beta^\uparrow}+\frac{2\Lambda_{\alpha^\uparrow}^2}{r_\perp^2}\phi_{\alpha^\uparrow}\phi_{\beta^\uparrow}\right)+\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(+)}\left(L_{\alpha^\downarrow}\phi_{\beta^\downarrow}+\phi_{\alpha^\downarrow}L_{\beta^\downarrow}+2\partial_z\phi_{\alpha^\downarrow}\partial_z\phi_{\beta^\downarrow}+2\partial_{r_\perp}\phi_{\alpha^\downarrow}\partial_{r_\perp}\phi_{\beta^\downarrow}+\frac{2\Lambda_{\alpha^\downarrow}^2}{r_\perp^2}\phi_{\alpha^\downarrow}\phi_{\beta^\downarrow}\right)\right].}
$$

The explicit angular terms cancel between the Laplacian amplitudes and the gradient product, consistent with

$$
\nabla^2\rho(z,r_\perp)=\left(\partial_z^2+\partial_{r_\perp}^2+\frac{1}{r_\perp}\partial_{r_\perp}\right)\rho(z,r_\perp).
$$

### Current Density

The current density is

$$
\mathbf j(\mathbf r)=\frac{1}{2i}\sum_{\alpha,\beta}\left\{\rho_{\alpha\beta}^{++}\left[(\phi_\beta^+)^\dagger\nabla\phi_\alpha^+-(\nabla(\phi_\beta^+)^\dagger)\phi_\alpha^+\right]+\rho_{\alpha\beta}^{--}\left[(\phi_\beta^-)^\dagger\nabla\phi_\alpha^--(\nabla(\phi_\beta^-)^\dagger)\phi_\alpha^-\right]\right\}.
$$

All basis spinors in this expression are evaluated at $\mathbf r$. The cylindrical derivatives are $\nabla_r=\partial_{r_\perp}$, $\nabla_\varphi=r_\perp^{-1}\partial_\varphi$, and $\nabla_z=\partial_z$. Only contractions within the same representative spin group survive. Within each spin group of a block, $\Lambda_\alpha=\Lambda_\beta$, so the azimuthal phases cancel.

##### $j_r$

The radial contractions in the two branches are

$$
\frac{1}{2i}\left[(\phi_\beta^\pm)^\dagger\nabla_r\phi_\alpha^\pm-(\nabla_r(\phi_\beta^\pm)^\dagger)\phi_\alpha^\pm\right]=\frac{1}{4\pi i}\left(\phi_\beta\partial_{r_\perp}\phi_\alpha-\phi_\alpha\partial_{r_\perp}\phi_\beta\right)e^{\pm i(\Lambda_\alpha-\Lambda_\beta)\varphi}\delta_{\Sigma_\alpha,\Sigma_\beta}.
$$

After phase cancellation, both branches share the same kernel and combine through $R^{(+)}$:

$$
\boxed{j_r(z,r_\perp)=\frac{1}{4\pi i}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(+)}\left(\phi_{\beta^\uparrow}\partial_{r_\perp}\phi_{\alpha^\uparrow}-\phi_{\alpha^\uparrow}\partial_{r_\perp}\phi_{\beta^\uparrow}\right)+\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(+)}\left(\phi_{\beta^\downarrow}\partial_{r_\perp}\phi_{\alpha^\downarrow}-\phi_{\alpha^\downarrow}\partial_{r_\perp}\phi_{\beta^\downarrow}\right)\right].}
$$

##### $j_\varphi$

The angular derivatives give opposite signs in the two branches:

$$
\frac{1}{2i}\left[(\phi_\beta^\pm)^\dagger\nabla_\varphi\phi_\alpha^\pm-(\nabla_\varphi(\phi_\beta^\pm)^\dagger)\phi_\alpha^\pm\right]=\pm\frac{\Lambda_\alpha+\Lambda_\beta}{4\pi r_\perp}\phi_\beta\phi_\alpha e^{\pm i(\Lambda_\alpha-\Lambda_\beta)\varphi}\delta_{\Sigma_\alpha,\Sigma_\beta}.
$$

After phase cancellation, $\Lambda_\alpha+\Lambda_\beta=2\Lambda_\alpha$, and the branch contributions combine through $R^{(-)}$:

$$
\boxed{j_\varphi(z,r_\perp)=\frac{1}{2\pi r_\perp}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(-)}\Lambda_{\alpha^\uparrow}\phi_{\alpha^\uparrow}\phi_{\beta^\uparrow}+\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(-)}\Lambda_{\alpha^\downarrow}\phi_{\alpha^\downarrow}\phi_{\beta^\downarrow}\right].}
$$

##### $j_z$

The longitudinal contractions in the two branches are

$$
\frac{1}{2i}\left[(\phi_\beta^\pm)^\dagger\nabla_z\phi_\alpha^\pm-(\nabla_z(\phi_\beta^\pm)^\dagger)\phi_\alpha^\pm\right]=\frac{1}{4\pi i}\left(\phi_\beta\partial_z\phi_\alpha-\phi_\alpha\partial_z\phi_\beta\right)e^{\pm i(\Lambda_\alpha-\Lambda_\beta)\varphi}\delta_{\Sigma_\alpha,\Sigma_\beta}.
$$

After phase cancellation, both branches share the same kernel and combine through $R^{(+)}$:

$$
\boxed{j_z(z,r_\perp)=\frac{1}{4\pi i}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(+)}\left(\phi_{\beta^\uparrow}\partial_z\phi_{\alpha^\uparrow}-\phi_{\alpha^\uparrow}\partial_z\phi_{\beta^\uparrow}\right)+\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(+)}\left(\phi_{\beta^\downarrow}\partial_z\phi_{\alpha^\downarrow}-\phi_{\alpha^\downarrow}\partial_z\phi_{\beta^\downarrow}\right)\right].}
$$

All spatial amplitudes are real, evaluated at $(z,r_\perp)$, and carry positive-branch labels. Hermiticity makes all three components real. For real density matrices, $j_r=j_z=0$, while $j_\varphi$ may remain nonzero. The definition of $\mathbf j$ does not include a factor of $\hbar/m$.

### Spin Density

The spin density is

$$
\mathbf s(\mathbf r)=\langle\hat\psi^\dagger(\mathbf r)\boldsymbol\sigma\hat\psi(\mathbf r)\rangle=\sum_{\alpha,\beta}\left[\rho_{\alpha\beta}^{++}(\phi_\beta^+)^\dagger(\mathbf r)\boldsymbol\sigma\phi_\alpha^+(\mathbf r)+\rho_{\alpha\beta}^{--}(\phi_\beta^-)^\dagger(\mathbf r)\boldsymbol\sigma\phi_\alpha^-(\mathbf r)\right].
$$

In the local cylindrical frame, the Pauli matrices are

$$
\sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad \sigma_r=\begin{pmatrix}0&e^{-i\varphi}\\e^{i\varphi}&0\end{pmatrix},\qquad \sigma_\varphi=\frac{1}{i}\begin{pmatrix}0&e^{-i\varphi}\\-e^{i\varphi}&0\end{pmatrix}.
$$

All spatial amplitudes below are real and evaluated at $(z,r_\perp)$. Arrows denote representative spins in both branches, as defined in [Spin Groups](#spin-groups).

##### $s_z$

The contractions within each representative spin group are

$$
(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_z\phi_{\alpha^\uparrow}^\pm=\pm\frac{\phi_{\beta^\uparrow}\phi_{\alpha^\uparrow}}{2\pi}e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\uparrow})\varphi},
$$

$$
(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_z\phi_{\alpha^\downarrow}^\pm=\mp\frac{\phi_{\beta^\downarrow}\phi_{\alpha^\downarrow}}{2\pi}e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\downarrow})\varphi}.
$$

Cross-group contractions vanish. Within each spin group of a block, $\Lambda_\alpha=\Lambda_\beta$, so both phases equal one. The negative branch reverses the spin sign, giving $R^{(-)}$:

$$
\boxed{s_z(z,r_\perp)=\frac{1}{2\pi}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(-)}\phi_{\alpha^\uparrow}\phi_{\beta^\uparrow}-\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(-)}\phi_{\alpha^\downarrow}\phi_{\beta^\downarrow}\right].}
$$

##### $s_r$

The cross-group contractions in the two branches are

$$
(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_r\phi_{\alpha^\uparrow}^\pm=\frac{\phi_{\beta^\downarrow}\phi_{\alpha^\uparrow}}{2\pi}e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi},
$$

$$
(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_r\phi_{\alpha^\downarrow}^\pm=\frac{\phi_{\beta^\uparrow}\phi_{\alpha^\downarrow}}{2\pi}e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi}.
$$

Same-group contractions vanish. Within a block,

$$
\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}=-1,\qquad \Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}=1,
$$

so all phases equal one. Both branches retain the same sign and combine through $R^{(+)}$:

$$
\boxed{s_r(z,r_\perp)=\frac{1}{2\pi}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\downarrow}R_{\alpha^\uparrow,\beta^\downarrow}^{(+)}\phi_{\alpha^\uparrow}\phi_{\beta^\downarrow}+\sum_{\alpha^\downarrow,\beta^\uparrow}R_{\alpha^\downarrow,\beta^\uparrow}^{(+)}\phi_{\alpha^\downarrow}\phi_{\beta^\uparrow}\right].}
$$

##### $s_\varphi$

The cross-group contractions are

$$
(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_\varphi\phi_{\alpha^\uparrow}^\pm=\mp\frac{\phi_{\beta^\downarrow}\phi_{\alpha^\uparrow}}{2\pi i}e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi},
$$

$$
(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_\varphi\phi_{\alpha^\downarrow}^\pm=\pm\frac{\phi_{\beta^\uparrow}\phi_{\alpha^\downarrow}}{2\pi i}e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi}.
$$

Same-group contractions vanish. The phase relations are the same as for $s_r$, but the negative branch reverses the contraction sign, giving $R^{(-)}$:

$$
\boxed{s_\varphi(z,r_\perp)=\frac{1}{2\pi i}\sum_{\text{blocks}}\left[\sum_{\alpha^\downarrow,\beta^\uparrow}R_{\alpha^\downarrow,\beta^\uparrow}^{(-)}\phi_{\alpha^\downarrow}\phi_{\beta^\uparrow}-\sum_{\alpha^\uparrow,\beta^\downarrow}R_{\alpha^\uparrow,\beta^\downarrow}^{(-)}\phi_{\alpha^\uparrow}\phi_{\beta^\downarrow}\right].}
$$

Hermiticity makes all three components real; for real density matrices, $s_\varphi=0$. The physical spin-angular-momentum density is $\hbar\mathbf s/2$.

### Spin-Current Tensor

The first index denotes the flow direction; the second denotes the spin direction. Define

$$
J_{\mu\nu}(\mathbf r)=\frac{1}{2i}\sum_{\alpha,\beta}\left\{\rho_{\alpha\beta}^{++}\left[(\phi_\beta^+)^\dagger\sigma_\nu\nabla_\mu\phi_\alpha^+-(\nabla_\mu(\phi_\beta^+)^\dagger)\sigma_\nu\phi_\alpha^+\right]+\rho_{\alpha\beta}^{--}\left[(\phi_\beta^-)^\dagger\sigma_\nu\nabla_\mu\phi_\alpha^--(\nabla_\mu(\phi_\beta^-)^\dagger)\sigma_\nu\phi_\alpha^-\right]\right\},
$$

where

$$
\nabla_z=\partial_z,\qquad \nabla_r=\partial_{r_\perp},\qquad \nabla_\varphi=\frac{1}{r_\perp}\partial_\varphi.
$$

All basis spinors are evaluated at $\mathbf r$. Derivatives act on the basis spinors, not on the Pauli matrices. Below, the upper and lower signs refer to $\phi^+$ and $\phi^-$ with the same representative indices; arrows always denote representative spins.

##### $J_{\varphi z}$

The equal-spin contractions are

$$
\frac{(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_z\nabla_\varphi\phi_{\alpha^\uparrow}^\pm-(\nabla_\varphi(\phi_{\beta^\uparrow}^\pm)^\dagger)\sigma_z\phi_{\alpha^\uparrow}^\pm}{2i}=\frac{\Lambda_{\alpha^\uparrow}+\Lambda_{\beta^\uparrow}}{4\pi r_\perp}\phi_{\beta^\uparrow}\phi_{\alpha^\uparrow}e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\uparrow})\varphi}.
$$

$$
\frac{(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_z\nabla_\varphi\phi_{\alpha^\downarrow}^\pm-(\nabla_\varphi(\phi_{\beta^\downarrow}^\pm)^\dagger)\sigma_z\phi_{\alpha^\downarrow}^\pm}{2i}=-\frac{\Lambda_{\alpha^\downarrow}+\Lambda_{\beta^\downarrow}}{4\pi r_\perp}\phi_{\beta^\downarrow}\phi_{\alpha^\downarrow}e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\downarrow})\varphi}.
$$

Cross-group terms vanish. Within each spin group, $\Lambda_\alpha=\Lambda_\beta$, so the phases equal one. Both spin and orbital projections reverse in the partner branch, leaving their product unchanged. Using [Spin Groups](#spin-groups),

$$
\boxed{J_{\varphi z}(z,r_\perp)=\frac{1}{2\pi r_\perp}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(+)}\Lambda_{\alpha^\uparrow}\phi_{\alpha^\uparrow}\phi_{\beta^\uparrow}-\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(+)}\Lambda_{\alpha^\downarrow}\phi_{\alpha^\downarrow}\phi_{\beta^\downarrow}\right].}
$$

On the right-hand sides, all spatial amplitudes are evaluated at $(z,r_\perp)$ and carry positive-branch labels.

##### $J_{z\varphi}$

The opposite-spin contractions are

$$
\frac{(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_\varphi\nabla_z\phi_{\alpha^\uparrow}^\pm-(\nabla_z(\phi_{\beta^\downarrow}^\pm)^\dagger)\sigma_\varphi\phi_{\alpha^\uparrow}^\pm}{2i}=\pm\frac{1}{4\pi}\left[\phi_{\beta^\downarrow}\partial_z\phi_{\alpha^\uparrow}-(\partial_z\phi_{\beta^\downarrow})\phi_{\alpha^\uparrow}\right]e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi}.
$$

$$
\frac{(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_\varphi\nabla_z\phi_{\alpha^\downarrow}^\pm-(\nabla_z(\phi_{\beta^\uparrow}^\pm)^\dagger)\sigma_\varphi\phi_{\alpha^\downarrow}^\pm}{2i}=\mp\frac{1}{4\pi}\left[\phi_{\beta^\uparrow}\partial_z\phi_{\alpha^\downarrow}-(\partial_z\phi_{\beta^\uparrow})\phi_{\alpha^\downarrow}\right]e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi}.
$$

Equal-spin terms vanish. Within one positive-branch block,

$$
\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}=-1,\qquad \Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}=1,
$$

so both phases equal one. The partner branch reverses the $\sigma_\varphi$ sign while leaving the spatial derivatives unchanged. Therefore,

$$
\boxed{J_{z\varphi}(z,r_\perp)=\frac{1}{4\pi}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\downarrow}R_{\alpha^\uparrow,\beta^\downarrow}^{(-)}\left(\phi_{\beta^\downarrow}\partial_z\phi_{\alpha^\uparrow}-\phi_{\alpha^\uparrow}\partial_z\phi_{\beta^\downarrow}\right)-\sum_{\alpha^\downarrow,\beta^\uparrow}R_{\alpha^\downarrow,\beta^\uparrow}^{(-)}\left(\phi_{\beta^\uparrow}\partial_z\phi_{\alpha^\downarrow}-\phi_{\alpha^\downarrow}\partial_z\phi_{\beta^\uparrow}\right)\right].}
$$

On the right-hand sides, all spatial amplitudes are evaluated at $(z,r_\perp)$ and carry positive-branch labels.

##### $J_{r\varphi}$

The opposite-spin contractions are

$$
\frac{(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_\varphi\nabla_r\phi_{\alpha^\uparrow}^\pm-(\nabla_r(\phi_{\beta^\downarrow}^\pm)^\dagger)\sigma_\varphi\phi_{\alpha^\uparrow}^\pm}{2i}=\pm\frac{1}{4\pi}\left[\phi_{\beta^\downarrow}\partial_{r_\perp}\phi_{\alpha^\uparrow}-(\partial_{r_\perp}\phi_{\beta^\downarrow})\phi_{\alpha^\uparrow}\right]e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi}.
$$

$$
\frac{(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_\varphi\nabla_r\phi_{\alpha^\downarrow}^\pm-(\nabla_r(\phi_{\beta^\uparrow}^\pm)^\dagger)\sigma_\varphi\phi_{\alpha^\downarrow}^\pm}{2i}=\mp\frac{1}{4\pi}\left[\phi_{\beta^\uparrow}\partial_{r_\perp}\phi_{\alpha^\downarrow}-(\partial_{r_\perp}\phi_{\beta^\uparrow})\phi_{\alpha^\downarrow}\right]e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi}.
$$

Equal-spin terms vanish. Within one positive-branch block,

$$
\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}=-1,\qquad \Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}=1,
$$

so both phases equal one. The partner branch reverses the $\sigma_\varphi$ sign while leaving the spatial derivatives unchanged. Therefore,

$$
\boxed{J_{r\varphi}(z,r_\perp)=\frac{1}{4\pi}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\downarrow}R_{\alpha^\uparrow,\beta^\downarrow}^{(-)}\left(\phi_{\beta^\downarrow}\partial_{r_\perp}\phi_{\alpha^\uparrow}-\phi_{\alpha^\uparrow}\partial_{r_\perp}\phi_{\beta^\downarrow}\right)-\sum_{\alpha^\downarrow,\beta^\uparrow}R_{\alpha^\downarrow,\beta^\uparrow}^{(-)}\left(\phi_{\beta^\uparrow}\partial_{r_\perp}\phi_{\alpha^\downarrow}-\phi_{\alpha^\downarrow}\partial_{r_\perp}\phi_{\beta^\uparrow}\right)\right].}
$$

On the right-hand sides, all spatial amplitudes are evaluated at $(z,r_\perp)$ and carry positive-branch labels.

##### $J_{\varphi r}$

The opposite-spin contractions are

$$
\frac{(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_r\nabla_\varphi\phi_{\alpha^\uparrow}^\pm-(\nabla_\varphi(\phi_{\beta^\downarrow}^\pm)^\dagger)\sigma_r\phi_{\alpha^\uparrow}^\pm}{2i}=\pm\frac{\Lambda_{\alpha^\uparrow}+\Lambda_{\beta^\downarrow}}{4\pi r_\perp}\phi_{\beta^\downarrow}\phi_{\alpha^\uparrow}e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi}.
$$

$$
\frac{(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_r\nabla_\varphi\phi_{\alpha^\downarrow}^\pm-(\nabla_\varphi(\phi_{\beta^\uparrow}^\pm)^\dagger)\sigma_r\phi_{\alpha^\downarrow}^\pm}{2i}=\pm\frac{\Lambda_{\alpha^\downarrow}+\Lambda_{\beta^\uparrow}}{4\pi r_\perp}\phi_{\beta^\uparrow}\phi_{\alpha^\downarrow}e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi}.
$$

Equal-spin terms vanish. Within one positive-branch block,

$$
\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}=-1,\qquad \Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}=1,\qquad \Lambda_{\alpha^\uparrow}+\Lambda_{\beta^\downarrow}=\Lambda_{\alpha^\downarrow}+\Lambda_{\beta^\uparrow}=2\Omega.
$$

Thus both phases equal one. The partner branch reverses the orbital projections while preserving the $\sigma_r$ contraction. Therefore,

$$
\boxed{J_{\varphi r}(z,r_\perp)=\frac{1}{2\pi r_\perp}\sum_{\text{blocks}}\Omega\left[\sum_{\alpha^\uparrow,\beta^\downarrow}R_{\alpha^\uparrow,\beta^\downarrow}^{(-)}\phi_{\alpha^\uparrow}\phi_{\beta^\downarrow}+\sum_{\alpha^\downarrow,\beta^\uparrow}R_{\alpha^\downarrow,\beta^\uparrow}^{(-)}\phi_{\alpha^\downarrow}\phi_{\beta^\uparrow}\right].}
$$

Here $\Omega>0$ is the projection of the current positive-branch block. All spatial amplitudes are evaluated at $(z,r_\perp)$ and carry positive-branch labels.

##### $J_{rz}$

The equal-spin contractions are

$$
\frac{(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_z\nabla_r\phi_{\alpha^\uparrow}^\pm-(\nabla_r(\phi_{\beta^\uparrow}^\pm)^\dagger)\sigma_z\phi_{\alpha^\uparrow}^\pm}{2i}=\pm\frac{1}{4\pi i}\left[\phi_{\beta^\uparrow}\partial_{r_\perp}\phi_{\alpha^\uparrow}-(\partial_{r_\perp}\phi_{\beta^\uparrow})\phi_{\alpha^\uparrow}\right]e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\uparrow})\varphi}.
$$

$$
\frac{(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_z\nabla_r\phi_{\alpha^\downarrow}^\pm-(\nabla_r(\phi_{\beta^\downarrow}^\pm)^\dagger)\sigma_z\phi_{\alpha^\downarrow}^\pm}{2i}=\mp\frac{1}{4\pi i}\left[\phi_{\beta^\downarrow}\partial_{r_\perp}\phi_{\alpha^\downarrow}-(\partial_{r_\perp}\phi_{\beta^\downarrow})\phi_{\alpha^\downarrow}\right]e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\downarrow})\varphi}.
$$

Cross-group terms vanish. Within each spin group, $\Lambda_\alpha=\Lambda_\beta$, so both phases equal one. The partner branch reverses the $\sigma_z$ sign while leaving the spatial derivatives unchanged. Therefore,

$$
\boxed{J_{rz}(z,r_\perp)=\frac{1}{4\pi i}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(-)}\left(\phi_{\beta^\uparrow}\partial_{r_\perp}\phi_{\alpha^\uparrow}-\phi_{\alpha^\uparrow}\partial_{r_\perp}\phi_{\beta^\uparrow}\right)-\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(-)}\left(\phi_{\beta^\downarrow}\partial_{r_\perp}\phi_{\alpha^\downarrow}-\phi_{\alpha^\downarrow}\partial_{r_\perp}\phi_{\beta^\downarrow}\right)\right].}
$$

All spatial amplitudes are evaluated at $(z,r_\perp)$ and carry positive-branch labels. For real density matrices, the symmetric $R^{(-)}$ contracts with an antisymmetric derivative kernel, giving $J_{rz}=0$.

##### $J_{zr}$

The opposite-spin contractions are

$$
\frac{(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_r\nabla_z\phi_{\alpha^\uparrow}^\pm-(\nabla_z(\phi_{\beta^\downarrow}^\pm)^\dagger)\sigma_r\phi_{\alpha^\uparrow}^\pm}{2i}=\frac{1}{4\pi i}\left[\phi_{\beta^\downarrow}\partial_z\phi_{\alpha^\uparrow}-(\partial_z\phi_{\beta^\downarrow})\phi_{\alpha^\uparrow}\right]e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi}.
$$

$$
\frac{(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_r\nabla_z\phi_{\alpha^\downarrow}^\pm-(\nabla_z(\phi_{\beta^\uparrow}^\pm)^\dagger)\sigma_r\phi_{\alpha^\downarrow}^\pm}{2i}=\frac{1}{4\pi i}\left[\phi_{\beta^\uparrow}\partial_z\phi_{\alpha^\downarrow}-(\partial_z\phi_{\beta^\uparrow})\phi_{\alpha^\downarrow}\right]e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi}.
$$

Equal-spin terms vanish. Within one positive-branch block,

$$
\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}=-1,\qquad \Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}=1,
$$

so both phases equal one. The partner branch preserves the $\sigma_r$ contraction and the spatial derivatives. Therefore,

$$
\boxed{J_{zr}(z,r_\perp)=\frac{1}{4\pi i}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\downarrow}R_{\alpha^\uparrow,\beta^\downarrow}^{(+)}\left(\phi_{\beta^\downarrow}\partial_z\phi_{\alpha^\uparrow}-\phi_{\alpha^\uparrow}\partial_z\phi_{\beta^\downarrow}\right)+\sum_{\alpha^\downarrow,\beta^\uparrow}R_{\alpha^\downarrow,\beta^\uparrow}^{(+)}\left(\phi_{\beta^\uparrow}\partial_z\phi_{\alpha^\downarrow}-\phi_{\alpha^\downarrow}\partial_z\phi_{\beta^\uparrow}\right)\right].}
$$

All spatial amplitudes are evaluated at $(z,r_\perp)$ and carry positive-branch labels. For real density matrices, exchanging $\alpha$ and $\beta$ makes the two sums cancel, giving $J_{zr}=0$.

##### $J_{zz}$

The equal-spin contractions are

$$
\frac{(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_z\nabla_z\phi_{\alpha^\uparrow}^\pm-(\nabla_z(\phi_{\beta^\uparrow}^\pm)^\dagger)\sigma_z\phi_{\alpha^\uparrow}^\pm}{2i}=\pm\frac{1}{4\pi i}\left[\phi_{\beta^\uparrow}\partial_z\phi_{\alpha^\uparrow}-(\partial_z\phi_{\beta^\uparrow})\phi_{\alpha^\uparrow}\right]e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\uparrow})\varphi}.
$$

$$
\frac{(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_z\nabla_z\phi_{\alpha^\downarrow}^\pm-(\nabla_z(\phi_{\beta^\downarrow}^\pm)^\dagger)\sigma_z\phi_{\alpha^\downarrow}^\pm}{2i}=\mp\frac{1}{4\pi i}\left[\phi_{\beta^\downarrow}\partial_z\phi_{\alpha^\downarrow}-(\partial_z\phi_{\beta^\downarrow})\phi_{\alpha^\downarrow}\right]e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\downarrow})\varphi}.
$$

Cross-group terms vanish. Within each spin group, $\Lambda_\alpha=\Lambda_\beta$, so both phases equal one. The partner branch reverses the $\sigma_z$ sign while leaving the spatial derivatives unchanged. Therefore,

$$
\boxed{J_{zz}(z,r_\perp)=\frac{1}{4\pi i}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(-)}\left(\phi_{\beta^\uparrow}\partial_z\phi_{\alpha^\uparrow}-\phi_{\alpha^\uparrow}\partial_z\phi_{\beta^\uparrow}\right)-\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(-)}\left(\phi_{\beta^\downarrow}\partial_z\phi_{\alpha^\downarrow}-\phi_{\alpha^\downarrow}\partial_z\phi_{\beta^\downarrow}\right)\right].}
$$

All spatial amplitudes are evaluated at $(z,r_\perp)$ and carry positive-branch labels. For real density matrices, the symmetric $R^{(-)}$ contracts with an antisymmetric derivative kernel, giving $J_{zz}=0$.

##### $J_{rr}$

The opposite-spin contractions are

$$
\frac{(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_r\nabla_r\phi_{\alpha^\uparrow}^\pm-(\nabla_r(\phi_{\beta^\downarrow}^\pm)^\dagger)\sigma_r\phi_{\alpha^\uparrow}^\pm}{2i}=\frac{1}{4\pi i}\left[\phi_{\beta^\downarrow}\partial_{r_\perp}\phi_{\alpha^\uparrow}-(\partial_{r_\perp}\phi_{\beta^\downarrow})\phi_{\alpha^\uparrow}\right]e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi}.
$$

$$
\frac{(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_r\nabla_r\phi_{\alpha^\downarrow}^\pm-(\nabla_r(\phi_{\beta^\uparrow}^\pm)^\dagger)\sigma_r\phi_{\alpha^\downarrow}^\pm}{2i}=\frac{1}{4\pi i}\left[\phi_{\beta^\uparrow}\partial_{r_\perp}\phi_{\alpha^\downarrow}-(\partial_{r_\perp}\phi_{\beta^\uparrow})\phi_{\alpha^\downarrow}\right]e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi}.
$$

Equal-spin terms vanish. Within one positive-branch block,

$$
\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}=-1,\qquad \Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}=1,
$$

so both phases equal one. The partner branch preserves the $\sigma_r$ contraction and the spatial derivatives. Therefore,

$$
\boxed{J_{rr}(z,r_\perp)=\frac{1}{4\pi i}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\downarrow}R_{\alpha^\uparrow,\beta^\downarrow}^{(+)}\left(\phi_{\beta^\downarrow}\partial_{r_\perp}\phi_{\alpha^\uparrow}-\phi_{\alpha^\uparrow}\partial_{r_\perp}\phi_{\beta^\downarrow}\right)+\sum_{\alpha^\downarrow,\beta^\uparrow}R_{\alpha^\downarrow,\beta^\uparrow}^{(+)}\left(\phi_{\beta^\uparrow}\partial_{r_\perp}\phi_{\alpha^\downarrow}-\phi_{\alpha^\downarrow}\partial_{r_\perp}\phi_{\beta^\uparrow}\right)\right].}
$$

All spatial amplitudes are evaluated at $(z,r_\perp)$ and carry positive-branch labels. For real density matrices, exchanging $\alpha$ and $\beta$ makes the two sums cancel, giving $J_{rr}=0$.

##### $J_{\varphi\varphi}$

The opposite-spin contractions are

$$
\frac{(\phi_{\beta^\downarrow}^\pm)^\dagger\sigma_\varphi\nabla_\varphi\phi_{\alpha^\uparrow}^\pm-(\nabla_\varphi(\phi_{\beta^\downarrow}^\pm)^\dagger)\sigma_\varphi\phi_{\alpha^\uparrow}^\pm}{2i}=-\frac{\Lambda_{\alpha^\uparrow}+\Lambda_{\beta^\downarrow}}{4\pi i r_\perp}\phi_{\beta^\downarrow}\phi_{\alpha^\uparrow}e^{\pm i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi}.
$$

$$
\frac{(\phi_{\beta^\uparrow}^\pm)^\dagger\sigma_\varphi\nabla_\varphi\phi_{\alpha^\downarrow}^\pm-(\nabla_\varphi(\phi_{\beta^\uparrow}^\pm)^\dagger)\sigma_\varphi\phi_{\alpha^\downarrow}^\pm}{2i}=\frac{\Lambda_{\alpha^\downarrow}+\Lambda_{\beta^\uparrow}}{4\pi i r_\perp}\phi_{\beta^\uparrow}\phi_{\alpha^\downarrow}e^{\pm i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi}.
$$

Equal-spin terms vanish. Within one positive-branch block,

$$
\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}=-1,\qquad \Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}=1,\qquad \Lambda_{\alpha^\uparrow}+\Lambda_{\beta^\downarrow}=\Lambda_{\alpha^\downarrow}+\Lambda_{\beta^\uparrow}=2\Omega.
$$

Thus both phases equal one. The partner branch reverses both the orbital projections and the $\sigma_\varphi$ contraction, leaving their product unchanged. Therefore,

$$
\boxed{J_{\varphi\varphi}(z,r_\perp)=\frac{1}{2\pi i r_\perp}\sum_{\text{blocks}}\Omega\left[\sum_{\alpha^\downarrow,\beta^\uparrow}R_{\alpha^\downarrow,\beta^\uparrow}^{(+)}\phi_{\alpha^\downarrow}\phi_{\beta^\uparrow}-\sum_{\alpha^\uparrow,\beta^\downarrow}R_{\alpha^\uparrow,\beta^\downarrow}^{(+)}\phi_{\alpha^\uparrow}\phi_{\beta^\downarrow}\right].}
$$

Here $\Omega>0$ is the projection of the current positive-branch block. All spatial amplitudes are evaluated at $(z,r_\perp)$ and carry positive-branch labels. For real density matrices, exchanging $\alpha$ and $\beta$ makes the two sums cancel, giving $J_{\varphi\varphi}=0$.

### Spin-Current Vector

Define the spin-current vector by

$$
J_k(\mathbf r)=\sum_{\mu,\nu}\epsilon_{k\mu\nu}J_{\mu\nu}(\mathbf r),
$$

where $\epsilon_{k\mu\nu}$ is the Levi-Civita symbol in the right-handed cylindrical frame $(r,\varphi,z)$.

The components are

$$
\boxed{J_r=J_{\varphi z}-J_{z\varphi},\qquad J_\varphi=J_{zr}-J_{rz},\qquad J_z=J_{r\varphi}-J_{\varphi r}.}
$$

Thus,

$$
\mathbf J(\mathbf r)=J_r(z,r_\perp)\mathbf e_r+J_\varphi(z,r_\perp)\mathbf e_\varphi+J_z(z,r_\perp)\mathbf e_z.
$$

For real density matrices, $J_{zr}=J_{rz}=0$, hence $J_\varphi=0$.

### Divergence of the Spin-Current Vector

For Cartesian directions $k,\mu,\nu$,

$$
\nabla\cdot\mathbf J=\sum_{k,\mu,\nu}\epsilon_{k\mu\nu}\nabla_kJ_{\mu\nu}=\frac{1}{2i}\sum_{\text{blocks}}\sum_{\alpha,\beta}\sum_{k,\mu,\nu}\epsilon_{k\mu\nu}\left[\rho_{\alpha\beta}^{++}\nabla_k\left[(\phi_\beta^+)^\dagger\sigma_\nu\nabla_\mu\phi_\alpha^+-(\nabla_\mu(\phi_\beta^+)^\dagger)\sigma_\nu\phi_\alpha^+\right]+\rho_{\alpha\beta}^{--}\nabla_k\left[(\phi_\beta^-)^\dagger\sigma_\nu\nabla_\mu\phi_\alpha^--(\nabla_\mu(\phi_\beta^-)^\dagger)\sigma_\nu\phi_\alpha^-\right]\right].
$$

The second derivatives cancel under the antisymmetric contraction, giving

$$
\nabla\cdot\mathbf J=\frac{1}{2i}\sum_{\text{blocks}}\sum_{\alpha,\beta}\sum_{k,\mu,\nu}\epsilon_{k\mu\nu}\left[\rho_{\alpha\beta}^{++}\left[(\nabla_k(\phi_\beta^+)^\dagger)\sigma_\nu(\nabla_\mu\phi_\alpha^+)-(\nabla_\mu(\phi_\beta^+)^\dagger)\sigma_\nu(\nabla_k\phi_\alpha^+)\right]+\rho_{\alpha\beta}^{--}\left[(\nabla_k(\phi_\beta^-)^\dagger)\sigma_\nu(\nabla_\mu\phi_\alpha^-)-(\nabla_\mu(\phi_\beta^-)^\dagger)\sigma_\nu(\nabla_k\phi_\alpha^-)\right]\right].
$$

Using $\rho_{\alpha\beta}^{++*}=\rho_{\beta\alpha}^{++}$, $\rho_{\alpha\beta}^{--*}=\rho_{\beta\alpha}^{--}$, and $\sigma_\nu^\dagger=\sigma_\nu$, exchange $\alpha,\beta$ in each subtracted term:

$$
\nabla\cdot\mathbf J=\frac{1}{2i}\sum_{\text{blocks}}\sum_{k,\mu,\nu}\epsilon_{k\mu\nu}\left\{\sum_{\alpha,\beta}\rho_{\alpha\beta}^{++}(\nabla_k(\phi_\beta^+)^\dagger)\sigma_\nu(\nabla_\mu\phi_\alpha^+)-\left[\sum_{\alpha,\beta}\rho_{\alpha\beta}^{++}(\nabla_k(\phi_\beta^+)^\dagger)\sigma_\nu(\nabla_\mu\phi_\alpha^+)\right]^*+\sum_{\alpha,\beta}\rho_{\alpha\beta}^{--}(\nabla_k(\phi_\beta^-)^\dagger)\sigma_\nu(\nabla_\mu\phi_\alpha^-)-\left[\sum_{\alpha,\beta}\rho_{\alpha\beta}^{--}(\nabla_k(\phi_\beta^-)^\dagger)\sigma_\nu(\nabla_\mu\phi_\alpha^-)\right]^*\right\}.
$$

Thus,

$$
\nabla\cdot\mathbf J=\operatorname{Re}\left\{\frac{1}{i}\sum_{\text{blocks}}\sum_{\alpha,\beta}\sum_{k,\mu,\nu}\epsilon_{k\mu\nu}\left[\rho_{\alpha\beta}^{++}(\nabla_k(\phi_\beta^+)^\dagger)\sigma_\nu(\nabla_\mu\phi_\alpha^+)+\rho_{\alpha\beta}^{--}(\nabla_k(\phi_\beta^-)^\dagger)\sigma_\nu(\nabla_\mu\phi_\alpha^-)\right]\right\}.
$$

In the right-handed cylindrical frame $(r,\varphi,z)$, the six terms of each branch pair into three:

$$
\boxed{\nabla\cdot\mathbf J=2\operatorname{Re}\left\{\frac{1}{i}\sum_{\text{blocks}}\sum_{\alpha,\beta}\left[\rho_{\alpha\beta}^{++}\left[(\nabla_\varphi(\phi_\beta^+)^\dagger)\sigma_r(\nabla_z\phi_\alpha^+)+(\nabla_z(\phi_\beta^+)^\dagger)\sigma_\varphi(\nabla_r\phi_\alpha^+)+(\nabla_r(\phi_\beta^+)^\dagger)\sigma_z(\nabla_\varphi\phi_\alpha^+)\right]+\rho_{\alpha\beta}^{--}\left[(\nabla_\varphi(\phi_\beta^-)^\dagger)\sigma_r(\nabla_z\phi_\alpha^-)+(\nabla_z(\phi_\beta^-)^\dagger)\sigma_\varphi(\nabla_r\phi_\alpha^-)+(\nabla_r(\phi_\beta^-)^\dagger)\sigma_z(\nabla_\varphi\phi_\alpha^-)\right]\right]\right\}.}
$$

Here $\nabla_r=\partial_{r_\perp}$, $\nabla_\varphi=r_\perp^{-1}\partial_\varphi$, and $\nabla_z=\partial_z$. All derivatives act on the basis spinors. The two branches are independent; no assumption of real density matrices or time-reversal symmetry is required.

##### $(\nabla_\varphi\phi_\beta^\dagger)\sigma_r(\nabla_z\phi_\alpha)$

For positive-branch representatives of opposite spin, the contractions in the two branches are

$$
\frac{1}{i}(\nabla_\varphi(\phi_{\beta^\downarrow}^{+})^\dagger)\sigma_r(\nabla_z\phi_{\alpha^\uparrow}^{+})=-\frac{\Lambda_{\beta^\downarrow}}{2\pi r_\perp}\phi_{\beta^\downarrow}\partial_z\phi_{\alpha^\uparrow}e^{i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi},\qquad \frac{1}{i}(\nabla_\varphi(\phi_{\beta^\downarrow}^{-})^\dagger)\sigma_r(\nabla_z\phi_{\alpha^\uparrow}^{-})=\frac{\Lambda_{\beta^\downarrow}}{2\pi r_\perp}\phi_{\beta^\downarrow}\partial_z\phi_{\alpha^\uparrow}e^{-i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi}.
$$

$$
\frac{1}{i}(\nabla_\varphi(\phi_{\beta^\uparrow}^{+})^\dagger)\sigma_r(\nabla_z\phi_{\alpha^\downarrow}^{+})=-\frac{\Lambda_{\beta^\uparrow}}{2\pi r_\perp}\phi_{\beta^\uparrow}\partial_z\phi_{\alpha^\downarrow}e^{i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi},\qquad \frac{1}{i}(\nabla_\varphi(\phi_{\beta^\uparrow}^{-})^\dagger)\sigma_r(\nabla_z\phi_{\alpha^\downarrow}^{-})=\frac{\Lambda_{\beta^\uparrow}}{2\pi r_\perp}\phi_{\beta^\uparrow}\partial_z\phi_{\alpha^\downarrow}e^{-i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi}.
$$

Equal-spin representative pairs vanish. Within one block,

$$
\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}=-1,\qquad \Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}=1.
$$

All four phases equal one. Reversing the branch changes the sign of the angular derivative while preserving the $\sigma_r$ spin factor. The contributions therefore combine through $R^{(-)}$:

$$
\boxed{\frac{1}{i}\sum_{\text{blocks}}\sum_{\alpha,\beta}\left[\rho_{\alpha\beta}^{++}(\nabla_\varphi(\phi_\beta^+)^\dagger)\sigma_r(\nabla_z\phi_\alpha^+)+\rho_{\alpha\beta}^{--}(\nabla_\varphi(\phi_\beta^-)^\dagger)\sigma_r(\nabla_z\phi_\alpha^-)\right]=-\frac{1}{2\pi r_\perp}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\downarrow}R_{\alpha^\uparrow,\beta^\downarrow}^{(-)}\Lambda_{\beta^\downarrow}\phi_{\beta^\downarrow}\partial_z\phi_{\alpha^\uparrow}+\sum_{\alpha^\downarrow,\beta^\uparrow}R_{\alpha^\downarrow,\beta^\uparrow}^{(-)}\Lambda_{\beta^\uparrow}\phi_{\beta^\uparrow}\partial_z\phi_{\alpha^\downarrow}\right].}
$$

All spatial amplitudes are real, evaluated at $(z,r_\perp)$, and carry positive-branch labels. The total divergence takes twice the real part of this contribution.

##### $(\nabla_z\phi_\beta^\dagger)\sigma_\varphi(\nabla_r\phi_\alpha)$

For positive-branch representatives of opposite spin, the contractions in the two branches are

$$
\frac{1}{i}(\nabla_z(\phi_{\beta^\downarrow}^{+})^\dagger)\sigma_\varphi(\nabla_r\phi_{\alpha^\uparrow}^{+})=\frac{1}{2\pi}(\partial_z\phi_{\beta^\downarrow})(\partial_{r_\perp}\phi_{\alpha^\uparrow})e^{i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi},\qquad \frac{1}{i}(\nabla_z(\phi_{\beta^\downarrow}^{-})^\dagger)\sigma_\varphi(\nabla_r\phi_{\alpha^\uparrow}^{-})=-\frac{1}{2\pi}(\partial_z\phi_{\beta^\downarrow})(\partial_{r_\perp}\phi_{\alpha^\uparrow})e^{-i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}+1)\varphi}.
$$

$$
\frac{1}{i}(\nabla_z(\phi_{\beta^\uparrow}^{+})^\dagger)\sigma_\varphi(\nabla_r\phi_{\alpha^\downarrow}^{+})=-\frac{1}{2\pi}(\partial_z\phi_{\beta^\uparrow})(\partial_{r_\perp}\phi_{\alpha^\downarrow})e^{i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi},\qquad \frac{1}{i}(\nabla_z(\phi_{\beta^\uparrow}^{-})^\dagger)\sigma_\varphi(\nabla_r\phi_{\alpha^\downarrow}^{-})=\frac{1}{2\pi}(\partial_z\phi_{\beta^\uparrow})(\partial_{r_\perp}\phi_{\alpha^\downarrow})e^{-i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}-1)\varphi}.
$$

Equal-spin representative pairs vanish. Within one block,

$$
\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\downarrow}=-1,\qquad \Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\uparrow}=1.
$$

All four phases equal one. Reversing the branch changes the sign of the $\sigma_\varphi$ spin factor while preserving the spatial derivatives. The contributions therefore combine through $R^{(-)}$:

$$
\boxed{\frac{1}{i}\sum_{\text{blocks}}\sum_{\alpha,\beta}\left[\rho_{\alpha\beta}^{++}(\nabla_z(\phi_\beta^+)^\dagger)\sigma_\varphi(\nabla_r\phi_\alpha^+)+\rho_{\alpha\beta}^{--}(\nabla_z(\phi_\beta^-)^\dagger)\sigma_\varphi(\nabla_r\phi_\alpha^-)\right]=\frac{1}{2\pi}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\downarrow}R_{\alpha^\uparrow,\beta^\downarrow}^{(-)}(\partial_z\phi_{\beta^\downarrow})(\partial_{r_\perp}\phi_{\alpha^\uparrow})-\sum_{\alpha^\downarrow,\beta^\uparrow}R_{\alpha^\downarrow,\beta^\uparrow}^{(-)}(\partial_z\phi_{\beta^\uparrow})(\partial_{r_\perp}\phi_{\alpha^\downarrow})\right].}
$$

All spatial amplitudes are real, evaluated at $(z,r_\perp)$, and carry positive-branch labels. The total divergence takes twice the real part of this contribution.

##### $(\nabla_r\phi_\beta^\dagger)\sigma_z(\nabla_\varphi\phi_\alpha)$

For positive-branch representatives of equal spin, the contractions in the two branches are

$$
\frac{1}{i}(\nabla_r(\phi_{\beta^\uparrow}^{+})^\dagger)\sigma_z(\nabla_\varphi\phi_{\alpha^\uparrow}^{+})=\frac{\Lambda_{\alpha^\uparrow}}{2\pi r_\perp}(\partial_{r_\perp}\phi_{\beta^\uparrow})\phi_{\alpha^\uparrow}e^{i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\uparrow})\varphi},\qquad \frac{1}{i}(\nabla_r(\phi_{\beta^\uparrow}^{-})^\dagger)\sigma_z(\nabla_\varphi\phi_{\alpha^\uparrow}^{-})=\frac{\Lambda_{\alpha^\uparrow}}{2\pi r_\perp}(\partial_{r_\perp}\phi_{\beta^\uparrow})\phi_{\alpha^\uparrow}e^{-i(\Lambda_{\alpha^\uparrow}-\Lambda_{\beta^\uparrow})\varphi}.
$$

$$
\frac{1}{i}(\nabla_r(\phi_{\beta^\downarrow}^{+})^\dagger)\sigma_z(\nabla_\varphi\phi_{\alpha^\downarrow}^{+})=-\frac{\Lambda_{\alpha^\downarrow}}{2\pi r_\perp}(\partial_{r_\perp}\phi_{\beta^\downarrow})\phi_{\alpha^\downarrow}e^{i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\downarrow})\varphi},\qquad \frac{1}{i}(\nabla_r(\phi_{\beta^\downarrow}^{-})^\dagger)\sigma_z(\nabla_\varphi\phi_{\alpha^\downarrow}^{-})=-\frac{\Lambda_{\alpha^\downarrow}}{2\pi r_\perp}(\partial_{r_\perp}\phi_{\beta^\downarrow})\phi_{\alpha^\downarrow}e^{-i(\Lambda_{\alpha^\downarrow}-\Lambda_{\beta^\downarrow})\varphi}.
$$

Opposite-spin representative pairs vanish. Within each equal-spin group of a block, $\Lambda_\alpha=\Lambda_\beta$, so all four phases equal one. Reversing the branch changes the signs of both the angular derivative and the $\sigma_z$ spin factor. Their product is unchanged, and the contributions combine through $R^{(+)}$:

$$
\boxed{\frac{1}{i}\sum_{\text{blocks}}\sum_{\alpha,\beta}\left[\rho_{\alpha\beta}^{++}(\nabla_r(\phi_\beta^+)^\dagger)\sigma_z(\nabla_\varphi\phi_\alpha^+)+\rho_{\alpha\beta}^{--}(\nabla_r(\phi_\beta^-)^\dagger)\sigma_z(\nabla_\varphi\phi_\alpha^-)\right]=\frac{1}{2\pi r_\perp}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}R_{\alpha^\uparrow,\beta^\uparrow}^{(+)}\Lambda_{\alpha^\uparrow}(\partial_{r_\perp}\phi_{\beta^\uparrow})\phi_{\alpha^\uparrow}-\sum_{\alpha^\downarrow,\beta^\downarrow}R_{\alpha^\downarrow,\beta^\downarrow}^{(+)}\Lambda_{\alpha^\downarrow}(\partial_{r_\perp}\phi_{\beta^\downarrow})\phi_{\alpha^\downarrow}\right].}
$$

All spatial amplitudes are real, evaluated at $(z,r_\perp)$, and carry positive-branch labels. The total divergence takes twice the real part of this contribution.

### Pairing Density

Define the local pairing tensor by

$$
\kappa_{\sigma\sigma'}(\mathbf r)=\left\langle\hat\psi_{\sigma'}(\mathbf r)\hat\psi_\sigma(\mathbf r)\right\rangle.
$$

The spin indices remain independent. The basis spinors form tensor products, not inner products:

$$
\boldsymbol\kappa(\mathbf r)=\sum_{\sigma,\sigma'}\kappa_{\sigma\sigma'}(\mathbf r)\chi_\sigma\otimes\chi_{\sigma'}=\sum_{\alpha,\beta}\left[\kappa_{\alpha\beta}^{+-}\phi_\alpha^+(\mathbf r)\otimes\phi_\beta^-(\mathbf r)+\kappa_{\alpha\beta}^{-+}\phi_\alpha^-(\mathbf r)\otimes\phi_\beta^+(\mathbf r)\right].
$$

No complex conjugation is applied. The spatial amplitudes and angular factors multiply normally:

$$
\phi_\alpha^+(\mathbf r)\otimes\phi_\beta^-(\mathbf r)=\frac{\phi_\alpha\phi_\beta}{2\pi}e^{i(\Lambda_\alpha-\Lambda_\beta)\varphi}\chi_{\Sigma_\alpha}\otimes\chi_{-\Sigma_\beta},
$$

$$
\phi_\alpha^-(\mathbf r)\otimes\phi_\beta^+(\mathbf r)=\frac{\phi_\alpha\phi_\beta}{2\pi}e^{-i(\Lambda_\alpha-\Lambda_\beta)\varphi}\chi_{-\Sigma_\alpha}\otimes\chi_{\Sigma_\beta}.
$$

Fermionic antisymmetry at the same spatial point gives

$$
\kappa_{\uparrow\uparrow}=\kappa_{\downarrow\downarrow}=0,\qquad \kappa_{\downarrow\uparrow}=-\kappa_{\uparrow\downarrow}.
$$

Define the independent scalar component as $\kappa(\mathbf r)=\kappa_{\uparrow\downarrow}(\mathbf r)$. Then

$$
\boldsymbol\kappa(\mathbf r)=\kappa(\mathbf r)\left(\chi_\uparrow\otimes\chi_\downarrow-\chi_\downarrow\otimes\chi_\uparrow\right).
$$

The $+-$ contribution to this component uses two up-spin representative labels; the $-+$ contribution uses two down-spin representative labels. Within each group, $\Lambda_\alpha=\Lambda_\beta$, so both azimuthal phases equal one:

$$
\kappa(z,r_\perp)=\frac{1}{2\pi}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}\kappa_{\alpha^\uparrow\beta^\uparrow}^{+-}\phi_{\alpha^\uparrow}\phi_{\beta^\uparrow}+\sum_{\alpha^\downarrow,\beta^\downarrow}\kappa_{\alpha^\downarrow\beta^\downarrow}^{-+}\phi_{\alpha^\downarrow}\phi_{\beta^\downarrow}\right].
$$

Antisymmetry gives

$$
\kappa_{\alpha^\downarrow\beta^\downarrow}^{-+}=-\kappa_{\beta^\downarrow\alpha^\downarrow}^{+-}.
$$

Exchanging the summation indices $\alpha,\beta$ in the second term, and using the commutativity of the scalar amplitudes, yields

$$
\boxed{\kappa(z,r_\perp)=\frac{1}{2\pi}\sum_{\text{blocks}}\left[\sum_{\alpha^\uparrow,\beta^\uparrow}\kappa_{\alpha^\uparrow\beta^\uparrow}^{+-}\phi_{\alpha^\uparrow}\phi_{\beta^\uparrow}-\sum_{\alpha^\downarrow,\beta^\downarrow}\kappa_{\alpha^\downarrow\beta^\downarrow}^{+-}\phi_{\alpha^\downarrow}\phi_{\beta^\downarrow}\right].}
$$

All spatial amplitudes are real and evaluated at $(z,r_\perp)$. Arrows on block-matrix indices denote representative spins. The scalar pairing density can be complex; no real-part operation is applied.
