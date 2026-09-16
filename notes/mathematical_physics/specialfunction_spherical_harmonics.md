# Spherical Harmonics

The unit direction $\hat{\mathbf r}=\mathbf r/r$, solid-angle measure $\mathrm d\hat{\mathbf r}$, and orbital angular-momentum operators used below are defined in [Spherical Coordinates](coordinate_spherical.md).

### Scalar Spherical Harmonics

Let $l=0,1,\ldots$, $m=-l,\ldots,l$, and let $P_l(x)$ and $P_l^m(x)$ denote the [Legendre polynomial and associated Legendre function](specialfunction_associated_legendre_polynomials.md). In the Condon–Shortley convention, for $m\geq0$,

$$
\boxed{Y_{lm}(\hat{\mathbf r})=Y_{lm}(\theta,\phi)=\sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}\,P_l^m(\cos\theta)e^{im\phi}}.
$$

The spherical harmonics satisfy:

- Eigenvalue equations

  $$
  \hat{\mathbf L}^2Y_{lm}=\hbar^2l(l+1)Y_{lm},\qquad\hat L_zY_{lm}=\hbar mY_{lm}.
  $$

- Conjugation and parity

  $$
  Y_{l,-m}(\hat{\mathbf r})=(-1)^mY_{lm}^*(\hat{\mathbf r}),\qquad Y_{lm}(-\hat{\mathbf r})=(-1)^lY_{lm}(\hat{\mathbf r}).
  $$

- Orthonormality and completeness, with $\delta_{S^2}$ denoting the Dirac delta on the unit sphere with respect to $\mathrm d\hat{\mathbf r}$

  $$
  \int \mathrm d\hat{\mathbf r}\,Y_{lm}^*(\hat{\mathbf r})Y_{l'm'}(\hat{\mathbf r})=\delta_{ll'}\delta_{mm'},\qquad\sum_{l=0}^{\infty}\sum_{m=-l}^{l}Y_{lm}(\hat{\mathbf r})Y_{lm}^*(\hat{\mathbf r}')=\delta_{S^2}(\hat{\mathbf r},\hat{\mathbf r}').
  $$

- Ladder relations, with $\hat L_\pm=\hat L_x\pm i\hat L_y$

  $$
  \hat L_\pm Y_{lm}=\hbar\sqrt{l(l+1)-m(m\pm1)}\,Y_{l,m\pm1}.
  $$

- Addition theorem, where $\cos\gamma=\hat{\mathbf r}\cdot\hat{\mathbf r}'$

  $$
  \sum_{m=-l}^{l}Y_{lm}(\hat{\mathbf r})Y_{lm}^*(\hat{\mathbf r}')=\frac{2l+1}{4\pi}P_l(\cos\gamma).
  $$

- Gaunt integral

  $$
  \int \mathrm d\hat{\mathbf r}\,Y_{l_1m_1}Y_{l_2m_2}Y_{l_3m_3}=\sqrt{\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi}}\begin{pmatrix}l_1&l_2&l_3\\0&0&0\end{pmatrix}\begin{pmatrix}l_1&l_2&l_3\\m_1&m_2&m_3\end{pmatrix}.
  $$

  The matrices are Wigner $3j$ symbols. The integral vanishes unless $m_1+m_2+m_3=0$, $|l_1-l_2|\leq l_3\leq l_1+l_2$, and $l_1+l_2+l_3$ is even.

A scalar function can be expanded as

$$
f(r,\hat{\mathbf r})=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\frac{u_{lm}(r)}{r}Y_{lm}(\hat{\mathbf r}),\qquad u_{lm}(r)=r\int \mathrm d\hat{\mathbf r}\,Y_{lm}^*(\hat{\mathbf r})f(r,\hat{\mathbf r}).
$$

### Vector Spherical Harmonics

Let $\mathbf e_0$ and $\mathbf e_\pm$ denote the spherical vector basis,

$$
\mathbf e_0=\hat{\mathbf z},\qquad \mathbf e_\pm=\mp\frac{\hat{\mathbf x}\pm i\hat{\mathbf y}}{\sqrt2}.
$$

Coupling $Y_{l m_l}(\hat{\mathbf r})$ to this spin-$1$ basis defines

$$
\boxed{\mathbf Y_{jm}^{l}(\hat{\mathbf r})=\sum_{m_l=-l}^{l}\sum_{\mu=-1}^{1}\langle l,m_l;1,\mu|j,m\rangle Y_{l m_l}(\hat{\mathbf r})\mathbf e_\mu},
$$

where $\langle l,m_l;1,\mu|j,m\rangle$ is a Clebsch–Gordan coefficient, and $\mathbf e_\mu$ denotes $\mathbf e_0$, $\mathbf e_+$, or $\mathbf e_-$ according to $\mu=0,+1,-1$:

$$
j=|l-1|,\ldots,l+1,\qquad m=-j,\ldots,j.
$$

Let $\hat{\mathbf S}$ act on the vector index and define $\hat{\mathbf J}=\hat{\mathbf L}+\hat{\mathbf S}$. The vector spherical harmonics satisfy:

- Eigenvalue equations

  $$
  \hat{\mathbf J}^2\mathbf Y_{jm}^{l}=\hbar^2j(j+1)\mathbf Y_{jm}^{l},\qquad \hat J_z\mathbf Y_{jm}^{l}=\hbar m\mathbf Y_{jm}^{l}.
  $$

  $$
  \hat{\mathbf L}^2\mathbf Y_{jm}^{l}=\hbar^2l(l+1)\mathbf Y_{jm}^{l},\qquad \hat{\mathbf S}^2\mathbf Y_{jm}^{l}=2\hbar^2\mathbf Y_{jm}^{l}.
  $$

- Conjugation and parity

  $$
  \left[\mathbf Y_{jm}^{l}(\hat{\mathbf r})\right]^*=(-1)^{j+l+m+1}\mathbf Y_{j,-m}^{l}(\hat{\mathbf r}),\qquad \mathbf Y_{jm}^{l}(-\hat{\mathbf r})=(-1)^l\mathbf Y_{jm}^{l}(\hat{\mathbf r}).
  $$

  The second relation transforms the argument while keeping the Cartesian basis fixed. A polar vector acquires one additional minus sign under the parity operator.

- Orthonormality and completeness

  $$
  \int \mathrm d\hat{\mathbf r}\,\left[\mathbf Y_{j'm'}^{l'}(\hat{\mathbf r})\right]^*\cdot\mathbf Y_{jm}^{l}(\hat{\mathbf r})=\delta_{jj'}\delta_{ll'}\delta_{mm'}.
  $$

  $$
  \sum_{l=0}^{\infty}\sum_{j=|l-1|}^{l+1}\sum_{m=-j}^{j}\left[\mathbf Y_{jm}^{l}(\hat{\mathbf r})\right]_a\left[\mathbf Y_{jm}^{l}(\hat{\mathbf r}')\right]_b^*=\delta_{ab}\delta_{S^2}(\hat{\mathbf r},\hat{\mathbf r}').
  $$

- Ladder relations, with $\hat J_\pm=\hat J_x\pm i\hat J_y$

  $$
  \hat J_\pm\mathbf Y_{jm}^{l}=\hbar\sqrt{j(j+1)-m(m\pm1)}\,\mathbf Y_{j,m\pm1}^{l}.
  $$

### Tangential Vector Spherical Harmonics

For $j\geq1$, define the radial and tangential harmonics by

$$
\mathbf Y_{jm}^{(r)}=\hat{\mathbf r}Y_{jm},\qquad \boldsymbol\Psi_{jm}=\frac{\nabla_{\hat{\mathbf r}}Y_{jm}}{\sqrt{j(j+1)}},\qquad \boldsymbol\Phi_{jm}=\frac{\hat{\mathbf r}\times\nabla_{\hat{\mathbf r}}Y_{jm}}{\sqrt{j(j+1)}}.
$$

Here $\nabla_{\hat{\mathbf r}}$ is the surface gradient on the unit sphere, and $\Delta_{\hat{\mathbf r}}=\nabla_{\hat{\mathbf r}}^2$. All angular functions are evaluated at $\hat{\mathbf r}$. The radial harmonic also exists for $j=0$.

- Tangentiality and cross products

  $$
  \hat{\mathbf r}\cdot\boldsymbol\Psi_{jm}=\hat{\mathbf r}\cdot\boldsymbol\Phi_{jm}=0,\qquad \hat{\mathbf r}\times\boldsymbol\Psi_{jm}=\boldsymbol\Phi_{jm},\qquad \hat{\mathbf r}\times\boldsymbol\Phi_{jm}=-\boldsymbol\Psi_{jm}.
  $$

- Conjugation

  $$
  \left[\mathbf Y_{jm}^{(r)}\right]^*=(-1)^m\mathbf Y_{j,-m}^{(r)},\qquad \boldsymbol\Psi_{jm}^{*}=(-1)^m\boldsymbol\Psi_{j,-m},\qquad \boldsymbol\Phi_{jm}^{*}=(-1)^m\boldsymbol\Phi_{j,-m}.
  $$

- Orthonormality

  $$
  \int \mathrm d\hat{\mathbf r}\,\boldsymbol\Psi_{jm}^{*}\cdot\boldsymbol\Psi_{j'm'}=\int \mathrm d\hat{\mathbf r}\,\boldsymbol\Phi_{jm}^{*}\cdot\boldsymbol\Phi_{j'm'}=\delta_{jj'}\delta_{mm'},\qquad \int \mathrm d\hat{\mathbf r}\,\boldsymbol\Psi_{jm}^{*}\cdot\boldsymbol\Phi_{j'm'}=0.
  $$

  The radial harmonics are orthogonal to both tangential families and satisfy

  $$
  \int \mathrm d\hat{\mathbf r}\,\left[\mathbf Y_{jm}^{(r)}\right]^*\cdot\mathbf Y_{j'm'}^{(r)}=\delta_{jj'}\delta_{mm'}.
  $$

- Surface divergence and curl

  For a tangential vector field $\mathbf V$, define the scalar surface curl by

  $$
  \operatorname{curl}_{\hat{\mathbf r}}\mathbf V\equiv-\nabla_{\hat{\mathbf r}}\cdot(\hat{\mathbf r}\times\mathbf V).
  $$

  Using $\Delta_{\hat{\mathbf r}}Y_{jm}=-j(j+1)Y_{jm}$ gives

  $$
  \boxed{\nabla_{\hat{\mathbf r}}\cdot\boldsymbol\Psi_{jm}=-\sqrt{j(j+1)}Y_{jm},\qquad \nabla_{\hat{\mathbf r}}\cdot\boldsymbol\Phi_{jm}=0,}
  $$

  $$
  \boxed{\operatorname{curl}_{\hat{\mathbf r}}\boldsymbol\Psi_{jm}=0,\qquad \operatorname{curl}_{\hat{\mathbf r}}\boldsymbol\Phi_{jm}=-\sqrt{j(j+1)}Y_{jm}.}
  $$

  Thus $\boldsymbol\Psi_{jm}$ is the gradient-type, curl-free harmonic, while $\boldsymbol\Phi_{jm}$ is the rotational-type, divergence-free harmonic.

- Completeness and tangential expansion

  The two families form a complete orthonormal basis for square-integrable tangential vector fields:

  $$
  \mathbf F_{\mathrm t}(r,\hat{\mathbf r})=\sum_{j=1}^{\infty}\sum_{m=-j}^{j}\left[F_{jm}^{(\Psi)}(r)\boldsymbol\Psi_{jm}(\hat{\mathbf r})+F_{jm}^{(\Phi)}(r)\boldsymbol\Phi_{jm}(\hat{\mathbf r})\right],
  $$

  $$
  F_{jm}^{(\Psi)}(r)=\int \mathrm d\hat{\mathbf r}\,\boldsymbol\Psi_{jm}^{*}\cdot\mathbf F_{\mathrm t},\qquad F_{jm}^{(\Phi)}(r)=\int \mathrm d\hat{\mathbf r}\,\boldsymbol\Phi_{jm}^{*}\cdot\mathbf F_{\mathrm t}.
  $$

  For sufficiently smooth fields, spherical integration by parts gives

  $$
  F_{jm}^{(\Psi)}(r)=-\frac1{\sqrt{j(j+1)}}\int \mathrm d\hat{\mathbf r}\,Y_{jm}^{*}\nabla_{\hat{\mathbf r}}\cdot\mathbf F_{\mathrm t},\qquad F_{jm}^{(\Phi)}(r)=-\frac1{\sqrt{j(j+1)}}\int \mathrm d\hat{\mathbf r}\,Y_{jm}^{*}\operatorname{curl}_{\hat{\mathbf r}}\mathbf F_{\mathrm t}.
  $$

- Relation to the coupled basis

  For fixed $(j,m)$, the full vector field has equivalent coupled and differential basis expansions:

  $$
  \mathbf F_{jm}(r,\hat{\mathbf r})=\sum_{l=|j-1|}^{j+1}F_{jlm}(r)\mathbf Y_{jm}^{l}(\hat{\mathbf r})=F_{jm}^{(r)}(r)\mathbf Y_{jm}^{(r)}(\hat{\mathbf r})+F_{jm}^{(\Psi)}(r)\boldsymbol\Psi_{jm}(\hat{\mathbf r})+F_{jm}^{(\Phi)}(r)\boldsymbol\Phi_{jm}(\hat{\mathbf r}).
  $$

  For $j\geq1$, the two bases span the same three-dimensional space and are related by a unitary change of basis. For $j=0$, only $\mathbf Y_{00}^{(r)}$ remains.

### Spinor Spherical Harmonics

Let $\hat I_2$ denote the two-dimensional identity operator. The Pauli matrices are

$$
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad \sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad \sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$

$$
\hat S_a=\frac{\hbar}{2}\sigma_a,\qquad a\in\{x,y,z\}.
$$

Their $\sigma_z$ eigenvectors are

$$
\chi_{1/2}\equiv\chi_+=\begin{pmatrix}1\\0\end{pmatrix},\qquad \chi_{-1/2}\equiv\chi_-=\begin{pmatrix}0\\1\end{pmatrix},\qquad \sigma_z\chi_\pm=\pm\chi_\pm.
$$

The Pauli algebra is collected in [Spin-$1/2$ in Matrix Form](../quantum_mechanics/angular_momentum.md#spin-12-in-matrix-form). On two-component functions,

$$
\hat{\mathbf J}=\hat{\mathbf L}\otimes\hat I_2+\hat{\mathbf S}.
$$

Coupling $Y_{lm_l}(\hat{\mathbf r})$ to the spin basis defines

$$
\boxed{\mathcal Y_{jm}^{l}(\hat{\mathbf r})=\sum_{m_l=-l}^{l}\sum_{\mu=\pm1/2}\left\langle l,m_l;\frac12,\mu\middle|j,m\right\rangle Y_{lm_l}(\hat{\mathbf r})\chi_\mu}.
$$

Here $j$ and $m$ are the total angular-momentum quantum numbers:

$$
j=l\pm\frac12,\qquad m=-j,-j+1,\ldots,j,
$$

with only $j=\frac12$ allowed for $l=0$. Since $m=m_l+\mu$, the two components are explicitly

$$
\mathcal Y_{jm}^{l}=\begin{pmatrix}\left\langle l,m-\frac12;\frac12,\frac12\middle|j,m\right\rangle Y_{l,m-1/2}\\\left\langle l,m+\frac12;\frac12,-\frac12\middle|j,m\right\rangle Y_{l,m+1/2}\end{pmatrix}.
$$

In the Condon–Shortley convention, for $j=l+\frac12$,

$$
\mathcal Y_{l+1/2,m}^{l}=\begin{pmatrix}\sqrt{\frac{l+m+1/2}{2l+1}}\,Y_{l,m-1/2}\\\sqrt{\frac{l-m+1/2}{2l+1}}\,Y_{l,m+1/2}\end{pmatrix}.
$$

For $j=l-\frac12$ with $l\geq1$,

$$
\mathcal Y_{l-1/2,m}^{l}=\begin{pmatrix}-\sqrt{\frac{l-m+1/2}{2l+1}}\,Y_{l,m-1/2}\\\sqrt{\frac{l+m+1/2}{2l+1}}\,Y_{l,m+1/2}\end{pmatrix}.
$$

The spinor spherical harmonics satisfy:

- Eigenvalue equations

  $$
  \hat{\mathbf J}^2\mathcal Y_{jm}^{l}=\hbar^2j(j+1)\mathcal Y_{jm}^{l},\qquad \hat J_z\mathcal Y_{jm}^{l}=\hbar m\mathcal Y_{jm}^{l}.
  $$

  $$
  \hat{\mathbf L}^2\mathcal Y_{jm}^{l}=\hbar^2l(l+1)\mathcal Y_{jm}^{l},\qquad \hat{\mathbf S}^2\mathcal Y_{jm}^{l}=\frac34\hbar^2\mathcal Y_{jm}^{l}.
  $$

- Complex conjugation and parity

  $$
  \left[\mathcal Y_{jm}^{l}(\hat{\mathbf r})\right]^*=(-1)^{l+j-m}i\sigma_y\mathcal Y_{j,-m}^{l}(\hat{\mathbf r}),\qquad \mathcal Y_{jm}^{l}(-\hat{\mathbf r})=(-1)^l\mathcal Y_{jm}^{l}(\hat{\mathbf r}).
  $$

- Orthonormality and completeness

  $$
  \int \mathrm d\hat{\mathbf r}\,\left[\mathcal Y_{j'm'}^{l'}(\hat{\mathbf r})\right]^\dagger\mathcal Y_{jm}^{l}(\hat{\mathbf r})=\delta_{jj'}\delta_{ll'}\delta_{mm'}.
  $$

  $$
  \sum_{l=0}^{\infty}\sum_{j=|l-1/2|}^{l+1/2}\sum_{m=-j}^{j}\mathcal Y_{jm}^{l}(\hat{\mathbf r})\left[\mathcal Y_{jm}^{l}(\hat{\mathbf r}')\right]^\dagger=\hat I_2\delta_{S^2}(\hat{\mathbf r},\hat{\mathbf r}').
  $$

- Ladder relations

  $$
  \hat J_\pm\mathcal Y_{jm}^{l}=\hbar\sqrt{j(j+1)-m(m\pm1)}\,\mathcal Y_{j,m\pm1}^{l}.
  $$

A two-component field can be expanded as

$$
\psi(r,\hat{\mathbf r})=\sum_{l=0}^{\infty}\sum_{j=|l-1/2|}^{l+1/2}\sum_{m=-j}^{j}\frac{u_{jlm}(r)}{r}\mathcal Y_{jm}^{l}(\hat{\mathbf r}),\qquad u_{jlm}(r)=r\int \mathrm d\hat{\mathbf r}\,\left[\mathcal Y_{jm}^{l}(\hat{\mathbf r})\right]^\dagger\psi(r,\hat{\mathbf r}).
$$
