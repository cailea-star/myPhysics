# R-Matrix Boundary Response

The R-matrix maps the internal interval $x\in[x_{\mathrm{in}},x_{\mathrm{out}}]$ to its boundary response. Define

$$
h_\mu=\frac{\hbar^2}{2\mu},\qquad T=-h_\mu\frac{d^2}{dx^2},\qquad H=T+V,
$$

where $h_\mu$ is the kinetic coefficient for the reduced mass $\mu$, and $V$ is Hermitian in channel space.

### Hermiticity on a Finite Interval

Since $V$ is Hermitian, the finite-interval residual comes from the kinetic term. Integration by parts gives

$$
\langle f|T|g\rangle=h_\mu\int_{x_{\mathrm{in}}}^{x_{\mathrm{out}}}f'^*(x)g'(x)\,dx-h_\mu\left[f^*(x)g'(x)\right]_{x_{\mathrm{in}}}^{x_{\mathrm{out}}}.
$$

Define the Wronskian at $x$ by

$$
W(x;f,g)=f^*(x)g'(x)-f'^*(x)g(x).
$$

The potential contribution vanishes:

$$
\langle f|V|g\rangle-\langle g|V|f\rangle^*=0.
$$

Therefore,

$$
\boxed{\langle f|H|g\rangle-\langle g|H|f\rangle^*=-h_\mu W(x_{\mathrm{out}};f,g)+h_\mu W(x_{\mathrm{in}};f,g)}.
$$

Thus, $H$ is not Hermitian on the finite interval when the boundary Wronskians are nonzero.

### Bloch Surface Operator

Assign the outward-normal signs

$$
n_{\mathrm{out}}=+1,\qquad n_{\mathrm{in}}=-1.
$$

For each endpoint $s\in\{x_{\mathrm{in}},x_{\mathrm{out}}\}$, define the Bloch surface operator

$$
\mathcal L_s(x)=n_s h_\mu\delta(x-s)\frac{d}{dx}.
$$

Its Hermitian residual is

$$
\langle f|\mathcal L_s|g\rangle-\langle g|\mathcal L_s|f\rangle^*=n_s h_\mu W(s;f,g).
$$

Summing over both endpoints cancels the boundary residual of $H$. Define

$$
\boxed{H_{\mathrm B}=H+\mathcal L_{\mathrm{in}}+\mathcal L_{\mathrm{out}}}.
$$

The kinetic contribution becomes

$$
\langle f|T+\mathcal L_{\mathrm{in}}+\mathcal L_{\mathrm{out}}|g\rangle=h_\mu\int_{x_{\mathrm{in}}}^{x_{\mathrm{out}}}f'^*(x)g'(x)\,dx.
$$

Thus, $H_{\mathrm B}$ is Hermitian on the finite interval without imposing vanishing boundary Wronskians.

### R-Matrix Response

The Schrödinger equation becomes

$$
(H-E)u=0\quad\Longrightarrow\quad(H_{\mathrm B}-E)u=(\mathcal L_{\mathrm{in}}+\mathcal L_{\mathrm{out}})u.
$$

Define the internal Green operator

$$
G(E)=(H_{\mathrm B}-E)^{-1},\qquad u=G(E)(\mathcal L_{\mathrm{in}}+\mathcal L_{\mathrm{out}})u.
$$

Let $\alpha,\beta$ label channels and $i,j$ label the non-orthogonal basis functions $B_i(x)$. The pairs $(\alpha,i)$ and $(\beta,j)$ are composite channel–basis indices.

Define the basis overlap matrix

$$
M_{ij}=\int_{x_{\mathrm{in}}}^{x_{\mathrm{out}}}B_i(x)B_j(x)\,dx.
$$

In the channel–basis space,

$$
\mathcal M_{\alpha i,\beta j}=\delta_{\alpha\beta}M_{ij},\qquad G_{\alpha i,\beta j}(E)=\left[H_{\mathrm B}-E\mathcal M\right]^{-1}_{\alpha i,\beta j}.
$$

The coordinate-space Green function is

$$
G_{\alpha\beta}(x,x';E)=\sum_{i,j}B_i(x)G_{\alpha i,\beta j}(E)B_j(x').
$$

The surface delta functions give

$$
u_\alpha(x)=\sum_{s\in\{x_{\mathrm{in}},x_{\mathrm{out}}\}}\sum_\beta n_s h_\mu G_{\alpha\beta}(x,s;E)u'_\beta(s).
$$

Define the R-matrix response kernel

$$
\boxed{R_{\alpha\beta}(x,s;E)=h_\mu G_{\alpha\beta}(x,s;E)}.
$$

Therefore,

$$
\boxed{u_\alpha(x)=\sum_{s\in\{x_{\mathrm{in}},x_{\mathrm{out}}\}}\sum_\beta n_sR_{\alpha\beta}(x,s;E)u'_\beta(s)}.
$$

### Inner-Boundary Self-Energy

Impose the logarithmic derivative directly at the inner boundary:

$$
u'(x_{\mathrm{in}})=\Lambda_{\mathrm{in}}u(x_{\mathrm{in}}).
$$

Since $n_{\mathrm{in}}=-1$, the inner Bloch source becomes

$$
\mathcal L_{\mathrm{in}}(x)u(x)=-h_\mu\Lambda_{\mathrm{in}}\delta(x-x_{\mathrm{in}})u(x_{\mathrm{in}}).
$$

Define the inner-boundary self-energy

$$
\boxed{\Sigma_{\mathrm{in}}(x)=-h_\mu\Lambda_{\mathrm{in}}\delta(x-x_{\mathrm{in}})}.
$$

Moving this term to the left leaves only the outer-boundary source:

$$
(H_{\mathrm B}-\Sigma_{\mathrm{in}}-E)u=\mathcal L_{\mathrm{out}}u.
$$

Define the modified Green operator and R-matrix response:

$$
G_\Sigma(E)=(H_{\mathrm B}-\Sigma_{\mathrm{in}}-E)^{-1},\qquad R_{\alpha\beta}(x,x_{\mathrm{out}};E)=h_\mu G_{\Sigma,\alpha\beta}(x,x_{\mathrm{out}};E).
$$

Therefore,

$$
\boxed{u_\alpha(x)=\sum_\beta R_{\alpha\beta}(x,x_{\mathrm{out}};E)u'_\beta(x_{\mathrm{out}})}.
$$

A complex absorbing choice of $\Lambda_{\mathrm{in}}$ makes the effective internal Hamiltonian non-Hermitian and, for flux-normalized open channels,

$$
\hat S^\dagger\hat S\leq I.
$$

### Outer-Boundary S-Matrix Matching

Collect the channel solutions into $\hat u$. At the outer boundary,

$$
\hat R(E)=\left[R_{\alpha\beta}(x_{\mathrm{out}},x_{\mathrm{out}};E)\right],\qquad \hat u(x_{\mathrm{out}})=\hat R(E)\hat u'(x_{\mathrm{out}}).
$$

The partial-wave convention relating this boundary form to $S_l=e^{2i\delta_l}$ and the scattering amplitude is given in [Phase-Shift Calculation](../../notes/quantum_mechanics/scattering_theory.md#phase-shift-calculation).

Let $\hat I(x)$ and $\hat O(x)$ denote the incoming- and outgoing-wave matrices:

$$
\hat u(x)=\hat I(x)-\hat O(x)\hat S,\qquad \hat u'(x)=\hat I'(x)-\hat O'(x)\hat S.
$$

Matching at $x_{\mathrm{out}}$ gives

$$
\hat I-\hat R\hat I'=\left(\hat O-\hat R\hat O'\right)\hat S.
$$

Therefore,

$$
\boxed{\hat S=\left(\hat O-\hat R\hat O'\right)^{-1}\left(\hat I-\hat R\hat I'\right)},
$$

where all external-wave matrices are evaluated at $x_{\mathrm{out}}$.

### Inverse R Matrix by Schur Complement

Partition $G_\Sigma^{-1}$ into the internal block $I$ and the outer-boundary block $B$:

$$
G_\Sigma^{-1}=\begin{pmatrix}(G_\Sigma^{-1})_{II}&(G_\Sigma^{-1})_{IB}\\(G_\Sigma^{-1})_{BI}&(G_\Sigma^{-1})_{BB}\end{pmatrix}.
$$

For B-splines that retain the outer endpoint,

$$
R=h_\mu(G_\Sigma)_{BB}.
$$

The block-inverse identity gives

$$
(G_\Sigma)_{BB}^{-1}=(G_\Sigma^{-1})_{BB}-(G_\Sigma^{-1})_{BI}\left[(G_\Sigma^{-1})_{II}\right]^{-1}(G_\Sigma^{-1})_{IB}.
$$

Therefore,

$$
\boxed{R^{-1}=h_\mu^{-1}\left\{(G_\Sigma^{-1})_{BB}-(G_\Sigma^{-1})_{BI}\left[(G_\Sigma^{-1})_{II}\right]^{-1}(G_\Sigma^{-1})_{IB}\right\}}.
$$

The code uses $h_\mu=1$. `calc_Rinv_matrix` solves the internal-block linear system without explicitly constructing $\left[(G_\Sigma^{-1})_{II}\right]^{-1}$.
