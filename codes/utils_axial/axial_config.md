# Axial Harmonic-Oscillator Configuration

An axial harmonic-oscillator configuration specifies the coordinate scales, finite single-particle space, symmetry blocks, and quadrature orders. Only states with positive $\Omega$ are stored explicitly; their negative-$\Omega$ partners follow from time reversal.

### Coordinate Scales

The coordinate conventions follow [Axial Harmonic-Oscillator Basis](axial_basis.md). Let $z$ and $r$ denote the axial and radial coordinates. The oscillator lengths $b_z$ and $b_r$ define

$$
\zeta=\frac{z}{b_z},\qquad \eta=\left(\frac{r}{b_r}\right)^2.
$$

The corresponding measures are

$$
dz=b_z\,d\zeta,\qquad r\,dr=\frac{b_r^2}{2}\,d\eta.
$$

Let $b_0$ be the mean oscillator length and $\beta_{20}$ the quadrupole deformation. The deformed lengths and volume conservation are

$$
b_z=b_0\exp\left(\sqrt{\frac{5}{16\pi}}\beta_{20}\right),\qquad b_r=b_0\exp\left(-\frac12\sqrt{\frac{5}{16\pi}}\beta_{20}\right),\qquad b_r^2b_z=b_0^3.
$$

The inverse transformation is

$$
b_0=(b_r^2b_z)^{1/3},\qquad \beta_{20}=\frac{8}{3}\sqrt{\frac{\pi}{5}}\log\frac{b_z}{b_r}.
$$

### Single-Particle Labels

The harmonic-oscillator single-particle states follow [Axial Harmonic-Oscillator Hamiltonian](axial_basis.md#axial-harmonic-oscillator-hamiltonian). Let $n_z$ and $n_r$ be nonnegative axial and radial quantum numbers. Let $\Lambda$ denote the $z$ projection of orbital angular momentum and $\Sigma=\pm\frac12$ the $z$ projection of spin. A state is labeled by

$$
\alpha\equiv(n_z,n_r,\Lambda,\Sigma).
$$

Its derived quantum numbers are

$$
N_\alpha=n_{z,\alpha}+2n_{r,\alpha}+\Lambda_\alpha,\qquad \Omega_\alpha=\Lambda_\alpha+\Sigma_\alpha,\qquad \pi_\alpha=(-1)^{n_{z,\alpha}+\Lambda_\alpha}.
$$

To avoid half-integer comparisons, introduce

$$
\Omega_2=2\Omega=2\Lambda\pm1,\qquad \Lambda_\uparrow=\frac{\Omega_2-1}{2},\qquad \Lambda_\downarrow=\frac{\Omega_2+1}{2}.
$$

Only positive odd values $\Omega_2\in\{1,3,5,\ldots\}$ are retained explicitly. The labels are enumerated in the order

$$
\Omega\rightarrow\pi\rightarrow n_r\rightarrow n_z\rightarrow\Lambda_\pm\rightarrow\Sigma_\pm.
$$

The time-reversed partner is

$$
\bar\alpha=(n_{z,\alpha},n_{r,\alpha},-\Lambda_\alpha,-\Sigma_\alpha),\qquad \Omega_{\bar\alpha}=-\Omega_\alpha,\qquad \pi_{\bar\alpha}=\pi_\alpha.
$$

### Symmetry Blocks

Axial symmetry conserves $\Omega$. Without reflection symmetry, the single-particle space decomposes as

$$
\mathcal B_\Omega=\{\alpha:\Omega_\alpha=\Omega\}.
$$

When reflection symmetry is imposed, parity is also conserved:

$$
\mathcal B_{\Omega,\pi}=\{\alpha:\Omega_\alpha=\Omega,\ \pi_\alpha=\pi\}.
$$

The positive-$\Omega$ blocks contain one representative of each time-reversal pair. Reflection symmetry also permits subsequent coordinate-space quadrature on the half-space $z\geq0$.

### Energy Cutoff

The single-particle spectrum follows [Axial Harmonic-Oscillator Hamiltonian](axial_basis.md#axial-harmonic-oscillator-hamiltonian). Let $m$ be the particle mass. The axial and transverse oscillator energies are

$$
\hbar\omega_z=\frac{\hbar^2}{mb_z^2},\qquad \hbar\omega_r=\frac{\hbar^2}{mb_r^2}.
$$

Define the transverse oscillator quantum number $\nu$ and the single-particle energy by

$$
\nu=2n_r+\Lambda,\qquad E(n_z,\nu)=\hbar\omega_z\left(n_z+\frac12\right)+\hbar\omega_r(\nu+1).
$$

For fixed $\nu$, the allowed radial quantum numbers and their multiplicity are

$$
n_r=0,\ldots,\left\lfloor\frac{\nu}{2}\right\rfloor,\qquad d_\nu=\left\lfloor\frac{\nu}{2}\right\rfloor+1.
$$

Let $N_{\rm shell}$ be the largest retained spherical oscillator shell. The corresponding number of spatial states is

$$
N_{\rm HO}=\sum_{N=0}^{N_{\rm shell}}\frac{(N+1)(N+2)}{2}=\frac{(N_{\rm shell}+1)(N_{\rm shell}+2)(N_{\rm shell}+3)}{6}.
$$

Candidate pairs $(n_z,\nu)$ are sorted by energy. Their multiplicities $d_\nu$ are accumulated until the count exceeds $N_{\rm HO}/2$; the current energy defines $E_{\rm cut}$. Every valid positive-$\Omega$ state with $E\leq E_{\rm cut}$ is retained. In the spherical limit, this reproduces the ordering by $N=n_z+\nu$.

### Quadrature Orders

Let $n_z^{\max}$, $n_r^{\max}$, and $\Lambda^{\max}$ be the largest quantum numbers among the retained labels. An $N$-point Gauss-Hermite or Gauss-Laguerre rule integrates weighted polynomials through degree $2N-1$.

Let $\deg_v$ denote the polynomial degree of the interaction potential. For products of four basis functions, the maximal axial and transverse degrees are

$$
p_z=4n_z^{\max}+\deg_v,\qquad p_r=4n_r^{\max}+2\Lambda^{\max}+\deg_v.
$$

The quadrature orders satisfy

$$
2N_z-1\geq p_z,\qquad 2N_r-1\geq p_r.
$$

Using $\deg_v=14$ gives

$$
\boxed{N_z=2n_z^{\max}+8,\qquad N_r=2n_r^{\max}+\Lambda^{\max}+8.}
$$
