### Axial Harmonic-Oscillator Configuration

An axial harmonic-oscillator configuration specifies the coordinate scales, finite single-particle space, symmetry blocks, and quadrature orders. Only states with positive $\Omega$ are stored explicitly; their negative-$\Omega$ partners follow from time reversal.

##### Coordinate Scales

Let $z$ and $r_\perp$ denote the axial and transverse coordinates. The oscillator lengths $b_z$ and $b_\perp$ define

$$
\zeta=\frac{z}{b_z},\qquad \eta=\left(\frac{r_\perp}{b_\perp}\right)^2.
$$

The corresponding measures are

$$
dz=b_z\,d\zeta,\qquad r_\perp dr_\perp=\frac{b_\perp^2}{2}\,d\eta.
$$

Let $b_0$ be the mean oscillator length and $\beta_{20}$ the quadrupole deformation. Define

$$
\alpha=\sqrt{\frac{5}{16\pi}},\qquad q=\exp(3\alpha\beta_{20}).
$$

The deformed lengths and volume conservation are

$$
b_z=b_0q^{1/3},\qquad b_\perp=b_0q^{-1/6},\qquad b_\perp^2b_z=b_0^3.
$$

The inverse transformation is

$$
b_0=(b_\perp^2b_z)^{1/3},\qquad \beta_{20}=\frac{2}{3\alpha}\log\frac{b_z}{b_\perp}.
$$

##### Single-Particle Labels

Let $n_z$ and $n_r$ be nonnegative axial and radial quantum numbers. Let $\Lambda$ denote the $z$ projection of orbital angular momentum and $\Sigma=\pm\frac12$ the $z$ projection of spin.

Define

$$
N=n_z+2n_r+\Lambda,\qquad \Omega=\Lambda+\Sigma,\qquad \pi=(-1)^{n_z+\Lambda}.
$$

A single-particle state is labeled by

$$
\xi=(N,n_z,n_r,\Lambda,\Omega,\Sigma,\pi).
$$

To avoid half-integer comparisons, introduce

$$
\Omega_2=2\Omega=2\Lambda\pm1,\qquad \Lambda_\uparrow=\frac{\Omega_2-1}{2},\qquad \Lambda_\downarrow=\frac{\Omega_2+1}{2}.
$$

Only positive odd values $\Omega_2\in\{1,3,5,\ldots\}$ are retained explicitly. The labels are enumerated in the order

$$
\Omega\rightarrow\pi\rightarrow n_r\rightarrow n_z\rightarrow\Lambda_\pm\rightarrow\Sigma_\pm.
$$

Time reversal maps

$$
(\Lambda,\Sigma,\Omega)\rightarrow(-\Lambda,-\Sigma,-\Omega)
$$

without changing $N$, $n_z$, $n_r$, or $\pi$.

##### Symmetry Blocks

Axial symmetry conserves $\Omega$. Without reflection symmetry, the single-particle space decomposes as

$$
\mathcal B_\Omega=\{\xi:\Omega_\xi=\Omega\}.
$$

When reflection symmetry is imposed, parity is also conserved:

$$
\mathcal B_{\Omega,\pi}=\{\xi:\Omega_\xi=\Omega,\ \pi_\xi=\pi\}.
$$

The positive-$\Omega$ blocks contain one representative of each time-reversal pair. Reflection symmetry also permits subsequent coordinate-space quadrature on the half-space $z\geq0$.

##### Energy Cutoff

Let $m$ be the particle mass. The axial and transverse oscillator energies are

$$
\hbar\omega_z=\frac{\hbar^2}{mb_z^2},\qquad \hbar\omega_\perp=\frac{\hbar^2}{mb_\perp^2}.
$$

Define the transverse oscillator quantum number $n_\perp$ and the single-particle energy by

$$
n_\perp=2n_r+\Lambda,\qquad E(n_z,n_\perp)=\hbar\omega_z\left(n_z+\frac12\right)+\hbar\omega_\perp(n_\perp+1).
$$

For fixed $n_\perp$, the allowed radial quantum numbers and their multiplicity are

$$
n_r=0,\ldots,\left\lfloor\frac{n_\perp}{2}\right\rfloor,\qquad d_{n_\perp}=\left\lfloor\frac{n_\perp}{2}\right\rfloor+1.
$$

Let $N_{\rm shell}$ be the largest retained spherical oscillator shell. The corresponding number of spatial states is

$$
N_{\rm HO}=\sum_{N=0}^{N_{\rm shell}}\frac{(N+1)(N+2)}{2}=\frac{(N_{\rm shell}+1)(N_{\rm shell}+2)(N_{\rm shell}+3)}{6}.
$$

Candidate pairs $(n_z,n_\perp)$ are sorted by energy. Their multiplicities $d_{n_\perp}$ are accumulated until the count exceeds $N_{\rm HO}/2$; the current energy defines $E_{\rm cut}$. Every valid positive-$\Omega$ state with $E\leq E_{\rm cut}$ is retained. In the spherical limit, this reproduces the ordering by $N=n_z+n_\perp$.

##### Quadrature Orders

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
