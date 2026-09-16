# Spherical Harmonic-Oscillator Setting

### Coordinate Scale

Let $r\geq0$ be the radial coordinate, $M$ the particle mass, and $\omega$ the oscillator angular frequency. The oscillator length $b$ defines the dimensionless radial variable $\eta$:

$$
b=\sqrt{\frac{\hbar}{M\omega}},\qquad \eta=\left(\frac{r}{b}\right)^2.
$$

The radial coordinate and volume measure become

$$
r=b\sqrt{\eta},\qquad r^2\,dr=\frac{b^3}{2}\eta^{1/2}\,d\eta.
$$

### Single-Particle Labels

Let $n=0,1,\ldots$ be the radial quantum number and $l=0,1,\ldots$ the orbital angular momentum. For spin $s=\frac12$, the total angular momentum $j$ and its projection $m$ satisfy

$$
j\in\left\{\left|l-\frac12\right|,l+\frac12\right\},\qquad m=-j,-j+1,\ldots,j.
$$

A single-particle state is labeled by

$$
\alpha\equiv(n,l,j,m).
$$

Its oscillator shell quantum number and parity are

$$
N_\alpha=2n_\alpha+l_\alpha,\qquad \pi_\alpha=(-1)^{l_\alpha}=(-1)^{N_\alpha}.
$$

An orbit $(n,l,j)$ contains $2j+1$ magnetic substates. Labels refer to one nucleon species; neutron and proton spaces are treated separately.

### Single-Particle Energies

The spherical harmonic-oscillator energy of a state $\alpha$ is

$$
E_\alpha=\hbar\omega\left(2n_\alpha+l_\alpha+\frac32\right)=\frac{\hbar^2}{Mb^2}\left(N_\alpha+\frac32\right).
$$

The energy depends only on $N_\alpha$; states with different $j$ or $m$ within the same shell are degenerate.

Let $N_{\mathrm{shell}}\geq0$ be the largest retained oscillator shell. The energy cutoff $E_{\mathrm{cut}}$ defines the finite single-particle space:

$$
E_{\mathrm{cut}}=\hbar\omega\left(N_{\mathrm{shell}}+\frac32\right),\qquad \mathcal B=\{\alpha:E_\alpha\leq E_{\mathrm{cut}}\}=\{\alpha:N_\alpha\leq N_{\mathrm{shell}}\}.
$$

### Shell Enumeration

Let $N\in\{N_1,N_2,\ldots\}$ denote the selected oscillator shells. The states within each shell form

$$
\mathcal B_N=\{\alpha:2n_\alpha+l_\alpha=N\}.
$$

Within each shell, enumerate the total angular momentum and its projection:

$$
j=\frac12,\frac32,\ldots,N+\frac12,\qquad m=-j,-j+1,\ldots,j.
$$

For each $(N,j)$, the orbital and radial quantum numbers satisfy

$$
l=j\pm\frac12,\qquad 0\leq l\leq N,\qquad (N-l)\bmod 2=0,\qquad n=\frac{N-l}{2}.
$$

The condition $(N-l)\bmod 2=0$ fixes $(n,l)$ for each $(N,j)$; then enumerate $m=-j,\ldots,j$.

$$
\boxed{N\rightarrow j(n,l)\rightarrow m.}
$$

Each shell has parity $\pi_N$ and dimension $d_N$, including all magnetic substates for one nucleon species:

$$
\pi_N=(-1)^N,\qquad d_N=\sum_{j=1/2}^{N+1/2}(2j+1)=(N+1)(N+2).
$$

### Time-Reversal Phases

With $\hat T=-i\sigma_yK$ and real Condon–Shortley Clebsch–Gordan coefficients, the spinor spherical harmonics are

$$
\mathcal Y^l_{jm}=\sum_{m_l,m_s}\langle l,m_l;\tfrac12,m_s|j,m\rangle Y_{lm_l}\chi_{m_s}.
$$

Using

$$
Y_{lm_l}^{*}=(-1)^{m_l}Y_{l,-m_l},\qquad \hat T\chi_{m_s}=(-1)^{1/2-m_s}\chi_{-m_s},
$$

the time-reversed spinor spherical harmonic becomes

$$
\hat T\mathcal Y^l_{jm}=(-1)^{1/2-m}\sum_{m_l,m_s}\langle l,m_l;\tfrac12,m_s|j,m\rangle Y_{l,-m_l}\chi_{-m_s}.
$$

The Clebsch–Gordan reversal identity is

$$
\langle l,m_l;\tfrac12,m_s|j,m\rangle=(-1)^{l+1/2-j}\langle l,-m_l;\tfrac12,-m_s|j,-m\rangle.
$$

Since $1-2j$ is even, the resulting phase is

$$
\boxed{\hat T\mathcal Y^l_{jm}=(-1)^{l+1-j-m}\mathcal Y^l_{j,-m}=(-1)^{l+j-m}\mathcal Y^l_{j,-m}.}
$$

For $\alpha=(n,l,j,m)$ and $\bar\alpha=(n,l,j,-m)$, the radial functions are real, so

$$
\hat T\phi_\alpha=\eta_\alpha\phi_{\bar\alpha},\qquad \eta_\alpha=(-1)^{l+j-m},\qquad \eta_\alpha\eta_{\bar\alpha}=-1.
$$

The phases $\eta_\alpha$ follow the single-particle ordering within each block. Under time-reversal reduction, only the $m>0$ representatives are retained.

### Quadrature Orders

Let $n^{\max}$ and $l^{\max}$ be the largest quantum numbers among the retained labels. Let $\deg_v$ denote the polynomial degree of the radial potential in $\eta=(r/b)^2$.

For products of four harmonic-oscillator radial functions at the same radius, the Gaussian factor is $e^{-2\eta}$. Introduce $x=2\eta$ and use generalized Gauss-Laguerre quadrature with weight

$$
w(x)=x^{1/2}e^{-x}.
$$

For parity-even products, the remaining polynomial degree is bounded by

$$
p_r=4n^{\max}+2l^{\max}+\deg_v.
$$

An $N_r$-point rule integrates weighted polynomials through degree $2N_r-1$. The quadrature order therefore satisfies

$$
2N_r-1\geq p_r.
$$

Using $\deg_v=14$ gives

$$
\boxed{N_r=2n^{\max}+l^{\max}+8.}
$$

Nonpolynomial interactions require numerical convergence checks.
