# Spherical Harmonic-Oscillator Rotations

### Rotation of Spherical Harmonic-Oscillator States

The eigenstates of the spherical harmonic-oscillator Hamiltonian are

$$
\hat H=\frac{\hat{\boldsymbol p}^{\,2}}{2M}+\frac12M\omega^2r^2,\qquad \hat H\phi_{nljm}=E_{nl}\phi_{nljm}.
$$

Following [Spherical Harmonic-Oscillator Basis](spherical_basis.md), the solutions separate into radial functions and spinor spherical harmonics:

$$
\boxed{\phi_{nljm}(\boldsymbol r)=\phi_{nl}(r)\mathcal Y^l_{jm}(\theta,\varphi),\qquad E_{nl}=\hbar\omega\left(2n+l+\frac32\right).}
$$

For a rotation $\hat R(\Omega)$, where $\Omega$ denotes the Euler angles, the spinor spherical harmonics transform as

$$
\boxed{\hat R(\Omega)\mathcal Y^l_{jm}=\sum_{m'=-j}^{j}\mathcal Y^l_{jm'}D^j_{m'm}(\Omega).}
$$

The radial coordinate is unchanged, so

$$
\hat R(\Omega)\phi_{nljm}=\sum_{m'=-j}^{j}\phi_{nljm'}D^j_{m'm}(\Omega).
$$

Thus, rotations preserve $n,l,j$ and mix only the magnetic substates within each multiplet.

### Rotation Matrix Elements

The rotation operators about the $z$ and $y$ axes are

$$
\hat R_z(\theta)=e^{-i\theta\hat J_z/\hbar},\qquad \hat R_y(\beta)=e^{-i\beta\hat J_y/\hbar}.
$$

Since $\hat J_z\phi_{nljm}=\hbar m\phi_{nljm}$, the $z$ rotation is diagonal:

$$
\boxed{\langle nljm|\hat R_z(\theta)|n'l'j'm'\rangle=\delta_{nn'}\delta_{ll'}\delta_{jj'}\delta_{mm'}e^{-im\theta}.}
$$

The $y$ rotation mixes magnetic substates through the real Wigner small-$d$ matrix:

$$
\boxed{\langle nljm|\hat R_y(\beta)|n'l'j'm'\rangle=\delta_{nn'}\delta_{ll'}\delta_{jj'}d^j_{mm'}(\beta),\qquad d^j_{mm'}(\beta)=\langle jm|e^{-i\beta\hat J_y/\hbar}|jm'\rangle.}
$$

For the Euler-angle convention $\Omega=(\alpha,\beta,\gamma)$,

$$
\hat R(\Omega)=\hat R_z(\alpha)\hat R_y(\beta)\hat R_z(\gamma),\qquad D^j_{mm'}(\Omega)=e^{-im\alpha}d^j_{mm'}(\beta)e^{-im'\gamma}.
$$
