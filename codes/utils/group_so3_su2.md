# Spin Representation

Theory: [SO(3) and SU(2)](../../notes/group_theory/so3_su2.md).

For fixed angular momentum $J=0,\tfrac12,1,\ldots$, use the ordered basis $|J,-J\rangle,\ldots,|J,J\rangle$. Angles are in radians.

### Angular-Momentum Generators

The $z$ component is diagonal in this basis:

$$
\langle JM'|\hat J_z|JM\rangle=\hbar M\,\delta_{M'M}.
$$

The ladder operators connect adjacent magnetic substates:

$$
\langle JM'|\hat J_+|JM\rangle=\hbar\sqrt{(J-M)(J+M+1)}\,\delta_{M',M+1}.
$$

$$
\langle JM'|\hat J_-|JM\rangle=\hbar\sqrt{(J+M)(J-M+1)}\,\delta_{M',M-1}.
$$

Their definitions relate them to the Cartesian components:

$$
\hat J_+=\hat J_x+i\hat J_y,\qquad \hat J_-=\hat J_x-i\hat J_y.
$$

The remaining components therefore follow from

$$
\hat J_x=\frac{\hat J_++\hat J_-}{2},\qquad \hat J_y=\frac{\hat J_+-\hat J_-}{2i}.
$$

### Rotations

The $z$-$y$-$z$ Euler-angle convention defines the rotation operator:

$$
\boxed{\hat R(\alpha,\beta,\gamma)=e^{-i\alpha\hat J_z/\hbar}e^{-i\beta\hat J_y/\hbar}e^{-i\gamma\hat J_z/\hbar}.}
$$

Let $J_y$ denote the matrix of $\hat J_y$. Its orthonormal eigenvectors form $Q$, with eigenvalues $\lambda_a$:

$$
J_y=Q\operatorname{diag}(\lambda_a)Q^\dagger.
$$

This decomposition gives the small $d$ matrix without repeated diagonalization:

$$
d^J(\beta)=Q\operatorname{diag}(e^{-i\beta\lambda_a/\hbar})Q^\dagger.
$$

The diagonal $z$ rotations supply the remaining phase factors:

$$
\boxed{\langle JM'|\hat R(\alpha,\beta,\gamma)|JM\rangle=e^{-iM'\alpha}d^J_{M'M}(\beta)e^{-iM\gamma}.}
$$
