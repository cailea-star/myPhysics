# Multipole Interactions in Second Quantization

### Separable Decomposition of Two-Body Interactions

Consider a central two-body interaction depending only on the interparticle distance:

$$
\hat V=\sum_{i<j}v(|\mathbf r_i-\mathbf r_j|)=\frac12\sum_{i\ne j}v(|\mathbf r_i-\mathbf r_j|).
$$

Its multipole expansion is

$$
v(|\mathbf r_i-\mathbf r_j|)=\sum_{\lambda,\mu}\frac{4\pi}{2\lambda+1}v_\lambda(r_i,r_j)Y_{\lambda\mu}^*(\hat{\mathbf r}_i)Y_{\lambda\mu}(\hat{\mathbf r}_j).
$$

Make the separable approximation with real functions $f_\lambda$ and real coefficients $\chi_\lambda$:

$$
\frac{4\pi}{2\lambda+1}v_\lambda(r_i,r_j)\approx\chi_\lambda f_\lambda(r_i)f_\lambda(r_j).
$$

Define the single-particle multipole function and the corresponding many-body operator by

$$
q_{\lambda\mu}(\mathbf r)=f_\lambda(r)Y_{\lambda\mu}(\hat{\mathbf r}),\qquad \hat Q_{\lambda\mu}=\sum_i q_{\lambda\mu}(\mathbf r_i)=\sum_{\alpha,\beta}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\hat c_\alpha^\dagger\hat c_\beta.
$$

Here $\hat q_{\lambda\mu}$ is the single-particle multiplication operator associated with $q_{\lambda\mu}(\mathbf r)$. Substitution gives

$$
\hat V\approx\frac12\sum_{\lambda,\mu}\chi_\lambda\sum_{i\ne j}q_{\lambda\mu}^*(\mathbf r_i)q_{\lambda\mu}(\mathbf r_j).
$$

Separating the $i=j$ contribution,

$$
\hat Q_{\lambda\mu}^\dagger\hat Q_{\lambda\mu}=\sum_{i\ne j}q_{\lambda\mu}^*(\mathbf r_i)q_{\lambda\mu}(\mathbf r_j)+\sum_i|q_{\lambda\mu}(\mathbf r_i)|^2.
$$

Therefore,

$$
\boxed{\hat V\approx\frac12\sum_{\lambda,\mu}\chi_\lambda\left[\hat Q_{\lambda\mu}^\dagger\hat Q_{\lambda\mu}-\sum_i|q_{\lambda\mu}(\mathbf r_i)|^2\right]}.
$$

Using the spherical-harmonic addition theorem,

$$
\sum_{\mu=-\lambda}^{\lambda}|q_{\lambda\mu}(\mathbf r)|^2=\frac{2\lambda+1}{4\pi}f_\lambda(r)^2.
$$

Define the multipole interaction and its one-body correction by

$$
\hat H_\lambda=\frac{\chi_\lambda}{2}\sum_{\mu=-\lambda}^{\lambda}\hat Q_{\lambda\mu}^\dagger\hat Q_{\lambda\mu},\qquad \hat h_{\mathrm{corr}}=-\sum_\lambda\frac{\chi_\lambda(2\lambda+1)}{8\pi}f_\lambda(r)^2.
$$

Then

$$
\boxed{\hat V\approx\sum_\lambda\hat H_\lambda+\sum_{\alpha,\beta}\langle\alpha|\hat h_{\mathrm{corr}}|\beta\rangle\hat c_\alpha^\dagger\hat c_\beta}.
$$

The one-body correction removes the $i=j$ contribution contained in $\hat Q^\dagger\hat Q$. A general radial kernel requires a sum over separable modes; the approximation above retains one mode per multipole channel.

### Multipole Interactions and Matrix Elements

For the separable approximation and its one-body correction, see [Separable Decomposition of Two-Body Interactions](#separable-decomposition-of-two-body-interactions).

Let $\hat q_{\lambda\mu}$ be a single-particle Hermitian spherical tensor of integer rank $\lambda$, with real coupling strength $\chi_\lambda$. Define

$$
\hat H_\lambda=\frac{\chi_\lambda}{2}\sum_{\mu=-\lambda}^{\lambda}\hat Q_{\lambda\mu}^\dagger\hat Q_{\lambda\mu},\qquad \hat Q_{\lambda\mu}=\sum_{\alpha,\beta}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\hat c_\alpha^\dagger\hat c_\beta,\qquad \hat Q_{\lambda\mu}^\dagger=(-1)^\mu\hat Q_{\lambda,-\mu}.
$$

Preserving the operator order gives

$$
\hat H_\lambda=\frac{\chi_\lambda}{2}\sum_{\mu,\alpha,\beta,\gamma,\delta}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle^*\langle\gamma|\hat q_{\lambda\mu}|\delta\rangle\hat c_\beta^\dagger\hat c_\alpha\hat c_\gamma^\dagger\hat c_\delta.
$$

Let $a$ label an orbital, including $j_a$ but excluding its magnetic quantum number:

$$
\alpha=\{a,m_a\},\qquad \sum_\alpha=\sum_a\sum_{|m_a|\leq j_a}.
$$

Using the convention in [Wigner Symbols](../mathematical_physics/wigner_symbols.md#reduced-matrix-elements),

$$
\langle a m_a|\hat q_{\lambda\mu}|b m_b\rangle=(-1)^{j_a-m_a}\begin{pmatrix}j_a&\lambda&j_b\\-m_a&\mu&m_b\end{pmatrix}\langle a\Vert\hat q_\lambda\Vert b\rangle=\frac{(-1)^{2j_a-j_b-m_b}}{\sqrt{2\lambda+1}}\langle j_a m_a;j_b,-m_b|\lambda\mu\rangle\langle a\Vert\hat q_\lambda\Vert b\rangle.
$$

To express the particle–hole operators as spherical tensors, define

$$
\tilde c_{b m_b}\equiv(-1)^{j_b+m_b}\hat c_{b,-m_b}.
$$

Their angular-momentum coupling is

$$
[\hat c_a^\dagger\otimes\tilde c_b]_{\lambda\mu}\equiv\sum_{m_a,m_b}\langle j_a m_a;j_b m_b|\lambda\mu\rangle\hat c_{a m_a}^\dagger\tilde c_{b m_b}=\sum_{m_a,m_b}(-1)^{j_b-m_b}\langle j_a m_a;j_b,-m_b|\lambda\mu\rangle\hat c_{a m_a}^\dagger\hat c_{b m_b}.
$$

For integer $\lambda$, nonzero reduced matrix elements require $j_a-j_b\in\mathbb Z$. Substitution therefore gives

$$
\boxed{\hat Q_{\lambda\mu}=\sum_{\alpha,\beta}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\hat c_\alpha^\dagger\hat c_\beta=\sum_{a,b}\frac{\langle a\Vert\hat q_\lambda\Vert b\rangle}{\sqrt{2\lambda+1}}[\hat c_a^\dagger\otimes\tilde c_b]_{\lambda\mu}}.
$$

Let $|\Phi^I_{M,\kappa}\rangle$ denote normalized angular-momentum eigenstates in a fixed-particle-number space, with $\kappa$ distinguishing states of the same $I$. The matrix element between these states is

$$
\langle\Phi^{I_1}_{M_1,\kappa_1}|\hat Q_{\lambda\mu}|\Phi^{I_2}_{M_2,\kappa_2}\rangle=\sum_{\alpha,\beta}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\langle\Phi^{I_1}_{M_1,\kappa_1}|\hat c_\alpha^\dagger\hat c_\beta|\Phi^{I_2}_{M_2,\kappa_2}\rangle.
$$

The reduced matrix element is

$$
\boxed{\langle\Phi^{I_1}_{\kappa_1}\Vert\hat Q_\lambda\Vert\Phi^{I_2}_{\kappa_2}\rangle=\sum_{a,b}\frac{\langle a\Vert\hat q_\lambda\Vert b\rangle}{\sqrt{2\lambda+1}}\langle\Phi^{I_1}_{\kappa_1}\Vert[\hat c_a^\dagger\otimes\tilde c_b]_\lambda\Vert\Phi^{I_2}_{\kappa_2}\rangle}.
$$

The interaction matrix element is

$$
\langle\Phi^{I_1}_{M_1,\kappa_1}|\hat H_\lambda|\Phi^{I_2}_{M_2,\kappa_2}\rangle=\frac{\chi_\lambda}{2}\sum_{\mu,\alpha,\beta,\gamma,\delta}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle^*\langle\gamma|\hat q_{\lambda\mu}|\delta\rangle\langle\Phi^{I_1}_{M_1,\kappa_1}|\hat c_\beta^\dagger\hat c_\alpha\hat c_\gamma^\dagger\hat c_\delta|\Phi^{I_2}_{M_2,\kappa_2}\rangle.
$$

For its reduced form, define the scalar operator

$$
\hat X^{(0)}_{\lambda;ab,cd}\equiv\sum_{\mu=-\lambda}^{\lambda}\left([\hat c_a^\dagger\otimes\tilde c_b]_{\lambda\mu}\right)^\dagger[\hat c_c^\dagger\otimes\tilde c_d]_{\lambda\mu}=\sum_{\mu,m_a,m_b,m_c,m_d}(-1)^{j_b-m_b+j_d-m_d}\langle j_a m_a;j_b,-m_b|\lambda\mu\rangle\langle j_c m_c;j_d,-m_d|\lambda\mu\rangle\hat c_{b m_b}^\dagger\hat c_{a m_a}\hat c_{c m_c}^\dagger\hat c_{d m_d}.
$$

Then

$$
\boxed{\langle\Phi^{I_1}_{\kappa_1}\Vert\hat H_\lambda\Vert\Phi^{I_2}_{\kappa_2}\rangle=\frac{\chi_\lambda}{2(2\lambda+1)}\sum_{a,b,c,d}\langle a\Vert\hat q_\lambda\Vert b\rangle^*\langle c\Vert\hat q_\lambda\Vert d\rangle\langle\Phi^{I_1}_{\kappa_1}\Vert\hat X^{(0)}_{\lambda;ab,cd}\Vert\Phi^{I_2}_{\kappa_2}\rangle}.
$$

### Multipole Pairing Interactions and Matrix Elements

Let $\hat q_{\lambda\mu}$ be the single-particle spherical tensor defined in [Multipole Interactions and Matrix Elements](#multipole-interactions-and-matrix-elements). For real coupling strength $G_\lambda$, define

$$
\hat H_\lambda^{\mathrm{pair}}=-G_\lambda\sum_{\mu=-\lambda}^{\lambda}\hat P_{\lambda\mu}^\dagger\hat P_{\lambda\mu},\qquad \hat P_{\lambda\mu}^\dagger=\frac12\sum_{\alpha,\beta}(-1)^{j_b-m_b}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\hat c_{a m_a}^\dagger\hat c_{b,-m_b}^\dagger.
$$

Preserving the operator order gives

$$
\hat H_\lambda^{\mathrm{pair}}=-\frac{G_\lambda}{4}\sum_{\mu,\alpha,\beta,\gamma,\delta}(-1)^{j_b-m_b+j_d-m_d}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\langle\gamma|\hat q_{\lambda\mu}|\delta\rangle^*\hat c_{a m_a}^\dagger\hat c_{b,-m_b}^\dagger\hat c_{d,-m_d}\hat c_{c m_c}.
$$

Let $a$ label an orbital, including $j_a$ but excluding its magnetic quantum number:

$$
\alpha=\{a,m_a\},\qquad \sum_\alpha=\sum_a\sum_{|m_a|\leq j_a}.
$$

Using the convention in [Wigner Symbols](../mathematical_physics/wigner_symbols.md#reduced-matrix-elements),

$$
\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle=(-1)^{j_a-m_a}\begin{pmatrix}j_a&\lambda&j_b\\-m_a&\mu&m_b\end{pmatrix}\langle a\Vert\hat q_\lambda\Vert b\rangle=\frac{(-1)^{2j_a-j_b-m_b}}{\sqrt{2\lambda+1}}\langle j_a m_a;j_b,-m_b|\lambda\mu\rangle\langle a\Vert\hat q_\lambda\Vert b\rangle.
$$

The coupled pair-creation operator is

$$
[\hat c_a^\dagger\otimes\hat c_b^\dagger]_{\lambda\mu}\equiv\sum_{m_a,m_b}\langle j_a m_a;j_b,-m_b|\lambda\mu\rangle\hat c_{a m_a}^\dagger\hat c_{b,-m_b}^\dagger.
$$

Taking the Hermitian adjoint reverses the operator order:

$$
\left([\hat c_a^\dagger\otimes\hat c_b^\dagger]_{\lambda\mu}\right)^\dagger=\sum_{m_a,m_b}\langle j_a m_a;j_b,-m_b|\lambda\mu\rangle\hat c_{b,-m_b}\hat c_{a m_a}.
$$

For integer $\lambda$, substitution gives

$$
\boxed{\hat P_{\lambda\mu}^\dagger=\frac12\sum_{\alpha,\beta}(-1)^{j_b-m_b}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\hat c_{a m_a}^\dagger\hat c_{b,-m_b}^\dagger=\sum_{a,b}\frac{\langle a\Vert\hat q_\lambda\Vert b\rangle}{2\sqrt{2\lambda+1}}[\hat c_a^\dagger\otimes\hat c_b^\dagger]_{\lambda\mu}}.
$$

Let $|\Phi^I_{M,\kappa}\rangle$ denote normalized angular-momentum eigenstates, with particle numbers $N_1,N_2$ for the bra and ket. Particle-number labels are suppressed. The matrix element is

$$
\langle\Phi^{I_1}_{M_1,\kappa_1}|\hat P_{\lambda\mu}^\dagger|\Phi^{I_2}_{M_2,\kappa_2}\rangle=\frac12\sum_{\alpha,\beta}(-1)^{j_b-m_b}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\langle\Phi^{I_1}_{M_1,\kappa_1}|\hat c_{a m_a}^\dagger\hat c_{b,-m_b}^\dagger|\Phi^{I_2}_{M_2,\kappa_2}\rangle.
$$

The reduced matrix element is

$$
\boxed{\langle\Phi^{I_1}_{\kappa_1}\Vert\hat P_\lambda^\dagger\Vert\Phi^{I_2}_{\kappa_2}\rangle=\sum_{a,b}\frac{\langle a\Vert\hat q_\lambda\Vert b\rangle}{2\sqrt{2\lambda+1}}\langle\Phi^{I_1}_{\kappa_1}\Vert[\hat c_a^\dagger\otimes\hat c_b^\dagger]_\lambda\Vert\Phi^{I_2}_{\kappa_2}\rangle}.
$$

These pair-addition matrix elements vanish unless $N_1=N_2+2$.

The interaction matrix element is

$$
\langle\Phi^{I_1}_{M_1,\kappa_1}|\hat H_\lambda^{\mathrm{pair}}|\Phi^{I_2}_{M_2,\kappa_2}\rangle=-\frac{G_\lambda}{4}\sum_{\mu,\alpha,\beta,\gamma,\delta}(-1)^{j_b-m_b+j_d-m_d}\langle\alpha|\hat q_{\lambda\mu}|\beta\rangle\langle\gamma|\hat q_{\lambda\mu}|\delta\rangle^*\langle\Phi^{I_1}_{M_1,\kappa_1}|\hat c_{a m_a}^\dagger\hat c_{b,-m_b}^\dagger\hat c_{d,-m_d}\hat c_{c m_c}|\Phi^{I_2}_{M_2,\kappa_2}\rangle.
$$

For its reduced form, define the scalar operator

$$
\hat X^{(0),\mathrm{pair}}_{\lambda;ab,cd}\equiv\sum_{\mu=-\lambda}^{\lambda}[\hat c_a^\dagger\otimes\hat c_b^\dagger]_{\lambda\mu}\left([\hat c_c^\dagger\otimes\hat c_d^\dagger]_{\lambda\mu}\right)^\dagger=\sum_{\mu,m_a,m_b,m_c,m_d}\langle j_a m_a;j_b,-m_b|\lambda\mu\rangle\langle j_c m_c;j_d,-m_d|\lambda\mu\rangle\hat c_{a m_a}^\dagger\hat c_{b,-m_b}^\dagger\hat c_{d,-m_d}\hat c_{c m_c}.
$$

Then

$$
\boxed{\langle\Phi^{I_1}_{\kappa_1}\Vert\hat H_\lambda^{\mathrm{pair}}\Vert\Phi^{I_2}_{\kappa_2}\rangle=-\frac{G_\lambda}{4(2\lambda+1)}\sum_{a,b,c,d}\langle a\Vert\hat q_\lambda\Vert b\rangle\langle c\Vert\hat q_\lambda\Vert d\rangle^*\langle\Phi^{I_1}_{\kappa_1}\Vert\hat X^{(0),\mathrm{pair}}_{\lambda;ab,cd}\Vert\Phi^{I_2}_{\kappa_2}\rangle}.
$$
