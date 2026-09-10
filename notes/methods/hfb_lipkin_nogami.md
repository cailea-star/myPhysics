# Lipkin–Nogami Method

### Formulation

Use the conventions of the [Hartree–Fock–Bogoliubov Method](../../codes/utils/hfb.md). Consider one particle species and a normalized zero-temperature HFB vacuum $|\Phi\rangle$, with $\langle\hat O\rangle\equiv\langle\Phi|\hat O|\Phi\rangle$.

The particle-number operator and its expectation value are

$$
\hat N=\sum_\alpha\hat c_\alpha^*\hat c_\alpha,\qquad \langle\hat N\rangle=\sum_\alpha\rho_{\alpha\alpha}=\operatorname{Tr}\rho.
$$

Using the fermionic anticommutation relations,

$$
\hat N\hat N=\sum_\alpha\sum_\beta\hat c_\alpha^*\hat c_\alpha\hat c_\beta^*\hat c_\beta=\hat N+\sum_\alpha\sum_\beta\hat c_\alpha^*\hat c_\beta^*\hat c_\beta\hat c_\alpha.
$$

Wick contraction gives

$$
\langle\hat c_\alpha^*\hat c_\beta^*\hat c_\beta\hat c_\alpha\rangle=\rho_{\alpha\alpha}\rho_{\beta\beta}-\rho_{\beta\alpha}\rho_{\alpha\beta}+\kappa_{\alpha\beta}^*\kappa_{\alpha\beta}.
$$

Therefore,

$$
\langle\hat N^2\rangle=\langle\hat N\rangle+\sum_\alpha\sum_\beta\left(\rho_{\alpha\alpha}\rho_{\beta\beta}-\rho_{\beta\alpha}\rho_{\alpha\beta}+\kappa_{\alpha\beta}^*\kappa_{\alpha\beta}\right)=\operatorname{Tr}\rho+(\operatorname{Tr}\rho)^2-\operatorname{Tr}(\rho\rho)+\operatorname{Tr}(\kappa^\dagger\kappa).
$$

The pure HFB identity $\rho-\rho^2=\kappa\kappa^\dagger$ yields the particle-number variance

$$
\boxed{\sigma_N^2\equiv\langle\hat N^2\rangle-\langle\hat N\rangle^2=2\operatorname{Tr}[\rho(1-\rho)]}.
$$

All traces cover the complete single-particle space; time-reversal symmetry is not assumed.

For a target particle number $N_0$, the Lipkin–Nogami prescription imposes $\langle\hat N\rangle=N_0$ and introduces a second-order coefficient $\lambda_2$. The corrected energy is

$$
\boxed{E_{\mathrm{LN}}=\langle\hat H\rangle-\lambda_2\sigma_N^2,\qquad \Delta E_{\mathrm{LN}}=-\lambda_2\sigma_N^2}.
$$

The coefficient $\lambda_2$ is held fixed during variation and updated between self-consistent iterations. With the chemical potential $\lambda$ enforcing the average particle number,

$$
\left.\delta_\Phi\!\left[\langle\hat H\rangle-\lambda_2\sigma_N^2-\lambda(\langle\hat N\rangle-N_0)\right]\right|_{\lambda_2\ \mathrm{fixed}}=0.
$$

This is the LN approximation to variation after particle-number projection. [The Lipkin–Nogami method](https://www.fuw.edu.pl/~dobaczew/vapnp23w/node7.html)

Using $\sigma_N^2=2\operatorname{Tr}[\rho(1-\rho)]$, the correction can be represented in the particle-hole field. For the original HFB fields $h,\Delta$ and single-particle identity $I$,

$$
\boxed{h_{\mathrm{LN}}=h+4\lambda_2\rho-2\lambda_2 I,\qquad \Delta_{\mathrm{LN}}=\Delta}.
$$

The HFB equation uses $h_{\mathrm{LN}}-\lambda I$, with the constant shift $-2\lambda_2 I$ retained in the field. [HFBTHO LN prescription](https://www.fuw.edu.pl/~dobaczew/hfbtho16w/node12.html)

### Effective Pairing Strength

For one particle species, consider a zero-temperature HFB vacuum. Let $(i,\bar i)$ denote canonical pairs, with each pair counted once, and choose their phases such that

$$
\rho_{ii}=\rho_{\bar i\bar i}=v_i^2,\qquad \kappa_{i\bar i}=u_iv_i,\qquad u_i,v_i\geq0,\qquad u_i^2+v_i^2=1.
$$

Define the canonical pairing gap $\Delta_i\equiv-\operatorname{Re}\Delta_{i\bar i}$. In the conventions of this note, the pairing energy and occupation-weighted average gap are

$$
E_{\mathrm{pair}}=\frac12\operatorname{Re}\operatorname{Tr}(\Delta\kappa^\dagger)=-\sum_i\Delta_i u_iv_i,\qquad \bar\Delta=\frac{\sum_i\Delta_i v_i^2}{\sum_i v_i^2}.
$$

The seniority-pairing approximation introduces an auxiliary effective strength,

$$
\boxed{G_{\mathrm{eff}}=-\frac{\bar\Delta^2}{E_{\mathrm{pair}}}},
$$

with $G_{\mathrm{eff}}>0$ for attractive pairing, $E_{\mathrm{pair}}<0$. This prescription estimates the LN curvature without replacing the pairing interaction used to construct $\Delta$. [Effective-strength prescription](https://www.fuw.edu.pl/~dobaczew/vapnp23w/node7.html)

Define the canonical sums

$$
S_{31}=\sum_i u_i^3v_i,\qquad S_{13}=\sum_i u_iv_i^3,\qquad S_{22}=\sum_i u_i^2v_i^2,\qquad S_{44}=\sum_i u_i^4v_i^4.
$$

The resulting second-order coefficient is

$$
\boxed{\lambda_2=\frac{G_{\mathrm{eff}}}{4}\frac{S_{31}S_{13}-S_{44}}{S_{22}^2-S_{44}}}.
$$

This is an effective seniority-pairing estimate of the general LN coefficient. It requires nonzero pairing energy and denominator; it is not a finite-temperature prescription. [HFBTHO LN prescription](https://www.fuw.edu.pl/~dobaczew/hfbtho16w/node12.html)
