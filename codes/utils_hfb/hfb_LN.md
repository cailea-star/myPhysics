# Lipkin–Nogami Method

### Formulation

Use the conventions of the [Hartree–Fock–Bogoliubov Method](hfb_full.md). Consider one particle species and a normalized zero-temperature HFB vacuum $|\Phi\rangle$, with $\langle\hat O\rangle\equiv\langle\Phi|\hat O|\Phi\rangle$.

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

For one particle species, define the pairing energy and an occupation-weighted average pairing gap by

$$
E_{\mathrm{pair}}=\frac12\operatorname{Re}\operatorname{Tr}(\Delta\kappa^\dagger),\qquad \bar\Delta=\frac{\sum_\alpha \Delta_\alpha n_\alpha}{\sum_\alpha n_\alpha}.
$$

Here $n_\alpha$ denotes a canonical occupation, and $\Delta_\alpha$ the corresponding pairing gap, with positive gaps for attractive pairing. The effective seniority-pairing strength is

$$
\boxed{G_{\mathrm{eff}}=-\frac{\bar\Delta^2}{E_{\mathrm{pair}}}}.
$$

This auxiliary strength estimates $\lambda_2$; it does not replace the interaction generating $\Delta$.

##### Excluding a Blocked Pair

Within the canonical blocking approximation, let $k$ denote the singly occupied pair. The remaining pairs have amplitudes $u_i,v_i$, with $u_i^2+v_i^2=1$. Count each unblocked pair once:

$$
S_{31}^{(k)}=\sum_{i\ne k}u_i^3v_i,\qquad S_{13}^{(k)}=\sum_{i\ne k}u_iv_i^3,\qquad S_{22}^{(k)}=\sum_{i\ne k}u_i^2v_i^2,\qquad S_{44}^{(k)}=\sum_{i\ne k}u_i^4v_i^4.
$$

The blocked-pair-exclusion estimate is

$$
\boxed{\lambda_2^{(k)}=\frac{G_{\mathrm{eff}}^{(k)}}{4}\frac{S_{31}^{(k)}S_{13}^{(k)}-S_{44}^{(k)}}{\left(S_{22}^{(k)}\right)^2-S_{44}^{(k)}}}.
$$

The pairing energy and average gap defining $G_{\mathrm{eff}}^{(k)}$ must follow the same blocked-state prescription. General quasiparticle blocking does not necessarily identify a single canonical pair $k$.

##### EFA Occupation-Based Estimate

For the EFA density, diagonalize the normal density matrix:

$$
\rho_{\mathrm{EFA}}W=W\,\operatorname{diag}(n_\alpha),\qquad 0\le n_\alpha\le1.
$$

Define auxiliary amplitudes

$$
\widetilde v_\alpha=\sqrt{n_\alpha},\qquad \widetilde u_\alpha=\sqrt{1-n_\alpha}.
$$

These amplitudes parameterize the normal occupations; they do not imply $\kappa_{\alpha\bar\alpha}=\widetilde u_\alpha\widetilde v_\alpha$.

With $\alpha$ covering the complete single-particle space, count each Kramers pair once through

$$
S_{31}^{\mathrm{EFA}}=\frac12\sum_\alpha\widetilde u_\alpha^3\widetilde v_\alpha,\qquad S_{13}^{\mathrm{EFA}}=\frac12\sum_\alpha\widetilde u_\alpha\widetilde v_\alpha^3,\qquad S_{22}^{\mathrm{EFA}}=\frac12\sum_\alpha\widetilde u_\alpha^2\widetilde v_\alpha^2,\qquad S_{44}^{\mathrm{EFA}}=\frac12\sum_\alpha\widetilde u_\alpha^4\widetilde v_\alpha^4.
$$

The occupation-based estimate is

$$
\boxed{\lambda_2^{\mathrm{EFA}}=\frac{G_{\mathrm{eff}}^{\mathrm{EFA}}}{4}\frac{S_{31}^{\mathrm{EFA}}S_{13}^{\mathrm{EFA}}-S_{44}^{\mathrm{EFA}}}{\left(S_{22}^{\mathrm{EFA}}\right)^2-S_{44}^{\mathrm{EFA}}}}.
$$

Here $G_{\mathrm{eff}}^{\mathrm{EFA}}$ uses the EFA pairing energy and average gap. This is a seniority-pairing estimate based on EFA occupations, not a derivation of the LN coefficient from EFA many-body moments.

Both estimates require nonzero pairing energy and denominator. The external model selects the prescription and supplies $\lambda_2$; the HFB solver only applies the supplied coefficient.
