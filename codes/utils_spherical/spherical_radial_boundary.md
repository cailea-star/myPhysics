# Spherical Radial Boundary Conditions

Define the reduced-mass kinetic coefficient

$$
h_\mu\equiv\frac{\hbar^2}{2\mu}.
$$

The reduced radial equation is

$$
\left[-h_\mu\frac{d^2}{dr^2}+V_l(r)\right]u_l(r)=Eu_l(r),
$$

where $u_l(r)=rR_l(r)$, $l$ is the orbital angular momentum, and $V_l(r)$ includes the centrifugal potential. A radial boundary condition supplies

$$
\mathcal B(r;l,h_\mu,Z,E)=\left(u_l(r),u_l'(r)\right),
$$

where $Z$ is the charge product entering $V_C(r)=Ze^2/r$. The overall normalization of the boundary pair is arbitrary.

### Regular Origin Boundary

For a potential less singular than $r^{-2}$ near the origin, the regular solution satisfies

$$
u_l(r)\propto r^{l+1},\qquad u_l'(r)\propto(l+1)r^l.
$$

Thus, at a small inner radius $r_{\min}>0$,

$$
\boxed{\mathcal B_{\mathrm{reg}}(r_{\min})=A\left(r_{\min}^{l+1},(l+1)r_{\min}^l\right),}
$$

where $A\neq0$ is an arbitrary scale.

### Bound-State Decaying Boundary

For $E<0$ and a short-range potential satisfying $V_l(r)\to0$, define

$$
\kappa=\sqrt{\frac{-2\mu E}{\hbar^2}}=\sqrt{\frac{-E}{h_\mu}}.
$$

The asymptotically decaying solution is

$$
u_l(r)\propto e^{-\kappa r},\qquad u_l'(r)=-\kappa u_l(r).
$$

Therefore, at a sufficiently large outer radius $r_{\max}$,

$$
\boxed{\mathcal B_{\mathrm{dec}}(r_{\max})=A\left(e^{-\kappa r_{\max}},-\kappa e^{-\kappa r_{\max}}\right).}
$$

This form assumes a short-range asymptotic potential; a nonzero Coulomb tail requires the corresponding Coulomb bound-state asymptotic function.

### Coulomb Incoming and Outgoing Boundaries

For $E>0$, define the wave number, dimensionless radius, and Sommerfeld parameter by

$$
k=\sqrt{\frac{2\mu E}{\hbar^2}}=\sqrt{\frac{E}{h_\mu}},\qquad \rho=kr,\qquad \eta=\frac{\mu Ze^2}{\hbar^2k}=\frac{Ze^2}{2h_\mu k}.
$$

Let $F_l(\eta,\rho)$ and $G_l(\eta,\rho)$ be the regular and irregular Coulomb functions. The Coulomb--Hankel functions are

$$
H_l^{(\pm)}(\eta,\rho)=G_l(\eta,\rho)\pm iF_l(\eta,\rho).
$$

Their radial derivatives satisfy

$$
\frac{dH_l^{(\pm)}}{dr}=k\left(\frac{\partial G_l}{\partial\rho}\pm i\frac{\partial F_l}{\partial\rho}\right).
$$

Hence,

$$
\boxed{\mathcal B_{\pm}(r)=\left(H_l^{(\pm)}(\eta,kr),k\frac{\partial H_l^{(\pm)}(\eta,kr)}{\partial\rho}\right).}
$$

With this convention, $H_l^{(+)}$ is outgoing and $H_l^{(-)}$ is incoming. These boundary pairs provide the initial values for inward or outward radial propagation.
