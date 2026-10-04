# Electromagnetic Radiation

## Source Equations for Electromagnetic Potentials

Let $\rho(\mathbf r,t)$ and $\mathbf j(\mathbf r,t)$ be the charge and current densities in vacuum, and let $\epsilon_0$ and $\mu_0$ be the vacuum permittivity and permeability. The [Maxwell equations](electromagnetic_maxwell_equations.md#maxwell-equations-in-vacuum) with sources are

$$
\nabla\cdot\mathbf E=\frac{\rho}{\epsilon_0},\qquad \nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t},
$$

$$
\nabla\cdot\mathbf B=0,\qquad \nabla\times\mathbf B=\mu_0\mathbf j+\mu_0\epsilon_0\frac{\partial\mathbf E}{\partial t}.
$$

Introduce the scalar and vector potentials:

$$
\mathbf E=-\nabla\Phi-\frac{\partial\mathbf A}{\partial t},\qquad \mathbf B=\nabla\times\mathbf A.
$$

Substitution into the two sourced equations gives

$$
\nabla^2\Phi+\frac{\partial}{\partial t}(\nabla\cdot\mathbf A)=-\frac{\rho}{\epsilon_0},
$$

$$
\left(\nabla^2-\mu_0\epsilon_0\frac{\partial^2}{\partial t^2}\right)\mathbf A-\nabla\left(\nabla\cdot\mathbf A+\mu_0\epsilon_0\frac{\partial\Phi}{\partial t}\right)=-\mu_0\mathbf j.
$$

Define the speed of light in vacuum by

$$
c^{-2}=\mu_0\epsilon_0.
$$

## Retarded Solutions in the Lorenz Gauge

Under a gauge transformation,

$$
\mathbf A'=\mathbf A+\nabla\chi,\qquad \Phi'=\Phi-\frac{\partial\chi}{\partial t}.
$$

Therefore,

$$
\nabla\cdot\mathbf A'+\frac{1}{c^2}\frac{\partial\Phi'}{\partial t}=\nabla\cdot\mathbf A+\frac{1}{c^2}\frac{\partial\Phi}{\partial t}+\left(\nabla^2-\frac{1}{c^2}\frac{\partial^2}{\partial t^2}\right)\chi.
$$

Choose $\chi$ so that the right-hand side vanishes. This imposes the Lorenz gauge:

$$
\boxed{\nabla\cdot\mathbf A+\frac{1}{c^2}\frac{\partial\Phi}{\partial t}=0},
$$

and the coupled equations from the preceding section become

$$
\boxed{\left(\nabla^2-\frac{1}{c^2}\frac{\partial^2}{\partial t^2}\right)\Phi=-\frac{\rho}{\epsilon_0},\qquad \left(\nabla^2-\frac{1}{c^2}\frac{\partial^2}{\partial t^2}\right)\mathbf A=-\mu_0\mathbf j.}
$$

Both potentials satisfy the same form of sourced wave equation. Define the retarded Green function by

$$
\left(\nabla^2-\frac{1}{c^2}\frac{\partial^2}{\partial t^2}\right)G_{\mathrm{ret}}(\mathbf r,t;\mathbf r',t')=-\delta^{(3)}(\mathbf r-\mathbf r')\delta(t-t'),
$$

where the derivatives act on the observation point $(\mathbf r,t)$. Let $\mathbf k$ be the spatial Fourier wave vector and $k=|\mathbf k|$. Expand directly in space and time:

$$
G_{\mathrm{ret}}(\mathbf r,t;\mathbf r',t')=\int\frac{\mathrm d^3\mathbf k\,\mathrm d\omega}{(2\pi)^4}\,G(\mathbf k,\omega)e^{i\mathbf k\cdot(\mathbf r-\mathbf r')-i\omega(t-t')}.
$$

Substitution into the Green-function equation gives the algebraic relation

$$
\left(k^2-\frac{\omega^2}{c^2}\right)G(\mathbf k,\omega)=1.
$$

The denominator has poles at $k^2=\frac{\omega^2}{c^2}$. Let $0^+$ denote a positive infinitesimal with units of inverse length. Select the retarded prescription:

$$
G(\mathbf k,\omega)=\frac{1}{k^2-\left(\frac{\omega}{c}+i0^+\right)^2}.
$$

First integrate over the direction of $\mathbf k$:

$$
\int\frac{\mathrm d^3\mathbf k}{(2\pi)^3}\frac{e^{i\mathbf k\cdot(\mathbf r-\mathbf r')}}{k^2-\left(\frac{\omega}{c}+i0^+\right)^2}=\frac{1}{2\pi^2|\mathbf r-\mathbf r'|}\int_0^\infty\frac{k\sin\!\left(k|\mathbf r-\mathbf r'|\right)}{k^2-\left(\frac{\omega}{c}+i0^+\right)^2}\,\mathrm dk.
$$

Extending the radial integral into the complex plane, the selected pole yields an outgoing spherical wave:

$$
\int\frac{\mathrm d^3\mathbf k}{(2\pi)^3}\frac{e^{i\mathbf k\cdot(\mathbf r-\mathbf r')}}{k^2-\left(\frac{\omega}{c}+i0^+\right)^2}=\frac{e^{i\frac{\omega}{c}|\mathbf r-\mathbf r'|}}{4\pi|\mathbf r-\mathbf r'|}.
$$

The inverse frequency transform then gives

$$
G_{\mathrm{ret}}(\mathbf r,t;\mathbf r',t')=\frac{1}{4\pi|\mathbf r-\mathbf r'|}\int_{-\infty}^{\infty}\frac{\mathrm d\omega}{2\pi}e^{-i\omega\left(t-t'-\frac{|\mathbf r-\mathbf r'|}{c}\right)}=\frac{\delta\!\left(t-t'-\frac{|\mathbf r-\mathbf r'|}{c}\right)}{4\pi|\mathbf r-\mathbf r'|}.
$$

Convolve this Green function with the sources:

$$
\Phi(\mathbf r,t)=\frac{1}{\epsilon_0}\int\mathrm d^3r'\int\mathrm dt'\,G_{\mathrm{ret}}(\mathbf r,t;\mathbf r',t')\rho(\mathbf r',t'),\qquad \mathbf A(\mathbf r,t)=\mu_0\int\mathrm d^3r'\int\mathrm dt'\,G_{\mathrm{ret}}(\mathbf r,t;\mathbf r',t')\mathbf j(\mathbf r',t').
$$

The delta function fixes the source time in the $t'$ integral:

$$
\int\mathrm dt'\,\delta\!\left(t-t'-\frac{|\mathbf r-\mathbf r'|}{c}\right)\rho(\mathbf r',t')=\rho\!\left(\mathbf r',t-\frac{|\mathbf r-\mathbf r'|}{c}\right).
$$

Therefore, the retarded potentials are

$$
\boxed{\Phi(\mathbf r,t)=\frac{1}{4\pi\epsilon_0}\int\frac{\rho\!\left(\mathbf r',t-\frac{|\mathbf r-\mathbf r'|}{c}\right)}{|\mathbf r-\mathbf r'|}\,\mathrm d^3r',\qquad \mathbf A(\mathbf r,t)=\frac{\mu_0}{4\pi}\int\frac{\mathbf j\!\left(\mathbf r',t-\frac{|\mathbf r-\mathbf r'|}{c}\right)}{|\mathbf r-\mathbf r'|}\,\mathrm d^3r'.}
$$

Only source events on the past light cone of $(\mathbf r,t)$ contribute to the retarded potentials.

## Radiation Zones and Multipole Expansion

Let the isolated source occupy a region of size $a$. For an observation point outside the source, define $r=|\mathbf r|$ and $\hat{\mathbf r}=\mathbf r/r$. At a characteristic angular frequency $\omega$, let

$$
k=\frac{\omega}{c},\qquad \lambda=\frac{2\pi}{k}.
$$

For $r\gg a$, expand the retarded source about $\mathbf r'=0$. When $ka\ll1$, the lowest multipoles dominate. For $F=\rho$ or any component of $\mathbf j$,

$$
\frac{F\!\left(\mathbf r',t-\frac{|\mathbf r-\mathbf r'|}{c}\right)}{|\mathbf r-\mathbf r'|}=\frac{F\!\left(\mathbf r',t-\frac{r}{c}\right)}{r}-\mathbf r'\cdot\nabla\left[\frac{F\!\left(\mathbf r',t-\frac{r}{c}\right)}{r}\right]+\cdots.
$$

Define the total charge and electric dipole moment:

$$
Q(t)=\int\rho(\mathbf r',t)\,\mathrm d^3r',\qquad \mathbf p(t)=\int\mathbf r'\rho(\mathbf r',t)\,\mathrm d^3r'.
$$

For a localized, isolated source, the continuity equation gives

$$
\frac{\mathrm dQ}{\mathrm dt}=-\int\nabla'\cdot\mathbf j(\mathbf r',t)\,\mathrm d^3r'=0,\qquad \frac{\mathrm d\mathbf p}{\mathrm dt}=-\int\mathbf r'\bigl(\nabla'\cdot\mathbf j(\mathbf r',t)\bigr)\,\mathrm d^3r'=\int\mathbf j(\mathbf r',t)\,\mathrm d^3r'.
$$

Thus the first terms of the retarded potentials are

$$
\Phi(\mathbf r,t)=\frac{1}{4\pi\epsilon_0}\left[\frac{Q}{r}-\nabla\cdot\frac{\mathbf p\!\left(t-\frac{r}{c}\right)}{r}+\cdots\right],\qquad \mathbf A(\mathbf r,t)=\frac{\mu_0}{4\pi}\left[\frac{\dot{\mathbf p}\!\left(t-\frac{r}{c}\right)}{r}+\cdots\right].
$$

The conserved monopole charge produces no radiation; the electric dipole is the first possible radiating term. For any source moment $f(t)$,

$$
\nabla\left[\frac{f\!\left(t-\frac{r}{c}\right)}{r}\right]=-\hat{\mathbf r}\left[\frac{f\!\left(t-\frac{r}{c}\right)}{r^2}+\frac{\dot f\!\left(t-\frac{r}{c}\right)}{cr}\right].
$$

For $f(t)\propto e^{-i\omega t}$, this becomes

$$
\nabla\left[\frac{f\!\left(t-\frac{r}{c}\right)}{r}\right]=\left(i\frac{\omega}{c}r-1\right)\frac{f\!\left(t-\frac{r}{c}\right)}{r^2}\hat{\mathbf r}=(ikr-1)\frac{f\!\left(t-\frac{r}{c}\right)}{r^2}\hat{\mathbf r}.
$$

Outside the source, the near zone has $kr\ll1$; the radiation zone has $kr\gg1$. In the radiation zone, terms proportional to $1/r$ dominate over terms proportional to $1/r^2$.
