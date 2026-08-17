### Shooting Method

Consider the boundary-value problem

$$
y''(x)=F(x;E)y(x),\qquad x\in[x_{\min},x_{\max}],
$$

where $E$ is the unknown spectral parameter. The inner and outer boundary conditions are

$$
\mathcal B_{\mathrm{in}}(E)=\bigl(y_{\mathrm{in}}(x_{\min}),y_{\mathrm{in}}'(x_{\min})\bigr),\qquad \mathcal B_{\mathrm{out}}(E)=\bigl(y_{\mathrm{out}}(x_{\max}),y_{\mathrm{out}}'(x_{\max})\bigr).
$$

Specific boundary conditions are supplied externally; spherical radial boundary conditions are given in [Spherical Radial Boundary Conditions](../utils_spherical/spherical_radial_boundary.md).

##### Bidirectional Propagation

Choose a matching point $x_m\in(x_{\min},x_{\max})$. Propagate the inner solution $y_{\mathrm{in}}$ from $x_{\min}$ to $x_m$, and the outer solution $y_{\mathrm{out}}$ from $x_{\max}$ to $x_m$.

Both solutions may be propagated using the [Numerov Method](ivp_numerov.md) or a Runge--Kutta method.

##### Matching Condition

At $x_m$, define

$$
W_m(E)=y_{\mathrm{in}}'(x_m)y_{\mathrm{out}}(x_m)-y_{\mathrm{in}}(x_m)y_{\mathrm{out}}'(x_m).
$$

The logarithmic-derivative residual is

$$
R(E)=\frac{y_{\mathrm{in}}'(x_m)}{y_{\mathrm{in}}(x_m)}-\frac{y_{\mathrm{out}}'(x_m)}{y_{\mathrm{out}}(x_m)}=\frac{W_m(E)}{y_{\mathrm{in}}(x_m)y_{\mathrm{out}}(x_m)}.
$$

The eigenenergy satisfies

$$
\boxed{R(E_n)=0.}
$$

The matching point must therefore avoid the nodes of both solutions.

##### Eigenenergy Search

For a real Sturm--Liouville bound-state problem, let $N(E)$ denote the number of nodes. The $n$th state satisfies

$$
N(E_n)=n.
$$

Within the selected energy interval $[E_{\mathrm{lo}},E_{\mathrm{up}}]$, the numerical eigenenergy is

$$
\boxed{E_n=\operatorname*{arg\,min}_{E\in[E_{\mathrm{lo}},E_{\mathrm{up}}]}|R(E)|.}
$$

For complex boundary conditions, node counting is generally unavailable, and the search interval must be specified independently.

##### Wave-Function Matching and Normalization

Scale the outer solution by

$$
A(E_n)=\frac{y_{\mathrm{in}}(x_m)}{y_{\mathrm{out}}(x_m)}.
$$

The matched solution is

$$
y(x)=\begin{cases}y_{\mathrm{in}}(x),&x\leq x_m,\\A(E_n)y_{\mathrm{out}}(x),&x>x_m.\end{cases}
$$

The scaling ensures continuity of the function value, while $R(E_n)=0$ ensures continuity of the first derivative. Finally, normalize the solution:

$$
\mathrm{Norm}=\sqrt{\int_{x_{\min}}^{x_{\max}}|y(x)|^2\,dx},\qquad y(x)\rightarrow\frac{y(x)}{\mathrm{Norm}}.
$$
