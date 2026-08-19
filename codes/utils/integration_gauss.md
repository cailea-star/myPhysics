# Gaussian Quadrature

### General Theory

##### Weighted Integral

Let $f(x)$ be the integrand and $W(x)>0$ the weight function on $[a,b]$. The weighted integral is

$$
I[f]=\int_a^b f(x)W(x)\,dx.
$$

Let $N_x$ be the number of quadrature nodes. The goal is to find $\{x_i\}_{i=1}^{N_x}$ and weights $\{w_i\}_{i=1}^{N_x}$ such that

$$
Q[f]=\sum_{i=1}^{N_x} f(x_i)w_i
$$

is exact for polynomials of the highest possible degree.

##### Interpolation Quadrature

For distinct nodes $\{x_i\}_{i=1}^{N_x}$, let $\delta_{ij}$ be the Kronecker delta. The Lagrange basis polynomial is

$$
l_i(x)=\prod_{j\ne i}\frac{x-x_j}{x_i-x_j},\qquad l_i(x_j)=\delta_{ij}.
$$

The interpolation polynomial $L(x)$ has degree $\deg L\leq N_x-1$ and satisfies

$$
L(x)=\sum_{i=1}^{N_x} f(x_i)l_i(x),\qquad L(x_i)=f(x_i).
$$

Since $f(x)$ and $L(x)$ agree at every node,

$$
Q[f]=Q[L]=\int_a^b L(x)W(x)\,dx=\sum_{i=1}^{N_x} f(x_i)\int_a^b l_i(x)W(x)\,dx=\sum_{i=1}^{N_x} f(x_i)w_i.
$$

Therefore, the quadrature weights are

$$
w_i=\int_a^b l_i(x)W(x)\,dx.
$$

If $\deg f\leq N_x-1$, then $L(x)=f(x)$ and

$$
Q[f]=I[f].
$$

Thus, interpolation quadrature is exact for polynomials of degree at most $N_x-1$.

##### Orthogonal Polynomial Space

Let $p_i(x)$ be a polynomial of degree $i$, and let $h_i>0$ be its squared norm. The orthogonal polynomials associated with $W(x)$ satisfy

$$
\int_a^b p_i(x)p_j(x)W(x)\,dx=h_i\delta_{ij}.
$$

##### Gaussian Exactness

Choose the nodes $\{x_i\}_{i=1}^{N_x}$ as the roots of $p_{N_x}(x)$:

$$
p_{N_x}(x_i)=0.
$$

For any polynomial $f(x)$ with $\deg f\leq 2N_x-1$, let $q(x)$ and $r(x)$ be the quotient and remainder after division by $p_{N_x}(x)$:

$$
f(x)=q(x)p_{N_x}(x)+r(x),\qquad \deg q\leq N_x-1,\qquad \deg r\leq N_x-1.
$$

Since $q(x)$ belongs to the space spanned by $\{p_0,\ldots,p_{N_x-1}\}$,

$$
\int_a^b q(x)p_{N_x}(x)W(x)\,dx=0.
$$

Because $\deg r\leq N_x-1$, interpolation quadrature is exact for $r(x)$. Therefore,

$$
I[f]=I[r]=Q[r].
$$

Since $p_{N_x}(x_i)=0$, $f(x_i)=r(x_i)$ at every node and $Q[f]=Q[r]$. Consequently,

$$
\boxed{I[f]=Q[f],\qquad \deg f\leq 2N_x-1.}
$$

##### Nodes and Weights

1. Choose the interval $[a,b]$ and weight $W(x)$.
2. Construct the corresponding orthogonal polynomial $p_{N_x}(x)$.
3. Use the roots of $p_{N_x}(x)$ as the nodes $x_i$.
4. Compute the weights from the Lagrange basis $l_i(x)$.

### Examples

##### Gauss-Hermite

For nonnegative integers $n$ and $k$, the physicists' [Hermite polynomials](../../notes/mathematical_physics/specialfunction_hermite_polynomials.md) satisfy

$$
\int_{-\infty}^{\infty}H_n(x)H_k(x)e^{-x^2}\,dx=2^n n!\sqrt{\pi}\,\delta_{nk},\qquad N_n=\frac{1}{\sqrt{2^n n!\sqrt{\pi}}}.
$$

##### Gauss-Laguerre

For nonnegative integers $n$ and $k$, and $\alpha>-1$, the [generalized Laguerre polynomials](../../notes/mathematical_physics/specialfunction_generalized_laguerre_polynomials.md) satisfy

$$
\int_0^\infty L_n^{(\alpha)}(x)L_k^{(\alpha)}(x)x^\alpha e^{-x}\,dx=\frac{\Gamma(n+\alpha+1)}{n!}\delta_{nk},\qquad N_n^\alpha=\sqrt{\frac{n!}{\Gamma(n+\alpha+1)}}.
$$

##### Gauss-Legendre

For integers $l,k\geq m\geq0$ with fixed $m$, the [associated Legendre functions](../../notes/mathematical_physics/specialfunction_associated_legendre_polynomials.md) satisfy

$$
\int_{-1}^1P_l^m(x)P_k^m(x)\,dx=\frac{2}{2l+1}\frac{(l+m)!}{(l-m)!}\delta_{lk},\qquad N_l^m=\sqrt{\frac{2l+1}{2}\frac{(l-m)!}{(l+m)!}}.
$$
