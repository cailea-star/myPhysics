# Numerical Integration

### Interpolation Quadrature

##### Theoretical Basis

For distinct nodes $\{x_i\}_{i=1}^{N_x}$, let $\delta_{ij}$ be the Kronecker delta. The Lagrange basis polynomial is

$$
l_i(x)=\prod_{j\ne i}\frac{x-x_j}{x_i-x_j},\qquad l_i(x_j)=\delta_{ij}.
$$

The interpolation polynomial $L(x)$ has degree $\deg L\leq N_x-1$ and satisfies

$$
L(x)=\sum_{i=1}^{N_x}f(x_i)l_i(x),\qquad L(x_i)=f(x_i).
$$

Numerical integration replaces $f(x)$ with $L(x)$:

$$
Q[f]=Q[L]=\int_a^bL(x)\,dx=\sum_{i=1}^{N_x}f(x_i)\int_a^bl_i(x)\,dx=\sum_{i=1}^{N_x}f(x_i)w_i.
$$

Therefore, the quadrature weights are

$$
w_i=\int_a^bl_i(x)\,dx.
$$

If $\deg f\leq N_x-1$, then $L(x)=f(x)$ and

$$
Q[f]=I[f].
$$

Thus, interpolation quadrature is exact for polynomials of degree at most $N_x-1$.

##### Examples

**Trapezoidal Rule.** Let $f_i=f(x_i)$ and $h_i=x_{i+1}-x_i$. The linear interpolation polynomial on $[x_i,x_{i+1}]$ is

$$
L_i^{(1)}(x)=f_i\frac{x_{i+1}-x}{h_i}+f_{i+1}\frac{x-x_i}{h_i}.
$$

Its integral is

$$
\int_{x_i}^{x_{i+1}}L_i^{(1)}(x)\,dx=\frac{h_i}{2}(f_i+f_{i+1}).
$$

Summing over all intervals gives the composite trapezoidal rule:

$$
\boxed{Q_T[f]=\sum_{i=0}^{N_x-2}\frac{x_{i+1}-x_i}{2}(f_i+f_{i+1})}.
$$

**Simpson's Rule.** Let $x_0,x_1,x_2$ be equally spaced nodes with $h=x_1-x_0=x_2-x_1$. The quadratic interpolation polynomial is

$$
L^{(2)}(x)=f_0\frac{(x-x_1)(x-x_2)}{2h^2}-f_1\frac{(x-x_0)(x-x_2)}{h^2}+f_2\frac{(x-x_0)(x-x_1)}{2h^2}.
$$

Its integral is

$$
\int_{x_0}^{x_2}L^{(2)}(x)\,dx=\frac{h}{3}(f_0+4f_1+f_2).
$$

For an odd number of points, let $h_k=x_{2k+1}-x_{2k}=x_{2k+2}-x_{2k+1}$. Summing over all three-point panels gives the composite Simpson rule:

$$
\boxed{Q_S[f]=\sum_{k=0}^{(N_x-3)/2}\frac{h_k}{3}(f_{2k}+4f_{2k+1}+f_{2k+2})}.
$$

For a uniform grid with spacing $h$ and sufficiently smooth $f(x)$,

$$
I[f]-Q_T[f]=\mathcal O(h^2),\qquad I[f]-Q_S[f]=\mathcal O(h^4).
$$
