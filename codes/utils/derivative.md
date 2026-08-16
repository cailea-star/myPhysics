### Finite-Difference Derivatives

Let $f(x)$ be a smooth function, $\Delta x$ a displacement from $x$, and $f^{(k)}(x)$ its $k$-th derivative. The Taylor expansion is

$$
f(x+\Delta x)=f(x)+f'(x)\Delta x+\frac{f''(x)}{2!}(\Delta x)^2+\cdots+\frac{f^{(k)}(x)}{k!}(\Delta x)^k+\cdots.
$$

For a uniform grid $x_i=x_{\min}+ih$ with $h\neq0$, define $f_i=f(x_i)$. The forward and backward approximations are

$$
f'(x_i)=\frac{f_{i+1}-f_i}{h}+\mathcal{O}(h),\qquad f'(x_i)=\frac{f_i-f_{i-1}}{h}+\mathcal{O}(h).
$$

The three-point and five-point central approximations of the first derivative are

$$
f'(x_i)=\frac{f_{i+1}-f_{i-1}}{2h}+\mathcal{O}(h^2)=\frac{8(f_{i+1}-f_{i-1})-(f_{i+2}-f_{i-2})}{12h}+\mathcal{O}(h^4).
$$

The corresponding central approximations of the second derivative are

$$
f''(x_i)=\frac{f_{i+1}-2f_i+f_{i-1}}{h^2}+\mathcal{O}(h^2)=\frac{16(f_{i+1}+f_{i-1})-(f_{i+2}+f_{i-2})-30f_i}{12h^2}+\mathcal{O}(h^4).
$$
