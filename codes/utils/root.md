# Scalar Root Finding

For a real-valued function $f:\mathbb R\to\mathbb R$, scalar root finding seeks $x^*$ such that

$$
f(x^*)=0.
$$

Bracketing methods preserve an interval containing a root, while local iterative methods construct successive approximations $\{x_h\}$. Convergence is commonly monitored through

$$
|f(x_h)|\leq\epsilon_f,\qquad |x_{h+1}-x_h|\leq\epsilon_x.
$$

### Bisection Method

Let $f(x)$ be continuous on $[a_0,b_0]$ and satisfy

$$
f(a_0)f(b_0)\leq0.
$$

At iteration $h$, define the midpoint

$$
c_h=\frac{a_h+b_h}{2}.
$$

Retain the half-interval containing a root:

$$
[a_{h+1},b_{h+1}]=\begin{cases}\left[a_h,c_h\right],&f(a_h)f(c_h)\leq0,\\\left[c_h,b_h\right],&f(a_h)f(c_h)>0.\end{cases}
$$

After $h$ iterations,

$$
b_h-a_h=\frac{b_0-a_0}{2^h},\qquad |c_h-x^*|\leq\frac{b_0-a_0}{2^{h+1}}.
$$

Thus, the bisection method converges linearly with convergence factor $1/2$.

### Newton's Method

Let $x^*$ be a simple root with $f'(x^*)\neq0$. The first-order Taylor expansion about $x_h$ gives

$$
f(x^*)\approx f(x_h)+f'(x_h)(x^*-x_h)=0.
$$

Therefore,

$$
\boxed{x_{h+1}=x_h-\frac{f(x_h)}{f'(x_h)}}.
$$

When the derivative is unavailable, use the forward difference

$$
f'(x_h)=\frac{f(x_h+\Delta x)-f(x_h)}{\Delta x}+\mathcal O(\Delta x).
$$

Let $e_h=x_h-x^*$. For an initial value sufficiently close to the root and an exact derivative,

$$
|e_{h+1}|=\mathcal O(|e_h|^2).
$$

Finite-difference error limits quadratic convergence.

Set $x_0=(x_{\mathrm{lo},0}+x_{\mathrm{up},0})/2$ and require

$$
f(x_{\mathrm{lo},0})f(x_{\mathrm{up},0})\leq0.
$$

After evaluating $f(x_h)$, retain the subinterval containing the root:

$$
[x_{\mathrm{lo},h+1},x_{\mathrm{up},h+1}]=\begin{cases}[x_{\mathrm{lo},h},x_h],&f(x_{\mathrm{lo},h})f(x_h)\leq0,\\ [x_h,x_{\mathrm{up},h}],&f(x_{\mathrm{lo},h})f(x_h)>0.\end{cases}
$$

Accept the Newton candidate only inside the updated interval; otherwise, use its midpoint.

### Brent's Method

Let $f(x)$ be continuous on $[a_h,b_h]$ and satisfy

$$
f(a_h)f(b_h)\leq0.
$$

The secant candidate is

$$
x_{\mathrm{sec}}=b_h-f(b_h)\frac{b_h-a_h}{f(b_h)-f(a_h)}.
$$

Given three points $\{x_i\}_{i=0}^2$ with distinct function values $f_i=f(x_i)$, inverse quadratic interpolation gives

$$
x_{\mathrm{IQI}}=\sum_{i=0}^2x_i\prod_{j\ne i}\frac{-f_j}{f_i-f_j}.
$$

Brent's method prefers an interpolation candidate. If the candidate lies outside the bracketing interval or fails to reduce it sufficiently, the method uses the bisection point

$$
x_{\mathrm{bis}}=\frac{a_h+b_h}{2}.
$$

Thus, Brent's method retains the robustness of bisection while using interpolation to obtain superlinear convergence near a simple root.
