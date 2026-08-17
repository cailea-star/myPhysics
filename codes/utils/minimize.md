# Scalar Minimization

For a real-valued function $f:\mathbb R\to\mathbb R$, scalar minimization seeks

$$
x^*=\underset{x\in[a,b]}{\operatorname{argmin}}\,f(x).
$$

For an interior minimum of a differentiable function,

$$
f'(x^*)=0,\qquad f''(x^*)>0.
$$

### Golden-Section Search

Assume that $f(x)$ is unimodal on $[a_0,b_0]$. Define

$$
\tau=\frac{\sqrt5-1}{2},\qquad c_h=b_h-\tau(b_h-a_h),\qquad d_h=a_h+\tau(b_h-a_h).
$$

If $f(c_h)<f(d_h)$, retain $\left[a_h,d_h\right]$; otherwise retain $\left[c_h,b_h\right]$. Since $\tau^2=1-\tau$, one interior point is reused at the next iteration.

After $h$ iterations,

$$
b_h-a_h=\tau^h(b_0-a_0).
$$

Thus, the golden-section search converges linearly with convergence factor $\tau$.

### Brent's Method

Let three points satisfy

$$
a_h<x_h<b_h,\qquad f(x_h)<f(a_h),\qquad f(x_h)<f(b_h).
$$

Construct the quadratic interpolation polynomial $p(x)$ through the three points:

$$
p(a_h)=f(a_h),\qquad p(x_h)=f(x_h),\qquad p(b_h)=f(b_h).
$$

The parabolic candidate satisfies

$$
p'(x_{\mathrm{par}})=0.
$$

Brent's method accepts $x_{\mathrm{par}}$ when it lies inside the bracketing interval and provides sufficient progress; otherwise it performs a golden-section step.

Thus, Brent's method retains the robustness of golden-section search while using parabolic interpolation to accelerate convergence near a smooth minimum.
