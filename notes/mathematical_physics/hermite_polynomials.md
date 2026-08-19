# Hermite Polynomials

Let $n\in\mathbb N_0=\{0,1,2,\ldots\}$. The physicists' Hermite polynomial $H_n(x)$ is the degree-$n$ solution of

$$
\boxed{H_n''(x)-2xH_n'(x)+2nH_n(x)=0}.
$$

Its normalization is fixed by the generating function below.

### Generating Function

With $H_0(x)=1$, the Hermite equation gives the raising relation

$$
H_{n+1}(x)=\left(2x-\frac{d}{dx}\right)H_n(x).
$$

Let $t$ be an auxiliary variable. For the exponential generating function $G(x,t)$,

$$
\boxed{G(x,t)=\sum_{n=0}^{\infty}H_n(x)\frac{t^n}{n!}}.
$$

The raising relation gives

$$
\frac{\partial G}{\partial t}=\sum_{n=0}^{\infty}H_{n+1}(x)\frac{t^n}{n!}=\left(2x-\frac{\partial}{\partial x}\right)G\quad\Longrightarrow\quad\frac{\partial G}{\partial t}+\frac{\partial G}{\partial x}=2xG.
$$

With the initial condition $G(x,0)=1$, the solution is

$$
\boxed{G(x,t)=e^{2xt-t^2}}.
$$

### Rodrigues Formula

Rewrite the generating function as

$$
G(x,t)=e^{x^2}e^{-(x-t)^2}.
$$

The Taylor expansion with respect to $t$ gives

$$
e^{-(x-t)^2}=\sum_{n=0}^{\infty}\frac{(-t)^n}{n!}\frac{d^n}{dx^n}e^{-x^2}.
$$

Comparing the coefficients of $t^n$ gives

$$
\boxed{H_n(x)=(-1)^ne^{x^2}\frac{d^n}{dx^n}e^{-x^2}}.
$$

### Recurrence and Derivatives

The generating function satisfies

$$
\frac{\partial G}{\partial t}=(2x-2t)G,\qquad \frac{\partial G}{\partial x}=2tG.
$$

Comparing the coefficients of $t^n$ gives

$$
\boxed{H_{n+1}(x)=2xH_n(x)-2nH_{n-1}(x),\qquad H_n'(x)=2nH_{n-1}(x)}.
$$

Equivalently,

$$
xH_n(x)=\frac12H_{n+1}(x)+nH_{n-1}(x),\qquad H_n''(x)=2xH_n'(x)-2nH_n(x).
$$

### Parity

Since $G(-x,t)=G(x,-t)$,

$$
H_n(-x)=(-1)^nH_n(x).
$$

At the origin,

$$
H_{2k+1}(0)=0,\qquad H_{2k}(0)=(-1)^k\frac{(2k)!}{k!},\qquad k\in\mathbb N_0.
$$

### Orthogonality and Normalization

For $m<n$, the Rodrigues formula and integration by parts give

$$
\int_{-\infty}^{\infty}e^{-x^2}H_m(x)H_n(x)\,dx=\int_{-\infty}^{\infty}e^{-x^2}\frac{d^nH_m(x)}{dx^n}\,dx=0.
$$

Let $\delta_{mn}$ denote the Kronecker delta. The complete orthogonality relation is

$$
\boxed{\int_{-\infty}^{\infty}e^{-x^2}H_m(x)H_n(x)\,dx=2^nn!\sqrt\pi\,\delta_{mn}}.
$$

The normalization coefficient is

$$
N_n=\frac{1}{\sqrt{2^nn!\sqrt\pi}}.
$$

For $f\in L^2(\mathbb R,e^{-x^2}dx)$, the Hermite expansion is

$$
f(x)=\sum_{n=0}^{\infty}c_nH_n(x),\qquad c_n=N_n^2\int_{-\infty}^{\infty}e^{-x^2}H_n(x)f(x)\,dx.
$$
