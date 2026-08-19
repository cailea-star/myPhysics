# Generalized Laguerre Polynomials

Let $n\in\mathbb N_0=\{0,1,2,\ldots\}$ and $\alpha>-1$. The generalized Laguerre polynomial $L_n^{(\alpha)}(x)$ is the degree-$n$ solution of

$$
\boxed{x\frac{d^2}{dx^2}L_n^{(\alpha)}(x)+(\alpha+1-x)\frac{d}{dx}L_n^{(\alpha)}(x)+nL_n^{(\alpha)}(x)=0}.
$$

For $\alpha\in\mathbb N_0$, this family includes the associated Laguerre polynomials used in radial problems.

### Generating Function

With $L_0^{(\alpha)}(x)=1$, the Laguerre equation gives the raising relation

$$
(n+1)L_{n+1}^{(\alpha)}(x)=(n+\alpha+1-x)L_n^{(\alpha)}(x)+x\frac{d}{dx}L_n^{(\alpha)}(x).
$$

Let $t$ be an auxiliary variable. For the ordinary generating function $G_\alpha(x,t)$,

$$
\boxed{G_\alpha(x,t)=\sum_{n=0}^{\infty}L_n^{(\alpha)}(x)t^n}.
$$

Termwise differentiation and reindexing give

$$
\frac{\partial G_\alpha}{\partial t}=\sum_{n=1}^{\infty}nL_n^{(\alpha)}(x)t^{n-1}=\sum_{n=0}^{\infty}(n+1)L_{n+1}^{(\alpha)}(x)t^n.
$$

Using the raising relation,

$$
\frac{\partial G_\alpha}{\partial t}=\sum_{n=0}^{\infty}\left[(n+\alpha+1-x)L_n^{(\alpha)}(x)+x\frac{d}{dx}L_n^{(\alpha)}(x)\right]t^n=t\frac{\partial G_\alpha}{\partial t}+(\alpha+1-x)G_\alpha+x\frac{\partial G_\alpha}{\partial x}.
$$

Therefore,

$$
\boxed{(1-t)\frac{\partial G_\alpha}{\partial t}-x\frac{\partial G_\alpha}{\partial x}=(\alpha+1-x)G_\alpha}.
$$

With the initial condition $G_\alpha(x,0)=1$, the solution is

$$
\boxed{G_\alpha(x,t)=\frac{1}{(1-t)^{\alpha+1}}\exp\left(-\frac{xt}{1-t}\right)}.
$$

### Rodrigues Formula

Expanding the generating function in powers of $t$ gives

$$
L_n^{(\alpha)}(x)=\sum_{k=0}^{n}\frac{(-x)^k}{k!}\frac{\Gamma(n+\alpha+1)}{\Gamma(k+\alpha+1)(n-k)!}.
$$

Applying the Leibniz rule gives

$$
\boxed{L_n^{(\alpha)}(x)=\frac{x^{-\alpha}e^x}{n!}\frac{d^n}{dx^n}\left(e^{-x}x^{n+\alpha}\right)}.
$$

### Recurrence and Derivatives

The generating function gives

$$
\boxed{(n+1)L_{n+1}^{(\alpha)}(x)=(2n+\alpha+1-x)L_n^{(\alpha)}(x)-(n+\alpha)L_{n-1}^{(\alpha)}(x),\qquad \frac{d}{dx}L_n^{(\alpha)}(x)=-L_{n-1}^{(\alpha+1)}(x)}.
$$

Equivalently,

$$
x\frac{d}{dx}L_n^{(\alpha)}(x)=nL_n^{(\alpha)}(x)-(n+\alpha)L_{n-1}^{(\alpha)}(x),\qquad x\frac{d^2}{dx^2}L_n^{(\alpha)}(x)=(x-\alpha-1)\frac{d}{dx}L_n^{(\alpha)}(x)-nL_n^{(\alpha)}(x).
$$

### Orthogonality and Normalization

For $m<n$, the Rodrigues formula and integration by parts give

$$
\int_0^\infty x^\alpha e^{-x}L_m^{(\alpha)}(x)L_n^{(\alpha)}(x)\,dx=\frac{(-1)^n}{n!}\int_0^\infty e^{-x}x^{n+\alpha}\frac{d^n}{dx^n}L_m^{(\alpha)}(x)\,dx=0.
$$

Let $\delta_{mn}$ denote the Kronecker delta. The complete orthogonality relation is

$$
\boxed{\int_0^\infty x^\alpha e^{-x}L_m^{(\alpha)}(x)L_n^{(\alpha)}(x)\,dx=\frac{\Gamma(n+\alpha+1)}{n!}\delta_{mn}}.
$$

The normalization coefficient is

$$
N_n^{(\alpha)}=\sqrt{\frac{n!}{\Gamma(n+\alpha+1)}}.
$$

For $f\in L^2([0,\infty),x^\alpha e^{-x}dx)$,

$$
f(x)=\sum_{n=0}^{\infty}c_nL_n^{(\alpha)}(x),\qquad c_n=\left[N_n^{(\alpha)}\right]^2\int_0^\infty x^\alpha e^{-x}L_n^{(\alpha)}(x)f(x)\,dx.
$$

### Value at the Origin

At $x=0$, the generating function becomes

$$
G_\alpha(0,t)=(1-t)^{-\alpha-1}.
$$

Therefore,

$$
L_n^{(\alpha)}(0)=\frac{\Gamma(n+\alpha+1)}{n!\Gamma(\alpha+1)}.
$$
