# Associated Legendre Polynomials

Let $l\in\mathbb N_0=\{0,1,2,\ldots\}$. The Legendre polynomial $P_l(x)$ is the degree-$l$ solution on $[-1,1]$ of

$$
\boxed{(1-x^2)\frac{d^2}{dx^2}P_l(x)-2x\frac{d}{dx}P_l(x)+l(l+1)P_l(x)=0}.
$$

Let $m\in\{0,1,\ldots,l\}$. In the Condon–Shortley convention, define the associated Legendre function by

$$
\boxed{P_l^m(x)=(-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}P_l(x)}.
$$

It satisfies

$$
\boxed{(1-x^2)\frac{d^2}{dx^2}P_l^m(x)-2x\frac{d}{dx}P_l^m(x)+\left[l(l+1)-\frac{m^2}{1-x^2}\right]P_l^m(x)=0}.
$$

### Generating Function

With $P_0(x)=1$ and $P_{-1}(x)=0$, the Legendre equation gives the raising relation

$$
(l+1)P_{l+1}(x)=(2l+1)xP_l(x)-lP_{l-1}(x).
$$

Let $t$ be an auxiliary variable. For the ordinary generating function $G(x,t)$,

$$
\boxed{G(x,t)=\sum_{l=0}^{\infty}P_l(x)t^l}.
$$

Termwise differentiation and reindexing give

$$
\frac{\partial G}{\partial t}=\sum_{l=1}^{\infty}lP_l(x)t^{l-1}=\sum_{l=0}^{\infty}(l+1)P_{l+1}(x)t^l.
$$

Using the raising relation,

$$
\frac{\partial G}{\partial t}=\sum_{l=0}^{\infty}\left[(2l+1)xP_l(x)-lP_{l-1}(x)\right]t^l=2xt\frac{\partial G}{\partial t}+xG-t^2\frac{\partial G}{\partial t}-tG.
$$

Therefore,

$$
\boxed{(1-2xt+t^2)\frac{\partial G}{\partial t}=(x-t)G}.
$$

With the initial condition $G(x,0)=1$, the solution is

$$
\boxed{G(x,t)=\frac{1}{\sqrt{1-2xt+t^2}}}.
$$

Using the definition of $P_l^m(x)$ gives

$$
\boxed{G_m(x,t)=\sum_{l=m}^{\infty}P_l^m(x)t^l=(-1)^m(1-x^2)^{m/2}\frac{\partial^mG}{\partial x^m}}.
$$

### Rodrigues Formula

The ordinary Legendre polynomial satisfies

$$
\boxed{P_l(x)=\frac{1}{2^ll!}\frac{d^l}{dx^l}(x^2-1)^l}.
$$

Substitution into the definition of $P_l^m(x)$ gives

$$
\boxed{P_l^m(x)=\frac{(-1)^m}{2^ll!}(1-x^2)^{m/2}\frac{d^{l+m}}{dx^{l+m}}(x^2-1)^l}.
$$

### Recurrence and Derivatives

For fixed $m$, the three-term recurrence and first derivative are

$$
\boxed{(l-m+1)P_{l+1}^m(x)=(2l+1)xP_l^m(x)-(l+m)P_{l-1}^m(x)}.
$$

$$
\boxed{(x^2-1)\frac{d}{dx}P_l^m(x)=lxP_l^m(x)-(l+m)P_{l-1}^m(x)}.
$$

Equivalently,

$$
(1-x^2)\frac{d^2}{dx^2}P_l^m(x)=2x\frac{d}{dx}P_l^m(x)-\left[l(l+1)-\frac{m^2}{1-x^2}\right]P_l^m(x).
$$

### Orthogonality and Normalization

For fixed $m$,

$$
\boxed{\int_{-1}^{1}P_l^m(x)P_k^m(x)\,dx=\frac{2}{2l+1}\frac{(l+m)!}{(l-m)!}\delta_{lk}}.
$$

The normalization coefficient is

$$
N_l^m=\sqrt{\frac{2l+1}{2}\frac{(l-m)!}{(l+m)!}}.
$$

For fixed $m$ and $f\in L^2([-1,1])$,

$$
f(x)=\sum_{l=m}^{\infty}c_lP_l^m(x),\qquad c_l=(N_l^m)^2\int_{-1}^{1}P_l^m(x)f(x)\,dx.
$$

### Parity

The Rodrigues formula gives

$$
P_l^m(-x)=(-1)^{l+m}P_l^m(x).
$$

At the origin,

$$
P_l^m(0)=0\quad\text{for odd }(l+m),\qquad P_l^m(0)=(-1)^{(l+m)/2}\frac{(l+m-1)!!}{(l-m)!!}\quad\text{for even }(l+m).
$$
