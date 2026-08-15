# Cubic Clamped B-Spline Basis

### Clamped Knot Vector

Let $[x_{\min},x_{\max}]$ be a finite interval and $\{x_i\}_{i=0}^{N_x-1}$ its $N_x$ distinct partition points:

$$
x_{\min}=x_0<x_1<\cdots<x_{N_x-1}=x_{\max}.
$$

Let $p$ be the polynomial degree. The clamped knot vector $\{t_\alpha\}_{\alpha=0}^{N_t-1}$ repeats each endpoint $p+1$ times and includes each interior partition point once:

$$
t_0=\cdots=t_p=x_{\min},\qquad t_{p+i}=x_i\ (i=1,\ldots,N_x-2),\qquad t_{p+N_x-1}=\cdots=t_{2p+N_x-1}=x_{\max}.
$$

Let $N_t$ be the number of knots and $N_B$ the number of degree-$p$ B-spline basis functions. Then

$$
N_t=N_x+2p,\qquad N_B=N_t-p-1=N_x+p-1.
$$

For cubic B-splines,

$$
p=3,\qquad \mathrm{order}=p+1=4,\qquad N_B=N_x+2.
$$

### B-Spline Basis

Let $\mathbf{1}_A(x)$ denote the indicator function of a set $A$. The degree-zero B-splines and their support intervals are

$$
B_i^{(0)}(x)=\mathbf{1}_{[t_i,t_{i+1})}(x),\qquad \operatorname{supp}B_i^{(0)}=[t_i,t_{i+1}].
$$

For $q=1,\ldots,p$, the Cox–de Boor recursion gives

$$
B_i^{(q)}(x)=\frac{x-t_i}{t_{i+q}-t_i}B_i^{(q-1)}(x)+\frac{t_{i+q+1}-x}{t_{i+q+1}-t_{i+1}}B_{i+1}^{(q-1)}(x),
$$

where a term with a zero denominator is defined as zero. The degree-$p$ basis functions are

$$
B_i(x)\equiv B_i^{(p)}(x),\qquad i=0,\ldots,N_B-1.
$$

With the explicit right-endpoint convention, the endpoint basis functions satisfy

$$
B_0(x_{\min})=1,\qquad B_{N_B-1}(x_{\max})=1.
$$

Each basis function has local support:

$$
\operatorname{supp}B_i=[t_i,t_{i+p+1}].
$$

Hence, at most $p+1$ basis functions are nonzero at any position $x$; for cubic B-splines, at most four are nonzero.

### Basis Expansion

Let $\{c_i\}_{i=0}^{N_B-1}$ be the expansion coefficients. A one-dimensional function is represented as

$$
u(x)=\sum_{i=0}^{N_B-1}c_iB_i(x).
$$

For derivative order $q=1,\ldots,p$,

$$
u^{(q)}(x)=\sum_{i=0}^{N_B-1}c_i\frac{d^qB_i(x)}{dx^q}.
$$

The clamped endpoint properties give

$$
u(x_{\min})=c_0,\qquad u(x_{\max})=c_{N_B-1}.
$$

Therefore, the homogeneous Dirichlet conditions $u(x_{\min})=u(x_{\max})=0$ are imposed by

$$
c_0=c_{N_B-1}=0,\qquad u(x)=\sum_{i=1}^{N_B-2}c_iB_i(x).
$$

### Matrix Elements in a Non-Orthogonal Basis

Let $O(x)$ represent a local multiplication operator and let $B_i'(x)=dB_i(x)/dx$. Define

$$
M_{ij}[O]=\int_{x_{\min}}^{x_{\max}}B_i(x)O(x)B_j(x)\,dx,\qquad K_{ij}[O]=\int_{x_{\min}}^{x_{\max}}B_i'(x)O(x)B_j'(x)\,dx.
$$

Because the B-spline basis is generally non-orthogonal, its overlap matrix is

$$
S_{ij}=M_{ij}[1]=\int_{x_{\min}}^{x_{\max}}B_i(x)B_j(x)\,dx.
$$

Let $T$ and $V$ be the kinetic and potential matrices, respectively. The Hamiltonian matrix is

$$
H=T+V.
$$

Let $E$ be an eigenvalue and $\boldsymbol{c}$ the corresponding coefficient vector. The finite-interval eigenvalue problem is

$$
\boxed{H\boldsymbol{c}=ES\boldsymbol{c}.}
$$

Under homogeneous Dirichlet boundary conditions, all matrices are restricted to the interior basis indices $i,j=1,\ldots,N_B-2$.
