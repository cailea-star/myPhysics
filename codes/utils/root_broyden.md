# Broyden Method

For the nonlinear system

$$
\vec f(\vec x)=0,\qquad \vec f:\mathbb R^n\to\mathbb R^n,
$$

Newton's method is

$$
\vec x_{h+1}=\vec x_h-\hat J_{\vec f}(\vec x_h)^{-1}\vec f(\vec x_h),
$$

Broyden's method approximates the Jacobian $\hat J_{\vec f}$ by $\hat B_h$, or its inverse by $\hat H_h$.

Define the successive differences and impose the secant equation:

$$
\Delta\vec x_h=\vec x_{h+1}-\vec x_h,\qquad \Delta\vec f_h=\vec f(\vec x_{h+1})-\vec f(\vec x_h),\qquad \hat B_{h+1}\Delta\vec x_h=\Delta\vec f_h.
$$

### Derivation of the Broyden Formula

The matrix $\hat B_{h+1}$ must satisfy the secant equation while remaining close to $\hat B_h$. Broyden's method requires $\hat B_{h+1}$ and $\hat B_h$ to act identically for every direction $\vec q\in\mathbb R^n$ orthogonal to $\Delta\vec x_h$:

$$
\Delta\hat B_h\equiv\hat B_{h+1}-\hat B_h,\qquad \Delta\hat B_h\vec q=0,\qquad \Delta\vec x_h^T\vec q=0\quad\Longrightarrow\quad\Delta\hat B_h=\vec u\Delta\vec x_h^T.
$$

Substitution into the secant equation gives

$$
\hat B_{h+1}\Delta\vec x_h=(\hat B_h+\Delta\hat B_h)\Delta\vec x_h=\Delta\vec f_h,\qquad \vec u=\frac{\Delta\vec f_h-\hat B_h\Delta\vec x_h}{\|\Delta\vec x_h\|^2}.
$$

Hence,

$$
\boxed{\hat B_{h+1}=\hat B_h+\frac{(\Delta\vec f_h-\hat B_h\Delta\vec x_h)\Delta\vec x_h^T}{\|\Delta\vec x_h\|^2}}.
$$

For $\hat H_h=\hat B_h^{-1}$, the corresponding inverse update is

$$
\boxed{\hat H_{h+1}=\hat H_h+\frac{(\Delta\vec x_h-\hat H_h\Delta\vec f_h)\Delta\vec x_h^T\hat H_h}{\Delta\vec x_h^T\hat H_h\Delta\vec f_h}}.
$$

### Broyden Iteration Procedure

1. Initialize $\vec x_0$, $\vec f_0=\vec f(\vec x_0)$, and $\hat H_0=-\alpha\hat I$.
2. Compute $\Delta\vec x_h=-\hat H_h\vec f_h$.
3. Compute $\vec x_{h+1}=\vec x_h+\Delta\vec x_h$ and $\vec f_{h+1}=\vec f(\vec x_{h+1})$.
4. Compute $\Delta\vec f_h=\vec f_{h+1}-\vec f_h$ and update $\hat H_{h+1}$ using the inverse Broyden formula.
5. Repeat Steps 2–4 until convergence.

# Modified Broyden Method

For large systems, storing the full inverse-Jacobian approximation is expensive. The modified Broyden method stores only a short history

$$
\mathrm{history}=\{(\Delta\vec x_h,\Delta\vec f_h)\}_{h=1}^{N_h},\qquad N_h\ll n,
$$

where $N_h$ is the history length and $n$ is the vector dimension.

### Derivation of the Modified Broyden Formula

Let $\hat H_0=-\alpha\hat I$, where $\alpha>0$. The matrix $\hat H_h$ must satisfy all stored secant equations while minimizing its correction from $\hat H_0$:

$$
\hat H_h\Delta\vec f_i=\Delta\vec x_i,\qquad i=1,\ldots,N_h,\qquad \hat C_h\equiv\hat H_h-\hat H_0,\qquad \|\hat C_h\|_F\to\min.
$$

Collect the historical differences into

$$
\Delta\hat F=\big[\Delta\vec f_1,\ldots,\Delta\vec f_{N_h}\big],\qquad \Delta\hat X=\big[\Delta\vec x_1,\ldots,\Delta\vec x_{N_h}\big],\qquad \hat H_h\Delta\hat F=\Delta\hat X.
$$

Assuming that the columns of $\Delta\hat F$ are linearly independent, the minimum-norm correction is

$$
\hat C_h\Delta\hat F=\Delta\hat X+\alpha\Delta\hat F\quad\Longrightarrow\quad\hat C_h=(\Delta\hat X+\alpha\Delta\hat F)(\Delta\hat F^T\Delta\hat F)^{-1}\Delta\hat F^T.
$$

Therefore,

$$
\boxed{\hat H_h=-\alpha\hat I+(\Delta\hat X+\alpha\Delta\hat F)(\Delta\hat F^T\Delta\hat F)^{-1}\Delta\hat F^T}.
$$

The quasi-Newton equation $\hat B_h\Delta\vec x_h=-\vec f_h$ gives

$$
\Delta\vec x_h=-\hat H_h\vec f_h=\alpha\vec f_h-(\Delta\hat X+\alpha\Delta\hat F)(\Delta\hat F^T\Delta\hat F)^{-1}\Delta\hat F^T\vec f_h.
$$

Define

$$
\boxed{\vec\gamma=(\Delta\hat F^T\Delta\hat F)^{-1}\Delta\hat F^T\vec f_h\qquad \vec x_{h+1}=\vec x_h+\alpha\vec f_h-(\Delta\hat X+\alpha\Delta\hat F)\vec\gamma}.
$$

### Modified Broyden Iteration Procedure

1. Initialize $\vec x_0$, $\vec f_0=\vec f(\vec x_0)$, and $\mathrm{history}=\varnothing$.
2. Compute the initial linear update $\vec x_1=\vec x_0+\alpha\vec f_0$ and $\vec f_1=\vec f(\vec x_1)$.
3. Store $(\Delta\vec x_0,\Delta\vec f_0)=(\vec x_1-\vec x_0,\vec f_1-\vec f_0)$ in the history.
4. Construct $\Delta\hat X$ and $\Delta\hat F$ from the history, then compute $\vec\gamma$ and $\vec x_{h+1}$ using the modified Broyden formula.
5. Compute $\vec f_{h+1}$, store $(\Delta\vec x_h,\Delta\vec f_h)$, discard the oldest pair when the history exceeds $N_h$, and repeat Steps 4–5 until convergence.

# Fixed Points for Large-Scale Iterations

A fixed point $\vec x^*$ of the map $G:\mathbb R^n\to\mathbb R^n$ satisfies

$$
\vec x^*=G(\vec x^*),\qquad \vec r(\vec x)\equiv G(\vec x)-\vec x,\qquad \vec r(\vec x^*)=0.
$$

Thus, the Broyden method is applied with $\vec f=\vec r$. Define the error $\vec e=\vec x-\vec x^*$. Near the fixed point,

$$
G(\vec x)\approx\vec x^*+\mathbf J_G\vec e,\qquad \mathbf J_G\equiv\left.\frac{\partial G}{\partial\vec x}\right|_{\vec x^*},\qquad \vec r(\vec x)\approx(\mathbf J_G-\mathbf I)\vec e=-\mathbf A\vec e,\qquad \mathbf A\equiv\mathbf I-\mathbf J_G.
$$

### Fixed-Point Convergence Criterion

The fixed-point error satisfies

$$
\vec e_{h+1}=G(\vec x_h)-G(\vec x^*)\approx\mathbf J_G\vec e_h,\qquad \vec e_h\approx\mathbf J_G^h\vec e_0.
$$

Therefore,

$$
\boxed{\rho(\mathbf J_G)\equiv\max_i|\lambda_i(\mathbf J_G)|<1\quad\Longrightarrow\quad\vec e_h\to0}.
$$

The eigenvalues of $\mathbf A=\mathbf I-\mathbf J_G$ satisfy

$$
\mu_i(\mathbf A)=1-\lambda_i(\mathbf J_G),\qquad \rho(\mathbf J_G)<1\quad\Longrightarrow\quad\rho(\mathbf A)<2.
$$

The converse does not generally hold. The equivalent condition expressed through $\mathbf A$ is

$$
\boxed{\rho(\mathbf I-\mathbf A)<1}.
$$

For the exact Newton correction $\Delta\vec x_{\mathrm{Newton}}=-\vec e$,

$$
\vec r^T\Delta\vec x_{\mathrm{Newton}}=(-\mathbf A\vec e)^T(-\vec e)=\vec e^T\mathbf A\vec e.
$$

Thus, $\vec r^T\Delta\vec x_{\mathrm{Newton}}>0$ when the symmetric part of $\mathbf A$ is positive definite. This motivates an alignment safeguard, while fixed-point convergence remains governed by $\rho(\mathbf I-\mathbf A)<1$.

For the relaxed linear iteration $\vec x_{h+1}=\vec x_h+\alpha\vec r_h$,

$$
\vec e_{h+1}\approx(\mathbf I-\alpha\mathbf A)\vec e_h,\qquad \boxed{\rho(\mathbf I-\alpha\mathbf A)<1}.
$$
