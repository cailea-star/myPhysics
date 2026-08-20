# Bayesian Target

For parameters $\mathbf x\in\mathbb R^d$ and data $\mathcal D$, the joint distribution has two factorizations:

$$
p(\mathbf x,\mathcal D)=p(\mathbf x\mid\mathcal D)p(\mathcal D)=p(\mathcal D\mid\mathbf x)p(\mathbf x).
$$

Therefore,

$$
\boxed{p(\mathbf x\mid\mathcal D)=\frac{p(\mathcal D\mid\mathbf x)p(\mathbf x)}{p(\mathcal D)}}.
$$

Here, $p(\mathcal D\mid\mathbf x)$ is the likelihood, $p(\mathbf x)$ is the prior, and $p(\mathcal D)$ is the evidence. Since the evidence is independent of $\mathbf x$,

$$
\pi(\mathbf x)=p(\mathbf x\mid\mathcal D)\propto p(\mathcal D\mid\mathbf x)p(\mathbf x),\qquad \ell(\mathbf x)=\log p(\mathcal D\mid\mathbf x)+\log p(\mathbf x).
$$

The parameter domain is $\mathbf x\in[\mathbf a,\mathbf b]$.

### Random-Walk Metropolis Method

Let $h=0,1,\ldots$ denote the Markov-chain step, $i=1,\ldots,d$ the parameter coordinate, and $x_{h,i}$ the $i$-th component of the state $\mathbf x_h$.

At step $h$, generate the candidate state coordinatewise:

$$
x_{h,i}'=R_i(x_{h,i}+s_i\xi_{h,i}),\qquad s_i=b_i-a_i,\qquad \xi_{h,i}\sim\mathcal N(0,1).
$$

Here, $s_i$ is the random-walk scale and $R_i$ is the reflecting boundary map. For this symmetric proposal, the acceptance probability is

$$
\alpha_h=\min\left(1,\frac{\pi(\mathbf x_h')}{\pi(\mathbf x_h)}\right)=\min\left(1,\exp[\ell(\mathbf x_h')-\ell(\mathbf x_h)]\right).
$$

Thus,

$$
\boxed{\mathbf x_{h+1}=\begin{cases}\mathbf x_h',&\text{with probability }\alpha_h,\\\mathbf x_h,&\text{with probability }1-\alpha_h.\end{cases}}
$$

Equivalently, for $u_h\sim U(0,1)$, accept $\mathbf x_h'$ when

$$
\log u_h<\ell(\mathbf x_h')-\ell(\mathbf x_h).
$$

### Reflecting Boundary

For coordinate $i$, define the interval width $L_i=b_i-a_i$. When $L_i>0$,

$$
y_i=(x_i-a_i)\bmod 2L_i,\qquad y_i\in[0,2L_i).
$$

The reflection map is

$$
\boxed{R_i(x_i)=\begin{cases}a_i+y_i,&y_i\leq L_i,\\a_i+2L_i-y_i,&y_i>L_i.\end{cases}}
$$

When $L_i\leq0$, the coordinate is fixed at $R_i(x_i)=a_i$.

### Sample Statistics

Rejected candidates retain the current state, so repeated states remain part of the Markov chain. Let the recorded samples be $\{\mathbf x^{(k)}\}_{k=1}^{N}$.

Before adding $\mathbf x^{(n+1)}$,

$$
\bar{\mathbf x}_n=\frac{1}{n}\sum_{k=1}^{n}\mathbf x^{(k)},\qquad \boldsymbol\delta_n=\mathbf x^{(n+1)}-\bar{\mathbf x}_n.
$$

The updated mean is

$$
\boxed{\bar{\mathbf x}_{n+1}=\frac{1}{n+1}\sum_{k=1}^{n+1}\mathbf x^{(k)}=\frac{n\bar{\mathbf x}_n+\mathbf x^{(n+1)}}{n+1}=\bar{\mathbf x}_n+\frac{\boldsymbol\delta_n}{n+1}}.
$$

Define the scatter matrix $\mathbf S_n$ and the unbiased sample covariance $\mathbf C_n$:

$$
\mathbf S_n=\sum_{k=1}^{n}(\mathbf x^{(k)}-\bar{\mathbf x}_n)(\mathbf x^{(k)}-\bar{\mathbf x}_n)^T,\qquad \mathbf C_n=\frac{\mathbf S_n}{n-1}.
$$

After adding the new sample, start from

$$
\mathbf S_{n+1}=\sum_{k=1}^{n}(\mathbf x^{(k)}-\bar{\mathbf x}_{n+1})(\mathbf x^{(k)}-\bar{\mathbf x}_{n+1})^T+(\mathbf x^{(n+1)}-\bar{\mathbf x}_{n+1})(\mathbf x^{(n+1)}-\bar{\mathbf x}_{n+1})^T.
$$

The mean update gives

$$
\mathbf x^{(k)}-\bar{\mathbf x}_{n+1}=\mathbf x^{(k)}-\bar{\mathbf x}_n-\frac{\boldsymbol\delta_n}{n+1},\qquad \mathbf x^{(n+1)}-\bar{\mathbf x}_{n+1}=\frac{n}{n+1}\boldsymbol\delta_n.
$$

Substitution yields

$$
\mathbf S_{n+1}=\sum_{k=1}^{n}\left(\mathbf x^{(k)}-\bar{\mathbf x}_n-\frac{\boldsymbol\delta_n}{n+1}\right)\left(\mathbf x^{(k)}-\bar{\mathbf x}_n-\frac{\boldsymbol\delta_n}{n+1}\right)^T+\frac{n^2}{(n+1)^2}\boldsymbol\delta_n\boldsymbol\delta_n^T.
$$

Since

$$
\sum_{k=1}^{n}(\mathbf x^{(k)}-\bar{\mathbf x}_n)=\mathbf0,
$$

the cross terms vanish and

$$
\mathbf S_{n+1}=\mathbf S_n+\frac{n}{(n+1)^2}\boldsymbol\delta_n\boldsymbol\delta_n^T+\frac{n^2}{(n+1)^2}\boldsymbol\delta_n\boldsymbol\delta_n^T=\mathbf S_n+\frac{n}{n+1}\boldsymbol\delta_n\boldsymbol\delta_n^T.
$$

Finally, using $\mathbf S_n=(n-1)\mathbf C_n$ and $\mathbf S_{n+1}=n\mathbf C_{n+1}$,

$$
\boxed{\mathbf C_{n+1}=\frac{n-1}{n}\mathbf C_n+\frac{\boldsymbol\delta_n\boldsymbol\delta_n^T}{n+1},\qquad \mathbf C_1=\mathbf0}.
$$
