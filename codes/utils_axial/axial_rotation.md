# Axial Harmonic-Oscillator Rotations

The oscillator lengths $b_{r_\perp},b_z>0$ follow [Axial Harmonic-Oscillator Basis](axial_basis.md).

### Cartesian Harmonic-Oscillator Rotations

For $n\in\mathbb N_0$ and oscillator length $b>0$, the normalized one-dimensional basis function is

$$
\phi_n(x;b)=\frac{N_n}{\sqrt b}H_n\!\left(\frac{x}{b}\right)e^{-x^2/(2b^2)},\qquad N_n=\frac1{\sqrt{2^nn!\sqrt\pi}}.
$$

The Cartesian basis is

$$
\boxed{\phi_{n_x,n_y,n_z}(x,y,z)=\phi_{n_x}(x;b_{r_\perp})\phi_{n_y}(y;b_{r_\perp})\phi_{n_z}(z;b_z)}.
$$

For an auxiliary variable $t$, the [Hermite generating function](../../notes/mathematical_physics/specialfunction_hermite_polynomials.md#generating-function) gives

$$
\sum_{n=0}^{\infty}H_n(x/b)\frac{t^n}{n!}=e^{2xt/b-t^2}.
$$

Rescaling $t\to t/\sqrt2$ and including the Gaussian normalization gives

$$
\mathcal G_b(x;t)\equiv\sum_{n=0}^{\infty}\frac{t^n}{\sqrt{n!}}\phi_n(x;b)=\frac1{\pi^{1/4}\sqrt b}\exp\left[-\frac{x^2}{2b^2}+\sqrt2\frac{xt}{b}-\frac{t^2}{2}\right].
$$

Let $\boldsymbol r=(x,y,z)^T$ and $\boldsymbol t=(t_x,t_y,t_z)^T$. The three-dimensional generating function is

$$
\mathcal G(\boldsymbol r;\boldsymbol t)=\sum_{n_x,n_y,n_z=0}^{\infty}\frac{t_x^{n_x}t_y^{n_y}t_z^{n_z}}{\sqrt{n_x!n_y!n_z!}}\phi_{n_x,n_y,n_z}(x,y,z).
$$

Multiplying the three one-dimensional generating functions gives

$$
\boxed{\mathcal G(\boldsymbol r;\boldsymbol t)=\frac1{\pi^{3/4}b_{r_\perp}\sqrt{b_z}}\exp\left[-\frac{x^2+y^2}{2b_{r_\perp}^2}-\frac{z^2}{2b_z^2}+\sqrt2\left(\frac{xt_x+yt_y}{b_{r_\perp}}+\frac{zt_z}{b_z}\right)-\frac{t_x^2+t_y^2+t_z^2}{2}\right]}.
$$

A spatial rotation $R$ acts on wavefunctions through inverse coordinates:

$$
[\hat R\phi](\boldsymbol r)=\phi(R^{-1}\boldsymbol r).
$$

For a rotation by $\beta$ about the $y$ axis,

$$
R_y(\beta)\begin{pmatrix}x\\y\\z\end{pmatrix}=\begin{pmatrix}\cos\beta&0&\sin\beta\\0&1&0\\-\sin\beta&0&\cos\beta\end{pmatrix}\begin{pmatrix}x\\y\\z\end{pmatrix}=\begin{pmatrix}x\cos\beta+z\sin\beta\\y\\-x\sin\beta+z\cos\beta\end{pmatrix}.
$$

Therefore,

$$
[\hat R_y(\beta)\mathcal G](\boldsymbol r;\boldsymbol t)=\mathcal G(R_y^{-1}(\beta)\boldsymbol r;\boldsymbol t).
$$

For independent auxiliary vectors $\boldsymbol t_1$ and $\boldsymbol t_2$, define the generating overlap

$$
\boxed{\mathcal K(\boldsymbol t_1,\boldsymbol t_2;\beta)=\int_{\mathbb R^3}\mathrm d^3\boldsymbol r\,\mathcal G(\boldsymbol r;\boldsymbol t_1)\mathcal G(R_y^{-1}(\beta)\boldsymbol r;\boldsymbol t_2)}.
$$

The $y$ integral separates:

$$
\int_{-\infty}^{\infty}\mathrm dy\,\mathcal G_{b_{r_\perp}}(y;t_{1y})\mathcal G_{b_{r_\perp}}(y;t_{2y})=e^{t_{1y}t_{2y}}.
$$

The remaining $xz$ integral is Gaussian and can be evaluated by completing the square. Define the orbital rotation matrix elements by

$$
R^{\mathrm{Cartesian}}_{(n_x,n_y,n_z),(n_x',n_y',n_z')}(\beta)=\langle n_x,n_y,n_z|\hat R_y(\beta)|n_x',n_y',n_z'\rangle.
$$

Expanding both generating functions gives

$$
\mathcal K(\boldsymbol t_1,\boldsymbol t_2;\beta)=\sum_{n_x,n_y,n_z,n_x',n_y',n_z'=0}^{\infty}\frac{t_{1x}^{n_x}t_{1y}^{n_y}t_{1z}^{n_z}t_{2x}^{n_x'}t_{2y}^{n_y'}t_{2z}^{n_z'}}{\sqrt{n_x!n_y!n_z!n_x'!n_y'!n_z'!}}R^{\mathrm{Cartesian}}_{(n_x,n_y,n_z),(n_x',n_y',n_z')}(\beta).
$$

Matching powers determines the rotation matrix elements. Equivalently,

$$
\boxed{R^{\mathrm{Cartesian}}_{(n_x,n_y,n_z),(n_x',n_y',n_z')}(\beta)=\left.\frac{\partial_{t_{1x}}^{n_x}\partial_{t_{1y}}^{n_y}\partial_{t_{1z}}^{n_z}\partial_{t_{2x}}^{n_x'}\partial_{t_{2y}}^{n_y'}\partial_{t_{2z}}^{n_z'}\mathcal K(\boldsymbol t_1,\boldsymbol t_2;\beta)}{\sqrt{n_x!n_y!n_z!n_x'!n_y'!n_z'!}}\right|_{\boldsymbol t_1=\boldsymbol t_2=\boldsymbol0}}.
$$

The factor $e^{t_{1y}t_{2y}}$ enforces $n_y=n_y'$. At $\beta=0$,

$$
\mathcal K(\boldsymbol t_1,\boldsymbol t_2;0)=e^{\boldsymbol t_1^T\boldsymbol t_2},\qquad R^{\mathrm{Cartesian}}_{(n_x,n_y,n_z),(n_x',n_y',n_z')}(0)=\delta_{n_xn_x'}\delta_{n_yn_y'}\delta_{n_zn_z'}.
$$

### Cartesian Harmonic-Oscillator Rotation Kernel

The generating overlap follows from the product
$$
\mathcal G(\boldsymbol r;\boldsymbol t_1)\mathcal G(R_y^{-1}(\beta)\boldsymbol r;\boldsymbol t_2)=\frac1{\pi^{3/2}b_{r_\perp}^2b_z}\exp\left[-\frac{x^2+2y^2+(x\cos\beta-z\sin\beta)^2}{2b_{r_\perp}^2}-\frac{z^2+(x\sin\beta+z\cos\beta)^2}{2b_z^2}+\sqrt2\left(\frac{xt_{1x}+y(t_{1y}+t_{2y})+(x\cos\beta-z\sin\beta)t_{2x}}{b_{r_\perp}}+\frac{zt_{1z}+(x\sin\beta+z\cos\beta)t_{2z}}{b_z}\right)-\frac{\boldsymbol t_1^T\boldsymbol t_1+\boldsymbol t_2^T\boldsymbol t_2}{2}\right].
$$

Integrating over $y$ produces $e^{t_{1y}t_{2y}}$. Collect the remaining coordinate terms into
$$
M=\begin{pmatrix}\dfrac{1+\cos^2\beta}{b_{r_\perp}^2}+\dfrac{\sin^2\beta}{b_z^2}&\sin\beta\cos\beta\left(\dfrac1{b_z^2}-\dfrac1{b_{r_\perp}^2}\right)\\\sin\beta\cos\beta\left(\dfrac1{b_z^2}-\dfrac1{b_{r_\perp}^2}\right)&\dfrac{\sin^2\beta}{b_{r_\perp}^2}+\dfrac{1+\cos^2\beta}{b_z^2}\end{pmatrix},
$$
$$
B_1=\begin{pmatrix}b_{r_\perp}^{-1}&0\\0&b_z^{-1}\end{pmatrix},\qquad B_2=\begin{pmatrix}b_{r_\perp}^{-1}\cos\beta&b_z^{-1}\sin\beta\\-b_{r_\perp}^{-1}\sin\beta&b_z^{-1}\cos\beta\end{pmatrix}.
$$

The remaining integral is
$$
\mathcal K=\frac{\exp\left[t_{1y}t_{2y}-\frac12\sum_{a=1}^2(t_{ax}^2+t_{az}^2)\right]}{\pi b_{r_\perp}b_z}\int_{\mathbb R^2}\mathrm dx\,\mathrm dz\,\exp\left[-\frac12\begin{pmatrix}x&z\end{pmatrix}M\begin{pmatrix}x\\z\end{pmatrix}+\sqrt2\begin{pmatrix}x&z\end{pmatrix}\sum_{a=1}^2B_a\begin{pmatrix}t_{ax}\\t_{az}\end{pmatrix}\right].
$$

Since $M$ is positive definite, completing the square gives
$$
\mathcal K=\frac{2}{b_{r_\perp}b_z\sqrt{\det M}}\exp\left[t_{1y}t_{2y}-\frac12\sum_{a=1}^2(t_{ax}^2+t_{az}^2)+\sum_{a,b=1}^2\begin{pmatrix}t_{ax}&t_{az}\end{pmatrix}B_a^TM^{-1}B_b\begin{pmatrix}t_{bx}\\t_{bz}\end{pmatrix}\right].
$$

Define the ground-state overlap $K_0$ and the $2\times2$ blocks $H_{ab}$:
$$
K_0(\beta)=\frac{2}{b_{r_\perp}b_z\sqrt{\det M}},\qquad H_{ab}(\beta)=2B_a^TM^{-1}B_b-\delta_{ab}I_2,\qquad a,b\in\{1,2\}.
$$

Collect the $xz$ source components and the corresponding blocks into
$$
\boldsymbol T=(t_{1x},t_{1z},t_{2x},t_{2z})^T,\qquad H(\beta)=\begin{pmatrix}H_{11}(\beta)&H_{12}(\beta)\\H_{21}(\beta)&H_{22}(\beta)\end{pmatrix}.
$$

Writing $H_{a\mu,b\nu}\equiv(H_{ab})_{\mu\nu}$, with $\mu,\nu\in\{x,z\}$, the generating overlap becomes
$$
\boxed{\mathcal K(\boldsymbol t_1,\boldsymbol t_2;\beta)=K_0(\beta)\exp\left[t_{1y}t_{2y}+\frac12\boldsymbol T^TH(\beta)\boldsymbol T\right]}.
$$

Define the normalized $xz$ coefficients by
$$
e^{-t_{1y}t_{2y}}\mathcal K=\sum_{n_{x1},n_{z1},n_{x2},n_{z2}=0}^{\infty}\frac{d_{n_{x1},n_{z1},n_{x2},n_{z2}}(\beta)}{\sqrt{n_{x1}!n_{z1}!n_{x2}!n_{z2}!}}t_{1x}^{n_{x1}}t_{1z}^{n_{z1}}t_{2x}^{n_{x2}}t_{2z}^{n_{z2}}.
$$

Differentiating with respect to a source component gives
$$
\frac{\partial\mathcal K}{\partial t_{a\mu}}=\sum_{b=1}^2\sum_{\nu\in\{x,z\}}H_{a\mu,b\nu}t_{b\nu}\mathcal K.
$$

For example, comparing powers after differentiation with respect to $t_{1x}$ gives
$$
\boxed{\sqrt{n_{x1}+1}\,d_{n_{x1}+1,n_{z1},n_{x2},n_{z2}}=H_{1x,1x}\sqrt{n_{x1}}\,d_{n_{x1}-1,n_{z1},n_{x2},n_{z2}}+H_{1x,1z}\sqrt{n_{z1}}\,d_{n_{x1},n_{z1}-1,n_{x2},n_{z2}}+H_{1x,2x}\sqrt{n_{x2}}\,d_{n_{x1},n_{z1},n_{x2}-1,n_{z2}}+H_{1x,2z}\sqrt{n_{z2}}\,d_{n_{x1},n_{z1},n_{x2},n_{z2}-1}}.
$$

The other three derivatives give analogous recurrences. Terms with negative indices vanish, and
$$
d_{0,0,0,0}=K_0.
$$

The full Cartesian matrix element is
$$
\boxed{R^{\mathrm{Cartesian}}_{(n_{x1},n_{y1},n_{z1}),(n_{x2},n_{y2},n_{z2})}(\beta)=\delta_{n_{y1}n_{y2}}d_{n_{x1},n_{z1},n_{x2},n_{z2}}(\beta)}.
$$

At $\beta=0$,
$$
K_0=1,\qquad H_{11}=H_{22}=0,\qquad H_{12}=H_{21}=I_2,\qquad d_{n_{x1},n_{z1},n_{x2},n_{z2}}(0)=\delta_{n_{x1}n_{x2}}\delta_{n_{z1}n_{z2}}.
$$

### Cartesian-to-Cylindrical Transformation

Using the radial functions in [Axial Harmonic-Oscillator Basis](axial_basis.md), define

$$
\langle r_\perp,\varphi|n_r,\Lambda\rangle=\phi_{n_r}^{|\Lambda|}(r_\perp)\frac{e^{i\Lambda\varphi}}{\sqrt{2\pi}},\qquad n_r\in\mathbb N_0,\quad\Lambda\in\mathbb Z.
$$

The transverse Cartesian creation operators and circular combinations are

$$
\hat a_x^\dagger=\frac{x/b_{r_\perp}-b_{r_\perp}\partial_x}{\sqrt2},\qquad \hat a_y^\dagger=\frac{y/b_{r_\perp}-b_{r_\perp}\partial_y}{\sqrt2},\qquad \hat a_\pm^\dagger=\frac{\hat a_x^\dagger\pm i\hat a_y^\dagger}{\sqrt2}.
$$

Since $\hat L_z=\hbar(\hat a_+^\dagger\hat a_+-\hat a_-^\dagger\hat a_-)$, the circular occupations are

$$
n_+=n_r+\frac{|\Lambda|+\Lambda}{2},\qquad n_-=n_r+\frac{|\Lambda|-\Lambda}{2}.
$$

For the transverse ground state $|0,0\rangle$, matching the radial Laguerre convention gives

$$
|n_r,\Lambda\rangle=\frac{(-1)^{n_r}(\hat a_+^\dagger)^{n_+}(\hat a_-^\dagger)^{n_-}}{\sqrt{n_+!n_-!}}|0,0\rangle=\frac{(-1)^{n_r}(\hat a_x^\dagger+i\hat a_y^\dagger)^{n_+}(\hat a_x^\dagger-i\hat a_y^\dagger)^{n_-}}{\sqrt{2^{n_++n_-}n_+!n_-!}}|0,0\rangle.
$$

Expanding both binomials and extracting the phase $(-i)^{n_y}$ defines real coefficients $C_{n_xn_y;n_r\Lambda}$:

$$
\boxed{|n_r,\Lambda\rangle=\sum_{\substack{n_x,n_y\ge0\\n_x+n_y=2n_r+|\Lambda|}}(-i)^{n_y}C_{n_xn_y;n_r\Lambda}|n_x,n_y\rangle},
$$

$$
\boxed{C_{n_xn_y;n_r\Lambda}=(-1)^{n_r}\sqrt{\frac{n_x!n_y!}{2^{n_++n_-}n_+!n_-!}}\sum_{\substack{0\le k\le n_+,\;0\le l\le n_-\\k+l=n_x}}\binom{n_+}{k}\binom{n_-}{l}(-1)^{n_+-k}\in\mathbb R}.
$$

The $z$-direction is unchanged, so the same coefficients transform the full basis:

$$
|n_z,n_r,\Lambda\rangle=\sum_{\substack{n_x,n_y\ge0\\n_x+n_y=2n_r+|\Lambda|}}(-i)^{n_y}C_{n_xn_y;n_r\Lambda}|n_x,n_y,n_z\rangle.
$$

### Axial-Basis Rotation Matrix

For $a=(n_{za},n_{ra},\Lambda_a,\Sigma_a)$, define

$$
|\phi_a\rangle=|n_{za},n_{ra},\Lambda_a\rangle|\Sigma_a\rangle,\qquad \Sigma_a=\pm\frac12,\qquad \Omega_a=\Lambda_a+\Sigma_a.
$$

Rotation about $y$ enforces $n_y=n_y'$, so the expansion phases cancel:

$$
i^{n_y}(-i)^{n_y'}=1\qquad(n_y=n_y').
$$

The orbital rotation therefore uses only the real coefficients:

$$
R^{\mathrm{orbital}}_{ab}(\beta)=\sum_{\substack{n_x+n_y=2n_{ra}+|\Lambda_a|\\n_x'+n_y'=2n_{rb}+|\Lambda_b|}}C_{n_xn_y;n_{ra}\Lambda_a}R^{\mathrm{Cartesian}}_{(n_x,n_y,n_{za}),(n_x',n_y',n_{zb})}(\beta)C_{n_x'n_y';n_{rb}\Lambda_b},
$$

where all summation indices are nonnegative.

Since $\hat J_y=\hat L_y+\hat S_y$, orbital and spin rotations factorize:

$$
\boxed{d_{ab}(\beta)\equiv\langle\phi_a|e^{-i\beta\hat J_y/\hbar}|\phi_b\rangle=R^{\mathrm{orbital}}_{ab}(\beta)d^{1/2}_{\Sigma_a\Sigma_b}(\beta)},\qquad d^{1/2}(\beta)=\begin{pmatrix}\cos(\beta/2)&-\sin(\beta/2)\\\sin(\beta/2)&\cos(\beta/2)\end{pmatrix}.
$$

The spin ordering is $(+\frac12,-\frac12)$. For Euler angles $(\alpha,\beta,\gamma)$,

$$
\boxed{D_{ab}(\alpha,\beta,\gamma)\equiv\langle\phi_a|e^{-i\alpha\hat J_z/\hbar}e^{-i\beta\hat J_y/\hbar}e^{-i\gamma\hat J_z/\hbar}|\phi_b\rangle=e^{-i\alpha\Omega_a}d_{ab}(\beta)e^{-i\gamma\Omega_b}}.
$$

Include both signs of $\Omega_a$. A truncated deformed basis generally gives a nonunitary restricted rotation matrix.

### Rotation in the Simplex Basis

Define the simplex operator using parity $\hat\Pi$:
$$
\hat{\mathcal S}_y\equiv\hat\Pi e^{-i\pi\hat J_y/\hbar},\qquad \hat{\mathcal S}_y^2=-I.
$$

For each state $|\phi_a\rangle=|n_{za},n_{ra},\Lambda_a\rangle|\Sigma_a\rangle$ with $\Omega_a>0$, define the reversed-label partner without an additional phase:
$$
|\phi_{\bar a}\rangle=|n_{za},n_{ra},-\Lambda_a\rangle|-\Sigma_a\rangle.
$$

The spatial reflection $y\to-y$ and spin rotation give
$$
\hat{\mathcal S}_y|\phi_a\rangle=2\Sigma_a|\phi_{\bar a}\rangle,\qquad \hat{\mathcal S}_y|\phi_{\bar a}\rangle=-2\Sigma_a|\phi_a\rangle.
$$

Hence the simplex eigenstates are
$$
\boxed{|\phi_{a,s}\rangle=\frac{|\phi_a\rangle-is(2\Sigma_a)|\phi_{\bar a}\rangle}{\sqrt2},\qquad \hat{\mathcal S}_y|\phi_{a,s}\rangle=is|\phi_{a,s}\rangle,\qquad s=\pm1}.
$$

For $n$ positive-$\Omega$ states, define $E_{ab}=2\Sigma_a\delta_{ab}$. The basis transformation is
$$
T=\frac1{\sqrt2}\begin{pmatrix}I_n&I_n\\-iE&iE\end{pmatrix},\qquad T^\dagger T=I_{2n}.
$$

Since $[\hat{\mathcal S}_y,e^{-i\beta\hat J_y/\hbar}]=0$, the rotation separates into simplex blocks:
$$
\boxed{T^\dagger d(\beta)T=\begin{pmatrix}d_+(\beta)&0\\0&d_-(\beta)\end{pmatrix},\qquad d_-(\beta)=d_+(\beta)^*}.
$$

The conjugation relation follows from the real axial rotation matrix and the chosen phases. Gauge rotation multiplies both blocks by $e^{i\varphi}$.

This decomposition requires a simplex-invariant retained space; it does not imply unitarity of the truncated rotation matrix.
