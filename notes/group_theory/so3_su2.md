### SO(3)

##### Definition and Group Space of SO(3)

Let $\underline x\in\mathbb R^N$.

- Preserve the origin: homogeneous linear transformations

$$
\underline x'=R\underline x.
$$

- Preserve distances: orthogonal transformations

$$
O(N)=\{R\in\mathbb R^{N\times N}:R^{\mathsf T}R=I_N\}.
$$

- Preserve orientation: special orthogonal transformations

$$
SO(N)=\{R\in O(N):\det R=1\}.
$$

Thus,

$$
\boxed{SO(3)=\{R\in\mathbb R^{3\times3}:R^{\mathsf T}R=I_3,\ \det R=1\}}.
$$

- Group parameters

Write a rotation through $\omega$ about $\widehat{\boldsymbol n}$ as

$$
R(\widehat{\boldsymbol n},\omega)=R(\boldsymbol\omega),\qquad \boldsymbol\omega=\omega\widehat{\boldsymbol n},\qquad 0\leq\omega\leq\pi.
$$

The parameter space is a solid ball of radius $\pi$, with antipodal boundary points identified:

$$
R(\widehat{\boldsymbol n},\pi)=R(-\widehat{\boldsymbol n},\pi).
$$

##### Coordinate-Space Representation

The natural coordinate-space representation of $SO(3)$ is

$$
\boxed{D^{(X)}(R)=R,\qquad \underline x'=R\underline x}.
$$

Henceforth, use $R$ directly for the coordinate-space representation matrix.

- Rotation about the $x$-axis

$$
x'=x,\qquad y'=y\cos\omega-z\sin\omega,\qquad z'=y\sin\omega+z\cos\omega.
$$

$$
\boxed{R_x(\omega)=\begin{pmatrix}1&0&0\\0&\cos\omega&-\sin\omega\\0&\sin\omega&\cos\omega\end{pmatrix}\underset{|\omega|\ll1}{=}\begin{pmatrix}1&0&0\\0&1&-\omega\\0&\omega&1\end{pmatrix}+O(\omega^2)=I_3+\omega\begin{pmatrix}0&0&0\\0&0&-1\\0&1&0\end{pmatrix}+O(\omega^2)}.
$$

- Rotation about the $y$-axis

$$
x'=x\cos\omega+z\sin\omega,\qquad y'=y,\qquad z'=-x\sin\omega+z\cos\omega.
$$

$$
\boxed{R_y(\omega)=\begin{pmatrix}\cos\omega&0&\sin\omega\\0&1&0\\-\sin\omega&0&\cos\omega\end{pmatrix}\underset{|\omega|\ll1}{=}\begin{pmatrix}1&0&\omega\\0&1&0\\-\omega&0&1\end{pmatrix}+O(\omega^2)=I_3+\omega\begin{pmatrix}0&0&1\\0&0&0\\-1&0&0\end{pmatrix}+O(\omega^2)}.
$$

- Rotation about the $z$-axis

$$
x'=x\cos\omega-y\sin\omega,\qquad y'=x\sin\omega+y\cos\omega,\qquad z'=z.
$$

$$
\boxed{R_z(\omega)=\begin{pmatrix}\cos\omega&-\sin\omega&0\\\sin\omega&\cos\omega&0\\0&0&1\end{pmatrix}\underset{|\omega|\ll1}{=}\begin{pmatrix}1&-\omega&0\\\omega&1&0\\0&0&1\end{pmatrix}+O(\omega^2)=I_3+\omega\begin{pmatrix}0&-1&0\\1&0&0\\0&0&0\end{pmatrix}+O(\omega^2)}.
$$

Comparing with

$$
R_a(\omega)=I_3-i\omega T_a^{(X)}+O(\omega^2),
$$

gives

$$
\boxed{T_x^{(X)}=\begin{pmatrix}0&0&0\\0&0&-i\\0&i&0\end{pmatrix},\qquad T_y^{(X)}=\begin{pmatrix}0&0&i\\0&0&0\\-i&0&0\end{pmatrix},\qquad T_z^{(X)}=\begin{pmatrix}0&-i&0\\i&0&0\\0&0&0\end{pmatrix}}.
$$

##### Infinitesimal Differential Operators

The general construction is given in [Lie Groups](lie_groups.md#infinitesimal-operators-on-functions). For a scalar function $\psi(\underline x)$,

$$
\boxed{(P_R\psi)(\underline x)=\psi(R^{-1}\underline x)}.
$$

For an infinitesimal rotation,

$$
P_{R(\delta\boldsymbol\omega)}=I-i\sum_{a=x,y,z}\delta\omega_a\widehat T_a+O\!\left(\|\delta\boldsymbol\omega\|^2\right),
$$

where

$$
\widehat T_a=-\sum_{i,j=1}^{3}[T_a^{(X)}]_{ji}x_i\frac{\partial}{\partial x_j}=i\sum_{i,j=1}^{3}\epsilon_{aji}x_i\frac{\partial}{\partial x_j}=-i\sum_{i,j=1}^{3}\epsilon_{aij}x_i\frac{\partial}{\partial x_j}.
$$

Substituting $T_x^{(X)},T_y^{(X)},T_z^{(X)}$ gives

$$
\boxed{\widehat{\boldsymbol L}\equiv\widehat{\boldsymbol T}=-i\boldsymbol x\times\boldsymbol\nabla}.
$$

Explicitly,

$$
\widehat L_x=-i\left(y\frac{\partial}{\partial z}-z\frac{\partial}{\partial y}\right),\qquad \widehat L_y=-i\left(z\frac{\partial}{\partial x}-x\frac{\partial}{\partial z}\right),\qquad \widehat L_z=-i\left(x\frac{\partial}{\partial y}-y\frac{\partial}{\partial x}\right).
$$

### SU(2)

##### Definition and Group Space of SU(2)

Let $\underline z\in\mathbb C^N$.

- Preserve the origin: homogeneous linear transformations

$$
\underline z'=U\underline z.
$$

- Preserve the Hermitian inner product: unitary transformations

$$
U(N)=\{U\in\mathbb C^{N\times N}:U^\dagger U=I_N\}.
$$

- Preserve the unit determinant: special unitary transformations

$$
SU(N)=\{U\in U(N):\det U=1\}.
$$

Thus,

$$
\boxed{SU(2)=\{U\in\mathbb C^{2\times2}:U^\dagger U=I_2,\ \det U=1\}}.
$$

Every $U\in SU(2)$ has the form

$$
U=\begin{pmatrix}a&b\\-b^*&a^*\end{pmatrix},\qquad |a|^2+|b|^2=1.
$$

- Group parameters

Write

$$
U=U(\widehat{\boldsymbol n},\omega)=U(\boldsymbol\omega),\qquad \boldsymbol\omega=\omega\widehat{\boldsymbol n},\qquad 0\leq\omega\leq2\pi.
$$

The parameter space is a solid ball of radius $2\pi$, with all boundary points identified:

$$
U(\widehat{\boldsymbol n},2\pi)=-I_2.
$$

Hence, the group space of $SU(2)$ is topologically $S^3$.

##### Fundamental Representation

The fundamental representation acts on $\underline z\in\mathbb C^2$:

$$
\boxed{D^{(1/2)}(U)=U,\qquad \underline z'=U\underline z}.
$$

For $\widehat{\boldsymbol n}=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$, define

$$
\widehat{\boldsymbol n}\cdot\boldsymbol\sigma\equiv\begin{pmatrix}\cos\theta&\sin\theta e^{-i\phi}\\\sin\theta e^{i\phi}&-\cos\theta\end{pmatrix}.
$$

Its Cartesian components are

$$
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad \sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad \sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$

Since $(\widehat{\boldsymbol n}\cdot\boldsymbol\sigma)^2=I_2$,

$$
\boxed{U(\widehat{\boldsymbol n},\omega)=\exp\left(-i\omega\frac{\widehat{\boldsymbol n}\cdot\boldsymbol\sigma}{2}\right)=I_2\cos\frac{\omega}{2}-i(\widehat{\boldsymbol n}\cdot\boldsymbol\sigma)\sin\frac{\omega}{2}}.
$$

For $|\omega|\ll1$,

$$
\boxed{U(\widehat{\boldsymbol n},\omega)=I_2-i\omega\frac{\widehat{\boldsymbol n}\cdot\boldsymbol\sigma}{2}+O(\omega^2),\qquad \boldsymbol T^{(1/2)}=\frac{\boldsymbol\sigma}{2}}.
$$

##### Homomorphism from SU(2) to SO(3)

For $\underline x=(x_1,x_2,x_3)\in\mathbb R^3$, define the traceless Hermitian matrix

$$
X(\underline x)=\boldsymbol\sigma\cdot\underline x=x_1\sigma_x+x_2\sigma_y+x_3\sigma_z=\begin{pmatrix}x_3&x_1-ix_2\\x_1+ix_2&-x_3\end{pmatrix}.
$$

Conversely,

$$
x_a=\frac12\operatorname{tr}(X\sigma_a),\qquad \det X=-|\underline x|^2.
$$

For $U\in SU(2)$, define

$$
X'=UXU^{-1}=\boldsymbol\sigma\cdot\underline x'.
$$

Conjugation preserves Hermiticity, trace, and determinant:

$$
X'^\dagger=X',\qquad \operatorname{tr}X'=0,\qquad \det X'=\det X.
$$

Therefore,

$$
|\underline x'|^2=|\underline x|^2,\qquad \underline x'=R(U)\underline x.
$$

Thus $R(U)\in O(3)$. Since $SU(2)$ is connected and $R(I_2)=I_3$,

$$
\det R(U)=1,\qquad R(U)\in SO(3).
$$

The representation matrix is determined by

$$
U\sigma_aU^{-1}=\sum_b\sigma_bR_{ba}(U),\qquad R_{ba}(U)=\frac12\operatorname{tr}\!\left(\sigma_bU\sigma_aU^{-1}\right).
$$

Conjugation preserves group multiplication, while $U$ and $-U$ induce the same rotation:

$$
R(UV)=R(U)R(V),\qquad R(U)=R(-U),\qquad \ker R=\{\pm I_2\}.
$$

Therefore,

$$
\boxed{SU(2)/\{\pm I_2\}\simeq SO(3)}.
$$

Thus, $SU(2)$ is the double cover of $SO(3)$.

### Inequivalent Irreducible Representations

##### Axis-Angle and Euler-Angle Parameterizations

- Axis-angle parameters

$$
\widehat{\boldsymbol n}=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta),\qquad R=R(\widehat{\boldsymbol n},\omega),\qquad U=U(\widehat{\boldsymbol n},\omega).
$$

For $SO(3)$,

$$
0\leq\theta\leq\pi,\qquad 0\leq\phi<2\pi,\qquad 0\leq\omega\leq\pi.
$$

For $SU(2)$,

$$
0\leq\theta\leq\pi,\qquad 0\leq\phi<2\pi,\qquad 0\leq\omega\leq2\pi.
$$

- Euler-angle parameters

Using the $z$-$y$-$z$ convention,

$$
\boxed{R(\alpha,\beta,\gamma)=R_z(\alpha)R_y(\beta)R_z(\gamma),\qquad U(\alpha,\beta,\gamma)=U_z(\alpha)U_y(\beta)U_z(\gamma)}.
$$

For $SO(3)$,

$$
0\leq\alpha<2\pi,\qquad 0\leq\beta\leq\pi,\qquad 0\leq\gamma<2\pi.
$$

For $SU(2)$,

$$
0\leq\alpha<2\pi,\qquad 0\leq\beta\leq\pi,\qquad 0\leq\gamma<4\pi.
$$

In the fundamental representation,

$$
U(\alpha,\beta,\gamma)=\begin{pmatrix}e^{-i(\alpha+\gamma)/2}\cos\frac{\beta}{2}&-e^{-i(\alpha-\gamma)/2}\sin\frac{\beta}{2}\\e^{i(\alpha-\gamma)/2}\sin\frac{\beta}{2}&e^{i(\alpha+\gamma)/2}\cos\frac{\beta}{2}\end{pmatrix}.
$$

The two parameterizations describe the same group element:

$$
R(\widehat{\boldsymbol n},\omega)=R(\alpha,\beta,\gamma),\qquad U(\widehat{\boldsymbol n},\omega)=U(\alpha,\beta,\gamma).
$$

Comparing the matrix elements gives

$$
\cos\frac{\omega}{2}=\cos\frac{\beta}{2}\cos\frac{\alpha+\gamma}{2},\qquad \sin\frac{\omega}{2}\cos\theta=\cos\frac{\beta}{2}\sin\frac{\alpha+\gamma}{2}.
$$

$$
\sin\frac{\omega}{2}\sin\theta=\sin\frac{\beta}{2},\qquad \phi=\frac{\pi+\alpha-\gamma}{2}\pmod{2\pi}.
$$

At $\beta=0$ or $\beta=\pi$, the Euler angles are not unique.

##### Irreducible Representations of SU(2)

Let $\underline z=(z_1,z_2)^{\mathsf T}\in\mathbb C^2$. The homogeneous polynomials of total degree $2j$ form the space

$$
V_j=\operatorname{span}\{\psi_\mu^{(j)}:\mu=-j,-j+1,\ldots,j\},\qquad j=0,\frac12,1,\frac32,\ldots.
$$

Choose the basis

$$
\psi_\mu^{(j)}(z_1,z_2)=\frac{(-1)^{j-\mu}z_1^{j-\mu}z_2^{j+\mu}}{\sqrt{(j+\mu)!(j-\mu)!}},\qquad (j-\mu)+(j+\mu)=2j.
$$

Hence,

$$
\dim V_j=2j+1.
$$

For $U\in SU(2)$, let

$$
\underline z'=U^{-1}\underline z.
$$

Since $z_1'$ and $z_2'$ are linear combinations of $z_1$ and $z_2$, $\psi_\mu^{(j)}(\underline z')$ remains homogeneous of total degree $2j$. Thus, $V_j$ is invariant and

$$
\boxed{P_U\psi_\mu^{(j)}(\underline z)=\psi_\mu^{(j)}(\underline z')=\sum_{\nu=-j}^{j}\psi_\nu^{(j)}(\underline z)D_{\nu\mu}^{(j)}(U)},\qquad \underline z'=U^{-1}\underline z.
$$

The coefficients $D_{\nu\mu}^{(j)}(U)$ form the $(2j+1)$-dimensional representation $D^{(j)}$.

In the basis $|j,\mu\rangle$,

$$
\boxed{T_z^{(j)}|j,\mu\rangle=\mu|j,\mu\rangle,\qquad \left(\boldsymbol T^{(j)}\right)^2|j,\mu\rangle=j(j+1)|j,\mu\rangle}.
$$

$$
\boxed{T_\pm^{(j)}=T_x^{(j)}\pm iT_y^{(j)},\qquad T_\pm^{(j)}|j,\mu\rangle=\sqrt{j(j+1)-\mu(\mu\pm1)}\,|j,\mu\pm1\rangle}.
$$

For a rotation through $\omega$,

$$
\chi_j(\omega)=\sum_{\mu=-j}^{j}e^{-i\mu\omega}=\frac{\sin\left[\left(j+\frac12\right)\omega\right]}{\sin\frac{\omega}{2}},\qquad \frac{1}{\pi}\int_0^{2\pi}\chi_i(\omega)^*\chi_j(\omega)\sin^2\frac{\omega}{2}\,d\omega=\delta_{ij}.
$$

The completeness of $\sin[(j+\frac12)\omega]$ shows that these exhaust all inequivalent irreducible representations:

$$
D^{(j)},\qquad j=0,\frac12,1,\frac32,\ldots,\qquad d_j=2j+1.
$$

In particular,

$$
D^{(0)}(U)=1,\qquad D^{(1/2)}(U)=U,\qquad D^{(1)}(U)\simeq R(U).
$$

##### Irreducible Representations of SO(3)

Since $R(U)=R(-U)$,

$$
D^{(j)}(-U)=(-1)^{2j}D^{(j)}(U).
$$

- $j=0,1,2,\ldots$: single-valued irreducible representations of $SO(3)$;
- $j=\frac12,\frac32,\ldots$: double-valued representations of $SO(3)$.

Thus,

$$
\boxed{D^{(j)},\qquad j=0,1,2,\ldots,\qquad d_j=2j+1}.
$$

##### Wigner D Functions

Using the basis order $|j,j\rangle,\ldots,|j,-j\rangle$, define

$$
T_y^{(j)}=\frac{T_+^{(j)}-T_-^{(j)}}{2i},\qquad d^{(j)}(\beta)=e^{-i\beta T_y^{(j)}}.
$$

For the $z$-$y$-$z$ Euler-angle parameterization,

$$
D_{\mu\nu}^{(j)}(\alpha,\beta,\gamma)=\left[e^{-i\alpha T_z^{(j)}}e^{-i\beta T_y^{(j)}}e^{-i\gamma T_z^{(j)}}\right]_{\mu\nu}=e^{-i\mu\alpha}d_{\mu\nu}^{(j)}(\beta)e^{-i\nu\gamma},\qquad \mu,\nu=-j,\ldots,j.
$$

- Composition and unitarity

$$
D^{(j)}(R_1R_2)=D^{(j)}(R_1)D^{(j)}(R_2),\qquad D^{(j)}(R^{-1})=D^{(j)}(R)^\dagger.
$$

- Complex conjugation

$$
D_{\mu\nu}^{(j)}(R)^*=(-1)^{\mu-\nu}D_{-\mu,-\nu}^{(j)}(R).
$$

- Orthogonality

For the normalized Haar measure $dR$,

$$
\int D_{\mu\nu}^{(j)}(R)^*D_{\mu'\nu'}^{(j')}(R)\,dR=\frac{1}{2j+1}\delta_{jj'}\delta_{\mu\mu'}\delta_{\nu\nu'}.
$$

- Spherical harmonics

For integer $\ell$,

$$
\mathcal H_\ell=\operatorname{span}\{Y_{\ell m}:m=-\ell,\ldots,\ell\},\qquad \dim\mathcal H_\ell=2\ell+1.
$$

$$
P_RY_{\ell m}(\widehat{\boldsymbol x})=Y_{\ell m}(R^{-1}\widehat{\boldsymbol x})=\sum_{m'=-\ell}^{\ell}Y_{\ell m'}(\widehat{\boldsymbol x})D_{m'm}^{(\ell)}(R).
$$

In the present convention,

$$
\boxed{Y_{\ell m}(\theta,\phi)=\sqrt{\frac{2\ell+1}{4\pi}}\,D_{m0}^{(\ell)}(\phi,\theta,0)^*}.
$$

Half-integer $j$ gives representations of $SU(2)$; integer $j=\ell$ gives representations of $SO(3)$.


- Low-dimensional cases

$$
T_y^{(0)}=0,\qquad d^{(0)}(\beta)=1.
$$

$$
T_y^{(1/2)}=\frac12\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad d^{(1/2)}(\beta)=\begin{pmatrix}\cos\frac{\beta}{2}&-\sin\frac{\beta}{2}\\\sin\frac{\beta}{2}&\cos\frac{\beta}{2}\end{pmatrix}.
$$

$$
T_y^{(1)}=\frac{1}{\sqrt2}\begin{pmatrix}0&-i&0\\i&0&-i\\0&i&0\end{pmatrix},\qquad d^{(1)}(\beta)=\begin{pmatrix}\cos^2\frac{\beta}{2}&-\frac{\sin\beta}{\sqrt2}&\sin^2\frac{\beta}{2}\\\frac{\sin\beta}{\sqrt2}&\cos\beta&-\frac{\sin\beta}{\sqrt2}\\\sin^2\frac{\beta}{2}&\frac{\sin\beta}{\sqrt2}&\cos^2\frac{\beta}{2}\end{pmatrix}.
$$

### Scalars, Vectors, Tensors, and Spinors

For a field defined on $\mathbb R^3$:

- $R$: acts on spatial coordinates;
- $P_R$: transforms the field function under the coordinate rotation;
- $Q_R$: transforms scalar, vector, or tensor components;
- $O_R=Q_RP_R$: performs the complete field transformation.

For spinor fields, use $U\in SU(2)$ with $R=R(U)$, and replace $Q_R,O_R$ by $Q_U,O_U$.

##### Scalar Fields

A scalar has one component and is unchanged by $Q_R$:

$$
Q_R=I,\qquad D^{(0)}(R)=1,\qquad O_R=Q_RP_R=P_R.
$$

Let $\underline x'$ be the rotated point and $\psi'$ the transformed scalar field:

$$
\boxed{\underline x'=R\underline x,\qquad \psi'=O_R\psi=P_R\psi,\qquad \psi'(\underline x')=\psi(\underline x)}.
$$

Equivalently, at the same coordinate $\underline x$,

$$
\boxed{\psi'(\underline x)=(P_R\psi)(\underline x)=\psi(R^{-1}\underline x)}.
$$

##### Vector and Tensor Fields

- Vector fields

A vector field $\boldsymbol V(\underline x)$ transforms under $D^{(1)}(R)=R$:

$$
Q_RV_a(\underline x)=\sum_{b=1}^{3}R_{ab}V_b(\underline x),\qquad P_RV_a(\underline x)=V_a(R^{-1}\underline x).
$$

With $O_R=Q_RP_R$,

$$
\boxed{\underline x'=R\underline x,\qquad \boldsymbol V'=O_R\boldsymbol V,\qquad V_a'(\underline x')=\sum_{b=1}^{3}R_{ab}V_b(\underline x)}.
$$

Equivalently,

$$
O_RV_a(\underline x)=\sum_{b=1}^{3}R_{ab}V_b(R^{-1}\underline x).
$$

- Tensor fields

A rank-$n$ Cartesian tensor transforms each index as a vector:

$$
Q_RT_{a_1\cdots a_n}(\underline x)=\sum_{b_1,\ldots,b_n=1}^{3}R_{a_1b_1}\cdots R_{a_nb_n}T_{b_1\cdots b_n}(\underline x).
$$

Thus,

$$
\boxed{T_{a_1\cdots a_n}'(\underline x')=\sum_{b_1,\ldots,b_n=1}^{3}R_{a_1b_1}\cdots R_{a_nb_n}T_{b_1\cdots b_n}(\underline x)}.
$$

Equivalently,

$$
O_RT_{a_1\cdots a_n}(\underline x)=\sum_{b_1,\ldots,b_n=1}^{3}R_{a_1b_1}\cdots R_{a_nb_n}T_{b_1\cdots b_n}(R^{-1}\underline x).
$$

The component space carries

$$
D^{(1)}(R)^{\otimes n}=R^{\otimes n}.
$$

Scalars and vectors are rank-$0$ and rank-$1$ tensors, respectively.

##### Spinor Fields

A spin-$s$ field has $2s+1$ components:

$$
\boldsymbol\psi^{(s)}(\underline x)=\left(\psi_s^{(s)}(\underline x),\ldots,\psi_{-s}^{(s)}(\underline x)\right)^{\mathsf T},\qquad s=0,\frac12,1,\ldots.
$$

Their actions are

$$
P_R\psi_\mu^{(s)}(\underline x)=\psi_\mu^{(s)}(R^{-1}\underline x),\qquad Q_U\psi_\mu^{(s)}(\underline x)=\sum_{\nu=-s}^{s}D_{\mu\nu}^{(s)}(U)\psi_\nu^{(s)}(\underline x).
$$

Thus,

$$
\boxed{\underline x'=R(U)\underline x,\qquad \boldsymbol\psi'^{(s)}=O_U\boldsymbol\psi^{(s)},\qquad \psi_\mu'^{(s)}(\underline x')=\sum_{\nu=-s}^{s}D_{\mu\nu}^{(s)}(U)\psi_\nu^{(s)}(\underline x)}.
$$

Equivalently,

$$
\boxed{O_U\psi_\mu^{(s)}(\underline x)=\sum_{\nu=-s}^{s}D_{\mu\nu}^{(s)}(U)\psi_\nu^{(s)}(R(U)^{-1}\underline x)}.
$$

Cartesian tensors contain only integer-spin representations. Since $U$ and $-U$ induce the same spatial rotation,

$$
D^{(s)}(-U)=(-1)^{2s}D^{(s)}(U).
$$

Hence,

- integer $s$: single-valued representations of $SO(3)$;
- half-integer $s$: spinor representations of $SU(2)$, not Cartesian tensors.

For an infinitesimal rotation,

- $\boldsymbol L=-i\boldsymbol x\times\boldsymbol\nabla$: differential operators acting on the field function under coordinate rotations;
- $\boldsymbol S^{(s)}$: generator matrices of $D^{(s)}$ acting on the spinor components.

Therefore,

$$
P_R=I-i\delta\boldsymbol\omega\cdot\boldsymbol L,\qquad Q_U=I-i\delta\boldsymbol\omega\cdot\boldsymbol S^{(s)}.
$$

$$
\boxed{O_U=I-i\delta\boldsymbol\omega\cdot\boldsymbol J,\qquad \boldsymbol J=\boldsymbol L+\boldsymbol S^{(s)}}.
$$

### Irreducible Tensor Operators and the Wigner–Eckart Theorem

##### Irreducible Tensor Operators

Here hats distinguish quantum operators, and $\hbar$ is restored explicitly.

Let $\hat{\mathbf J}$ be the angular-momentum operator defined in [Angular Momentum](../quantum_mechanics/angular_momentum.md), with $\hat J_\pm=\hat J_x\pm i\hat J_y$.

An irreducible tensor operator of rank $k$ consists of $2k+1$ components

$$
\hat T_q^{(k)},\qquad q=-k,\ldots,k.
$$

Under $U\in SU(2)$,

$$
O_U\hat T_q^{(k)}O_U^{-1}=\sum_{q'=-k}^{k}\hat T_{q'}^{(k)}D_{q'q}^{(k)}(U).
$$

Equivalently,

$$
\boxed{[\hat J_z,\hat T_q^{(k)}]=\hbar q\hat T_q^{(k)},\qquad [\hat J_\pm,\hat T_q^{(k)}]=\hbar\sqrt{k(k+1)-q(q\pm1)}\,\hat T_{q\pm1}^{(k)}}.
$$

A rank-$0$ tensor is a rotational scalar:

$$
[\hat J_i,\hat T_0^{(0)}]=0,\qquad i\in\{x,y,z\}.
$$

For a vector operator $\hat{\mathbf A}=(\hat A_x,\hat A_y,\hat A_z)$, the rank-$1$ spherical components are

$$
\hat A_0^{(1)}=\hat A_z,\qquad \hat A_{\pm1}^{(1)}=\mp\frac{\hat A_x\pm i\hat A_y}{\sqrt2}.
$$

A Hermitian spherical tensor satisfies

$$
\left(\hat T_q^{(k)}\right)^\dagger=(-1)^q\hat T_{-q}^{(k)}.
$$

The tensor product of ranks $k_1$ and $k_2$ can be coupled to rank $k$:

$$
\boxed{\left[\hat A^{(k_1)}\otimes\hat B^{(k_2)}\right]_q^{(k)}=\sum_{q_1,q_2}\langle k_1q_1;k_2q_2|kq\rangle\hat A_{q_1}^{(k_1)}\hat B_{q_2}^{(k_2)}}.
$$

The allowed ranks are

$$
|k_1-k_2|\leq k\leq k_1+k_2,\qquad q=q_1+q_2.
$$

##### Wigner–Eckart Theorem

Define

$$
\hat j\equiv\sqrt{2j+1}.
$$

Let $|\alpha_i j_i m_i\rangle$, with $i\in\{1,2\}$, be angular-momentum eigenstates, where $\alpha_i$ denotes all quantum numbers other than $j_i,m_i$.

With the Condon–Shortley convention,

$$
\boxed{\langle\alpha_1j_1m_1|\hat T_q^{(k)}|\alpha_2j_2m_2\rangle=(-1)^{j_1-m_1}\begin{pmatrix}j_1&k&j_2\\-m_1&q&m_2\end{pmatrix}\langle\alpha_1j_1\|\hat T^{(k)}\|\alpha_2j_2\rangle=\frac{\langle j_2m_2;kq|j_1m_1\rangle}{\hat j_1}\langle\alpha_1j_1\|\hat T^{(k)}\|\alpha_2j_2\rangle}.
$$

The reduced matrix element is independent of $m_1,m_2$, and $q$:

$$
\langle\alpha_1j_1\|\hat T^{(k)}\|\alpha_2j_2\rangle.
$$

It can be recovered from the full matrix elements:

$$
\langle\alpha_1j_1\|\hat T^{(k)}\|\alpha_2j_2\rangle=\frac{1}{\hat j_1}\sum_{m_1,m_2,q}\langle j_2m_2;kq|j_1m_1\rangle\langle\alpha_1j_1m_1|\hat T_q^{(k)}|\alpha_2j_2m_2\rangle.
$$

The angular-momentum selection rules are

$$
m_1=m_2+q,\qquad |j_2-k|\leq j_1\leq j_2+k,\qquad j_1+j_2+k\in\mathbb Z.
$$

For a scalar operator,

$$
\langle\alpha_1j_1m_1|\hat T_0^{(0)}|\alpha_2j_2m_2\rangle=\frac{\delta_{j_1j_2}\delta_{m_1m_2}}{\hat j_1}\langle\alpha_1j_1\|\hat T^{(0)}\|\alpha_2j_2\rangle.
$$
