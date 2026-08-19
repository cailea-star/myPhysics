### Log-Derivative Propagation

Consider the fundamental matrix equation

$$
\hat\Psi''(x)=\hat F(x)\hat\Psi(x),\qquad \hat Y(x)=\hat\Psi'(x)\hat\Psi^{-1}(x).
$$

A general differential-equation solver defines the step propagator $\hat P_i$ by

$$
\begin{pmatrix}\hat\Psi_{i+1}\\ \hat\Psi'_{i+1}\end{pmatrix}=\hat P_i\begin{pmatrix}\hat\Psi_i\\ \hat\Psi'_i\end{pmatrix}.
$$

At $x_i$, apply the QR decomposition

$$
\hat\Psi_i=\hat Q(i)\hat R(i),\qquad \hat Q^\dagger(i)\hat Q(i)=\mathbb I,
$$

where $\hat Q(i)$ is an orthonormal solution basis and $\hat R(i)$ contains its scale and column mixing. Holding $\hat R(i)$ fixed during the step, define

$$
\hat\Psi_j(i)\equiv\hat\Psi_j[\hat R(i)]^{-1},\qquad \hat\Psi'_j(i)\equiv\hat\Psi'_j[\hat R(i)]^{-1},
$$

so that $\hat\Psi_i(i)=\hat Q(i)$. The normalized basis is propagated by

$$
\begin{pmatrix}\hat\Psi_{i+1}(i)\\ \hat\Psi'_{i+1}(i)\end{pmatrix}=\hat P_i\begin{pmatrix}\hat\Psi_i(i)\\ \hat\Psi'_i(i)\end{pmatrix}.
$$

The original solution satisfies

$$
\hat\Psi_{i+1}=\hat\Psi_{i+1}(i)\hat R(i),\qquad \hat\Psi'_{i+1}=\hat\Psi'_{i+1}(i)\hat R(i).
$$

Therefore,

$$
\boxed{\hat Y_{i+1}=\hat\Psi'_{i+1}\hat\Psi_{i+1}^{-1}=\hat\Psi'_{i+1}(i)[\hat\Psi_{i+1}(i)]^{-1}.}
$$

The ill-conditioned scale $\hat R(i)$ is removed at every step, while only the orthonormal basis $\hat Q(i)$ and its transformed derivative are propagated.
