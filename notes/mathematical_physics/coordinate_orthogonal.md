# Orthogonal Curvilinear Coordinates

### Coordinate Geometry

Let $(q_1,q_2,q_3)$ be a right-handed orthogonal coordinate system defined by

$$
x=x(q_1,q_2,q_3),\qquad y=y(q_1,q_2,q_3),\qquad z=z(q_1,q_2,q_3).
$$

The position vector is

$$
\mathbf R(q_1,q_2,q_3)=x\hat{\mathbf x}+y\hat{\mathbf y}+z\hat{\mathbf z}.
$$

Write the coordinate-line tangent vectors as

$$
\frac{\partial\mathbf R}{\partial q_i}=\frac{\partial x}{\partial q_i}\hat{\mathbf x}+\frac{\partial y}{\partial q_i}\hat{\mathbf y}+\frac{\partial z}{\partial q_i}\hat{\mathbf z}=h_i\hat{\mathbf e}_i,
$$

$$
\hat{\mathbf e}_i\cdot\hat{\mathbf e}_j=\delta_{ij},\qquad\hat{\mathbf e}_1\times\hat{\mathbf e}_2=\hat{\mathbf e}_3.
$$

Then

$$
\boxed{d\mathbf R=\sum_{i=1}^{3}\frac{\partial\mathbf R}{\partial q_i}\,dq_i=\sum_{i=1}^{3}h_i\hat{\mathbf e}_i\,dq_i,\qquad ds^2=\sum_{i=1}^{3}h_i^2\,dq_i^2}.
$$

For cyclic permutations $(i,j,k)$,

$$
d\mathbf S_i=\hat{\mathbf e}_i\,h_jh_k\,dq_jdq_k,\qquad dV=h_1h_2h_3\,dq_1dq_2dq_3.
$$

### Differential Operators

The gradient operator is

$$
\nabla=\hat{\mathbf e}_1\frac{1}{h_1}\frac{\partial}{\partial q_1}+\hat{\mathbf e}_2\frac{1}{h_2}\frac{\partial}{\partial q_2}+\hat{\mathbf e}_3\frac{1}{h_3}\frac{\partial}{\partial q_3}.
$$

Using $\hat{\mathbf e}_i\cdot$ to select the component of the field to its right, the divergence operator is

$$
\nabla\cdot=\frac{1}{h_1h_2h_3}\left[\frac{\partial}{\partial q_1}\left(h_2h_3\hat{\mathbf e}_1\cdot\right)+\frac{\partial}{\partial q_2}\left(h_3h_1\hat{\mathbf e}_2\cdot\right)+\frac{\partial}{\partial q_3}\left(h_1h_2\hat{\mathbf e}_3\cdot\right)\right].
$$

The curl operator is

$$
\nabla\times=\frac{\hat{\mathbf e}_1}{h_2h_3}\left[\frac{\partial}{\partial q_2}\left(h_3\hat{\mathbf e}_3\cdot\right)-\frac{\partial}{\partial q_3}\left(h_2\hat{\mathbf e}_2\cdot\right)\right]+\frac{\hat{\mathbf e}_2}{h_3h_1}\left[\frac{\partial}{\partial q_3}\left(h_1\hat{\mathbf e}_1\cdot\right)-\frac{\partial}{\partial q_1}\left(h_3\hat{\mathbf e}_3\cdot\right)\right]+\frac{\hat{\mathbf e}_3}{h_1h_2}\left[\frac{\partial}{\partial q_1}\left(h_2\hat{\mathbf e}_2\cdot\right)-\frac{\partial}{\partial q_2}\left(h_1\hat{\mathbf e}_1\cdot\right)\right].
$$

The scalar Laplacian is

$$
\boxed{\nabla^2=\frac{1}{h_1h_2h_3}\left[\frac{\partial}{\partial q_1}\left(\frac{h_2h_3}{h_1}\frac{\partial}{\partial q_1}\right)+\frac{\partial}{\partial q_2}\left(\frac{h_3h_1}{h_2}\frac{\partial}{\partial q_2}\right)+\frac{\partial}{\partial q_3}\left(\frac{h_1h_2}{h_3}\frac{\partial}{\partial q_3}\right)\right]}.
$$
