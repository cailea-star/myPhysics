### Numerov Method

Consider the coupled second-order differential equation

$$
\vec y''(x)=\hat F(x)\vec y(x),
$$

where $\vec y(x)$ is the channel vector and $\hat F(x)$ is the channel-coupling matrix. On the uniform grid $x_n=x_0+nh$, define $\vec y_n=\vec y(x_n)$ and $\hat F_n=\hat F(x_n)$.

The general Taylor expansion is

$$
\vec y(x_n+\Delta x)=\vec y_n+\vec y_n'\Delta x+\frac{\vec y_n''}{2!}(\Delta x)^2+\frac{\vec y_n^{(3)}}{3!}(\Delta x)^3+\frac{\vec y_n^{(4)}}{4!}(\Delta x)^4+\cdots.
$$

Setting $\Delta x=\pm h$ and adding the two expansions gives

$$
\vec y_{n+1}-2\vec y_n+\vec y_{n-1}=h^2\vec y_n''+\frac{h^4}{12}\vec y_n^{(4)}+\mathcal O(h^6).
$$

From $\vec y''=\hat F\vec y$,

$$
\vec y_n^{(4)}=(\hat F\vec y)_n''=\frac{\hat F_{n+1}\vec y_{n+1}-2\hat F_n\vec y_n+\hat F_{n-1}\vec y_{n-1}}{h^2}+\mathcal O(h^2).
$$

Therefore,

$$
\vec y_{n+1}-2\vec y_n+\vec y_{n-1}=\frac{h^2}{12}\left(\hat F_{n+1}\vec y_{n+1}+10\hat F_n\vec y_n+\hat F_{n-1}\vec y_{n-1}\right)+\mathcal O(h^6).
$$

Hence,

$$
\boxed{\left(\mathbb I-\frac{h^2}{12}\hat F_{n+1}\right)\vec y_{n+1}=\left(2\mathbb I+\frac{5h^2}{6}\hat F_n\right)\vec y_n-\left(\mathbb I-\frac{h^2}{12}\hat F_{n-1}\right)\vec y_{n-1}.}
$$

The recurrence requires two initial values, $\vec y_0$ and $\vec y_1$. Multiple independent solutions may be arranged as matrix columns and propagated by the same recurrence.
