# Runge-Kutta Method

For the first-order differential equation $\vec y'=\vec f(t,\vec y)$, the fourth-order Runge-Kutta iteration is

$$
\vec{y}_{n+1} = \vec{y}_n + \frac{h}{6}(\vec{k}_1 + 2\vec{k}_2 + 2\vec{k}_3 + \vec{k}_4)
$$

where
$$
\vec{k}_1 = \vec{f}(t_n, \vec{y}_n), \quad \vec{k}_2 = \vec{f}(t_n + \frac{h}{2}, \vec{y}_n + \frac{h}{2}\vec{k}_1)
$$

$$
\vec{k}_3 = \vec{f}(t_n + \frac{h}{2}, \vec{y}_n + \frac{h}{2}\vec{k}_2), \quad \vec{k}_4 = \vec{f}(t_n + h, \vec{y}_n + h\vec{k}_3)
$$
