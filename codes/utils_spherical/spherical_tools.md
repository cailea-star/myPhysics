### Spherical Tools

The definitions and conventions of spherical harmonics and Fourier transforms are given in [Spherical Harmonics](../../notes/mathematical_physics/spherical_harmonics.md#scalar-spherical-harmonics) and [Fourier Transform](../../notes/mathematical_physics/fourier_transform.md#spherical-harmonic-decomposition).

##### Fourier–Bessel Transform

Let $f_{lm}(r)$ and $F_{lm}(k)$ be the radial components in coordinate and momentum space, respectively. Their transform pair is

$$
\boxed{F_{lm}(k)=4\pi\int_0^\infty r^2f_{lm}(r)j_l(kr)\,dr},\qquad \boxed{f_{lm}(r)=\frac{1}{2\pi^2}\int_0^\infty k^2F_{lm}(k)j_l(kr)\,dk}.
$$

##### Real Spherical-Harmonic Expansion

The spherical-harmonic expansion and its coefficients are

$$
f(\Omega)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}f_{lm}Y_{lm}(\Omega),\qquad f_{lm}=\int d\Omega\,Y_{lm}^{*}(\Omega)f(\Omega).
$$

The spherical harmonics satisfy

$$
Y_{l,-m}(\Omega)=(-1)^mY_{lm}^{*}(\Omega).
$$

For a real function $f(\Omega)$,

$$
f^{*}(\Omega)=f(\Omega)\ \Longrightarrow\ f_{l,-m}=(-1)^mf_{lm}^{*}.
$$

Combining the $\pm m$ components gives

$$
\boxed{f(\Omega)=\sum_{l=0}^{\infty}f_{l0}Y_{l0}(\Omega)+2\sum_{l=0}^{\infty}\sum_{m=1}^{l}\operatorname{Re}\!\left[f_{lm}Y_{lm}(\Omega)\right]}.
$$
