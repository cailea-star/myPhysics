# Fourier Transform

## Convention

We use the Fourier-transform convention

$$
\tilde f(\mathbf k)
=
\int d^3\mathbf r\,
f(\mathbf r)e^{-i\mathbf k\cdot\mathbf r},
\qquad
f(\mathbf r)
=
\frac{1}{(2\pi)^3}
\int d^3\mathbf k\,
\tilde f(\mathbf k)e^{i\mathbf k\cdot\mathbf r}.
$$

## Spherical Harmonic Decomposition

The spherical-coordinate conventions used below are defined in [Spherical Coordinates](coordinate_spherical.md).

Expand the coordinate-space function in spherical harmonics:

$$
f(\mathbf r)
=
\sum_{l,m}f_{lm}(r)Y_{lm}(\hat{\mathbf r}),
\qquad
f_{lm}(r)
=
\int d\Omega_{\mathbf r}\,
Y_{lm}^*(\hat{\mathbf r})f(\mathbf r).
$$

Using the plane-wave expansion

$$
e^{i\mathbf k\cdot\mathbf r}
=
4\pi
\sum_{l,m}
i^l\, j_l(kr)
Y_{lm}^*(\hat{\mathbf r})
Y_{lm}(\hat{\mathbf k}),
$$

the momentum-space function becomes

$$
\tilde f(\mathbf k)
=
\int d^3\mathbf r\,
f(\mathbf r)e^{-i\mathbf k\cdot\mathbf r}
$$

$$
=4\pi\sum_{l,m}
i^{-l}\,Y_{lm}(\hat{\mathbf k})
\int_0^\infty r^2\,dr\,
f_{lm}(r)j_l(kr)
$$

$$
=\sum_{l,m}
i^{-l}\,\tilde f_{lm}(k)Y_{lm}(\hat{\mathbf k}),
$$

Then we can expand the momentum-space function in spherical harmonics:

$$
\tilde f(\mathbf k)
=
\sum_{l,m}i^{-l}\tilde f_{lm}(k)Y_{lm}(\hat{\mathbf k}),
\qquad
\tilde f_{lm}(k)
=
4\pi
\int_0^\infty r^2\,dr\,
f_{lm}(r)j_l(kr).
$$

The inverse transform is

$$
f_{lm}(r)
=
\frac{1}{2\pi^2}
\int_0^\infty k^2\,dk\,
\tilde f_{lm}(k)j_l(kr).
$$

Thus, after the spherical harmonic decomposition, the three-dimensional Fourier transform reduces to a one-dimensional spherical Bessel transform for each $(l,m)$ component.
