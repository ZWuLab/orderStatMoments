# Lloyd's best linear unbiased estimators of location and scale

For a location-scale family with standard member \\Z\\, the order
statistics at the ranks of `mom` satisfy \\X\_{(r_i)} = \mu + \sigma
Z\_{(r_i)}\\, with \\E\[Z\_{(r_i)}\] = \alpha_i\\ and
\\\mathrm{Cov}(Z\_{(r_i)}, Z\_{(r_j)}) = B\_{ij}\\. Lloyd's estimators
are the generalized least squares estimators of \\(\mu, \sigma)\\ in
this linear model: \\\hat\mu = \bm w\_\mu'\bm X\\ and \\\hat\sigma = \bm
w\_\sigma'\bm X\\ with the weights of the Details. They are unbiased and
have the smallest variance among linear unbiased estimators from these
ranks.

## Usage

``` r
os_blue(mom)
```

## Arguments

- mom:

  An object of class `"orderStatMoments"` computed with `cov = TRUE` at
  two or more ranks. The weights are computed from its `standard`
  component, the moments of the standard member of the family (the
  member with the shape parameters of the call and no shift or
  multiplier), so they do not depend on the other parameters of the
  family or on the `os_shift` and `os_multiplier` it was computed with.

## Value

A list with `w_mu` and `w_sigma` (the weight vectors of \\\hat\mu\\ and
\\\hat\sigma\\, named by rank), `var_mu`, `var_sigma` and `cov` (the
variances and the covariance of the two estimators per unit of the
squared scale parameter, so that \\\mathrm{Var}(\hat\mu) = \sigma^2
\cdot\\ `var_mu`), `ranks`, `n` and `dist`. The names say which
parameter of the location-scale family each output belongs to, \\\mu\\
or \\\sigma\\. These are the unknown parameters that are estimated, as
distinct from the known constants `os_shift` and `os_multiplier` of
[`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md).

## Details

With \\\bm\Omega = \bm B^{-1}\\, \\a\_{11} = \mathbf 1'\bm\Omega\mathbf
1\\, \\a\_{1\alpha} = \mathbf 1'\bm\Omega\bm\alpha\\,
\\a\_{\alpha\alpha} = \bm\alpha'\bm\Omega\bm\alpha\\ and \\\Delta =
a\_{11} a\_{\alpha\alpha} - a\_{1\alpha}^2\\, \$\$\bm w\_\mu' =
(a\_{\alpha\alpha}\mathbf 1' - a\_{1\alpha}\bm\alpha')\bm\Omega/\Delta,
\qquad \bm w\_\sigma' = (a\_{11}\bm\alpha' - a\_{1\alpha}\mathbf
1')\bm\Omega/\Delta,\$\$ \$\$\mathrm{Var}(\hat\mu) = \sigma^2
a\_{\alpha\alpha}/\Delta, \qquad \mathrm{Var}(\hat\sigma) = \sigma^2
a\_{11}/\Delta, \qquad \mathrm{Cov}(\hat\mu, \hat\sigma) = -\sigma^2
a\_{1\alpha}/\Delta.\$\$

**Computation.** The formulas above define the estimators. The function
computes them by whitening the design \\\[\mathbf 1, \bm\alpha\]\\ with
the Cholesky factor of \\\bm B\\ and solving the least-squares problem
by QR, which keeps its accuracy when \\\mathbf 1\\ and \\\bm\alpha\\ are
nearly collinear, as for two close ranks at a large \\n\\, where
\\\Delta\\ itself loses its digits. It then checks the four unbiasedness
conditions, \\\bm w\_\mu'\mathbf 1 = 1\\, \\\bm w\_\mu'\bm\alpha = 0\\,
\\\bm w\_\sigma'\mathbf 1 = 0\\ and \\\bm w\_\sigma'\bm\alpha = 1\\, and
warns when one fails by more than \\10^{-6}\\.

**Conditions.** The standard member of `mom` is taken as the \\Z\\ of
the location-scale family \\\mu + \sigma Z\\. For a built-in family it
is the member with the shape parameters of the call (the degrees of
freedom of a \\t\\, the shape of a gamma or a Weibull, the `sdlog` of a
lognormal), and for a user's distribution it is the distribution that
[`os_dist()`](https://zwulab.github.io/orderStatMoments/reference/os_dist.md)
was given. The required moments exist, which
[`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
has checked; at least two distinct ranks are present, so that the design
\\\[\mathbf 1, \bm\alpha\]\\ has rank two; and \\\bm B\\ is positive
definite. The \\\sigma\\ that is estimated is the multiplier of the
standard member: the standard deviation for the normal, the scale of the
logistic, the Gumbel, the gamma or the Weibull, the length of the
support for the uniform, and \\\exp(\mathrm{meanlog})\\ for the
lognormal, with \\\mu\\ a shift (a threshold) of the family.

## References

Lloyd, E. H. (1952). Least-squares estimation of location and scale
parameters using order statistics. *Biometrika*, 39(1-2), 88-95.

## See also

[`os_moments`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md),
[`os_lin_moments`](https://zwulab.github.io/orderStatMoments/reference/os_lin_moments.md)

## Examples

``` r
m <- os_moments(9, ranks = c(1, 3, 5, 7, 9))
b <- os_blue(m)
b$w_mu                           # sums to one
#>        r1        r3        r5        r7        r9 
#> 0.1457682 0.2382168 0.2320301 0.2382168 0.1457682 
sqrt(b$var_mu)                   # SD of the location estimate per unit sigma
#> [1] 0.3375654
```
