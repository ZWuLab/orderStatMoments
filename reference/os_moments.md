# Exact moments of order statistics

Computes the means \\\alpha_i = E\[X\_{(r_i)}\]\\, the variances and the
covariance matrix \\B\_{ij} = \mathrm{Cov}(X\_{(r_i)}, X\_{(r_j)})\\ of
the order statistics \\X\_{(1)} \le \dots \le X\_{(n)}\\ of an i.i.d.
sample of size `n` from a continuous distribution, at the ranks `ranks`,
by numerical evaluation of the exact integrals, or from a closed form
where the family has one. The distribution is one of thirteen built-in
families, given as in R by its name and its parameters, or a
distribution that the user defines with
[`os_dist`](https://zwulab.github.io/orderStatMoments/reference/os_dist.md).

## Usage

``` r
os_moments(
  n,
  ranks = seq_len(n),
  dist = "norm",
  ...,
  os_shift = 0,
  os_multiplier = 1,
  cov = TRUE,
  method = c("auto", "integrate", "closed"),
  rel.tol = 1e-11,
  check = TRUE
)
```

## Arguments

- n:

  The sample size, an integer \\\ge 1\\.

- ranks:

  The ranks, strictly increasing integers in \\1, \dots, n\\. The
  default is every rank, which gives the full covariance matrix.

- dist:

  The distribution: the name of a built-in family, `"norm"`, `"unif"`,
  `"exp"`, `"logis"`, `"gumbel"`, `"t"`, `"cauchy"`, `"gamma"`,
  `"chisq"`, `"beta"`, `"lnorm"`, `"weibull"` or `"f"`, or an object
  made by
  [`os_dist`](https://zwulab.github.io/orderStatMoments/reference/os_dist.md).

- ...:

  The parameters of a built-in family, under the names and with the
  defaults of R's own distribution functions (the table in Details), for
  example `mean` and `sd` for `"norm"`, `shape` and `rate` for
  `"gamma"`, `df` for `"t"`. A parameter without a default in R must be
  given. Nothing is passed here when `dist` is an object made by
  [`os_dist()`](https://zwulab.github.io/orderStatMoments/reference/os_dist.md).

- os_shift, os_multiplier:

  The constants \\a\\ and \\b \> 0\\ of \\Y = a + bX\\. The function
  returns the moments of the order statistics of \\Y\\, where \\X\\ has
  the distribution that `dist` and its parameters give: the means become
  \\a + b\\\alpha\\ and the covariances \\b^2 B\\. They serve where the
  distribution has no parameter of its own for a shift or a multiplier:
  a \\t\\ distribution with a location and a scale, a shifted
  exponential, a three-parameter gamma, Weibull or lognormal, or a
  user's distribution.

- cov:

  `TRUE` for the full covariance matrix, `FALSE` for the means and the
  variances only, which need single integrals and no double integral.

- method:

  `"auto"` takes the closed form where the family has one and the
  quadrature elsewhere, `"integrate"` forces the quadrature, `"closed"`
  requires a closed form for every requested quantity and stops
  otherwise.

- rel.tol:

  The relative tolerance of the quadrature. The single integrals and the
  inner integral of a product moment run at `rel.tol`, the outer
  integral at `100 * rel.tol`.

- check:

  Whether to run the run-time checks of the Details and warn when a
  residual exceeds \\10^{-6}\\.

## Value

An object of class `"orderStatMoments"`, a list with the components `n`,
`ranks`, `alpha` (the means, named by rank), `var` (the variances), `B`
(the covariance matrix, `NULL` when `cov = FALSE`), `quad_err` (the
summed error estimates of the quadrature, see Details), `standard` (a
list with the `alpha`, `var` and `B` of the standard member of the
family, which
[`os_blue`](https://zwulab.github.io/orderStatMoments/reference/os_blue.md)
reads), `dist` (the `os_dist` object of the standard member), `par` (the
parameters of a built-in family, `NULL` for a user's distribution),
`os_shift`, `os_multiplier`, `affine` (the shift and the multiplier that
carry the standard member to \\Y\\, which combine the family's own
parameters with `os_shift` and `os_multiplier`), `cov`, `method`,
`route` (`"closed"`, `"integrate"` or `"mixed"`), `checks` (the
residuals of the run-time checks, `NA` where a check did not apply) and
`call`.

## Details

**The built-in families.** A family is given as in R: `dist` is the name
that R's functions carry after `p`, `d` and `q` (`"gamma"` for
`pgamma`), and its parameters are passed by R's own names. Every
argument of `os_moments()` after `dist` is given by its full name.

|  |  |  |
|----|----|----|
| `dist` | Parameters, with R's defaults | Closed forms returned by `method = "auto"` |
| `"norm"` | `mean = 0`, `sd = 1` | none |
| `"unif"` | `min = 0`, `max = 1` | every mean, variance and covariance |
| `"exp"` | `rate = 1` | every mean, variance and covariance |
| `"logis"` | `location = 0`, `scale = 1` | the means and the variances |
| `"gumbel"` | `location = 0`, `scale = 1` | the mean and the variance of the maximum |
| `"t"` | `df` | none |
| `"cauchy"` | `location = 0`, `scale = 1` | none |
| `"gamma"` | `shape`, and `rate = 1` or `scale` | none |
| `"chisq"` | `df` | none |
| `"beta"` | `shape1`, `shape2` | none |
| `"lnorm"` | `meanlog = 0`, `sdlog = 1` | none |
| `"weibull"` | `shape`, `scale = 1` | none |
| `"f"` | `df1`, `df2` | none |

The Gumbel is in the maximum form, \\F(x) = \exp(-e^{-(x -
\mathrm{location})/\mathrm{scale}})\\, and is the one family that R does
not provide. The noncentral versions of the \\t\\, the chi-square, the
beta and the \\F\\ (R's `ncp`) are not built in, and `df1` of the \\F\\
must be finite. The package computes the moments of the standard member
of each family, the member whose shape parameters are those given and
whose other parameters are at R's defaults, and applies the remaining
parameters afterwards, exactly: the mean and the standard deviation of
the normal, the ends of the uniform, the rate of the exponential, the
gamma's rate or scale, the Weibull's scale, and
\\\exp(\mathrm{meanlog})\\ of the lognormal are a shift and a multiplier
of the standard member. A distribution with a large mean or a small rate
therefore loses no digits.

**The integrals.** With cdf \\F\\, density \\f\\ and \\y \< z\\,
\$\$f\_{r:n}(x) = \frac{n!}{(r-1)!\\(n-r)!} F(x)^{r-1} \\1 -
F(x)\\^{n-r} f(x),\$\$ \$\$f\_{r,s:n}(y, z) =
\frac{n!}{(r-1)!\\(s-r-1)!\\(n-s)!} F(y)^{r-1} \\F(z) - F(y)\\^{s-r-1}
\\1 - F(z)\\^{n-s} f(y) f(z),\$\$ and the moments are \$\$\mu\_{r:n} =
\int x f\_{r:n}(x)\\dx, \qquad \mu^{(2)}\_{r:n} = \int x^2
f\_{r:n}(x)\\dx, \qquad \mu\_{r,s:n} = \int\\\\\int\_{y \< z} y z
f\_{r,s:n}(y, z)\\dy\\dz,\$\$ with \\B\_{ij} = \mu\_{r_i, r_j:n} -
\alpha_i \alpha_j\\. The integrals are evaluated by adaptive quadrature
([`integrate`](https://rdrr.io/r/stats/integrate.html)). The numerical
design has three parts.

*What is computed.* The integration is in standardized units of the
standard member, \\t = \\x - Q(1/2)\\/\tilde\sigma\\, where the spread
\\\tilde\sigma\\ is the interquartile range divided by that of the
standard normal (so \\\tilde\sigma = 1\\ for the normal). The absolute
tolerances and the cut points then mean the same for a distribution on
any scale. Each moment is taken about the center of its own order
statistic, \\c_r = Q(p_r)\\ with \\p_r = r/(n+1)\\, so that a variance
is not the small difference of two large moments: taken about the
median, the variance of the uniform minimum kept seven digits at \\n =
1000\\ and none at \\n = 10^5\\.

*How it is integrated.* The range is cut at \\c_r \pm h_r\\, with \\h_r
= d_r + 4/\sqrt{n+1}\\, \\d_r = 10 s_r / f(c_r)\\ and \\s_r =
\sqrt{p_r(1-p_r)/(n+2)}\\, so that the narrow peak of the order
statistic's density at large \\n\\ lies inside one piece whose ends the
adaptive rule sees. The second term is a minimum half-width. The first
moment, whose integrand changes sign, is cut at \\c_r\\ as well, and the
outer range of a product moment is cut at the points of both ranks. The
absolute tolerances are scaled to the width \\d_r\\ of the order
statistic. A piece that the integrator flags is accepted only when its
error estimate is negligible, and the computation stops otherwise (see
"What stops with an error").

*How the integrand is evaluated.* The integrand is assembled on the
logarithmic scale, from
[`lbeta()`](https://rdrr.io/r/base/Special.html), the log-tails
`p(x, log.p = TRUE)` and `p(x, lower.tail = FALSE, log.p = TRUE)` and
the log-density, so that the powers of the tails never underflow. The
difference \\F(z) - F(y)\\ is formed on the side where the two numbers
are small, so that nothing cancels.

**A density that is infinite at an endpoint.** A gamma, a chi-square, a
Weibull, a beta or an \\F\\ distribution with a shape below 1 has a
density that grows like \\e^{a - 1}\\ in the distance \\e\\ from a
finite endpoint of its support. Next to such an endpoint the package
integrates in the variable \\v = e^a\\, in which the integrand is
bounded. The distance from the endpoint is then the integration variable
itself, so an order statistic that sits far closer to the endpoint than
the centered variable \\t\\ can resolve is still placed. The range is
also cut around the peak of the order statistic in that variable. The
package's tests assert the result against exact values: the minimum of a
million observations from a beta(a, 1) at \\a = 0.5\\, 0.1 and 0.05,
within \\10^{-10}\\ for the mean and \\10^{-9}\\ for the variance in
relative terms, and the covariances next to the endpoint at \\n = 1001\\
within \\10^{-7}\\. For a built-in family \\a\\ is the family's own
shape. For a user's distribution whose density is infinite at a lower
endpoint, \\a\\ is estimated from the distribution function, and the
endpoint should be at 0, where the distance is represented without
rounding. The upper endpoint of the beta is treated in the same way
through the reflected variable \\1 - X\\. For a user's distribution with
an infinite density at an upper endpoint, the reflected variable is not
available: the order statistics next to that endpoint are computed less
accurately (about \\10^{-10}\\ for the maximum of a beta(1/2, 1/2)
supplied through
[`os_dist()`](https://zwulab.github.io/orderStatMoments/reference/os_dist.md)
at \\n = 100\\), and the computation stops once they concentrate against
the endpoint. Such a distribution is better supplied reflected, with the
endpoint at 0.

**The moment condition.** A moment is computed only if it exists: with
the tail indices \\(\nu_L, \nu_U)\\ of the distribution,
\\E\|X\_{(r)}\|^2 \< \infty\\ if and only if \\\nu_L r \> 2\\ and
\\\nu_U (n - r + 1) \> 2\\, and a rank that fails stops the computation
with an error. The indices of the built-in families are set by the
package: \\(\mathrm{df}, \mathrm{df})\\ for `"t"`, \\(1, 1)\\ for
`"cauchy"`, \\(\infty, \mathrm{df2}/2)\\ for `"f"`, and infinite for the
others, which always pass. The median of five Cauchy observations has a
variance, and the minimum and the maximum have no mean.

**The error estimate.** `quad_err` is the sum of the error estimates
that [`integrate()`](https://rdrr.io/r/stats/integrate.html) reports for
the single integrals and for the outer integrals of the product moments,
in the units of the standard member. It is an estimate, not a bound: the
inner integral's estimate is not carried into the outer one, and the
errors of the means are not propagated into the covariances. The
accuracy evidence is the agreement with the closed forms, the identities
and the printed tables, documented in the package's tests.

**Accuracy at large n.** With `method = "integrate"`, the means and the
variances of the uniform and the exponential at the minimum, the
quartiles, the median and the maximum agree with their closed forms to
within \\10^{-9}\\ in relative terms at \\n = 10^4\\, \\10^5\\ and
\\10^6\\, and so does the correlation of every pair of them, in absolute
terms. Up to the relative errors of the two variances, the relative
error of a covariance is the error of its correlation divided by the
correlation, so the covariance of two nearly independent order
statistics keeps fewer digits: the two extremes, whose correlation is of
order \\1/n\\, keep five to ten significant digits of their covariance
at those sizes, depending on the distribution and on \\n\\. The
package's tests assert the variances and the covariances of adjacent
ranks at \\10^{-9}\\ and the extremes' covariance at \\10^{-3}\\ beyond
\\n = 10^4\\. A distribution as heavy-tailed as the lognormal with
`sdlog = 2` is accurate to about \\10^{-9}\\ in relative terms rather
than \\10^{-11}\\, because the far tail of the inner integral of a
covariance is resolved at a fixed tolerance. These are the cases that
were measured, not a bound for every distribution.

**What has been tested.** The tests cover three groups of cases. (1) The
normal, uniform, exponential, logistic and Gumbel distributions: the
closed forms of the last four at sample sizes up to 1001, the identities
on full covariance matrices up to \\n = 100\\, the printed tables of the
literature, and the large samples of the paragraph above. (2) The other
eight families at ordinary parameter values: the minimum, the quartiles,
the median and the maximum at sample sizes up to 1001, at the ranks
whose variance exists, against an integral in the probability scale that
uses none of the package's integrands, and the identities on full
matrices at \\n = 10\\ and 25 for the distributions that have a
variance. (3) Cases that were chosen because they are hard. The first
kind is a density that is infinite at an endpoint. Gamma, chi-square,
Weibull, beta and \\F\\ distributions with \\a\\ from 0.8 down to 0.05
were run at sample sizes up to \\10^6 + 1\\. The beta(\\1/K\\, 1) and
beta(1, \\1/K\\) distributions, whose moments are exact finite products,
were run for \\K\\ from 2 to 40 (\\a\\ down to 0.025), at sample sizes
from \\10^3\\ to \\10^7\\ and at ranks from the extreme to \\0.999 n\\.
The beta(\\1/K\\, 0.8) distribution, whose density is infinite at both
endpoints, was run for \\K\\ from 20 to 100 next to its lower endpoint,
and its reflection, the beta(0.8, \\1/K\\) distribution, for \\K\\ of 20
and 30. The second kind is a variance that barely exists: a \\t\\
distribution with 1.005 to 3 degrees of freedom and an \\F\\
distribution with `df2` from 4.02 to 6, at the ranks where the variance
exists and at sample sizes up to \\10^6 + 1\\. In group (3) the package
was asked for the mean and the variance of 25898 order statistics. It
returned 25595 of them, each within \\10^{-9}\\ of its reference (a
variance in relative terms, a mean relative to the larger of the mean
and the standard deviation of its order statistic), and it stopped with
an error for the other 303. Of those, 206 have a variance that rounds to
0 in double precision, and 67 are beta distributions that are refused
when they are built (see "What stops with an error"). A variance below
the smallest normal double, \\2.2 \times 10^{-308}\\, is returned to the
precision of a subnormal number only, about three digits at
\\10^{-320}\\. Of 594 covariances next to an endpoint with an infinite
density, at sample sizes from 1000 to \\10^6\\, the package returned
590, each within \\10^{-9}\\ of its exact value in units of the
correlation (the error divided by the two standard deviations) and in
relative terms, and it stopped for four. The largest sample size in
these tests is \\10^6 + 1\\, apart from the beta(\\1/K\\, 1) and beta(1,
\\1/K\\) distributions, which were run up to \\n = 10^7\\.

**What stops with an error.** The quadrature follows the peak of each
order statistic's density and the tails of the distribution. It does not
follow a second mode: for a mixture of two narrow normals, the mode that
holds a tenth of the mass fell between the nodes, and the moments came
back as those of the other mode. The package therefore integrates the
density of every order statistic it computes by quadrature and stops
when the mass is not 1 within \\10^{-7}\\, naming the rank and the mass
it found. It also stops when a piece of the range cannot be resolved.
This happens for a moment that does not exist when the tail indices of a
user's distribution were not declared, for a tail so heavy that the mass
of the integrand lies far out (the covariances of a lognormal with
`sdlog = 3` or of a Weibull with `shape = 0.2`), for the covariance of
an extreme order statistic whose variance barely exists (the minimum of
20 observations from a \\t\\ distribution with 2.05 to 2.3 degrees of
freedom), and for a variance that comes out zero or negative. A beta
distribution whose second shape is below about 0.03 can be refused when
it is built, with a message that `p` and `q` disagree: its median lies
so close to 1 that the rounding of the median moves the distribution
function, and the two cannot be checked against each other there. Its
order statistics are one minus those of the beta distribution with the
two shapes exchanged, taken at the ranks \\n + 1 - r\\, which the
package computes. A wrong number is never returned knowingly. Extreme
cases cannot be checked exhaustively, and further development can
address the most difficult of them.

**The run-time checks.** For a symmetric distribution and a rank set
closed under \\r \mapsto n + 1 - r\\, the means are antisymmetric and
the second moments symmetric about the median, which the integrators do
not impose; the product moments are then computed for the pairs with
\\r + s \le n + 1\\ only and reflected. For a complete rank set, the
order statistics are a permutation of the sample, so \\\sum_r \alpha_r =
n E\[X\]\\, \\\sum_r \mu^{(2)}\_{r:n} = n E\[X^2\]\\ and \\\sum\_{r,s}
B\_{rs} = n\\\mathrm{Var}(X)\\ for any distribution with a finite second
moment, and for the normal (recognized by its functions `pnorm`, `dnorm`
and `qnorm`, not by its name) every row of \\B\\ sums to \\\sigma^2\\ as
well. The residuals are stored in `checks`, in the standardized units of
the integration, so that the threshold means the same on every scale,
and a residual above \\10^{-6}\\ raises a warning, as does a covariance
matrix with a negative eigenvalue beyond rounding. The checks detect
some failures of the quadrature and do not certify its accuracy: a
failure that every integral shares, such as a peak missed in the same
way at every rank, can satisfy them.

**Cost.** A single integral takes a fraction of a millisecond and a
double integral a few hundredths of a second. When the covariances are
not closed forms, a full matrix needs one double integral per pair,
\\n(n-1)/2\\ of them, halved for a symmetric family.

## References

Wu, Z. and Tong, T. (2026). orderStatMoments: Exact means, variances,
and covariances of order statistics in R. arXiv preprint
arXiv:2611.19980.
[doi:10.48550/arXiv.2611.19980](https://doi.org/10.48550/arXiv.2611.19980)
. The computation, its verification and the comparison with other
software.

Royston, J. P. (1982). Algorithm AS 177: Expected normal order
statistics (exact and approximate). *Applied Statistics*, 31(2),
161-165.

## See also

[`os_dist`](https://zwulab.github.io/orderStatMoments/reference/os_dist.md)

## Examples

``` r
m <- os_moments(5)            # all ranks of n = 5, standard normal
m$alpha
#>            r1            r2            r3            r4            r5 
#> -1.162964e+00 -4.950190e-01  5.413025e-29  4.950190e-01  1.162964e+00 
m$B
#>            r1        r2        r3        r4         r5
#> r1 0.44753407 0.2243310 0.1481477 0.1057720 0.07421527
#> r2 0.22433096 0.3115190 0.2084354 0.1499427 0.10577198
#> r3 0.14814773 0.2084354 0.2868337 0.2084354 0.14814773
#> r4 0.10577198 0.1499427 0.2084354 0.3115190 0.22433096
#> r5 0.07421527 0.1057720 0.1481477 0.2243310 0.44753407
os_moments(20, ranks = c(1, 5, 10, 20), dist = "exp", rate = 0.5)
#> Exact moments of order statistics: exp (rate = 0.5), n = 20, ranks 1 5 10 20
#> alpha (means):
#>        r1        r5       r10       r20 
#> 0.1000000 0.5590213 1.3375430 7.1954790 
#> B (covariance matrix):
#>       r1         r5        r10        r20
#> r1  0.01 0.01000000 0.01000000 0.01000000
#> r5  0.01 0.06289184 0.06289184 0.06289184
#> r10 0.01 0.06289184 0.18558210 0.18558210
#> r20 0.01 0.06289184 0.18558210 6.38465300
#> route: closed; quadrature error estimate: 0
#> checks: reflection n/a; identities n/a
os_moments(100, ranks = c(1, 100), dist = "gumbel", cov = FALSE)
#> Exact moments of order statistics: gumbel (location = 0, scale = 1), n = 100, ranks 1 100
#> alpha (means):
#>        r1      r100 
#> -1.618284  5.182386 
#> var (variances; cov = FALSE):
#>         r1       r100 
#> 0.05442057 1.64493400 
#> route: mixed; quadrature error estimate: 1.1e-13
#> checks: reflection n/a; identities n/a

# a gamma as in rgamma(), and a t distribution with a location and a scale
os_moments(10, ranks = c(1, 5, 10), dist = "gamma", shape = 2, rate = 3)
#> Exact moments of order statistics: gamma (shape = 2, rate = 3), n = 10, ranks 1 5 10
#> alpha (means):
#>        r1        r5       r10 
#> 0.1553405 0.5250292 1.5409860 
#> B (covariance matrix):
#>              r1          r5        r10
#> r1  0.008447579 0.004891628 0.00368668
#> r5  0.004891628 0.024771420 0.01888001
#> r10 0.003686680 0.018880010 0.24793220
#> route: integrate; quadrature error estimate: 1.11e-10
#> checks: reflection n/a; identities n/a
os_moments(7, dist = "t", df = 5, os_shift = 3, os_multiplier = 2, cov = FALSE)
#> Exact moments of order statistics: t (df = 5), os_shift = 3, os_multiplier = 2, n = 7, ranks 1 2 3 4 5 6 7
#> alpha (means):
#>        r1        r2        r3        r4        r5        r6        r7 
#> -0.436658  1.276920  2.224316  3.000000  3.775684  4.723080  6.436658 
#> var (variances; cov = FALSE):
#>        r1        r2        r3        r4        r5        r6        r7 
#> 4.7801480 1.5722180 1.1000210 0.9992729 1.1000210 1.5722180 4.7801480 
#> route: integrate; quadrature error estimate: 5.03e-12
#> checks: reflection 2.27e-17; identities sum_alpha 4.11e-17, sum_sq 3.55e-15

# the Cauchy: the median of five has a variance, the minimum has none
os_moments(5, ranks = 3, dist = "cauchy")
#> Exact moments of order statistics: cauchy (location = 0, scale = 1), n = 5, ranks 3
#> alpha (means):
#>            r3 
#> -1.285952e-17 
#> B (covariance matrix):
#>          r3
#> r3 1.221253
#> route: integrate; quadrature error estimate: 3.89e-14
#> checks: reflection 1.73e-17; identities n/a
try(os_moments(5, ranks = 1, dist = "cauchy"))
#> Error : os_moments(): the moment of order 2 of X_(1) does not exist: the lower tail index 1 times the rank 1 is not above 2
```
