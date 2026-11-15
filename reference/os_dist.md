# Define a distribution for the moments of its order statistics

Builds the distribution object that
[`os_moments`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
takes in its `dist` argument, from the distribution function, the
density and the quantile function of a continuous distribution that is
not among the built-in families. The three functions define the
distribution itself, with its parameters fixed inside them.
[`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
computes the moments of the order statistics of this distribution, and
its arguments `os_shift` and `os_multiplier` give those of \\a + bX\\.

## Usage

``` r
os_dist(
  p,
  d,
  q,
  support = c(-Inf, Inf),
  symmetric = FALSE,
  tail = c(Inf, Inf),
  name = "user",
  breaks = NULL
)
```

## Arguments

- p:

  The distribution function, with the signature of R's `p*` functions:
  `p(q, lower.tail = TRUE, log.p = FALSE)`, vectorized in `q`. The
  `log.p` and `lower.tail` arguments must be honored, because the
  integrands are evaluated on the log scale in both tails.

- d:

  The density, `d(x, log = FALSE)`, vectorized in `x`, with `log`
  honored.

- q:

  The quantile function, `q(p)`, vectorized in `p`.

- support:

  Numeric vector of length two, the endpoints of the support, finite or
  infinite. They are the limits of integration.

- symmetric:

  Logical. `TRUE` if the density is symmetric about its median, that is
  \\Q(1 - u) = 2 Q(1/2) - Q(u)\\. A symmetric distribution lets
  [`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
  halve the number of double integrals of a full matrix and check the
  reflection symmetry of its result.

- tail:

  Numeric vector of length two, the tail indices \\(\nu_L, \nu_U)\\: the
  density decays like \\\|x\|^{-\nu_L - 1}\\ in the lower tail and like
  \\x^{-\nu_U - 1}\\ in the upper tail. `Inf` for a tail that decays
  faster than any power, or for a bounded tail. The indices decide which
  moments exist (see Details).

- name:

  A label for printing.

- breaks:

  Numeric vector, or `NULL` for none: the points strictly inside the
  support at which the density is not smooth, that is, where it has a
  kink or a jump (the median of the Laplace distribution, the mode of a
  triangular distribution, the knots of a piecewise density). The
  density must be finite there, and it may be zero.
  [`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
  cuts every integral at these points (see Details). For a symmetric
  distribution the points must be symmetric about the median.

## Value

An object of class `"os_dist"`: a list with the components `name`, `p`,
`d`, `q`, `support`, `symmetric`, `tail`, `center` (the median
\\Q(1/2)\\), `closed` (the closed-form moments the distribution has,
none for a user's distribution) and `breaks` (sorted, `NULL` when there
are none).

## Details

The thirteen built-in families need no call of `os_dist()`: the normal,
uniform, exponential, logistic, Gumbel, Student \\t\\, Cauchy, gamma,
chi-square, beta, lognormal, Weibull and \\F\\ distributions are given
to
[`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
by name, with R's own parameters, and the package sets their supports
and tail indices.

**The moment condition.** For a distribution whose density decays like
\\\|x\|^{-\nu - 1}\\, the density of \\X\_{(r)}\\ decays like
\\\|x\|^{-\nu r - 1}\\ in the lower tail and like \\x^{-\nu (n -
r + 1) - 1}\\ in the upper tail, so \\E\|X\_{(r)}\|^k \< \infty\\ if and
only if \\\nu_L r \> k\\ and \\\nu_U (n - r + 1) \> k\\.
[`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
applies this condition with \\k = 2\\ to every requested rank and
refuses a moment that does not exist. The package cannot check the
declared indices. The default `tail = c(Inf, Inf)` declares that every
moment exists, so the indices of a heavy-tailed distribution must be
given.

**Accuracy requirement.** The log-tails `p(x, log.p = TRUE)` and
`p(x, lower.tail = FALSE, log.p = TRUE)` must keep relative precision
deep in both tails, as R's base distributions do. A function that
computes `log(p(x))` from an underflowing probability loses the tail,
and the integrals with it. `os_dist()` checks that the three functions
honor their arguments at the median and are vectorized, and that `d`
integrates to \\1/2\\ between the quartiles, where `p` rises by \\1/2\\.
[`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
checks that the log-tails and the log-density are finite at the median
and at the cut points of the smallest and the largest requested rank,
and that the density of every order statistic it integrates has mass 1.
None of these certifies accuracy.

**A distribution far from the origin.** The integrals are taken about
the median of the supplied distribution, in units of its spread (its
interquartile range divided by that of the standard normal), so the
scale on which the distribution is supplied does not matter: the
package's tests check a normal supplied with standard deviations from
\\10^{-3}\\ to \\10^6\\. Its position does: a distribution whose median
is large against its spread, say \\10^8\\ with unit spread, loses the
digits of \\x - Q(1/2)\\ in double precision and the computation stops
with an error. Such a distribution is supplied near the origin, and the
shift is passed to
[`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
as `os_shift`, which is applied exactly. This also keeps
[`os_blue()`](https://zwulab.github.io/orderStatMoments/reference/os_blue.md)
exact.

**A density that is infinite at an endpoint.** When the density is
infinite at a finite lower endpoint of the support, as for a power or a
gamma-type density with a shape below 1,
[`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
integrates next to that endpoint in the distance from it, raised to a
power that it estimates from the distribution function. The distance is
represented without rounding when the endpoint is 0, so such a
distribution is supplied with its endpoint at 0. An infinite density at
an upper endpoint has no such treatment for a user's distribution, and
the order statistics next to it are computed less accurately or refused.
Supplying the reflected variable, whose singular endpoint is the lower
one, brings that endpoint under the treatment above.

**A density that is not smooth inside the support.** The quadrature rule
assumes an integrand that is smooth inside each piece of the integration
range. A point where the density has a kink or a jump, such as the
median of the Laplace distribution, is therefore declared in `breaks`,
and
[`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
then cuts every single integral and both levels of every double integral
there, so that each piece is smooth. Without the declaration the
covariances lose accuracy: for the Laplace example below, at \\n = 5\\,
the sum of the covariance matrix misses its exact value by about
\\10^{-8}\\, and with `breaks = 0` by less than \\10^{-12}\\. The means
and the variances are affected less. The error estimate `quad_err` does
not show the loss, and the residual `sum_B` of the run-time checks,
which is available for a complete rank set, does. The package cannot
find such points itself. The built-in families have none.

## See also

[`os_moments`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)

## Examples

``` r
# A Pareto distribution with shape 3 on [1, Inf): the upper tail index is 3,
# so the variance of the maximum exists (3 * 1 > 2)
pareto <- function(shape) os_dist(
  p = function(q, lower.tail = TRUE, log.p = FALSE) {
    ls <- -shape * log(pmax(q, 1))                  # log S(x) = -shape * log(x)
    if (lower.tail) { if (log.p) log(-expm1(ls)) else -expm1(ls) }
    else if (log.p) ls else exp(ls)
  },
  d = function(x, log = FALSE) {
    ld <- ifelse(x >= 1, log(shape) - (shape + 1) * log(x), -Inf)
    if (log) ld else exp(ld)
  },
  q = function(p) (1 - p)^(-1 / shape),
  support = c(1, Inf), tail = c(Inf, shape), name = "pareto")
os_moments(5, dist = pareto(3), cov = FALSE)
#> Exact moments of order statistics: pareto, n = 5, ranks 1 2 3 4 5
#> alpha (means):
#>       r1       r2       r3       r4       r5 
#> 1.071429 1.168831 1.314935 1.577922 2.366883 
#> var (variances; cov = FALSE):
#>         r1         r2         r3         r4         r5 
#> 0.00588697 0.01844908 0.05116556 0.18049160 2.40885300 
#> route: integrate; quadrature error estimate: 3.14e-12
#> checks: reflection n/a; identities sum_alpha 4.44e-16, sum_sq 7.11e-15

# With shape 1.5 the maximum has no variance, and the request is refused
try(os_moments(5, ranks = 5, dist = pareto(1.5)))
#> Error : os_moments(): the moment of order 2 of X_(5) does not exist: the upper tail index 1.5 times n - r + 1 = 1 is not above 2
os_moments(5, ranks = 1:4, dist = pareto(1.5), cov = FALSE)$var
#>         r1         r2         r3         r4 
#> 0.03227542 0.12829478 0.51263572 3.91479400 

# The Laplace (double exponential) distribution, symmetric about 0, where its
# density is not smooth
laplace <- os_dist(
  p = function(q, lower.tail = TRUE, log.p = FALSE) {
    s <- if (lower.tail) q else -q                  # F(x) = S(-x) by symmetry
    lp <- ifelse(s <= 0, log(0.5) + s, log1p(-0.5 * exp(-abs(s))))
    if (log.p) lp else exp(lp)
  },
  d = function(x, log = FALSE) if (log) log(0.5) - abs(x) else 0.5 * exp(-abs(x)),
  q = function(p) ifelse(p <= 0.5, log(2 * p), -log(2 * (1 - p))),
  symmetric = TRUE, name = "laplace",
  breaks = 0)                                       # the density has a kink at 0
os_moments(5, dist = laplace, os_shift = 10, os_multiplier = 2)
#> Exact moments of order statistics: laplace, os_shift = 10, os_multiplier = 2, n = 5, ranks 1 2 3 4 5
#> alpha (means):
#>        r1        r2        r3        r4        r5 
#>  6.822917  8.854167 10.000000 11.145830 13.177080 
#> B (covariance matrix):
#>           r1       r2       r3       r4        r5
#> r1 5.8810030 1.930286 1.056111 0.797908 0.7388585
#> r2 1.9302860 2.009844 1.120694 0.857934 0.7979080
#> r3 1.0561110 1.120694 1.404722 1.120694 1.0561110
#> r4 0.7979080 0.857934 1.120694 2.009844 1.9302860
#> r5 0.7388585 0.797908 1.056111 1.930286 5.8810030
#> route: integrate; quadrature error estimate: 4.82e-10
#> checks: reflection 4.44e-16; identities sum_alpha 2.45e-16, sum_sq 3.55e-15, sum_B 3.55e-15
```
