# orderStatMoments

The means, variances and covariances of order statistics, the sorted
values of a sample, are what the standard error of a sample quantile,
the variance of a range, or the best linear estimator of a mean from a
five-number summary are built from. Printed tables of the covariances of
normal order statistics stop at n = 50 in the sources we found, and the
software we found gives the moments of single order statistics or
approximates the covariances of the normal by a series. orderStatMoments
computes the exact means, variances and covariances of order statistics,
at any ranks and any sample size, by numerical evaluation of the
defining integrals. Thirteen families are built in and are given as in
R, by name and with R’s own parameters: the normal, uniform,
exponential, logistic, Gumbel, Student t, Cauchy, gamma, chi-square,
beta, lognormal, Weibull and F distributions. Any other continuous
distribution enters through its distribution, density and quantile
functions. On top of the moments it gives those of sample quantiles and
of any linear combination of order statistics, and Lloyd’s best linear
unbiased estimators of location and scale.

The package has two requirements: the distribution must be continuous,
and a moment is computed only if it exists, which the package checks
rank by rank from the tail indices of the distribution. The integrals
are evaluated in standardized units, each moment about the center of its
own order statistic, with the integration range cut around the peak of
each order statistic’s density and the integrand on the log scale. So a
few ranks at a sample size in the thousands cost about as much as at n =
20, while a full covariance matrix grows with the square of n. Where a
family has closed forms, they are returned instead of the quadrature.

The package’s tests check the quadrature against four kinds of evidence.
The closed forms of the uniform and the exponential are reproduced at
the minimum, quartiles, median and maximum and all their pairs for n =
20, 100 and 1001 within 1e-10. The ten-decimal normal tables of Tietjen,
Kahaner and Beckman (1976) are reproduced within their stated accuracy
of 1e-7 on the entries that are printed correctly. Harter’s (1960)
moments of the normal range are reproduced within 2e-8. The identities
that hold on a complete set of ranks are satisfied within 1e-9. The
vignette states each check and its limits.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("ZWuLab/orderStatMoments", build_vignettes = TRUE)
```

## Example

``` r

library(orderStatMoments)

# the means and the covariance matrix of the minimum, the quartiles and the maximum
# of 25 standard normal observations
m <- os_moments(n = 25, ranks = c(1, 7, 13, 19, 25), dist = "norm")
m$alpha
m$B

# the same ranks for a gamma distribution, with the parameters of rgamma()
os_moments(25, ranks = c(1, 7, 13, 19, 25), dist = "gamma", shape = 2, rate = 3)

# a t distribution with 5 degrees of freedom, shifted by 10 and multiplied by 2
os_moments(25, ranks = c(1, 7, 13, 19, 25), dist = "t", df = 5,
           os_shift = 10, os_multiplier = 2)

# the standard deviation of the median of an even-sized normal sample
qw <- os_quantile_weights(p = 0.5, n = 10)
sqrt(os_lin_moments(qw$W, os_moments(10, ranks = qw$ranks))$var)

# Lloyd's BLUE of the mean and the standard deviation from a five-number summary
os_blue(m)
```

A distribution that is not built in enters through
[`os_dist()`](https://zwulab.github.io/orderStatMoments/reference/os_dist.md),
with its distribution, density and quantile functions and the tail
indices that decide which moments exist. Here is a Pareto distribution
with shape 3 on \[1, Inf):

``` r

pareto3 <- os_dist(
  p = function(q, lower.tail = TRUE, log.p = FALSE) {
    ls <- -3 * log(pmax(q, 1))                        # log S(x) = -3 log(x)
    if (lower.tail) { if (log.p) log(-expm1(ls)) else -expm1(ls) }
    else if (log.p) ls else exp(ls)
  },
  d = function(x, log = FALSE) {
    ld <- ifelse(x >= 1, log(3) - 4 * log(x), -Inf)
    if (log) ld else exp(ld)
  },
  q = function(p) (1 - p)^(-1 / 3),
  support = c(1, Inf), tail = c(Inf, 3), name = "pareto3"
)
os_moments(5, dist = pareto3)
```

The vignette walks through the families and their parameters, a
user-supplied distribution, the linear layer, what the result holds and
how the integrals are evaluated. It is at
<https://zwulab.github.io/orderStatMoments/> and, after the installation
above,
[`vignette("orderStatMoments")`](https://zwulab.github.io/orderStatMoments/articles/orderStatMoments.md).

## Citation

The paper that describes the package, its verification and the
comparison with other software is

Wu, Z. and Tong, T. (2026). orderStatMoments: Exact means, variances,
and covariances of order statistics in R. arXiv preprint
arXiv:2611.19980. <https://doi.org/10.48550/arXiv.2611.19980>

`citation("orderStatMoments")` prints this reference.

## References

Harter, H. L. (1960). Tables of range and studentized range. The Annals
of Mathematical Statistics, 31(4), 1122-1147.

Tietjen, G. L., Kahaner, D. K. and Beckman, R. J. (1976). Variances and
covariances of the normal order statistics for sample sizes 2 to 50.
Report LA-6482-MS, Los Alamos Scientific Laboratory.

## License

GPL (\>= 3).

> Developed with AI assistance; see
> [AI-USE.md](https://zwulab.github.io/orderStatMoments/AI-USE.md).
