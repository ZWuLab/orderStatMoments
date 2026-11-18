# orderStatMoments: Exact Moments of Order Statistics

Computes the exact means, variances and covariances of the order
statistics \\X\_{(1)} \le \dots \le X\_{(n)}\\ of an i.i.d. sample from
a continuous distribution, at any ranks and any sample size, by
numerical evaluation of the defining integrals, with closed forms where
the distribution has them. Thirteen families are built in and are given
as in R, by name and with R's own parameters: the normal, uniform,
exponential, logistic, Gumbel, Student \\t\\, Cauchy, gamma, chi-square,
beta, lognormal, Weibull and \\F\\ distributions. Any other continuous
distribution enters through its distribution, density and quantile
functions via
[`os_dist`](https://zwulab.github.io/orderStatMoments/reference/os_dist.md).

## Details

The main function is
[`os_moments`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md),
which returns the mean vector and the covariance matrix of the order
statistics at the requested ranks. On top of it,
[`os_quantile_weights`](https://zwulab.github.io/orderStatMoments/reference/os_quantile_weights.md),
[`os_lin_moments`](https://zwulab.github.io/orderStatMoments/reference/os_lin_moments.md)
and
[`os_blue`](https://zwulab.github.io/orderStatMoments/reference/os_blue.md)
give the moments of sample quantiles and of any linear combination of
order statistics, and Lloyd's best linear unbiased estimators of
location and scale. The vignette walks through all of it:
[`vignette("orderStatMoments")`](https://zwulab.github.io/orderStatMoments/articles/orderStatMoments.md).
The paper of Wu and Tong (2026) describes the computation, its
verification against closed forms, identities and printed tables, and
the comparison with other software; `citation("orderStatMoments")` gives
the reference.

## References

Wu, Z. and Tong, T. (2026). orderStatMoments: Exact means, variances,
and covariances of order statistics in R. arXiv preprint
arXiv:2611.19980.
[doi:10.48550/arXiv.2611.19980](https://doi.org/10.48550/arXiv.2611.19980)
.

## See also

Useful links:

- <https://github.com/ZWuLab/orderStatMoments>

- <https://zwulab.github.io/orderStatMoments/>

- Report bugs at <https://github.com/ZWuLab/orderStatMoments/issues>

## Author

**Maintainer**: Zheyang Wu <zheyangwu@wpi.edu>
