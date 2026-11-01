# orderStatMoments 0.2.0

* The package is renamed `orderStatMoments`, from `osmoments`, so that its name
  says what it computes. The functions keep their names (`os_moments()`,
  `os_dist()`, `os_quantile_weights()`, `os_lin_moments()`, `os_blue()`). The
  object that `os_moments()` returns now has class `"orderStatMoments"`
  (formerly `"osmoments"`), and its print method is renamed with it. Install
  with `remotes::install_github("ZWuLab/orderStatMoments")`. The former
  repository address redirects there.
* The vignette's Appendix A and the help page of `os_moments()` say "approach"
  where they said "device" for the four numerical approaches of the
  computation. The code is unchanged.

# osmoments 0.1.0

* First release.
* `os_moments()` computes the exact means, variances and covariances of order
  statistics at any ranks and any sample size, for the normal, uniform,
  exponential, logistic and Gumbel distributions and for a user-supplied
  distribution built with `os_dist()`. Closed forms are returned where the
  family has them (the uniform, the exponential, the logistic single moments,
  the Gumbel maximum) and the exact integrals are evaluated by adaptive
  quadrature elsewhere.
* `os_quantile_weights()` gives the ranks and weights of the sample quantiles
  under the nine definitions of `quantile()`, `os_lin_moments()` the moments
  of any linear combination of order statistics, and `os_blue()` Lloyd's best
  linear unbiased estimators of location and scale.
