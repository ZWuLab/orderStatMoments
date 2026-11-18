# orderStatMoments 0.3.1

* The paper that describes the package, Wu and Tong (2026), is now cited as
  arXiv:2611.19980 (<https://doi.org/10.48550/arXiv.2611.19980>), where it was
  cited as a manuscript. The new citation is given by
  `citation("orderStatMoments")`, the help pages of the package and of
  `os_moments()`, the README and the vignette.
* The DOIs in the reference list of the vignette are links to
  <https://doi.org/>.
* The reference index of the package's website names thirteen built-in
  families, where it said five.
* No function, argument or result changes: the computations are those of
  0.3.0.

# orderStatMoments 0.3.0

* Eight families are added, so that the built-in families are the twelve
  continuous distributions for which base R provides the `p`, `d` and `q`
  functions, and the Gumbel: `"t"`, `"cauchy"`, `"gamma"`,
  `"chisq"`, `"beta"`, `"lnorm"`, `"weibull"` and `"f"` join `"norm"`, `"unif"`,
  `"exp"`, `"logis"` and `"gumbel"`. The package sets their supports and their
  tail indices, so a moment that does not exist (the variance of the minimum of
  a Cauchy sample) is refused without a declaration by the user. The noncentral
  versions (R's `ncp`) are not built in.
* A built-in family is given as in R. Its parameters are passed to
  `os_moments()` by R's own names, with R's defaults:
  `os_moments(10, dist = "gamma", shape = 2, rate = 3)`,
  `os_moments(10, dist = "norm", mean = 3, sd = 2)`. The package splits them
  into the shape parameters, which define the standard member it integrates,
  and a shift and a multiplier, which it applies exactly afterwards.
* The arguments `location` and `scale` of `os_moments()` are replaced by
  `os_shift` and `os_multiplier`, the constants `a` and `b > 0` of `a + bX`,
  which apply to any distribution. Code that passed `location` or `scale` needs
  the family's own parameter (`mean` and `sd` for the normal, `rate` for the
  exponential, `min` and `max` for the uniform) or the two new arguments. The
  logistic keeps `location` and `scale`, which are its parameters in R, and
  so does the Gumbel, which R does not provide and whose parameters the
  package names in the same way. A call with a parameter that the family does not have stops
  with a message that names the right one. Every argument after `dist` is given
  by its full name.
* The object that `os_moments()` returns stores `par`, `os_shift`,
  `os_multiplier` and `affine` in place of `location` and `scale`.
* A density that is infinite at a finite endpoint of the support (a gamma, a
  chi-square, a Weibull, a beta or an F distribution with a shape below 1) is
  integrated next to that endpoint in a variable in which the integrand is
  bounded. Such distributions are computed where the previous version stopped
  with an error (a Weibull with shape 0.5 at the minimum of 101, a gamma with
  shape 0.1 at any sample size). An extreme order statistic that concentrates
  against the endpoint is computed in the distance from it: the tests assert
  the minimum of a million observations from a beta(a, 1), at a = 0.5, 0.1
  and 0.05, against its exact moments, within 1e-10 for the mean and 1e-9 for
  the variance in relative terms. The moments of a distribution
  whose density is finite at its endpoints, and for which no `breaks` are
  declared, are computed by the same calls as before and are unchanged.
* `os_dist()` takes a new argument `breaks`, the points inside the support at
  which the density of a user-supplied distribution has a kink or a jump (the
  median of the Laplace distribution, the mode of a triangular distribution).
  `os_moments()` cuts every integral at these points, because the quadrature
  rule assumes a smooth integrand inside each piece. For the Laplace
  distribution at `n = 5`, the covariances sum to their exact value within
  1e-12 with `breaks = 0`, and within about 1e-8 without it. The density may
  be zero at a declared point. The built-in families have no such point.
* The F family computes its quantiles through `qbeta()` in each tail, because
  R's `qf()` loses the digits of a small quantile (it returns 0 at `p = 1e-4`
  for `df1 = 0.5` and `df2 = 10`).
* The help page of `os_moments()` lists the cases on which the accuracy has
  been tested. They include cases chosen because they are hard: a density that
  is infinite at an endpoint with an index down to 0.025, sample sizes up to
  1e7, and variances that barely exist.
* The help pages, the vignette and the README are rewritten for the new
  interface. They present the numerical design in three parts (what is
  computed, how it is integrated, how the integrand is evaluated), and the
  example of a user-supplied distribution is a Pareto distribution.
* `os_blue()` returns the weights and the variance factors of the two estimators
  under names that say which parameter of the location-scale family they belong
  to: `w_mu` and `w_sigma` (formerly `w_location` and `w_scale`), and `var_mu`
  and `var_sigma` (formerly `var_location` and `var_scale`). `cov`, `ranks`, `n`
  and `dist` are unchanged. The former names read like the arguments `location`
  and `scale` of `os_moments()`, which are known constants that shift and scale
  the distribution, whereas `os_blue()` estimates the unknown location
  and scale parameters of the family. Code that reads the old names needs the new ones. The
  computation is unchanged.

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
