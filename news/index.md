# Changelog

## orderStatMoments 0.2.0

- The package is renamed `orderStatMoments`, from `osmoments`, so that
  its name says what it computes. The functions keep their names
  ([`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md),
  [`os_dist()`](https://zwulab.github.io/orderStatMoments/reference/os_dist.md),
  [`os_quantile_weights()`](https://zwulab.github.io/orderStatMoments/reference/os_quantile_weights.md),
  [`os_lin_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_lin_moments.md),
  [`os_blue()`](https://zwulab.github.io/orderStatMoments/reference/os_blue.md)).
  The object that
  [`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
  returns now has class `"orderStatMoments"` (formerly `"osmoments"`),
  and its print method is renamed with it. Install with
  `remotes::install_github("ZWuLab/orderStatMoments")`. The former
  repository address redirects there.
- The vignette’s Appendix A and the help page of
  [`os_moments()`](https://zwulab.github.io/orderStatMoments/reference/os_moments.md)
  say “approach” where they said “device” for the four numerical
  approaches of the computation. The code is unchanged.
