# ==============================================================================
# File: test-os_families.R
# Purpose: Tests of the built-in families given by R's own parameters: the
#          split of the parameters into the standard member and a shift and a
#          multiplier, the closed forms that exist for the new families (the
#          Weibull at any shape, the t with 2 degrees of freedom), the special
#          cases that coincide with a family that has closed forms, the
#          permutation identities against the exact moments of the parent, the
#          moment condition of the t, the Cauchy and the F, the densities that
#          are infinite at an endpoint (the gamma, the chi-square, the Weibull,
#          the F and the beta with a shape below 1), and the error messages.

# E[X_(r)^k] of a Weibull with shape c and scale 1: expanding F^(r-1) with
# F(x) = 1 - exp(-x^c) gives a finite alternating sum, exact for a small rank
weibull_moment <- function(r, n, k, c) {
  j <- 0:(r - 1)
  r * choose(n, r) * gamma(1 + k / c) * sum((-1)^j * choose(r - 1, j) * (n - r + 1 + j)^(-(1 + k / c)))
}
# The exact mean and variance of X_(r) for X = U^K, a beta(1/K, 1):
# E[X_(r)^k] = prod_{j < Kk} (r + j)/(n + 1 + j). The variance is mu^2 {E[X^2]/mu^2 - 1}, with
# log(E[X^2]/mu^2) = sum_{j < K} log1p(K/(r + j)) - log1p(K/(n + 1 + j)), a sum of small terms
# of one sign each. (Forming that logarithm as a difference of two sums of logs loses its
# digits when the variance is small against mu^2: at n = 1e7 the reference was then 2e-8 off.)
power_exact <- function(n, r, K) {
  j <- 0:(K - 1)
  mu <- exp(sum(log(r + j) - log(n + 1 + j)))
  c(mu, mu^2 * expm1(sum(log1p(K / (r + j)) - log1p(K / (n + 1 + j)))))
}
# The t with 2 degrees of freedom has Q(u) = (2u - 1)/sqrt(2u(1 - u)), so its single
# moments are beta functions; the second moment exists for 2 <= r <= n - 1
t2_mean <- function(r, n)
  (2 * beta(r + 0.5, n - r + 0.5) - beta(r - 0.5, n - r + 0.5)) / (sqrt(2) * beta(r, n - r + 1))
t2_m2 <- function(r, n)
  (4 * beta(r + 1, n - r) - 4 * beta(r, n - r) + beta(r - 1, n - r)) / (2 * beta(r, n - r + 1))

test_that("no parameter of a family has the name of an argument of os_moments()", {
  spec <- orderStatMoments:::.os_family_spec
  pars <- unique(unlist(lapply(spec, function(s) names(s$par))))
  args <- setdiff(names(formals(os_moments)), "...")
  expect_length(intersect(pars, args), 0L)
  # nor is one a leading part of an argument that stands before the dots, which R would match
  before <- args[seq_len(which(names(formals(os_moments)) == "...") - 1L)]
  expect_false(any(vapply(pars, function(p) any(startsWith(before, p)), logical(1))))
  expect_identical(names(spec), c("norm", "unif", "exp", "logis", "gumbel", "t", "cauchy",
                                  "gamma", "chisq", "beta", "lnorm", "weibull", "f"))
})

test_that("R's parameters are split into the standard member and a shift and a multiplier", {
  rk <- c(1, 3, 6)
  same <- function(m, z, a, b) {
    expect_equal(m$alpha, a + b * z$alpha, tolerance = 1e-14)
    expect_equal(m$B, b^2 * z$B, tolerance = 1e-14)
    expect_identical(m$standard, z$standard)
    expect_equal(unname(m$affine), c(a, b))
  }
  same(os_moments(6, rk, "norm", mean = 3, sd = 2), os_moments(6, rk, "norm"), 3, 2)
  same(os_moments(6, rk, "unif", min = -1, max = 4), os_moments(6, rk, "unif"), -1, 5)
  same(os_moments(6, rk, "exp", rate = 4), os_moments(6, rk, "exp"), 0, 1 / 4)
  same(os_moments(6, rk, "logis", location = 1, scale = 3), os_moments(6, rk, "logis"), 1, 3)
  same(os_moments(6, rk, "gumbel", location = 1, scale = 3), os_moments(6, rk, "gumbel"), 1, 3)
  same(os_moments(6, 3:4, "cauchy", location = 1, scale = 3), os_moments(6, 3:4, "cauchy"), 1, 3)
  same(os_moments(6, rk, "gamma", shape = 2, rate = 3), os_moments(6, rk, "gamma", shape = 2), 0, 1 / 3)
  same(os_moments(6, rk, "gamma", shape = 2, scale = 3), os_moments(6, rk, "gamma", shape = 2), 0, 3)
  same(os_moments(6, rk, "lnorm", meanlog = 1, sdlog = 0.5),
       os_moments(6, rk, "lnorm", sdlog = 0.5), 0, exp(1))
  same(os_moments(6, rk, "weibull", shape = 1.5, scale = 2),
       os_moments(6, rk, "weibull", shape = 1.5), 0, 2)
  # os_shift and os_multiplier act on top of the family's own parameters
  same(os_moments(6, rk, "gamma", shape = 2, rate = 3, os_shift = 5, os_multiplier = 2),
       os_moments(6, rk, "gamma", shape = 2), 5, 2 / 3)
  same(os_moments(6, rk, "t", df = 5, os_shift = 3, os_multiplier = 2),
       os_moments(6, rk, "t", df = 5), 3, 2)
  # a large mean is applied after the integration and costs no digits
  big <- os_moments(6, rk, "norm", mean = 1e15, sd = 2)
  expect_equal(diag(big$B), 4 * diag(os_moments(6, rk)$B), tolerance = 1e-14)
  expect_equal(big$par, c(mean = 1e15, sd = 2))
  expect_null(os_moments(3, dist = os_dist(pnorm, dnorm, qnorm, symmetric = TRUE))$par)
})

test_that("the built-in call agrees with R's own random variable through the object's fields", {
  m <- os_moments(4, dist = "gamma", shape = 2, rate = 3, os_shift = 1, os_multiplier = 2)
  expect_identical(m$dist$name, "gamma")
  expect_equal(m$par, c(shape = 2, rate = 3))
  expect_identical(c(m$os_shift, m$os_multiplier), c(1, 2))
  expect_output(print(m), "gamma \\(shape = 2, rate = 3\\), os_shift = 1, os_multiplier = 2, n = 4")
  # n = 1 is the distribution itself: mean 1 + 2 * 2/3, variance 4 * 2/9
  one <- os_moments(1, dist = "gamma", shape = 2, rate = 3, os_shift = 1, os_multiplier = 2)
  expect_equal(unname(one$alpha), 1 + 4 / 3, tolerance = 1e-12)
  expect_equal(unname(one$var), 8 / 9, tolerance = 1e-12)
})

test_that("special cases agree with the families that have closed forms", {
  rk <- c(1, 4, 9)
  e <- os_moments(9, rk, "exp", method = "closed")
  for (m in list(os_moments(9, rk, "gamma", shape = 1), os_moments(9, rk, "weibull", shape = 1))) {
    expect_equal(m$alpha, e$alpha, tolerance = 1e-11)
    expect_equal(m$B, e$B, tolerance = 1e-10)
  }
  u <- os_moments(9, rk, "unif", method = "closed")
  b <- os_moments(9, rk, "beta", shape1 = 1, shape2 = 1)
  expect_equal(b$alpha, u$alpha, tolerance = 1e-11)
  expect_equal(b$B, u$B, tolerance = 1e-10)
  # the chi-square with df degrees of freedom is the gamma with shape df/2 and scale 2
  ch <- os_moments(9, rk, "chisq", df = 3); g <- os_moments(9, rk, "gamma", shape = 1.5, scale = 2)
  expect_equal(ch$alpha, g$alpha, tolerance = 1e-11)
  expect_equal(ch$B, g$B, tolerance = 1e-10)
  # the Cauchy is the t with one degree of freedom
  expect_equal(os_moments(7, 3:5, "cauchy")$B, os_moments(7, 3:5, "t", df = 1)$B, tolerance = 1e-10)
})

test_that("the Weibull agrees with its closed form at shapes above and below 1", {
  for (c in c(0.3, 0.5, 2)) for (n in c(5, 101)) for (r in 1:2) {
    m <- os_moments(n, r, "weibull", shape = c, cov = FALSE)
    mu <- weibull_moment(r, n, 1, c)
    expect_lt(abs(m$alpha[[1]] / mu - 1), 1e-10)
    expect_lt(abs(m$var[[1]] / (weibull_moment(r, n, 2, c) - mu^2) - 1), 1e-9)
  }
  # the minimum of 1001 at shape 0.3 sits 1e-10 spreads from the endpoint, where the
  # centered coordinate cannot place it; the endpoint coordinate keeps its relative precision
  m <- os_moments(1001, 1, "weibull", shape = 0.3, cov = FALSE)
  expect_lt(abs(m$alpha[[1]] / weibull_moment(1, 1001, 1, 0.3) - 1), 1e-10)
  expect_lt(abs(m$var[[1]] / (weibull_moment(1, 1001, 2, 0.3) - weibull_moment(1, 1001, 1, 0.3)^2) - 1), 1e-9)
})

test_that("the t with 2 degrees of freedom agrees with its closed form", {
  for (n in c(5, 25, 101)) {
    rk <- sort(unique(c(2, n %/% 2 + 1, n - 1)))
    m <- os_moments(n, rk, "t", df = 2, cov = FALSE)
    expect_lt(max(abs(m$alpha - t2_mean(rk, n))), 1e-11)
    expect_lt(max(abs(m$var / (t2_m2(rk, n) - t2_mean(rk, n)^2) - 1)), 1e-9)
  }
})

test_that("the moment condition follows the tail indices of the t, the Cauchy and the F", {
  expect_error(os_moments(5, 1, "t", df = 2), "does not exist: the lower tail index 2 times the rank 1")
  expect_silent(os_moments(5, 2, "t", df = 2, cov = FALSE))
  expect_silent(os_moments(5, 1, "t", df = 2.5, cov = FALSE))
  expect_error(os_moments(5, 2, "cauchy"), "lower tail index 1 times the rank 2")
  expect_error(os_moments(5, 4, "cauchy"), "upper tail index 1")
  expect_silent(os_moments(5, 3, "cauchy"))
  # the F: the upper tail index is df2/2, so the maximum needs df2 > 4 and the next rank df2 > 2
  expect_error(os_moments(5, 5, "f", df1 = 3, df2 = 4), "upper tail index 2 times n - r \\+ 1 = 1")
  expect_silent(os_moments(5, 4, "f", df1 = 3, df2 = 4, cov = FALSE))
  expect_silent(os_moments(5, 5, "f", df1 = 3, df2 = 5, cov = FALSE))
  # the variance of the median of five Cauchy observations, against the integral in u
  v <- integrate(function(u) qcauchy(u)^2 * dbeta(u, 3, 3), 0, 1, rel.tol = 1e-12)$value
  expect_equal(unname(os_moments(5, 3, "cauchy")$var), v, tolerance = 1e-9)
})

test_that("the permutation identities hold for the new families against the exact parent moments", {
  n <- 8
  chk <- function(m, EX, VX, tol = 1e-9) {
    expect_lt(abs(sum(m$alpha) - n * EX) / (n * max(1, abs(EX))), tol)
    expect_lt(abs(sum(m$B) / (n * VX) - 1), tol)
    expect_gt(min(eigen(m$B, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
  chk(os_moments(n, dist = "t", df = 5), 0, 5 / 3)
  chk(os_moments(n, dist = "gamma", shape = 2.5), 2.5, 2.5)
  chk(os_moments(n, dist = "chisq", df = 4), 4, 8)
  chk(os_moments(n, dist = "beta", shape1 = 2, shape2 = 3), 2 / 5, 1 / 25)
  chk(os_moments(n, dist = "lnorm", sdlog = 1), exp(0.5), (exp(1) - 1) * exp(1))
  chk(os_moments(n, dist = "weibull", shape = 1.5), gamma(1 + 1 / 1.5),
      gamma(1 + 2 / 1.5) - gamma(1 + 1 / 1.5)^2)
  chk(os_moments(n, dist = "f", df1 = 4, df2 = 10), 10 / 8, 2 * 100 * 12 / (4 * 64 * 6))
})

test_that("a density that is infinite at an endpoint is integrated in the endpoint's coordinate", {
  n <- 8
  chk <- function(m, EX, VX, tol = 1e-10) {
    expect_lt(abs(sum(m$alpha) / (n * EX) - 1), tol)
    expect_lt(abs(sum(m$B) / (n * VX) - 1), tol)
  }
  chk(os_moments(n, dist = "gamma", shape = 0.5), 0.5, 0.5)
  chk(os_moments(n, dist = "gamma", shape = 0.1), 0.1, 0.1)
  chk(os_moments(n, dist = "chisq", df = 1), 1, 2)
  chk(os_moments(n, dist = "weibull", shape = 0.5), 2, 20)
  chk(os_moments(n, dist = "f", df1 = 1, df2 = 10), 10 / 8, 2 * 100 * 9 / (64 * 6))
  bm <- function(a, b) c(a / (a + b), a * b / ((a + b)^2 * (a + b + 1)))
  for (ab in list(c(0.5, 0.5), c(0.2, 0.2), c(0.3, 3), c(3, 0.3), c(0.2, 0.45))) {
    m <- os_moments(n, dist = "beta", shape1 = ab[1], shape2 = ab[2])
    chk(m, bm(ab[1], ab[2])[1], bm(ab[1], ab[2])[2])
  }
  # the beta with its shapes exchanged is the reflected variable 1 - X: the upper endpoint
  # of one is the lower endpoint of the other, and the two computations agree
  m1 <- os_moments(n, dist = "beta", shape1 = 3, shape2 = 0.3)
  m2 <- os_moments(n, dist = "beta", shape1 = 0.3, shape2 = 3)
  expect_equal(unname(m1$alpha), 1 - rev(unname(m2$alpha)), tolerance = 1e-12)
  expect_equal(unname(m1$B), unname(m2$B)[n:1, n:1], tolerance = 1e-9)
  # the arcsine distribution at n = 2: Cov(X_(1), X_(2)) = 4/pi^4 from X_(1) X_(2) = X_1 X_2
  a2 <- os_moments(2, dist = "beta", shape1 = 0.5, shape2 = 0.5)
  expect_lt(abs(a2$B[1, 2] - 4 / pi^4), 1e-12)
  # extremes that concentrate against the endpoint: the minimum and the maximum of 1001
  ex <- os_moments(1001, c(1, 1001), dist = "beta", shape1 = 0.5, shape2 = 0.5, cov = FALSE)
  expect_equal(1 - ex$alpha[[1]], ex$alpha[[2]], tolerance = 1e-14)   # the two ends mirror each other
  expect_equal(ex$var[[1]], ex$var[[2]], tolerance = 1e-9)
  ref <- integrate(function(v) qbeta(v^2, 0.5, 0.5) * 1001 * (1 - v^2)^1000 * 2 * v, 0, 1,
                   rel.tol = 1e-13)$value        # E[X_(1)] in u = v^2, where the integrand is smooth
  expect_lt(abs(ex$alpha[[1]] / ref - 1), 1e-9)
})

test_that("an order statistic closer to the endpoint than the centered coordinate resolves", {
  # A beta(0.1, 1) variable is U^10 with U uniform, so the minimum of n has
  # E[X_(1)^k] = prod_{j=1}^{10k} j/(n + j). At n = 101 its center lies 2e-19 spreads above 0
  # and equals the endpoint in the centered coordinate. The first engine with endpoint
  # coordinates placed the cut at the center by that rounded value, so the first piece
  # started at the center instead of at the endpoint, and it lost the scaling of the
  # tolerances: the mean came back 2.5e-7 off at n = 101, 5 percent off at n = 1001 and 44
  # percent off at n = 3001, without a word. The cuts now carry their distances from the
  # endpoint.
  for (n in c(101, 1001, 3001)) {
    mu <- prod((1:10) / (n + 1:10)); v <- prod((1:20) / (n + 1:20)) - mu^2
    z <- os_moments(n, 1, "beta", shape1 = 0.1, shape2 = 1, cov = FALSE)
    expect_lt(abs(z$alpha[[1]] / mu - 1), 1e-10)
    expect_lt(abs(z$var[[1]] / v - 1), 1e-9)
    # the maximum of the reflected beta(1, 0.1) is one minus that minimum
    zr <- os_moments(n, n, "beta", shape1 = 1, shape2 = 0.1, cov = FALSE)
    expect_lt(abs(zr$var[[1]] / v - 1), 1e-9)
  }
  # At a million observations the peak of the minimum is 1e-5 of the first piece in the
  # endpoint's variable, and for a small shape the integrand of a moment peaks far above
  # the density's own peak. Before the cuts of .os_end_cuts() the minimum of 1e6 from a
  # beta(1/2, 1) stopped at the mass check, and with two cuts only the variance for a
  # beta(0.05, 1) came back 48 percent off without a word.
  for (a in c(0.5, 0.1, 0.05)) {
    n <- 1e6; cc <- round(1 / a)
    mu <- exp(sum(log((1:cc) / (n + 1:cc))))
    v <- exp(sum(log((1:(2 * cc)) / (n + 1:(2 * cc))))) - mu^2
    z <- os_moments(n, 1, "beta", shape1 = a, shape2 = 1, cov = FALSE)
    expect_lt(abs(z$alpha[[1]] / mu - 1), 1e-10)
    expect_lt(abs(z$var[[1]] / v - 1), 1e-9)
    expect_lt(abs(os_moments(n, n, "beta", shape1 = 1, shape2 = a, cov = FALSE)$var[[1]] / v - 1), 1e-9)
  }
  # the density of X_(1) has mass 1 whether or not the probability of the center is passed,
  # and no piece of its range has its limits reversed
  ns <- asNamespace("orderStatMoments")
  for (sh in list(c(0.1, 1), c(1, 0.1))) {
    n <- 101L; r <- if (sh[1] < 1) 1L else n
    cd <- ns$.os_centered(ns$.os_family("beta", list(shape1 = sh[1], shape2 = sh[2]))$dist)
    ctrl <- ns$.os_ctrl(); u <- r / (n + 1)
    expect_lt(abs(ns$.os_single(r, n, 0, cd, ctrl, cd$q(u), u)$value - 1), 1e-9)
    expect_lt(abs(ns$.os_single(r, n, 0, cd, ctrl, cd$q(u))$value - 1), 1e-9)
    ds <- ns$.os_shift_dist(cd, cd$q(u), u)
    ends <- ns$.os_ends(cd, function(ev, dev, side) function(e) 1 + 0 * e, ds[["lo"]], ds[["hi"]])
    cuts <- rbind(ns$.os_cuts(cd, ns$.os_breaks(r, n, cd, ctrl$k)), ns$.os_cuts(cd, cd$q(u), ds),
                  ns$.os_end_cuts(r, n, cd, ctrl$k))
    pcs <- ns$.os_pieces(cd$lower, cd$upper, cuts, ends, NULL)
    expect_true(all(vapply(pcs, function(p) p$a <= p$b, logical(1))))
    # the pieces cover the support once: their lengths in t add up to its length
    len <- sum(vapply(pcs, function(p) if (p$side == "t") p$b - p$a else p$b - p$a, numeric(1)))
    expect_equal(len, cd$upper - cd$lower, tolerance = 1e-12)
    # an inner range from a point next to the upper endpoint is not reversed either
    if (!is.null(cd$hi)) {
      inner <- ns$.os_pieces(cd$upper - 1e-19, cd$upper, cuts, ends, list(side = "hi", e = 1e-19))
      expect_true(all(vapply(inner, function(p) p$side == "hi" && p$a <= p$b, logical(1))))
      expect_equal(sum(vapply(inner, function(p) p$b - p$a, numeric(1))), 1e-19)
    }
  }
})

test_that("the covariances next to an endpoint agree with their exact values for a beta(a, 1)", {
  # X = U^c with c = 1/a. For r < s, U_(r) = U_(s) W with W ~ beta(r, s - r) independent of
  # U_(s) ~ beta(s, n - s + 1), so X_(r) = X_(s) W^c and
  # Cov(X_(r), X_(s)) = E[W^c] Var(U_(s)^c), with E[B^k] = B(al + k, be)/B(al, be).
  bm <- function(al, be, k) exp(lbeta(al + k, be) - lbeta(al, be))
  for (case in list(list(n = 7L, rk = 1:7, tol = 1e-9), list(n = 1001L, rk = c(1L, 2L, 1000L, 1001L), tol = 1e-7))) {
    n <- case$n; rk <- case$rk; cc <- 10
    z <- os_moments(n, rk, "beta", shape1 = 0.1, shape2 = 1)
    zr <- os_moments(n, rev(n + 1L - rk), "beta", shape1 = 1, shape2 = 0.1)   # the reflection
    k <- length(rk)
    for (i in seq_len(k - 1L)) for (j in (i + 1L):k) {
      r <- rk[i]; s <- rk[j]
      ex <- bm(r, s - r, cc) * (bm(s, n - s + 1, 2 * cc) - bm(s, n - s + 1, cc)^2)
      expect_lt(abs(z$B[i, j] / ex - 1), case$tol)
      expect_lt(abs(zr$B[k + 1L - j, k + 1L - i] / ex - 1), case$tol)
    }
  }
})

test_that("a user's distribution with an infinite density at 0 takes the same route", {
  g <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                 pgamma(q, shape = 0.5, lower.tail = lower.tail, log.p = log.p),
               function(x, log = FALSE) dgamma(x, shape = 0.5, log = log),
               function(p) qgamma(p, shape = 0.5), support = c(0, Inf))
  # the built-in gamma uses its own index 1/2 and the user's object an estimate of it, so the
  # two agree to the accuracy of the quadrature, not to the last bit
  expect_equal(os_moments(6, dist = g)$B, os_moments(6, dist = "gamma", shape = 0.5)$B,
               tolerance = 1e-9)
  cd <- orderStatMoments:::.os_centered(g)
  expect_equal(cd$lo$pw, 2, tolerance = 1e-5)   # the density behaves like e^(1/2 - 1)
  expect_identical(orderStatMoments:::.os_centered(
    orderStatMoments:::.os_family("gamma", list(shape = 0.5))$dist)$lo$pw, 2)
  expect_null(cd$hi)
  expect_null(orderStatMoments:::.os_centered(orderStatMoments:::.os_family("exp")$dist)$lo)
})

test_that("the F family takes its quantiles through the beta quantile of each tail", {
  # R's qf() forms 1/q - 1 from an upper-tail beta quantile and loses the digits of a small
  # quantile: it returns 0 at p = 1e-4 for df1 = 0.5 and df2 = 10
  q <- orderStatMoments:::.qf_beta(0.5, 10)
  p <- c(1e-6, 1e-4, 1e-3, 0.3, 0.5, 0.7, 0.999)
  expect_lt(max(abs(pf(q(p), 0.5, 10) / p - 1)), 1e-12)
  expect_gt(q(1e-4), 0)
  p2 <- c(0.01, 0.3, 0.5, 0.7, 0.99)
  expect_equal(orderStatMoments:::.qf_beta(4, 10)(p2), qf(p2, 4, 10), tolerance = 1e-12)
  # the minimum of 100000 at df1 = 0.5 is computed, with a positive variance
  expect_gt(os_moments(100000, 1, "f", df1 = 0.5, df2 = 10, cov = FALSE)$var[[1]], 0)
  # a quantile function that returns the endpoint itself still gives a usable center
  g <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                 pgamma(q, shape = 0.5, lower.tail = lower.tail, log.p = log.p),
               function(x, log = FALSE) dgamma(x, shape = 0.5, log = log),
               function(p) ifelse(p < 0.01, 0, qgamma(p, shape = 0.5)), support = c(0, Inf))
  expect_equal(os_moments(1001, 1, g, cov = FALSE)$var,
               os_moments(1001, 1, "gamma", shape = 0.5, cov = FALSE)$var, tolerance = 1e-8)
})

test_that("the parameters are validated with messages that say what to give", {
  expect_error(os_moments(5, dist = "norm", scale = 2),
               "'scale' is not a parameter of the \"norm\" family, whose parameters are mean and sd.*os_multiplier")
  expect_error(os_moments(5, dist = "exp", location = 1), "os_shift")
  expect_error(os_moments(5, dist = "gamma"), "needs 'shape'")
  expect_error(os_moments(5, dist = "t"), "needs 'df'")
  expect_error(os_moments(5, dist = "beta", shape1 = 2), "needs 'shape2'")
  expect_error(os_moments(5, dist = "gamma", shape = 2, rate = 1, scale = 1), "not both")
  expect_error(os_moments(5, dist = "t", df = 5, ncp = 1), "noncentral")
  expect_error(os_moments(5, dist = "chisq", df = 3, ncp = 1), "noncentral")
  expect_error(os_moments(5, dist = "exp", covv = FALSE), "full name")
  expect_error(os_moments(5, dist = "gamma", shape = -1), "'shape'")
  expect_error(os_moments(5, dist = "gamma", shape = c(1, 2)), "one number")
  expect_error(os_moments(5, dist = "unif", min = 2, max = 1), "'max' must be above 'min'")
  expect_error(os_moments(5, dist = "lnorm", sdlog = 0), "'sdlog'")
  expect_error(os_moments(5, dist = "weibull", shape = 2, scale = Inf), "'scale'")
  expect_error(os_moments(5, dist = os_dist(pnorm, dnorm, qnorm), mean = 2), "os_dist\\(\\)")
  expect_error(os_moments(5, dist = "norm", os_multiplier = 0), "'os_multiplier'")
})

test_that("a piece that the standardized coordinate does not resolve is integrated in the distance", {
  # X = U^K is a beta(1/K, 1), so E[X_(r)^k] = prod_{j < Kk} (r + j)/(n + 1 + j) exactly. At
  # these ranks a piece of the range lies 1e-18 spreads from the endpoint 0 and is about as
  # long as it is far, so that the rule for a singular endpoint (the endpoint no farther than
  # the piece is long) did not apply, and in the standardized coordinate its two limits are
  # one rounding step apart, or equal. The mean came back 2.9 percent low and the variance
  # 14.5 percent low without an error (at n = 1e7: 0.12 and 0.05 percent off).
  for (a in list(c(1e6, 39200, 14), c(1e6, 45000, 15), c(1e6, 57800, 17), c(1e6, 72200, 19),
                 c(1e7, 78400, 14))) {
    o <- os_moments(a[1], ranks = a[2], dist = "beta", shape1 = 1 / a[3], shape2 = 1, cov = FALSE)
    e <- power_exact(a[1], a[2], a[3])
    expect_lt(abs(o$alpha[[1]] / e[1] - 1), 1e-10)
    expect_lt(abs(o$var[[1]] / e[2] - 1), 1e-9)
  }
})

test_that("the rounding of the argument x decides the coordinate of a piece, not that of t alone", {
  # 1 - U^40 is a beta(1, 1/40), and its order statistic of rank s is one minus the order
  # statistic of rank n + 1 - s of U^40. Its median is 1 - 9e-13 and its spread 7.5e-6, so
  # the argument x = median + spread * t is rounded 1.3e5 times more coarsely than t. A
  # piece that t resolves and x does not stayed in t, and the variance of rank 587180 of
  # n = 1e6 came back 39 percent low (rank 580001: 26 percent low) without an error.
  n <- 1e6; K <- 40
  for (s in c(587180, 580001)) {
    o <- os_moments(n, ranks = s, dist = "beta", shape1 = 1, shape2 = 1 / K, cov = FALSE)
    e <- power_exact(n, n + 1 - s, K)
    expect_lt(abs(o$alpha[[1]] - (1 - e[1])), 1e-14)
    expect_lt(abs(o$var[[1]] / e[2] - 1), 1e-9)
    # the reflection, through the lower endpoint's view
    lo <- os_moments(n, ranks = n + 1 - s, dist = "beta", shape1 = 1 / K, shape2 = 1, cov = FALSE)
    expect_lt(abs(lo$alpha[[1]] / e[1] - 1), 1e-10)
    expect_lt(abs(lo$var[[1]] / e[2] - 1), 1e-9)
  }
})

test_that("an order statistic is taken from the endpoint it is nearer to, not from the side of its rank", {
  # A beta(1/K, b) with b = 0.8 has an infinite density at both endpoints. Next to 0 its
  # distribution function is x^a/(a B(a, b)) {1 + O(x)} with a = 1/K, so X_(r) = (c U_(r))^K
  # up to a relative error of about 0.2 x, with c = a B(a, b): the mean and the variance are
  # c^K and c^(2K) times those of U_(r)^K (power_exact()), and the covariance is
  # c^(2K) E[W^K] Var(U_(s)^K) as in the test of the covariances above. The median of a
  # beta(1/66, 0.8) is 2e-20, so the ranks just above n/2 lie next to the lower endpoint.
  # They were taken from the upper endpoint, because their probability is above one half:
  # the means came back as 0 and, at n = 1e6, the variances and the covariance 21 to 28
  # percent too high, without an error and with a mass of 1.
  b <- 0.8
  near0 <- function(n, r, K)
    power_exact(n, r, K) * exp(c(K, 2 * K) * (lgamma(1 + 1 / K) + lgamma(b) - lgamma(1 / K + b)))
  bm <- function(al, be, k) exp(lbeta(al + k, be) - lbeta(al, be))
  n <- 1e6; K <- 66; rk <- c(500100, 510000)
  z <- os_moments(n, rk, "beta", shape1 = 1 / K, shape2 = b)
  e1 <- near0(n, rk[1], K); e2 <- near0(n, rk[2], K)
  expect_lt(max(abs(z$alpha / c(e1[1], e2[1]) - 1)), 1e-10)
  expect_lt(max(abs(z$var / c(e1[2], e2[2]) - 1)), 1e-9)
  expect_lt(abs(z$B[1, 2] / (bm(rk[1], rk[2] - rk[1], K) * e2[2]) - 1), 1e-8)
  # a pair farther apart, both above the median rank (the correlation is 0.83)
  rk <- c(510000, 600000)
  z <- os_moments(n, rk, "beta", shape1 = 1 / K, shape2 = b)
  e2 <- near0(n, rk[2], K)
  expect_lt(abs(z$B[1, 2] / (bm(rk[1], rk[2] - rk[1], K) * e2[2]) - 1), 1e-8)
  # the same at smaller sample sizes and other indices: the mean of each came back as 0, or
  # with four or five correct digits (K = 40, where it is 3e-12 and was assembled as 1 - (1 - x))
  for (a in list(c(1e6, 510000, 40), c(1e5, 51000, 50), c(1e4, 5100, 66), c(1e4, 6000, 100))) {
    o <- os_moments(a[1], ranks = a[2], dist = "beta", shape1 = 1 / a[3], shape2 = b, cov = FALSE)
    e <- near0(a[1], a[2], a[3])
    expect_lt(abs(o$alpha[[1]] / e[1] - 1), 1e-10)
    expect_lt(abs(o$var[[1]] / e[2] - 1), 1e-9)
  }
})

test_that("an order statistic above the median rank and next to the lower endpoint is computed", {
  # A beta(1/K, 1) has one view, at 0, and its rank 0.6 n lies 2e-7 above 0 for K = 30. The
  # cuts around the peak were made for the ranks below the median only, and the first of
  # these calls stopped at the mass check ("the density of X_(600000) integrates to
  # 5.4e-25"). The other two were computed before and guard the same path. The reflections
  # go through the upper endpoint's view.
  for (a in list(c(1e6, 600000, 30), c(1e6, 700000, 40), c(1e5, 80000, 40))) {
    e <- power_exact(a[1], a[2], a[3])
    lo <- os_moments(a[1], ranks = a[2], dist = "beta", shape1 = 1 / a[3], shape2 = 1, cov = FALSE)
    expect_lt(abs(lo$alpha[[1]] / e[1] - 1), 1e-10)
    expect_lt(abs(lo$var[[1]] / e[2] - 1), 1e-9)
    up <- os_moments(a[1], ranks = a[1] + 1 - a[2], dist = "beta", shape1 = 1, shape2 = 1 / a[3], cov = FALSE)
    expect_lt(abs(up$alpha[[1]] - (1 - e[1])), 1e-14)
    expect_lt(abs(up$var[[1]] / e[2] - 1), 1e-9)
  }
})

test_that("a ratio of tail probabilities that rounding puts above 1 raises no warning", {
  # The ranks 999991 and 999999 of a million from a beta(1, 1/3) are the reflections of the
  # ranks 2 and 10 of U^3. The value was right, and R warned "NaNs produced" from inside the
  # difference of two tail probabilities.
  n <- 1e6
  bm <- function(al, be, k) exp(lbeta(al + k, be) - lbeta(al, be))
  expect_no_warning(z <- os_moments(n, c(999991, 999999), "beta", shape1 = 1, shape2 = 1 / 3))
  ex <- bm(2, 8, 3) * (bm(10, n - 9, 6) - bm(10, n - 9, 3)^2)
  expect_lt(abs(z$B[1, 2] / ex - 1), 1e-9)
})

test_that("a covariance across the two endpoints of a beta(1/2, 1/2) agrees with its series", {
  # The beta(1/2, 1/2) (arcsine) distribution has Q(u) = sin^2(pi u / 2) = sum_k c_k u^(2k),
  # c_k = (-1)^(k + 1) pi^(2k) / {2 (2k)!}. For r < s, X_(r) = Q(U) and X_(s) = 1 - Q(V) with
  # U = U_(r) and V = 1 - U_(s), two parts of a Dirichlet(r, s - r, t) vector, t = n + 1 - s,
  # whose moments are ratios of rising factorials a^(k) = a (a + 1) ... (a + k - 1):
  #   Cov(X_(r), X_(s)) = sum_{i,j} c_i c_j {E[U^m] E[V^l] - E[U^m V^l]},  m = 2i, l = 2j,
  #   E[U^m] E[V^l] - E[U^m V^l] = r^(m) t^(l) / {(n+1)^(m) (n+1)^(l)}
  #                                 x [1 - prod_{k < l} (n + 1 + k)/(n + 1 + m + k)],
  # the bracket formed by expm1 of a sum of log1p terms, so that nothing cancels. Both order
  # statistics lie next to an endpoint where the density is infinite, one at each end, and
  # they are nearly independent: the covariance is accurate in units of the correlation and
  # keeps fewer digits in relative terms, as the help of os_moments() says of two extremes.
  arc_cov <- function(n, r, s, terms = 12L) {
    t <- n + 1 - s
    lrise <- function(a, k) sum(log(a + seq_len(k) - 1))          # log of the rising factorial
    ck <- function(k) (-1)^(k + 1) * pi^(2 * k) / (2 * factorial(2 * k))
    tot <- 0
    for (i in seq_len(terms)) for (j in seq_len(terms)) {
      m <- 2 * i; l <- 2 * j
      tot <- tot + ck(i) * ck(j) * exp(lrise(r, m) + lrise(t, l) - lrise(n + 1, m) - lrise(n + 1, l)) *
        (-expm1(-sum(log1p(m / (n + seq_len(l))))))
    }
    tot
  }
  # the series itself, at n = 2, where the covariance is E[X]^2 - alpha_1 alpha_2 exactly
  z <- os_moments(2, dist = "beta", shape1 = 0.5, shape2 = 0.5)
  expect_lt(abs(arc_cov(2, 1, 2, terms = 30L) / (0.25 - prod(z$alpha)) - 1), 1e-12)
  for (a in list(c(1e3, 1, 1e-9), c(1e6, 1, 1e-3), c(1e6, 2, 1e-3))) {
    n <- a[1]; rk <- c(a[2], n + 1 - a[2])
    z <- os_moments(n, rk, "beta", shape1 = 0.5, shape2 = 0.5)
    ex <- arc_cov(n, rk[1], rk[2])
    expect_lt(abs(z$B[1, 2] - ex) / sqrt(prod(z$var)), 1e-10)     # in units of the correlation
    expect_lt(abs(z$B[1, 2] / ex - 1), a[3])                      # in relative terms
  }
})
