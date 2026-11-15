# ==============================================================================
# File: test-os_engine.R
# Purpose: The package's normal moments against the reference engine, an
#          earlier implementation of the same method for the normal written
#          separately before the package, on 39 sample sizes from 5 to 1001 at
#          the five-number ranks (fixtures/01-frozen-grid.csv), the
#          user-distribution path against the built-in normal, and the
#          cancellation-free difference.

test_that("the normal at the five-number ranks agrees with the reference engine on the grid", {
  skip_on_cran()                     # about 15 s
  fx <- read.csv(test_path("fixtures", "01-frozen-grid.csv"))
  expect_equal(length(unique(fx$n)), 39L)
  worst_alpha <- 0; worst_B <- 0
  for (n in unique(fx$n)) {
    f <- fx[fx$n == n, ]
    rk <- f$rank_i[f$what == "alpha"]
    m <- os_moments(n, rk, "norm")
    worst_alpha <- max(worst_alpha, abs(m$alpha - f$value[f$what == "alpha"]))
    fb <- f[f$what == "B", ]
    worst_B <- max(worst_B, abs(m$B[cbind(fb$i, fb$j)] - fb$value))
  }
  expect_lt(worst_alpha, 1e-10)
  expect_lt(worst_B, 1e-10)
})

test_that("a normal supplied through os_dist() reproduces the built-in normal", {
  user <- os_dist(pnorm, dnorm, qnorm, symmetric = TRUE, name = "user-normal")
  a <- os_moments(20, c(1, 5, 10, 20), user)
  b <- os_moments(20, c(1, 5, 10, 20), "norm")
  expect_lt(max(abs(a$alpha - b$alpha)), 1e-12)
  expect_lt(max(abs(a$B - b$B)), 1e-12)
  expect_identical(a$route, "integrate")
})

test_that("the cancellation-free difference agrees with the direct difference where that is accurate", {
  cd <- orderStatMoments:::.os_centered(orderStatMoments:::.os_family("norm")$dist)
  y <- -1.2; z <- c(-0.9, -0.3, 0.4, 2.5)
  expect_equal(orderStatMoments:::.os_logdiff(y, z, cd), log(pnorm(z) - pnorm(y)), tolerance = 1e-13)
  y <- 0.7; z <- c(0.8, 1.5, 3)
  expect_equal(orderStatMoments:::.os_logdiff(y, z, cd), log(pnorm(z) - pnorm(y)), tolerance = 1e-13)
  # deep in the upper tail the direct difference is zero in doubles and the tail form is not
  expect_true(is.finite(orderStatMoments:::.os_logdiff(40, 41, cd)))
  expect_identical(orderStatMoments:::.os_logdiff(0.5, 0.5, cd), -Inf)
})

test_that("a divergent second moment stops even when the tail index is not declared", {
  # t with 1.9 degrees of freedom has no variance. Declared without its tail index, the
  # tail pieces of E[X^2] came back as -11.3 each, flagged probably divergent with an
  # error estimate below the acceptance threshold, and gave a variance of -19.
  t19 <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                   pt(q, 1.9, lower.tail = lower.tail, log.p = log.p),
                 function(x, log = FALSE) dt(x, 1.9, log = log),
                 function(p) qt(p, 1.9), symmetric = TRUE)
  expect_error(os_moments(1, dist = t19), "probably divergent")
})

test_that("a divergence flag on a bounded piece is spurious and the piece is accepted", {
  # the beta(1/2, 1/2) has an infinite density at both endpoints; integrate() flagged an
  # inner piece near 1 as probably divergent with an error estimate of 1e-11, and the
  # rule that stopped on any such flag refused the covariance of two observations, whose
  # exact value is 4/pi^4 (X_(1) X_(2) = X_1 X_2 and the beta moments)
  b <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                 pbeta(q, 0.5, 0.5, lower.tail = lower.tail, log.p = log.p),
               function(x, log = FALSE) dbeta(x, 0.5, 0.5, log = log),
               function(p) qbeta(p, 0.5, 0.5), support = c(0, 1), symmetric = TRUE)
  m <- os_moments(2, dist = b)
  expect_lt(abs(m$B[1, 2] - 4 / pi^4), 1e-9)
})

# ---- a density that is not smooth inside the support: os_dist(breaks = ) ----
# The exact values used below hold for any parent distribution: the order statistics are a
# permutation of the sample, so sum(alpha) = n E[X] and sum(B) = n Var(X), and at n = 2,
# X_(1) X_(2) = X_1 X_2, so Cov(X_(1), X_(2)) = E[X]^2 - alpha_1 alpha_2.

.laplace <- function(breaks = NULL) os_dist(
  p = function(q, lower.tail = TRUE, log.p = FALSE) {
    s <- if (lower.tail) q else -q
    lp <- ifelse(s <= 0, log(0.5) + s, log1p(-0.5 * exp(-abs(s))))
    if (log.p) lp else exp(lp)
  },
  d = function(x, log = FALSE) if (log) log(0.5) - abs(x) else 0.5 * exp(-abs(x)),
  q = function(p) ifelse(p <= 0.5, log(2 * p), -log(2 * (1 - p))),
  symmetric = TRUE, name = "laplace", breaks = breaks)

test_that("the Laplace distribution with its kink declared has exact covariances", {
  # n = 2: the maximum has expected value E|X_1 - X_2|/2 = 3/4, so the covariance is 9/16
  m2 <- os_moments(2, dist = .laplace(0))
  expect_equal(unname(m2$alpha), c(-0.75, 0.75), tolerance = 1e-13)
  expect_lt(abs(m2$B[1, 2] - 9 / 16), 1e-12)
  # n = 5: the covariances sum to n Var(X) = 10
  m5 <- os_moments(5, dist = .laplace(0))
  expect_lt(abs(sum(m5$B) - 10), 1e-11)
  # without the declaration the same sum is off by more than 1e-9 (1.7e-8 when written),
  # which is what the argument is for; the means and the variances agree far more closely
  u5 <- os_moments(5, dist = .laplace())
  expect_gt(abs(sum(u5$B) - 10), 1e-9)
  expect_lt(max(abs(u5$alpha - m5$alpha)), 1e-11)
  expect_lt(max(abs(u5$var - m5$var)), 1e-11)
})

test_that("a kink away from the median is cut in every integral (a triangular distribution)", {
  m <- 0.3                                    # the mode, where the density has its kink
  tri <- os_dist(
    p = function(q, lower.tail = TRUE, log.p = FALSE) {
      q <- pmin(pmax(q, 0), 1)
      lo <- q^2 / m; up <- (1 - q)^2 / (1 - m)        # F below the mode, 1 - F above it
      v <- if (lower.tail) ifelse(q <= m, lo, 1 - up) else ifelse(q <= m, 1 - lo, up)
      if (log.p) log(v) else v
    },
    d = function(x, log = FALSE) {
      v <- ifelse(x < 0 | x > 1, 0, ifelse(x <= m, 2 * x / m, 2 * (1 - x) / (1 - m)))
      if (log) log(v) else v
    },
    q = function(p) ifelse(p <= m, sqrt(p * m), 1 - sqrt((1 - p) * (1 - m))),
    support = c(0, 1), name = "triangular", breaks = m)
  mu <- (1 + m) / 3; v <- (1 + m^2 - m) / 18  # E[X] and Var(X) of the triangular on [0, 1]
  for (n in c(2, 5)) {
    o <- os_moments(n, dist = tri)
    expect_lt(abs(sum(o$alpha) - n * mu), 1e-12)
    expect_lt(abs(sum(o$B) - n * v), 1e-12)
  }
  o <- os_moments(2, dist = tri)
  expect_lt(abs(o$B[1, 2] - (mu^2 - o$alpha[[1]] * o$alpha[[2]])), 1e-12)
})

test_that("a jump of the density is cut, and stops the computation when it is not declared", {
  # the density 3/2 on [0, 1/2] and 1/2 on (1/2, 1]
  two <- function(breaks = NULL) os_dist(
    p = function(q, lower.tail = TRUE, log.p = FALSE) {
      q <- pmin(pmax(q, 0), 1)
      v <- if (lower.tail) ifelse(q <= 0.5, 1.5 * q, 0.75 + 0.5 * (q - 0.5))
           else ifelse(q <= 0.5, 1 - 1.5 * q, 0.5 * (1 - q))
      if (log.p) log(v) else v
    },
    d = function(x, log = FALSE) {
      v <- ifelse(x < 0 | x > 1, 0, ifelse(x <= 0.5, 1.5, 0.5))
      if (log) log(v) else v
    },
    q = function(p) ifelse(p <= 0.75, p / 1.5, 0.5 + (p - 0.75) / 0.5),
    support = c(0, 1), name = "two-level", breaks = breaks)
  mu <- 3 / 8; v <- 1.5 * (1 / 24) + 0.5 * (7 / 24) - mu^2    # E[X] = 3/8, E[X^2] = 5/24
  o <- os_moments(5, dist = two(0.5))
  expect_lt(abs(sum(o$alpha) - 5 * mu), 1e-12)
  expect_lt(abs(sum(o$B) - 5 * v), 1e-12)
  expect_error(os_moments(5, dist = two()), "breaks")
})

test_that("a declared point where the density is zero is accepted", {
  # f(x) = 3.2 |x - 1/4| on [0, 1]: a kink at 1/4 where the density is zero, so that its
  # log is -Inf there, which the probe of the log-density refused
  zero <- os_dist(
    p = function(q, lower.tail = TRUE, log.p = FALSE) {
      q <- pmin(pmax(q, 0), 1)
      v <- ifelse(q <= 0.25, 3.2 * (0.25 * q - q^2 / 2), 0.1 + 1.6 * (q - 0.25)^2)
      if (!lower.tail) v <- 1 - v
      if (log.p) log(v) else v
    },
    d = function(x, log = FALSE) {
      v <- ifelse(x < 0 | x > 1, 0, 3.2 * abs(x - 0.25))
      if (log) log(v) else v
    },
    q = function(p) ifelse(p <= 0.1, 0.25 - sqrt(pmax(0, 0.0625 - p / 1.6)),
                           0.25 + sqrt(pmax(0, (p - 0.1) / 1.6))),
    support = c(0, 1), breaks = 0.25)
  # E[X^k] = 3.2 {int_0^(1/4) x^k (1/4 - x) dx + int_(1/4)^1 x^k (x - 1/4) dx}
  mom <- function(k) {
    G <- function(x) x^(k + 2) / (k + 2) - 0.25 * x^(k + 1) / (k + 1)   # antiderivative of x^k (x - 1/4)
    3.2 * (-(G(0.25) - G(0)) + (G(1) - G(0.25)))
  }
  o <- os_moments(5, dist = zero)
  expect_lt(abs(sum(o$alpha) - 5 * mom(1)), 1e-12)
  expect_lt(abs(sum(o$B) - 5 * (mom(2) - mom(1)^2)), 1e-12)
})

test_that("a declared point that rounds to an endpoint in the standardized units is still a cut", {
  # X = U^10, where U has the density 2 on [0, 0.001] and (1 - 0.002)/0.999 above: the
  # density of X is infinite at 0 and jumps at 1e-30, which equals the endpoint 0 in the
  # standardized units. The cut is carried by its distance from the endpoint.
  pw <- 10; knot <- 0.001; h1 <- 2; h2 <- (1 - h1 * knot) / (1 - knot); p0 <- h1 * knot
  d <- os_dist(
    p = function(q, lower.tail = TRUE, log.p = FALSE) {
      u <- pmin(pmax(q, 0), 1)^(1 / pw)
      v <- if (lower.tail) ifelse(u <= knot, h1 * u, p0 + h2 * (u - knot))
           else ifelse(u <= knot, 1 - h1 * u, h2 * (1 - u))
      if (log.p) log(v) else v
    },
    d = function(x, log = FALSE) {
      u <- pmin(pmax(x, 0), 1)^(1 / pw)
      ld <- log(ifelse(u <= knot, h1, h2)) - log(pw) + (1 / pw - 1) * log(pmax(x, 0))
      ld[x < 0 | x > 1] <- -Inf
      if (log) ld else exp(ld)
    },
    q = function(p) ifelse(p <= p0, p / h1, knot + (p - p0) / h2)^pw,
    support = c(0, 1), breaks = knot^pw)
  cd <- .os_centered(d)
  expect_identical(unname(cd$kinks), cd$lower)          # the point is the endpoint in t
  expect_gt(cd$kink_rows[1, "lo"], 0)                   # and is apart from it by its distance
  # n = 2 against the exact parent moments E[X^k] = E[U^(10 k)]
  mom <- function(k) (h1 * knot^(pw * k + 1) + h2 * (1 - knot^(pw * k + 1))) / (pw * k + 1)
  o <- os_moments(2, dist = d)
  expect_lt(abs(sum(o$alpha) / (2 * mom(1)) - 1), 1e-12)
  expect_lt(abs(sum(o$B) / (2 * (mom(2) - mom(1)^2)) - 1), 1e-12)
  # and the covariance itself: X_(1) X_(2) = X_1 X_2, so Cov = E[X]^2 - alpha_1 alpha_2
  expect_lt(abs(o$B[1, 2] / (mom(1)^2 - o$alpha[[1]] * o$alpha[[2]]) - 1), 1e-11)
  # the minimum of n = 1000, which sits at the jump: E[X_(1)^k] = int 10k u^(10k - 1) {1 - F_U(u)}^n du,
  # an incomplete beta function on each side of the jump
  n <- 1000
  exmin <- function(k) pw * k * beta(pw * k, n + 1) *
    (pbeta(h1 * knot, pw * k, n + 1) / h1^(pw * k) +
       h2^n * pbeta(knot, pw * k, n + 1, lower.tail = FALSE))
  o <- os_moments(n, ranks = 1:2, dist = d)
  expect_lt(abs(o$alpha[[1]] / exmin(1) - 1), 1e-10)
  expect_lt(abs(o$var[[1]] / (exmin(2) - exmin(1)^2) - 1), 1e-9)
})

test_that("the points of a range are ordered and kept by their distances, not by a rounded t", {
  # Two cuts next to the lower endpoint -1 whose t, formed in two ways, differ by a rounding
  # step in the order opposite to their distances (the second review's fix left the order to
  # t wherever the t differ, which gave a piece with its limits reversed).
  t_far <- -1 + 2^-52; t_near <- -1 + 2^-53          # t_near < t_far
  cuts <- rbind(c(t = t_far, lo = 3.0e-18, hi = NA),   # the smaller distance has the larger t
                c(t = t_near, lo = 3.1e-18, hi = NA))
  ends <- list(lo = list(at = -1), hi = NULL, off = 0)
  pcs <- .os_pieces(-1, 2, cuts, ends, NULL)
  a <- vapply(pcs, function(p) p$a, 0); b <- vapply(pcs, function(p) p$b, 0)
  expect_identical(vapply(pcs, function(p) p$side, ""), rep("lo", 3L))
  expect_equal(a, c(0, 3.0e-18, 3.1e-18)); expect_equal(b, c(3.0e-18, 3.1e-18, 3))
  expect_true(all(a < b))
  # an inner range that starts at the distance 3.05e-18 keeps the cut above it only
  pcs <- .os_pieces(-1, 2, cuts, ends, list(side = "lo", e = 3.05e-18))
  expect_equal(c(pcs[[1]]$a, pcs[[1]]$b), c(3.05e-18, 3.1e-18))
  expect_length(pcs, 2L)
})

test_that("a declared zero of the density that is also an automatic cut point is accepted", {
  # the normal with its density set to 0 at one point, the lower cut point of the minimum of
  # five: the same distribution, and the probe of the cut points met a log-density of -Inf
  b <- .os_breaks(1, 5, .os_centered(os_dist(pnorm, dnorm, qnorm)), 10, kinks = FALSE)[1]
  d <- os_dist(pnorm, function(x, log = FALSE) {
    v <- dnorm(x, log = log)
    v[x == b] <- if (log) -Inf else 0
    v
  }, qnorm, breaks = b)
  expect_lt(max(abs(os_moments(5, dist = d)$B - os_moments(5)$B)), 1e-12)
})

test_that("next to an endpoint view the mass and the two moments are integrated over the same cuts", {
  # The mass check guards the pieces of the moments only if they are the same pieces. The
  # cuts that .os_single() hands to the integrator are recorded for the three powers.
  cd <- .os_centered(.os_resolve_dist("beta", list(shape1 = 1, shape2 = 1 / 40))$dist)
  ctrl <- .os_ctrl(1e-11)
  n <- 1e6; r <- 587180; p <- r / (n + 1)
  seen <- list()
  local_mocked_bindings(.os_integrate_split = function(f, lower, upper, breaks, ...) {
    seen[[length(seen) + 1L]] <<- breaks
    list(value = 1, abs.error = 0)
  })
  for (pow in 0:2) .os_single(r, n, pow, cd, ctrl, shift = cd$q(p), u = p)
  expect_true(is.matrix(seen[[1]]) && nrow(seen[[1]]) > 2L)
  expect_identical(seen[[1]], seen[[2]])
  expect_identical(seen[[2]], seen[[3]])
  # on the side without a view (the lower one here) the cut at the shift is the first moment's only
  seen <- list(); r <- 100; p <- r / (n + 1)
  for (pow in 0:2) .os_single(r, n, pow, cd, ctrl, shift = cd$q(p), u = p)
  expect_identical(seen[[1]], seen[[3]])
  expect_identical(nrow(seen[[2]]), nrow(seen[[1]]) + 1L)
})

test_that("the endpoint an order statistic is next to is decided by its distance, not by its rank", {
  side <- function(cd, p, ...) vapply(p, function(u) .os_side(cd, u, ...), "")
  p <- c(0.2, 0.49, 0.51, 0.8)
  beta_cd <- function(a, b) .os_centered(.os_resolve_dist("beta", list(shape1 = a, shape2 = b))$dist)
  # a view at both endpoints: the median of a beta(0.015, 0.8) is 1e-20, so every rank up
  # to 0.8 n is next to the lower endpoint, and the mirror image holds for a beta(0.8, 0.05)
  expect_identical(side(beta_cd(0.015, 0.8), p), rep("lo", 4L))
  expect_identical(side(beta_cd(0.8, 0.05), p), rep("hi", 4L))
  # a symmetric distribution: the former rule, the lower endpoint up to the median
  expect_identical(side(beta_cd(0.5, 0.5), c(p, 0.5)), c("lo", "lo", "hi", "hi", "lo"))
  # one view and an infinite other endpoint: the view at every rank
  expect_identical(side(.os_centered(.os_resolve_dist("gamma", list(shape = 0.5))$dist), c(0.01, 0.99)),
                   c("lo", "lo"))
  # one view and a finite other endpoint without one (a user's beta(1/2, 1/2)): next to the
  # lower endpoint below the median, next to no view above it, where the density at the
  # center is still taken through the one view there is
  ub <- .os_centered(os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                               pbeta(q, 0.5, 0.5, lower.tail = lower.tail, log.p = log.p),
                             function(x, log = FALSE) dbeta(x, 0.5, 0.5, log = log),
                             function(p) qbeta(p, 0.5, 0.5), support = c(0, 1)))
  expect_identical(side(ub, c(0.2, 0.99)), c("lo", ""))
  expect_identical(side(ub, 0.99, nearest = FALSE), "lo")
  # no view: nothing changes
  nrm <- .os_centered(os_dist(pnorm, dnorm, qnorm))
  expect_identical(c(.os_side(nrm, 0.3), .os_side(nrm, 0.3, nearest = FALSE)), c("", ""))
  # the width of the rank 510000 of a million from a beta(0.015, 0.8): taken through the
  # upper endpoint's view it was 0, which left the absolute tolerances at the scale one for
  # an order statistic whose standard deviation is 3e-21
  w <- .os_width(510000, 1e6, beta_cd(0.015, 0.8), 10)$w
  expect_gt(w, 0)
  expect_lt(w, 1e-10)
})
