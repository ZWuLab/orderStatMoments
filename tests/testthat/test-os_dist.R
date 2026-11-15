# ==============================================================================
# File: test-os_dist.R
# Purpose: The distribution object: the validation of os_dist(), the built-in
#          families' centers and symmetry flags, and the Gumbel functions'
#          log-tails.

test_that("os_dist() validates its arguments and the functions' behavior at the median", {
  expect_error(os_dist(1, dnorm, qnorm), "functions")
  expect_error(os_dist(pnorm, dnorm, qnorm, support = c(1, -1)), "support")
  expect_error(os_dist(pnorm, dnorm, qnorm, symmetric = NA), "symmetric")
  expect_error(os_dist(pnorm, dnorm, qnorm, tail = c(0, 1)), "tail")
  expect_error(os_dist(pnorm, dnorm, qnorm, name = 1), "name")
  expect_error(os_dist(pnorm, dnorm, qexp), "disagree")                     # p and q of different laws
  expect_error(os_dist(function(q, lower.tail = TRUE, log.p = FALSE) pnorm(q), dnorm, qnorm),
               "lower.tail")                                                # ignores both arguments
  expect_error(os_dist(function(q, lower.tail = TRUE, log.p = FALSE) pnorm(q, lower.tail = lower.tail),
                       dnorm, qnorm), "log.p")                              # ignores log.p
  expect_error(os_dist(pnorm, function(x, log = FALSE) dnorm(x), qnorm), "log = TRUE")
  expect_error(os_dist(pnorm, dnorm, qnorm, support = c(1, 2)), "outside")
  d <- os_dist(pnorm, dnorm, qnorm, symmetric = TRUE, name = "z")
  expect_s3_class(d, "os_dist")
  expect_equal(d$center, 0)
  expect_identical(d$tail, c(Inf, Inf))
  expect_null(d$closed$mean)
})

test_that("the built-in families have the right centers, supports and symmetry", {
  f <- function(name) orderStatMoments:::.os_family(name)$dist
  expect_equal(f("norm")$center, 0); expect_true(f("norm")$symmetric)
  expect_equal(f("unif")$center, 0.5); expect_true(f("unif")$symmetric)
  expect_identical(f("unif")$support, c(0, 1))
  expect_equal(f("exp")$center, log(2)); expect_false(f("exp")$symmetric)
  expect_identical(f("exp")$support, c(0, Inf))
  expect_equal(f("logis")$center, 0); expect_true(f("logis")$symmetric)
  expect_equal(f("gumbel")$center, -log(log(2))); expect_false(f("gumbel")$symmetric)
  expect_error(orderStatMoments:::.os_resolve_dist("pois"), "os_dist")
})

test_that("the Gumbel functions keep their log-tails finite and accurate", {
  p <- orderStatMoments:::.pgumbel; d <- orderStatMoments:::.dgumbel; q <- orderStatMoments:::.qgumbel
  x <- seq(-5, 700, by = 0.5)
  expect_true(all(is.finite(p(x, log.p = TRUE))))
  expect_true(all(is.finite(p(x, lower.tail = FALSE, log.p = TRUE))))
  # F + S = 1 in the center
  xc <- seq(-3, 3, by = 0.25)
  expect_lt(max(abs(p(xc) + p(xc, lower.tail = FALSE) - 1)), 1e-15)
  # the two branches of log S agree where both are accurate, and log S ~ -x far out
  xb <- c(29, 29.5, 29.9, 30.1, 30.5, 31)
  expect_lt(max(abs(p(xb, lower.tail = FALSE, log.p = TRUE) - (-xb + log1p(-exp(-xb) / 2)))), 1e-14)
  xf <- c(40, 100, 500)
  expect_lt(max(abs(p(xf, lower.tail = FALSE, log.p = TRUE) + xf)), 1e-15)
  # the density integrates to one and q inverts p
  expect_lt(abs(integrate(d, -Inf, Inf, rel.tol = 1e-12)$value - 1), 1e-10)
  expect_lt(max(abs(q(p(xc)) - xc)), 1e-12)
  expect_lt(max(abs(d(xc, log = TRUE) - log(d(xc)))), 1e-14)
})

test_that("a gamma supplied through os_dist() reproduces its exact rational moments at n = 5 and 10", {
  # The gamma is not a built-in family, so this exercises the user-distribution path end
  # to end against an external exact value. For the gamma of shape 2,
  # F(x) = 1 - (1 + x) exp(-x), so every moment of an order statistic is a finite sum of
  # rationals: expand F(x)^(k-1) and then (1 + x)^m by the binomial theorem and integrate
  # x^(s+1+i) exp(-c x) termwise, which gives (s+1+i)!/c^(s+2+i). The fixture holds those
  # sums evaluated in exact rational arithmetic (fixtures/README.md); the printed table
  # of Gupta (1960, Technometrics 2, Table I) is wrong at n = 10 from rank 4 on, so the
  # exact values and not the table are the reference here.
  ex <- read.csv(test_path("fixtures", "04-gamma2-exact.csv"))
  gamma2 <- os_dist(p = function(q, lower.tail = TRUE, log.p = FALSE)
                      pgamma(q, shape = 2, lower.tail = lower.tail, log.p = log.p),
                    d = function(x, log = FALSE) dgamma(x, shape = 2, log = log),
                    q = function(p) qgamma(p, shape = 2),
                    support = c(0, Inf), name = "gamma(shape = 2)")
  for (n in c(5L, 10L)) {
    m <- os_moments(n, dist = gamma2, cov = FALSE)
    e <- ex[ex$n == n, ]
    expect_identical(m$route, "integrate")
    expect_lt(max(abs(unname(m$alpha) - e$mean)), 1e-10)
    expect_lt(max(abs(unname(m$var) - e$var)), 1e-10)
  }
})

test_that("a function without the required arguments, or an NA, gives an error that says what to fix", {
  expect_error(os_dist(function(q) pnorm(q), dnorm, qnorm),
               "p\\(q, lower.tail = TRUE, log.p = FALSE\\)", class = "error")
  expect_error(os_dist(pnorm, function(x) dnorm(x), qnorm), "d\\(x, log = FALSE\\)")
  expect_error(os_dist(function(q, lower.tail = TRUE, log.p = FALSE) NA_real_, dnorm, qnorm),
               "disagree")
  # a gamma declared on the whole line: its log-tail is -Inf below zero
  g <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                 pgamma(q, 2, lower.tail = lower.tail, log.p = log.p),
               function(x, log = FALSE) dgamma(x, 2, log = log), function(p) qgamma(p, 2))
  expect_error(os_moments(5, dist = g), "support")
})

test_that("a density that disagrees with the distribution function is refused", {
  # 1.001 dnorm passed every check at the median and gave a variance of 1.001
  expect_error(os_dist(pnorm, function(x, log = FALSE) if (log) log(1.001) + dnorm(x, log = TRUE)
                                                        else 1.001 * dnorm(x), qnorm),
               "integrates to")
})

test_that("os_dist() validates the points where the density is not smooth", {
  expect_null(os_dist(pnorm, dnorm, qnorm)$breaks)                          # none by default
  expect_error(os_dist(pnorm, dnorm, qnorm, breaks = "a"), "breaks")
  expect_error(os_dist(pnorm, dnorm, qnorm, breaks = c(1, NA)), "finite numbers")
  expect_error(os_dist(pnorm, dnorm, qnorm, breaks = Inf), "finite numbers")
  expect_error(os_dist(pnorm, dnorm, qnorm, breaks = numeric(0)), "finite numbers")
  # an endpoint of the support is not a point inside it
  expect_error(os_dist(pexp, dexp, qexp, support = c(0, Inf), breaks = 0), "strictly inside")
  expect_error(os_dist(pexp, dexp, qexp, support = c(0, Inf), breaks = c(1, -1)), "strictly inside")
  # a point where the density is infinite is not a kink or a jump
  expect_error(os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                         pchisq(q + 1, 0.5, lower.tail = lower.tail, log.p = log.p),
                       function(x, log = FALSE) dchisq(x + 1, 0.5, log = log),
                       function(p) qchisq(p, 0.5) - 1, breaks = -1), "finite at 'breaks'")
  # a symmetric density has its non-smooth points symmetric about the median
  expect_error(os_dist(pnorm, dnorm, qnorm, symmetric = TRUE, breaks = 1), "symmetric about")
  expect_identical(os_dist(pnorm, dnorm, qnorm, symmetric = TRUE, breaks = c(1, -1))$breaks, c(-1, 1))
  # stored sorted, without repeats
  expect_identical(os_dist(pnorm, dnorm, qnorm, breaks = c(2, -1, 2))$breaks, c(-1, 2))
})

test_that("the density is integrated between the quartiles piece by piece between the declared points", {
  # 0.9 Uniform(0, 1) + 0.1 Uniform(0.399, 0.401): the narrow component lies between two
  # jumps, and one integral over the quartiles misses it (0.4 where F rises by 0.5)
  lo <- 0.399; hi <- 0.401; w <- hi - lo
  p <- function(q, lower.tail = TRUE, log.p = FALSE) {
    x <- pmin(pmax(q, 0), 1)
    v <- if (lower.tail) 0.9 * x + 0.1 * pmin(pmax((x - lo) / w, 0), 1)
         else 0.9 * (1 - x) + 0.1 * pmin(pmax((hi - x) / w, 0), 1)
    if (log.p) log(v) else v
  }
  d <- function(x, log = FALSE) {
    v <- 0.9 * (x >= 0 & x <= 1) + 0.1 / w * (x >= lo & x <= hi)
    if (log) log(v) else v
  }
  q <- function(p) ifelse(p < 0.9 * lo, p / 0.9,
                          ifelse(p > 0.9 * hi + 0.1, (p - 0.1) / 0.9, lo + (p - 0.9 * lo) / (0.9 + 0.1 / w)))
  expect_error(os_dist(p, d, q, support = c(0, 1)), "integrates to")       # the jumps not declared
  mix <- os_dist(p, d, q, support = c(0, 1), breaks = c(lo, hi))
  expect_identical(mix$breaks, c(lo, hi))
  # E[X] = 0.9/2 + 0.1 * 0.4 and E[X^2] = 0.9/3 + 0.1 (0.4^2 + w^2/12), and the order
  # statistics are a permutation of the sample
  mu <- 0.49; v <- 0.3 + 0.1 * (0.16 + w^2 / 12) - mu^2
  o <- os_moments(5, dist = mix)
  expect_lt(abs(sum(o$alpha) - 5 * mu), 1e-12)
  expect_lt(abs(sum(o$B) - 5 * v), 1e-12)
})

test_that("the symmetry of the declared points is judged on the distribution's own scale", {
  # a normal with a small standard deviation and one point a standard deviation above the
  # median, with no mirror point: refused at every scale (a tolerance with a floor of 1e-8
  # in absolute units accepted it at sd = 1e-9)
  for (sc in c(1, 1e-3, 1e-9, 1e6))
    expect_error(os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                           pnorm(q, sd = sc, lower.tail = lower.tail, log.p = log.p),
                         function(x, log = FALSE) dnorm(x, sd = sc, log = log),
                         function(p) qnorm(p, sd = sc), symmetric = TRUE, breaks = sc),
                 "symmetric about")
  # a median so large against the spread that the rounding of a point reaches the tolerance:
  # refused with the advice to recenter (a uniform on 2^50 +/- 2 with one point was accepted)
  cc <- 2^50
  expect_error(os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                         punif(q, cc - 2, cc + 2, lower.tail = lower.tail, log.p = log.p),
                       function(x, log = FALSE) dunif(x, cc - 2, cc + 2, log = log),
                       function(p) qunif(p, cc - 2, cc + 2), support = c(cc - 2, cc + 2),
                       symmetric = TRUE, breaks = cc + 1), "too large against the spread")
  # the mirror pair is accepted at every scale, also about a median far from 0
  for (sc in c(1, 1e-9, 1e6)) {
    o <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                   pnorm(q, 3 * sc, sc, lower.tail = lower.tail, log.p = log.p),
                 function(x, log = FALSE) dnorm(x, 3 * sc, sc, log = log),
                 function(p) qnorm(p, 3 * sc, sc), symmetric = TRUE, breaks = c(2, 4) * sc)
    expect_length(o$breaks, 2L)
  }
})

test_that("the density is integrated in log|x| between quartiles that lie orders of magnitude apart", {
  # X = Y^14, where Y has a piecewise constant density on [0, 1] with three jumps: the
  # density of X is infinite at 0, its quartiles lie a factor 1e7 apart, and the plain
  # integral between them returned 0.74 where the distribution function rises by 0.5
  K <- 14; edges <- c(0, 0.005, 0.01, 0.02, 1); h <- c(2, 3, 0.1)
  h <- c(h, (1 - sum(diff(edges)[1:3] * h)) / (1 - edges[4]))
  probs <- c(0, cumsum(diff(edges) * h)); probs[5] <- 1
  py <- function(y, lower.tail) {
    j <- findInterval(y, edges, all.inside = TRUE)
    if (lower.tail) probs[j] + h[j] * (y - edges[j]) else (1 - probs[j + 1]) + h[j] * (edges[j + 1] - y)
  }
  d <- os_dist(
    p = function(q, lower.tail = TRUE, log.p = FALSE) {
      v <- py(pmin(pmax(q, 0), 1)^(1 / K), lower.tail)
      if (log.p) log(v) else v
    },
    d = function(x, log = FALSE) {
      j <- findInterval(pmin(pmax(x, 0), 1)^(1 / K), edges, all.inside = TRUE)
      v <- log(h[j] / K) + (1 / K - 1) * log(pmax(x, 0))
      v[x < 0 | x > 1] <- -Inf
      if (log) v else exp(v)
    },
    q = function(p) {
      j <- findInterval(p, probs, all.inside = TRUE)
      (edges[j] + (p - probs[j]) / h[j])^K
    },
    support = c(0, 1), breaks = edges[2:4]^K)
  expect_s3_class(d, "os_dist")
  # E[X^k] = E[Y^(14 k)] = sum_j h_j (e_(j+1)^(m+1) - e_j^(m+1))/(m + 1), m = 14 k
  mom <- function(k) { m <- K * k; sum(h * diff(edges^(m + 1))) / (m + 1) }
  o <- os_moments(2, dist = d)
  expect_lt(abs(sum(o$alpha) / (2 * mom(1)) - 1), 1e-11)
  expect_lt(abs(sum(o$B) / (2 * (mom(2) - mom(1)^2)) - 1), 1e-11)
})
