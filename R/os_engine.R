# ==============================================================================
# File: os_engine.R
# Purpose: The internal integrators of orderStatMoments: the exact single and product
#          moments of order statistics by adaptive quadrature of the defining
#          integrals. The numerical design has three parts. What is computed:
#          integration in the standardized units t = {x - Q(1/2)}/spread, each
#          moment taken about the center of its own order statistic,
#          c_r = Q(r/(n+1)), so that E[X^2] - E[X]^2 does not cancel at a large
#          location or at a large n. How it is integrated: the range cut around
#          the peak of the order statistic's density, the tolerances scaled to
#          the order statistic, a flagged piece accepted only when negligible,
#          and the pieces next to an endpoint where the density is infinite
#          taken in the distance from that endpoint. How the integrand is
#          evaluated: on the logarithmic scale through log.p and lower.tail,
#          with log{F(z) - F(y)} formed on the side where the numbers are small.
#          Everything here works in the standardized units; os_moments() brings
#          the moments back and applies the shift and the multiplier.
#
# INTERNAL HELPERS:
#   - .os_ctrl()            The quadrature control constants from rel.tol
#   - .os_centered()        The log-tails, log-density and quantile function in centered units
#   - .os_end_view()        The same functions in the distance from a finite endpoint of the
#                           support where the density is infinite, with the endpoint's power
#   - .os_inf_density()     Whether the density is infinite at a point
#   - .os_ends()            The integrand in the endpoint coordinates of both ends
#   - .os_shift_dist()      The distances of a shift from the two endpoints
#   - .os_side()            The endpoint that an order statistic is next to, by its distance
#   - .os_cuts()            Cut points with their distances from the two endpoints
#   - .os_end_cuts()        The cuts around the peak of an order statistic next to such an endpoint
#   - .os_pieces()          The pieces of an integration range and the coordinate of each
#   - .os_breaks()          The cut points of rank r, c_r +/- h_r, and the declared points
#                           where the density is not smooth
#   - .os_width()           c_r and the delta-method half-width w_r of rank r
#   - .os_tolscale()        min(1, w): the scale of the absolute tolerances of rank r
#   - .os_integrate_split() integrate() over the pieces, values and error estimates summed
#   - .os_logdiff()         log{F(z) - F(y)} on the side where the numbers are small
#   - .os_xlog()            k * log(v) with the k = 0 guard against 0 * (-Inf)
#   - .os_single()          E[X_(r)^pow], pow = 1 or 2, one integral
#   - .os_prod()            E[X_(r) X_(s)], r < s, nested quadrature

#################### INTERNAL HELPER FUNCTIONS ####################

#' The quadrature control constants
#'
#' The constants of the reference engine, an earlier implementation of the same
#' method for the normal, written separately before the package: the single
#' integrals and the inner integral of
#' a product moment at rel.tol (1e-11 by default), the outer integral at 100
#' times rel.tol, so that the inner quadrature noise sits two orders below what
#' the outer adaptive rule chases; the absolute tolerances 1e-13 (single), 1e-15
#' (inner) and 1e-12 (outer); 400 subdivisions; a non-"OK" piece accepted when
#' its error estimate is below 1e-10; the half-width k = 10 of the split in
#' units of the large-sample standard deviation of the uniform order statistic.
#' @keywords internal
#' @noRd
.os_ctrl <- function(rel.tol = 1e-11) {
  list(rel.tol = rel.tol, rel.tol.outer = 100 * rel.tol,
       abs.tol.single = 1e-13, abs.tol.inner = 1e-15, abs.tol.outer = 1e-12,
       subdivisions = 400L, accept.err = 1e-10, k = 10)
}

#' The distribution in standardized units t = {x - Q(1/2)} / spread
#'
#' The standardized units. The spread is the standard member's interquartile range divided
#' by that of the standard normal, 2 qnorm(3/4) = 1.349, so that it is exactly 1
#' for the normal, for which t = x. Returns closures lF(t) = log F, lS(t) =
#' log(1 - F), ld(t) = log of the density of t (the density of x times the
#' spread), pf(t) = F and q(u) = {Q(u) - c} / spread, the support endpoints in
#' the same units, and the center c = Q(1/2) and the spread, which os_moments()
#' uses to bring the results back to the standard member's units. In these
#' units the median is 0, the quartiles are 1.349 apart (at -0.674 and 0.674 for
#' a symmetric parent, as for the normal), and the minimum half-width of the cut and the integrator's absolute tolerances, which
#' are fixed numbers, mean the same thing for every distribution, whatever the
#' scale on which it is supplied. (With the minimum half-width and the tolerances in absolute
#' units, a normal supplied with standard deviation 0.01 came back with variance
#' exactly zero, its whole mass between the nodes of a central piece 290 of its
#' standard deviations wide, and at 0.001 the second moments, of order 1e-6,
#' were integrated to seven digits only.)
#' @keywords internal
#' @noRd
.os_centered <- function(dist) {
  c0 <- dist$center
  p <- dist$p; d <- dist$d; q <- dist$q
  sp <- (q(0.75) - q(0.25)) / (2 * qnorm(0.75))
  lo <- .os_end_view(dist, sp, "lo"); hi <- .os_end_view(dist, sp, "hi")
  # the declared points where the density is not smooth, as rows (t, lo, hi) like the cuts
  # of .os_cuts(): the distance from an endpoint with a view is taken from the point and the
  # endpoint themselves, because a point next to the endpoint can round to it in t
  br <- as.numeric(dist$breaks)
  kink_rows <- cbind(t = (br - c0) / sp,
                     lo = if (is.null(lo)) rep(NA_real_, length(br)) else (br - dist$support[1]) / sp,
                     hi = if (is.null(hi)) rep(NA_real_, length(br)) else (dist$support[2] - br) / sp)
  list(lF = function(t) p(c0 + sp * t, lower.tail = TRUE, log.p = TRUE),
       lS = function(t) p(c0 + sp * t, lower.tail = FALSE, log.p = TRUE),
       ld = function(t) d(c0 + sp * t, log = TRUE) + log(sp),
       pf = function(t) p(c0 + sp * t),
       q  = function(u) (q(u) - c0) / sp,
       lower = (dist$support[1] - c0) / sp, upper = (dist$support[2] - c0) / sp,
       center = c0, spread = sp,
       kinks = kink_rows[, "t"], kink_rows = kink_rows,   # none: length 0 and no rows
       lo = lo, hi = hi)
}

#' The distribution in the distance from an endpoint where its density is infinite
#'
#' Why this exists: a gamma, a Weibull or a beta with a shape a below 1 has a
#' density that behaves like e^(a - 1) in the distance e from a finite endpoint
#' of its support, and so does the density of an order statistic near that
#' endpoint, with a times the rank in place of a. Two things then fail in the
#' standardized units centered at the median. First, the adaptive rule does not
#' resolve the algebraic singularity to the tolerance (a Weibull with shape 0.5
#' stopped with a failed integral for the minimum at n = 101, a gamma with shape
#' 0.1 at every n). Second, an order statistic that sits within 1e-13 spreads of
#' the endpoint cannot be placed at all, because t = (x - median)/spread carries
#' the distance from the endpoint as the small difference of two numbers of
#' order one. Both are repaired by integrating the pieces next to the endpoint in
#' the variable v = e^a, with e = v^(1/a) the distance from the endpoint in units
#' of the spread: the parent's distribution function is then nearly linear in v,
#' the integrand is bounded, and the distance is the integration variable itself,
#' so that it keeps its relative precision however small it is.
#'
#' Returns NULL when the endpoint is infinite or the density is finite there, and
#' then nothing changes. Otherwise returns the closures lF(e), lS(e), ld(e) and
#' pf(e) of .os_centered() as functions of the distance e, the function dq(u),
#' the distance of the quantile Q(u) from the endpoint, and the power pw = 1/a.
#' A built-in family carries its index a in `end` (the shape of a gamma or a
#' Weibull, half the degrees of freedom of a chi-square or of the numerator of
#' an F, the two shapes of a beta), and the integrand in v is then bounded. For
#' a user's distribution a is estimated as the slope of the log of the tail
#' probability against log e between e = 1e-9 and 1e-6 spreads. The estimate
#' is a heuristic: with a slightly wrong power the integrand in v keeps a weak
#' singularity v^(a/a_est - 1) (an exponent of -3e-8 for a beta(0.5, 0.2)), which
#' the adaptive rule resolves, and an estimate at or above 1 gives no view, so
#' that the endpoint is integrated in t as before the views existed. At the lower endpoint L the functions
#' are evaluated at x = L + spread * e, which is exact when L = 0, the case of
#' the gamma, the Weibull, the chi-square, the F and the beta. At the upper
#' endpoint U, x = U - spread * e loses the digits of e below the double epsilon
#' of U, so the view exists only when the object carries in `upper` the
#' functions p, d and q of the reflected variable U - X, which are exact (the
#' built-in beta, whose reflection is the beta with its two shapes exchanged).
#' Without them the upper endpoint is integrated in t like any other, as it was
#' before the views: a user's beta(1/2, 1/2) keeps the accuracy it had, about
#' 1e-10 for the maximum at n = 100, and stops at larger n.
#' @keywords internal
#' @noRd
.os_end_view <- function(dist, sp, side) {
  p <- dist$p; d <- dist$d; q <- dist$q
  E <- dist$support[if (side == "lo") 1L else 2L]
  if (!is.finite(E) || !.os_inf_density(d, E)) return(NULL)
  if (side == "hi" && is.null(dist$upper)) return(NULL)   # no exact functions of U - X
  lsp <- log(sp)
  if (side == "lo") {
    v <- list(lF = function(e) p(E + sp * e, lower.tail = TRUE, log.p = TRUE),
              lS = function(e) p(E + sp * e, lower.tail = FALSE, log.p = TRUE),
              ld = function(e) d(E + sp * e, log = TRUE) + lsp,
              pf = function(e) p(E + sp * e),
              dq = function(u) (q(u) - E) / sp)
    tail <- v$lF
  } else {
    up <- dist$upper                           # p, d and q of the reflected variable U - X
    v <- list(lF = function(e) up$p(sp * e, lower.tail = FALSE, log.p = TRUE),
              lS = function(e) up$p(sp * e, lower.tail = TRUE, log.p = TRUE),
              ld = function(e) up$d(sp * e, log = TRUE) + lsp,
              pf = function(e) up$p(sp * e, lower.tail = FALSE),
              dq = function(u) up$q(1 - u) / sp)
    tail <- v$lS
  }
  # the tail probability behaves like e^a: the index a is the family's own where the object
  # carries it in `end` (the built-in families), and is estimated otherwise
  a <- if (!is.null(dist$end)) dist$end[[side]] else (tail(1e-6) - tail(1e-9)) / log(1e3)
  if (!is.finite(a) || a <= 0 || a >= 1) return(NULL)
  v$pw <- 1 / a
  # A quantile function may return the endpoint itself for a small probability (R's qf()
  # returns 0 at p = 1e-4 for df1 = 0.5). The distance of Q(u) from the endpoint is then
  # taken from the power law log(tail probability) = l0 + a log(e) fitted above, which is
  # accurate enough for a center and a width.
  l0 <- tail(1e-6) - a * log(1e-6)
  dq <- v$dq
  v$dq <- function(u) {
    d <- dq(u)
    if (is.finite(d) && d > 0) d else exp((log(if (side == "lo") u else 1 - u) - l0) / a)
  }
  v
}

#' Whether a density is infinite at x (an error or a missing value counts as no)
#' @keywords internal
#' @noRd
.os_inf_density <- function(d, x) {
  v <- tryCatch(suppressWarnings(d(x)), error = function(e) NA_real_)
  is.numeric(v) && length(v) == 1L && !is.na(v) && v == Inf
}

#' The integrand in the endpoint coordinates, for .os_integrate_split()
#'
#' `make(ev, dev)` builds an integrand from a set of closures ev (lF, lS, ld, pf
#' in some coordinate) and the function dev that returns the deviation from the
#' shift in that coordinate. For each endpoint that has a view, the integrand is
#' built in the distance e from the endpoint, where the deviation t - shift is
#' e - dlo at the lower endpoint and dhi - e at the upper, with dlo and dhi the
#' distances of the shift from the two endpoints. Returns NULL when neither
#' endpoint has a view.
#' @keywords internal
#' @noRd
.os_ends <- function(cd, make, dlo, dhi) {
  if (is.null(cd$lo) && is.null(cd$hi)) return(NULL)
  list(lo = if (!is.null(cd$lo))
         list(f = make(cd$lo, function(e) e - dlo, "lo"), pw = cd$lo$pw, at = cd$lower),
       hi = if (!is.null(cd$hi))
         list(f = make(cd$hi, function(e) dhi - e, "hi"), pw = cd$hi$pw, at = cd$upper),
       off = cd$center / cd$spread)   # x/spread = t + off: the size of the argument x at t
}

#' The distances of a shift from the two endpoints
#'
#' When the shift is the quantile Q(u) of the centered distribution and u is
#' given, the distances come from the views' own dq(u), which do not pass
#' through the centered coordinate. Otherwise they are the plain differences.
#' @keywords internal
#' @noRd
.os_shift_dist <- function(cd, shift, u = NULL) {
  c(lo = if (is.null(cd$lo)) NA_real_ else if (is.null(u)) shift - cd$lower else cd$lo$dq(u),
    hi = if (is.null(cd$hi)) NA_real_ else if (is.null(u)) cd$upper - shift else cd$hi$dq(u))
}

#' The endpoint that an order statistic is next to
#'
#' The center Q(p) of the order statistic of probability p = r/(n+1) lies at the
#' distance dlo from the lower endpoint and dhi from the upper one, in units of
#' the spread. Each distance is taken from the endpoint's view where it has one
#' (the view's own dq(p), which keeps its relative precision however small it
#' is) and from the centered coordinate otherwise, and an infinite endpoint is
#' infinitely far. With `nearest = TRUE` the result is "lo" or "hi" when the
#' nearer endpoint has a view, and "" when it has none: the order statistic is
#' then next to an endpoint whose coordinate exists, and the cuts around its
#' peak (.os_end_cuts()) and the cut at its center for every power
#' (.os_single()) are made in that coordinate. With `nearest = FALSE` the result
#' is the nearer endpoint among those that have a view, "" when neither has
#' one: the view through which the density at the center (.os_width()) and the
#' mean (os_moments()) are formed. Why the distance and not the probability: a
#' beta with shapes 0.015 and 0.8 has a view at both endpoints and its median
#' lies 1e-20 from the lower one, so the rank 510000 of a million, whose
#' probability is above one half, is next to the lower endpoint. Chosen by
#' p <= 1/2, the upper endpoint's coordinate resolved neither the order
#' statistic nor its moments: the reflected argument rounded to 1, the width
#' came out as 0, the mean was assembled as 1 - 1 = 0, and the variance was 28
#' percent too high, with a mass of 1 to ten digits. A tie goes to the lower
#' endpoint, which is the former rule for a symmetric distribution.
#' @keywords internal
#' @noRd
.os_side <- function(cd, p, nearest = TRUE) {
  has_lo <- !is.null(cd$lo); has_hi <- !is.null(cd$hi)
  if (!has_lo && !has_hi) return("")
  dlo <- if (has_lo) cd$lo$dq(p) else cd$q(p) - cd$lower
  dhi <- if (has_hi) cd$hi$dq(p) else cd$upper - cd$q(p)
  if (dlo <= dhi) {
    if (has_lo) "lo" else if (nearest) "" else "hi"
  } else {
    if (has_hi) "hi" else if (nearest) "" else "lo"
  }
}

#' Cut points with their distances from the two endpoints
#'
#' A cut point is carried as a row (t, lo, hi): its value in standardized units
#' and its distances from the lower and the upper endpoint, NA without a view.
#' An order statistic can sit closer to an endpoint than the centered coordinate
#' resolves (the minimum of 1001 observations from a gamma with shape 0.1 lies
#' 2e-29 spreads above 0, and its center c_r equals the endpoint in t), so a cut
#' at c_r is ordered and placed by its distance, which is exact, and not by its
#' t. `d` gives those distances for the cuts that have them (.os_shift_dist());
#' the others are the plain differences.
#' @keywords internal
#' @noRd
.os_cuts <- function(cd, t, d = NULL) {
  if (!length(t)) return(matrix(numeric(0), 0L, 3L, dimnames = list(NULL, c("t", "lo", "hi"))))
  lo <- if (is.null(cd$lo)) rep(NA_real_, length(t)) else if (is.null(d)) t - cd$lower else rep(d[["lo"]], length(t))
  hi <- if (is.null(cd$hi)) rep(NA_real_, length(t)) else if (is.null(d)) cd$upper - t else rep(d[["hi"]], length(t))
  cbind(t = t, lo = lo, hi = hi)
}

#' The cuts around the peak of an order statistic next to an endpoint with a view
#'
#' Why this exists: the cut of .os_breaks() has a minimum half-width of
#' 4/sqrt(n + 1) spreads, which suits a density that is finite at the endpoint.
#' Next to an endpoint where the density is infinite, an extreme order statistic
#' is far narrower than that: in the variable v = e^a the distribution function
#' is nearly linear, F = C v, so that v at X_(r) is nearly a beta(r, n - r + 1)
#' variable divided by C, with center v_r and relative standard deviation
#' w = s_r/p_r, while the first piece runs to about (4/sqrt(n + 1))^a. At
#' n = 1e6 and a = 1/2 the peak is 1e-5 of the piece, the first nodes miss it,
#' and the density of the minimum of a gamma with shape 1/2 integrated to 3e-21
#' (the mass check stopped it), while for a beta(1/2, 1/2) the mass came out
#' right and the variance 2 percent off, without a word. The cuts here do in
#' v what .os_breaks() does in t: at v_r (1 - 10 w) and v_r (1 - 40 w) below the
#' center, where they are positive, and above it at v_r (1 + 10 w 4^j) for
#' j = 0, 1, 2, ..., a ladder of cuts whose spacing grows by the factor 4, up to
#' the upper cut of .os_breaks(). The ladder is needed because the integrand of
#' a moment is the density times e^k = v^(k/a), which for a small a has its own
#' peak far above that of the density: for a = 0.05 the integrand of a second
#' moment of the minimum peaks near 40 v_r, and with cuts at 10 and 40 standard
#' deviations only, the variance of the minimum of 1e6 observations from a
#' beta(0.05, 1) came back 48 percent off against its exact value. On a ladder
#' of ratio 4 a peak of any position is at least a few percent of its piece wide.
#' In the distance e = v^(1/a) the cuts are e_r (1 + k w)^(1/a). They are made
#' for the endpoint that the order statistic is next to (.os_side()), the upper
#' one with 1 - p_r in place of p_r. Returns the rows of .os_cuts(), or NULL when
#' the nearer endpoint has no view.
#' @keywords internal
#' @noRd
.os_end_cuts <- function(r, n, cd, k0) {
  p <- r / (n + 1)
  s <- sqrt(p * (1 - p) / (n + 2))
  side <- .os_side(cd, p)                     # the endpoint the order statistic is next to
  if (side == "") return(NULL)
  lower_side <- side == "lo"
  view <- cd[[side]]
  w <- s / (if (lower_side) p else 1 - p)     # the relative standard deviation in v
  k <- c(-40, -10, 10 * 4^(0:12))
  k <- k[1 + k * w > 0]
  e <- view$dq(p) * exp(log1p(k * w) * view$pw)   # e_r (1 + k w)^(1/a), the distances
  # no cut beyond the ordinary one on the far side of the peak, c_r + h_r or c_r - h_r
  cw <- .os_width(r, n, cd, k0); h <- cw$w + 4 / sqrt(n + 1)
  e <- e[is.finite(e) & e < (if (lower_side) cw$c + h - cd$lower else cd$upper - (cw$c - h))]
  if (!length(e)) return(NULL)
  t <- if (lower_side) cd$lower + e else cd$upper - e
  cbind(t = t,
        lo = if (lower_side) e else if (is.null(cd$lo)) NA_real_ else t - cd$lower,
        hi = if (!lower_side) e else if (is.null(cd$hi)) NA_real_ else cd$upper - t)
}

#' The pieces of an integration range, each with the coordinate it is integrated in
#'
#' Without endpoint views, `breaks` are cut points in standardized units: those
#' inside (lower, upper) are kept and sorted, and every piece is integrated in
#' t, exactly as before the views existed. With a view, `breaks` is the matrix
#' of .os_cuts(), and the points are ordered, and kept or dropped against the
#' limits of the range, by their distance from the endpoint they are nearer to,
#' which is exact where t is rounded. A finite
#' piece is integrated in the distance from an endpoint when that endpoint is no
#' farther from the piece than the piece is long, which is where the singularity
#' is felt, and a piece that is close to both endpoints in this sense is halved
#' first. A finite piece is also integrated in the distance from the nearer
#' endpoint when t does not resolve it: a node t is rounded to eps |t|, and the
#' argument x = median + spread t of the distribution's functions to eps |x|,
#' and when either rounding exceeds 1e-12 of the length of
#' the piece (the length taken from the distances), the piece goes to the
#' endpoint's coordinate, where the distance is the integration variable and
#' keeps its relative precision. Why this exists: with the first rule alone, a
#' piece 1.5e-18 spreads from the endpoint and 1.5e-18 long, one of the cuts of
#' .os_end_cuts() around the rank 39200 of n = 1e6 from a beta(1/14, 1), was
#' integrated in t, where its two limits are one rounding step apart, and the
#' mean came back 2.9 percent low and the variance 14.5 percent low, the mass
#' check passing because its own cuts differ. And with the rounding of t alone
#' in the rule, a piece of the first moment of rank 587180 of n = 1e6 from a
#' beta(1, 1/40) stayed in t, where a thousand points of it map to three values
#' of x, and the variance came back 39 percent low. The limits of the range itself are exact: the lower limit is the
#' support's endpoint (distance 0), or the point that `from` gives as a distance
#' from an endpoint, for the inner integral of a product moment whose outer
#' point is in endpoint coordinates. With from$side "lo" the first piece starts
#' at that distance, and with "hi" the whole range lies above the outer point
#' and is taken in the distance from the upper endpoint.
#' Returns a list of pieces list(side, a, b, A, B): the limits a < b in the
#' piece's coordinate and the limits A < B in t, which the error messages print.
#' @keywords internal
#' @noRd
.os_pieces <- function(lower, upper, breaks, ends, from) {
  if (is.null(ends)) {
    breaks <- sort(breaks[breaks > lower & breaks < upper])
    bounds <- c(lower, breaks, upper)
    return(lapply(seq_len(length(bounds) - 1L), function(i)
      list(side = "t", a = bounds[i], b = bounds[i + 1L], A = bounds[i], B = bounds[i + 1L])))
  }
  has_lo <- !is.null(ends$lo); has_hi <- !is.null(ends$hi)
  from_lo <- !is.null(from) && from$side == "lo"
  from_hi <- !is.null(from) && from$side == "hi"
  # the two limits of the range as points (t, lo, hi)
  lo_pt <- c(t = lower,
             lo = if (!has_lo) NA_real_ else if (from_lo) from$e else lower - ends$lo$at,
             hi = if (!has_hi) NA_real_ else if (from_hi) from$e else ends$hi$at - lower)
  up_pt <- c(t = upper, lo = if (has_lo) upper - ends$lo$at else NA_real_,
             hi = if (has_hi) ends$hi$at - upper else NA_real_)
  # A point is placed by its distance from an endpoint with a view, and not by its t: the
  # distance is exact for the cuts that carry one, while two ways of forming the t of one
  # point (lower + distance, or (x - median)/spread) can differ by a rounding step and put
  # two points in the wrong order. With a view at both endpoints, a point is placed by the
  # distance from the endpoint it is nearer to (`up_side`), and the points nearer to the
  # lower endpoint come first.
  up_side <- function(lo, hi) if (has_lo && has_hi) lo > hi else rep(!has_lo, length(hi))
  above <- function(b, p) {                   # the cuts strictly above the point p
    gb <- up_side(b[, "lo"], b[, "hi"]); gp <- up_side(p[["lo"]], p[["hi"]])
    ifelse(gb, !gp | b[, "hi"] < p[["hi"]], !gp & b[, "lo"] > p[["lo"]])
  }
  below <- function(b, p) {                   # the cuts strictly below the point p
    gb <- up_side(b[, "lo"], b[, "hi"]); gp <- up_side(p[["lo"]], p[["hi"]])
    ifelse(gb, gp & b[, "hi"] > p[["hi"]], gp | b[, "lo"] < p[["lo"]])
  }
  b <- breaks[above(breaks, lo_pt) & below(breaks, up_pt), , drop = FALSE]
  if (nrow(b)) {
    g <- up_side(b[, "lo"], b[, "hi"])
    key <- ifelse(g, -b[, "hi"], b[, "lo"])   # increasing along the range within each side
    o <- order(g, key)
    b <- b[o, , drop = FALSE]
    b <- b[!duplicated(cbind(g, key)[o, , drop = FALSE]), , drop = FALSE]
  }
  pts <- rbind(lo_pt, b, up_pt)
  eps <- .Machine$double.eps
  tres <- 1e-12       # the rounding of a node in t allowed against the length of its piece
  out <- list()
  add <- function(side, a, b, A, B) out[[length(out) + 1L]] <<- list(side = side, a = a, b = b, A = A, B = B)
  for (i in seq_len(nrow(pts) - 1L)) {
    A <- pts[i, "t"]; B <- pts[i + 1L, "t"]
    near_lo <- has_lo && is.finite(B) && !from_hi &&
      ((i == 1L && from_lo) || pts[i, "lo"] <= pts[i + 1L, "lo"] - pts[i, "lo"])
    near_hi <- has_hi && is.finite(A) &&
      (from_hi || pts[i + 1L, "hi"] <= pts[i, "hi"] - pts[i + 1L, "hi"])
    # A piece that t does not resolve goes to the endpoint it is nearer to. A node t of the
    # rule is rounded to eps * |t|, and the argument x = median + spread * t at which the
    # distribution's functions are evaluated is rounded to eps * |x|, which is
    # eps * |t + off| in units of the spread and is far larger than eps * |t| for a
    # distribution whose median is large against its spread (the median of a beta(1, 1/40)
    # is 1 - 9e-13 and its spread 7.5e-6, so off = 1.3e5). The piece is resolved when both
    # roundings are below tres of its length, which is taken from the distances, since
    # B - A may itself be rounded to nothing.
    if (!near_lo && !near_hi && is.finite(A) && is.finite(B)) {
      mt <- max(abs(A), abs(B), abs(A + ends$off), abs(B + ends$off))
      len_lo <- if (has_lo) pts[i + 1L, "lo"] - pts[i, "lo"] else NA_real_
      len_hi <- if (has_hi) pts[i, "hi"] - pts[i + 1L, "hi"] else NA_real_
      fine_lo <- has_lo && !from_hi && len_lo * tres < eps * mt && pts[i + 1L, "lo"] < mt
      fine_hi <- has_hi && len_hi * tres < eps * mt && pts[i, "hi"] < mt
      if (fine_lo && fine_hi) {               # nearer to which endpoint
        if (pts[i + 1L, "lo"] <= pts[i, "hi"]) fine_hi <- FALSE else fine_lo <- FALSE
      }
      near_lo <- fine_lo; near_hi <- fine_hi
    }
    if (near_lo && near_hi) {                  # close to both endpoints: halve the piece
      M <- (A + B) / 2
      add("lo", pts[i, "lo"], M - ends$lo$at, A, M)
      add("hi", pts[i + 1L, "hi"], ends$hi$at - M, M, B)
    } else if (near_lo) {
      add("lo", pts[i, "lo"], pts[i + 1L, "lo"], A, B)
    } else if (near_hi) {
      add("hi", pts[i + 1L, "hi"], pts[i, "hi"], A, B)
    } else {
      add("t", A, B, A, B)
    }
  }
  out
}

#' The cut points of rank r
#'
#' Why this exists: the density of X_(r) concentrates around c_r = Q(p_r),
#' p_r = r/(n+1), in a region whose width shrinks like n^{-1/2} for a central rank. For
#' an extreme rank the rate depends on the distribution (1/n for the uniform minimum,
#' no shrinking for the exponential maximum). Meanwhile the
#' adaptive rule starts from a fixed set of nodes on the whole (transformed)
#' range. When no initial node lands on the spike, the rule returns a value near
#' zero with a small claimed error (observed with the reference engine for
#' the pair (1, 251) at n = 1001). Splitting the range at c_r - h_r and c_r + h_r puts the spike
#' inside one finite piece whose endpoints the rule sees, and leaves monotone
#' tail pieces. The half-width is the reference engine's,
#'   h_r = 10 s_r / f(c_r) + 4 / sqrt(n + 1),
#' in the standardized units of .os_centered(), with s_r = sqrt{p_r (1 - p_r)/(n + 2)}
#' the standard deviation of the uniform order statistic U_(r), so that
#' s_r / f(c_r) is the large-sample standard deviation of X_(r) by the delta
#' method, generous by the factor 10 for the skewed densities of the extreme
#' ranks, with the minimum half-width 4/sqrt(n + 1) spreads of the distribution, which is
#' 4/sqrt(n + 1) in the normal's own units. The delta-method half-width
#' w_r = k s_r / f(c_r), returned by .os_width(), also sets the scale of the
#' absolute tolerances of the integrals of rank r (see .os_single()). A split at
#' c_r +/- w_r as well, when w_r is below the minimum half-width, was tried and dropped: it
#' drives the bisection onto an endpoint where the density is infinite (the
#' beta(1/2, 1/2)). The scaled tolerances recover the sharp extremes of a density that is
#' finite at its endpoints, and next to an endpoint where the density is infinite the
#' cuts of .os_end_cuts() are needed as well.
#' A rule in probability units, Q(p_r +/- 10 s_r), was
#' tried first and rejected: for an extreme rank it gives at most one break, so
#' the spike sits in a semi-infinite piece, and the rule then converges to a
#' value off by 1e-10 (the normal at n = 97) to 1e-7 (the logistic at n = 9),
#' which a tighter tolerance does not repair. The breaks are in centered units,
#' and a break outside the support is dropped. If f(c_r) is zero or not finite,
#' only the minimum half-width is used.
#'
#' The points that the distribution object declares in `breaks`, where its
#' density has a kink or a jump, are returned with the cut points of every
#' rank, because the rule assumes a smooth integrand inside a piece. Every
#' integral takes its cuts from here (the single integrals, and the inner and
#' the outer integral of a product moment), so every piece is smooth. Without
#' such a cut, the sum of the covariance matrix of the Laplace distribution at
#' n = 5 missed 2n by 1.7e-8 with an error estimate of 6e-10, the pair (1, 2)
#' being 4e-9 off, while the pairs with the median, whose own cut at c_3 = 0
#' falls on the kink, were right to 1e-12. `kinks = FALSE` leaves the declared points out,
#' for the callers that take them from the rows `kink_rows` of .os_centered() with their
#' own distances from an endpoint (the integrals with an endpoint view), and for the probe
#' of .os_tail_probe(), which evaluates them in the distribution's own units.
#' @keywords internal
#' @noRd
.os_breaks <- function(r, n, cd, k, kinks = TRUE) {
  cw <- .os_width(r, n, cd, k)
  h <- cw$w + 4 / sqrt(n + 1)               # the minimum half-width, in the standardized units
  b <- c(cw$c - h, cw$c + h)
  if (kinks) b <- unique(c(b, cd$kinks))    # and the points where the density is not smooth
  b[b > cd$lower & b < cd$upper]
}

#' The location c_r and the delta-method half-width w_r = k s_r / f(c_r) of rank r
#'
#' In standardized units; w_r is 0 where f(c_r) is zero or not finite. Besides
#' the split of .os_breaks(), w_r sets the scale of the absolute tolerances of
#' the integrals of rank r through .os_tolscale(): the moments of a sharp order
#' statistic are small (the variance of the minimum of 10^6 uniform
#' observations is 7e-12 in these units), and the fixed absolute tolerances of
#' the reference engine, set for the normal's scale one, were met before
#' such a moment was resolved (the minimum of 10^6 uniform or exponential
#' observations by quadrature was off by 1e-4 to 1e-3). Each absolute tolerance
#' and acceptance threshold is multiplied by min(1, w_r)^pow for a single
#' integral, min(1, w_s) for the inner integral of a product moment and
#' min(1, w_r) min(1, w_s) for its outer integral, so they are never looser than
#' the engine's constants and tighten with the order statistic; the normal, whose
#' w_r is of order one at n in the thousands, moves by 2e-14.
#' @keywords internal
#' @noRd
.os_width <- function(r, n, cd, k) {
  p <- r / (n + 1)
  s <- sqrt(p * (1 - p) / (n + 2))
  c0 <- cd$q(p)
  # the log-density at c_r, taken in the distance from an endpoint where the density is
  # infinite: an order statistic within the rounding of the centered coordinate of such an
  # endpoint has c0 equal to the endpoint itself, where the density is not finite. With a
  # view at both endpoints it is the view of the nearer one (.os_side()): from the far
  # endpoint the argument rounds to that endpoint's opposite, and the width came out as 0.
  side <- .os_side(cd, p, nearest = FALSE)
  ldc <- if (side == "") cd$ld(c0) else cd[[side]]$ld(cd[[side]]$dq(p))
  w <- k * s / exp(ldc)                     # k large-sample standard deviations of X_(r)
  list(c = c0, w = if (is.finite(w)) w else 0)
}

#' @keywords internal
#' @noRd
.os_tolscale <- function(w) if (w > 0) min(1, w) else 1

#' integrate() over the pieces cut by the breaks, values and estimates summed
#'
#' Pieces entirely in a tail are monotone-decaying and safe for the adaptive
#' rule. A piece whose convergence message is not "OK" (for example "roundoff
#' error was detected" on a near-zero tail piece, where the relative tolerance is
#' unattainable) is accepted if its reported error estimate is below
#' ctrl$accept.err, and stops the computation otherwise, naming the piece. A
#' piece reported as "probably divergent" is treated by where it lies. On a
#' finite piece the flag is spurious for a probability density, since any power
#' of t times an integrable density is integrable on a bounded interval, and it
#' arises at an integrable singularity of the density (the beta(1/2, 1/2) at its
#' endpoints), so the piece is accepted on its error estimate like any other.
#' On an infinite piece it is accepted only when its value is below
#' ctrl$accept.err as well: a negligible tail piece of a convergent integral can
#' carry the flag (a t parent with 1.5 degrees of freedom gave one of size
#' 1e-15), while a divergent second moment returns a finite value of ordinary
#' size with a small error estimate (a t parent with 1.9 degrees of freedom
#' declared without its tail index gave -11.3 on each tail piece). The error
#' names the piece and the usual causes.
#' @keywords internal
#' @noRd
.os_integrate_split <- function(f, lower, upper, breaks, rel.tol, abs.tol, ctrl, ctx,
                                accept = ctrl$accept.err, ends = NULL, from = NULL) {
  val <- 0; err <- 0
  for (pc in .os_pieces(lower, upper, breaks, ends, from)) {
    if (pc$side == "t") {
      g <- f; ga <- pc$a; gb <- pc$b
    } else {                                   # the distance e = v^pw from an endpoint
      fe <- ends[[pc$side]]$f; pw <- ends[[pc$side]]$pw
      g <- function(v) fe(v^pw) * pw * v^(pw - 1)
      ga <- max(pc$a, 0)^(1 / pw); gb <- max(pc$b, 0)^(1 / pw)
    }
    bounds_i <- unname(c(pc$A, pc$B))
    out <- tryCatch(integrate(g, ga, gb, rel.tol = rel.tol, abs.tol = abs.tol,
                              subdivisions = ctrl$subdivisions, stop.on.error = FALSE),
                    error = function(e)
      stop("os_moments(): the integral failed (", ctx, "), piece [", signif(bounds_i[1], 4), ", ",
           signif(bounds_i[2], 4), "] in standardized units: ", conditionMessage(e),
           ". A non-finite value of the integrand means that the density is infinite or ",
           "undefined at a point the quadrature reached: an order statistic within about ",
           "1e-13 of a finite endpoint of the support, relative to the spread, cannot be ",
           "placed in double precision", call. = FALSE))
    divergent <- out$message == "the integral is probably divergent"
    infinite_piece <- !is.finite(bounds_i[1]) || !is.finite(bounds_i[2])
    if (out$message != "OK" &&
        (out$abs.error > accept ||
         (divergent && infinite_piece && abs(out$value) > accept)))
      stop("os_moments(): the integral failed (", ctx, "), piece [",
           signif(bounds_i[1], 4), ", ", signif(bounds_i[2], 4), "] in standardized units: ",
           out$message, " (error estimate ", signif(out$abs.error, 3), "). ",
           if (infinite_piece)
             paste0("On a tail piece this usually means that the moment does not exist ",
                    "(declare the tail indices in os_dist()), or that the tail is so heavy ",
                    "that the mass of the integrand lies far out, as for a lognormal with ",
                    "log-SD 3, which the quadrature does not resolve")
           else
             paste0("On a bounded piece this points at a density that is zero, infinite ",
                    "or not smooth inside the piece (a kink or a jump of the density is ",
                    "declared in the 'breaks' argument of os_dist())"),
           call. = FALSE)
    val <- val + out$value; err <- err + out$abs.error
  }
  list(value = val, abs.error = err)
}

#' log{F(z) - F(y)} for y <= z without cancellation
#'
#' Why this exists: the joint density carries {F(z) - F(y)}^(s-r-1), and the
#' log-scale integrand takes (s-r-1) log of the difference, which therefore
#' needs relative precision even when y and z sit deep in the same tail, where
#' F(y) and F(z) are both near 0 or both near 1. The difference is formed on
#' the side where the numbers are small: below the median (z <= 0 in centered
#' units) as log F(z) + log{1 - F(y)/F(z)}, above it (y >= 0) as
#' log S(y) + log{1 - S(z)/S(y)} with S = 1 - F, and directly when the pair
#' straddles the median, where the difference is of ordinary size. A pair whose
#' two log-tails both underflow to -Inf gives NaN, which is set to -Inf: the
#' difference is zero in doubles. The ratio of the two tail probabilities is
#' capped at 1: for z within a rounding step of y it can come out just above 1,
#' and log1p() of a number below -1 would make R warn "NaNs produced" for a
#' value that is -Inf either way (the covariance of the ranks 999991 and 999999
#' of a million from a beta(1, 1/3) returned its value with that warning). y is a scalar (the outer variable) and z a
#' vector (the inner nodes). The caller may pass the values it already holds,
#' log F(y), log S(y), F(y) and log S(z), so that each is evaluated once per
#' node; the arithmetic is the same either way. When z is in the coordinate of
#' an endpoint view (.os_end_view()), `ev` is that view and `zt` holds the same
#' points in t, which only decide the side of the median, and the caller passes
#' the values at y from the view of y's own coordinate.
#' @keywords internal
#' @noRd
.os_logdiff <- function(y, z, cd, lFy = cd$lF(y), lSy = cd$lS(y), pfy = cd$pf(y),
                        lSz = NULL, ev = cd, zt = z) {
  out <- numeric(length(z))
  a <- zt <= 0
  b <- !a & (y >= 0)
  m <- !a & !b
  if (any(a)) { lFz <- ev$lF(z[a]); out[a] <- lFz + log1p(-pmin(1, exp(lFy - lFz))) }
  if (any(b)) out[b] <- lSy + log1p(-pmin(1, exp((if (is.null(lSz)) ev$lS(z[b]) else lSz[b]) - lSy)))
  if (any(m)) out[m] <- log(ev$pf(z[m]) - pfy)
  out[is.nan(out)] <- -Inf
  out
}

#' k * log(v), returning 0 when k = 0 (guards 0 * -Inf = NaN)
#' @keywords internal
#' @noRd
.os_xlog <- function(k, logv) if (k == 0) 0 else k * logv

#' E[(X_(r) - shift)^pow] in the engine's units by one exact integral
#'
#' Density: f_{r:n}(t) = C F(t)^(r-1) {1 - F(t)}^(n-r) f(t) with
#' C = n!/{(r-1)! (n-r)!}, assembled on the logarithmic scale: lbeta for the
#' coefficient (the difference of three lgamma values lost 1e-9 of a central
#' rank's covariance at n = 1e6), the log-tails for the powers, the log-density,
#' exponentiated once. With pow = 0 the integral is the mass of the density,
#' which os_moments() checks against 1.
#' The moment is taken about `shift`, which os_moments() sets to the order
#' statistic's own center c_r = Q(r/(n+1)): with the moments taken
#' about the median, a variance was the small difference of two moments of size
#' (c_r - c)^2, and its relative error grew like that size over the variance,
#' n^2 for the uniform minimum (7 digits left at n = 1000, none at
#' n = 1e5). About c_r the first moment is small and nothing cancels. For pow = 1,
#' the only power at which the integrand changes sign, the shift is added to the
#' cut points, so that (t - shift) times the density is one-signed on each
#' piece: a piece straddling the sign change has a near-zero value that the
#' relative tolerance cannot resolve, and integrate() then flags roundoff. For
#' pow = 0 and 2 the split would only drive the bisection onto a singular
#' endpoint. The absolute tolerance and the acceptance threshold are scaled by
#' the order statistic's width to the power pow (.os_width()). `u` is the
#' probability whose quantile the shift is, when it is one (r/(n+1) for c_r):
#' the endpoint coordinates then take the distance of the shift from an endpoint
#' from the quantile function directly (.os_shift_dist()).
#' @keywords internal
#' @noRd
.os_single <- function(r, n, pow, cd, ctrl, shift = 0, u = NULL) {
  lC <- -lbeta(r, n - r + 1)                 # log n!/{(r-1)! (n-r)!}, without the cancellation
                                             # of three large lgamma values at n = 1e6
  f <- function(t) {
    lw <- lC + .os_xlog(r - 1, cd$lF(t)) + .os_xlog(n - r, cd$lS(t)) + cd$ld(t)
    (t - shift)^pow * exp(lw)
  }
  sc <- .os_tolscale(.os_width(r, n, cd, ctrl$k)$w)^pow   # the moment's own scale, at most 1
  br <- .os_breaks(r, n, cd, ctrl$k)
  if (pow == 1) br <- c(br, shift)
  # the same integrand in the distance from an endpoint where the density is infinite
  # (NULL for every distribution whose density is finite at its endpoints)
  ds <- .os_shift_dist(cd, shift, u)
  ends <- .os_ends(cd, function(ev, dev, side) function(e) {
    lw <- lC + .os_xlog(r - 1, ev$lF(e)) + .os_xlog(n - r, ev$lS(e)) + ev$ld(e)
    dev(e)^pow * exp(lw)
  }, ds[["lo"]], ds[["hi"]])
  # The cuts with their distances: the cut at the shift is placed by its own distance from
  # the endpoint, and so are the declared points. When the endpoint that the order
  # statistic is next to has a view (.os_side(), as in .os_end_cuts()),
  # the cut at the shift is made for every power, so that the mass (pow = 0), the first and
  # the second moment are integrated over the same pieces in the same coordinates. The
  # mass check of os_moments() then guards the pieces of the moments themselves: with the
  # cut made for the first moment only, a piece that was integrated in a coordinate that
  # did not resolve it gave a wrong mean or variance while the mass, cut elsewhere, came
  # out as 1 (a beta(1/14, 1) and a beta(1, 1/40) at n = 1e6, 2.9 and 39 percent off
  # without an error). With the same cuts the mass of those cases is 3.7 and 1.2, and the
  # computation stops. Next to an endpoint without a view the cut stays with the first
  # moment: next to a singular endpoint that is integrated in t, it would add a piece that
  # the rule flags for the mass alone (a user's beta(1/2, 1/2), the maximum of 100). The
  # same cuts do not certify the moments: a mass of 1 was returned with a variance 28
  # percent off when the cuts were made from the wrong endpoint (.os_side()).
  if (!is.null(ends)) {
    own_view <- .os_side(cd, r / (n + 1)) != ""
    br <- rbind(.os_cuts(cd, .os_breaks(r, n, cd, ctrl$k, kinks = FALSE)),
                if (pow == 1 || own_view) .os_cuts(cd, shift, ds),
                .os_end_cuts(r, n, cd, ctrl$k), cd$kink_rows)
  }
  .os_integrate_split(f, cd$lower, cd$upper, br, ctrl$rel.tol, ctrl$abs.tol.single * sc, ctrl,
                      ctx = paste0("E[X^", pow, "], r = ", r, ", n = ", n),
                      accept = ctrl$accept.err * sc, ends = ends)
}

#' E[(X_(r) - sr)(X_(s) - ss)] in the engine's units, r < s, by nested exact quadrature
#'
#' Joint density for y < z:
#'   f(y, z) = C F(y)^(r-1) {F(z) - F(y)}^(s-r-1) {1 - F(z)}^(n-s) f(y) f(z),
#'   C = n!/{(r-1)! (s-r-1)! (n-s)!}.
#' The outer integral runs over y in the support and the inner over z in
#' (y, upper), split at the breaks of rank s. The outer range is split at the
#' breaks of both ranks: at those of rank r for the spike of y's own density, and
#' at those of rank s because the inner quadrature changes its piece structure
#' where y crosses them, a feature of the outer integrand that the adaptive rule
#' otherwise does not see. With the breaks of rank r alone, the pair (9, 20) at
#' n = 70 of the normal converged to a value 5.5e-10 off with an error estimate
#' of 1e-12, which the permutation identity caught and an integration in the
#' probability scale confirmed; with both sets of breaks the value is recovered
#' to 3e-14. The middle factor uses .os_logdiff() and is skipped when s = r + 1
#' (exponent 0). The inner error estimate is not carried into the outer one: the
#' returned estimate is the outer integral's own. The tails at the outer point y
#' are evaluated once per inner integral, and log S(z) once per inner node, which
#' saves about a tenth of the time and leaves every value unchanged. The moment
#' is taken about the two centers sr and ss of the order statistics (see
#' .os_single()), which are added to the cut points of their ranks, so
#' that the covariance E[(Y - c_r)(Z - c_s)] - m_r m_s is formed from small
#' quantities. The inner absolute tolerance is divided by max(1, |y - sr|),
#' the outer factor the inner value is multiplied by: at a fixed 1e-15 a
#' lognormal with log-SD 2, whose outer nodes reach y = 1e5, lost 3.5e-6 of a
#' covariance of 38.8 at n = 2, and the scaling leaves 5e-8 at no measurable
#' cost (the outer nodes of the ordinary families have |y| below 10). The
#' inner absolute tolerance and acceptance threshold are scaled by the width of
#' rank s, the outer ones by the widths of both ranks (.os_width()). `ur` and
#' `us` are the probabilities whose quantiles sr and ss are. When the density is
#' infinite at a finite endpoint, the pieces next to that endpoint, of the outer
#' and of the inner integral, are taken in the distance from it (.os_end_view()),
#' and the outer point is handed to the inner integral in that coordinate, so
#' that F(z) - F(y) and the deviations keep their precision next to the endpoint.
#' @keywords internal
#' @noRd
.os_prod <- function(r, s, n, cd, ctrl, sr = 0, ss = 0, ur = NULL, us = NULL) {
  lC <- -lbeta(r, n - r + 1) - lbeta(s - r, n - s + 1)   # log n!/{(r-1)! (s-r-1)! (n-s)!}
  m <- s - r - 1
  br_r <- c(.os_breaks(r, n, cd, ctrl$k), sr)
  br_s <- c(.os_breaks(s, n, cd, ctrl$k), ss)
  sc_r <- .os_tolscale(.os_width(r, n, cd, ctrl$k)$w)   # the two order statistics' widths
  sc_s <- .os_tolscale(.os_width(s, n, cd, ctrl$k)$w)
  ctx_in <- paste0("inner, r = ", r, ", s = ", s, ", n = ", n)
  inner <- function(y) {                      # y scalar: int_y^upper z f(y, z) dz
    lFy <- cd$lF(y); lSy <- cd$lS(y); pfy <- cd$pf(y)
    ly <- lC + .os_xlog(r - 1, lFy) + cd$ld(y)
    g <- function(z) {
      lSz <- cd$lS(z)
      lmid <- if (m == 0) 0 else m * .os_logdiff(y, z, cd, lFy, lSy, pfy, lSz)
      (z - ss) * exp(ly + lmid + .os_xlog(n - s, lSz) + cd$ld(z))
    }
    .os_integrate_split(g, y, cd$upper, br_s, ctrl$rel.tol,
                        ctrl$abs.tol.inner * sc_s / max(1, abs(y - sr)), ctrl, ctx = ctx_in,
                        accept = ctrl$accept.err * sc_s)$value
  }
  f_outer <- function(y) (y - sr) * vapply(y, inner, numeric(1))
  if (is.null(cd$lo) && is.null(cd$hi))       # the density is finite at both endpoints
    return(.os_integrate_split(f_outer, cd$lower, cd$upper, sort(unique(c(br_r, br_s))),
                               ctrl$rel.tol.outer, ctrl$abs.tol.outer * sc_r * sc_s, ctrl,
                               ctx = paste0("outer, r = ", r, ", s = ", s, ", n = ", n),
                               accept = ctrl$accept.err * sc_r * sc_s))

  # ---- with an endpoint view: the same double integral, each point in the coordinate ----
  # of its piece. A point is a pair (value, side): t in standardized units for side "t",
  # the distance from the lower or the upper endpoint for "lo" and "hi".
  dr <- .os_shift_dist(cd, sr, ur); ds <- .os_shift_dist(cd, ss, us)
  view <- list(t = cd, lo = cd$lo, hi = cd$hi)
  tval <- list(t = function(x) x, lo = function(e) cd$lower + e, hi = function(e) cd$upper - e)
  dev_r <- list(t = function(x) x - sr, lo = function(e) e - dr[["lo"]], hi = function(e) dr[["hi"]] - e)
  dev_s <- list(t = function(x) x - ss, lo = function(e) e - ds[["lo"]], hi = function(e) ds[["hi"]] - e)
  cuts_r <- rbind(.os_cuts(cd, .os_breaks(r, n, cd, ctrl$k, kinks = FALSE)), .os_cuts(cd, sr, dr),
                  .os_end_cuts(r, n, cd, ctrl$k), cd$kink_rows)
  cuts_s <- rbind(.os_cuts(cd, .os_breaks(s, n, cd, ctrl$k, kinks = FALSE)), .os_cuts(cd, ss, ds),
                  .os_end_cuts(s, n, cd, ctrl$k), cd$kink_rows)
  inner_at <- function(y, side) {             # the inner integral from the point (y, side)
    ev <- view[[side]]
    lFy <- ev$lF(y); lSy <- ev$lS(y); pfy <- ev$pf(y)
    yt <- tval[[side]](y)
    ly <- lC + .os_xlog(r - 1, lFy) + ev$ld(y)
    make <- function(evz, dev, zside) function(z) {
      lSz <- evz$lS(z)
      lmid <- if (m == 0) 0 else
        m * .os_logdiff(yt, z, cd, lFy, lSy, pfy, lSz, ev = evz, zt = tval[[zside]](z))
      dev(z) * exp(ly + lmid + .os_xlog(n - s, lSz) + evz$ld(z))
    }
    .os_integrate_split(make(cd, dev_s$t, "t"), yt, cd$upper, cuts_s, ctrl$rel.tol,
                        ctrl$abs.tol.inner * sc_s / max(1, abs(dev_r[[side]](y))), ctrl,
                        ctx = ctx_in, accept = ctrl$accept.err * sc_s,
                        ends = .os_ends(cd, make, ds[["lo"]], ds[["hi"]]),
                        from = if (side == "t") NULL else list(side = side, e = y))$value
  }
  outer_make <- function(ev, dev, side) function(y)
    dev(y) * vapply(y, function(p) inner_at(p, side), numeric(1))
  .os_integrate_split(outer_make(cd, dev_r$t, "t"), cd$lower, cd$upper,
                      rbind(cuts_r, cuts_s), ctrl$rel.tol.outer,
                      ctrl$abs.tol.outer * sc_r * sc_s, ctrl,
                      ctx = paste0("outer, r = ", r, ", s = ", s, ", n = ", n),
                      accept = ctrl$accept.err * sc_r * sc_s,
                      ends = .os_ends(cd, outer_make, dr[["lo"]], dr[["hi"]]))
}
