# ==============================================================================
# File: os_dist.R
# Purpose: The distribution object of orderStatMoments: the public constructor
#          os_dist() for a user-supplied distribution, the internal constructor
#          of the thirteen built-in families from R's own parameters, the
#          validator, and the internal Gumbel functions (maximum form) with
#          accurate log-tails.
#
# EXPORTED FUNCTIONS:
#   - os_dist()           Build a distribution object from (p, d, q)
#
# INTERNAL HELPERS:
#   - .os_family_spec     The built-in families with R's parameters and defaults
#   - .os_family()        A family and its parameters -> the standard member and the
#                         shift and multiplier that its parameters give
#   - .os_validate_dist() Check the functions, the support, the tails, the center
#   - .os_resolve_dist()  A family name with parameters, or an os_dist object, resolved
#   - .qf_beta()          The F quantile function through the beta quantile of each tail
#   - .pgumbel(), .dgumbel(), .qgumbel()  Gumbel (max form), F(x) = exp(-e^{-x})

#################### EXPORTED MAIN FUNCTIONS ####################

#' Define a distribution for the moments of its order statistics
#'
#' Builds the distribution object that \code{\link{os_moments}} takes in its
#' \code{dist} argument, from the distribution function, the density and the
#' quantile function of a continuous distribution that is not among the
#' built-in families. The three functions define the distribution itself, with
#' its parameters fixed inside them. \code{os_moments()} computes the moments
#' of the order statistics of this distribution, and its arguments
#' \code{os_shift} and \code{os_multiplier} give those of \eqn{a + bX}.
#'
#' The thirteen built-in families need no call of \code{os_dist()}: the normal,
#' uniform, exponential, logistic, Gumbel, Student \eqn{t}, Cauchy, gamma,
#' chi-square, beta, lognormal, Weibull and \eqn{F} distributions are given to
#' \code{os_moments()} by name, with R's own parameters, and the package sets
#' their supports and tail indices.
#'
#' @param p The distribution function, with the signature of R's
#'   \code{p*} functions: \code{p(q, lower.tail = TRUE, log.p = FALSE)},
#'   vectorized in \code{q}. The \code{log.p} and \code{lower.tail} arguments
#'   must be honored, because the integrands are evaluated on the log scale in
#'   both tails.
#' @param d The density, \code{d(x, log = FALSE)}, vectorized in \code{x}, with
#'   \code{log} honored.
#' @param q The quantile function, \code{q(p)}, vectorized in \code{p}.
#' @param support Numeric vector of length two, the endpoints of the support,
#'   finite or infinite. They are the limits of integration.
#' @param symmetric Logical. \code{TRUE} if the density is symmetric about its
#'   median, that is \eqn{Q(1 - u) = 2 Q(1/2) - Q(u)}. A symmetric distribution
#'   lets \code{os_moments()} halve the number of double integrals of a full
#'   matrix and check the reflection symmetry of its result.
#' @param tail Numeric vector of length two, the tail indices
#'   \eqn{(\nu_L, \nu_U)}: the density decays like \eqn{|x|^{-\nu_L - 1}} in the
#'   lower tail and like \eqn{x^{-\nu_U - 1}} in the upper tail. \code{Inf}
#'   for a tail that decays faster than any power, or for a bounded tail. The
#'   indices decide which moments exist (see Details).
#' @param name A label for printing.
#' @param breaks Numeric vector, or \code{NULL} for none: the points strictly
#'   inside the support at which the density is not smooth, that is, where it
#'   has a kink or a jump (the median of the Laplace distribution, the mode of
#'   a triangular distribution, the knots of a piecewise density). The density
#'   must be finite there, and it may be zero. \code{os_moments()} cuts every integral at these
#'   points (see Details). For a symmetric distribution the points must be
#'   symmetric about the median.
#'
#' @return An object of class \code{"os_dist"}: a list with the components
#'   \code{name}, \code{p}, \code{d}, \code{q}, \code{support},
#'   \code{symmetric}, \code{tail}, \code{center} (the median \eqn{Q(1/2)}),
#'   \code{closed} (the closed-form moments the distribution has, none for a
#'   user's distribution) and \code{breaks} (sorted, \code{NULL} when there are
#'   none).
#'
#' @details
#' \strong{The moment condition.} For a distribution whose density decays like
#' \eqn{|x|^{-\nu - 1}}, the density of \eqn{X_{(r)}} decays like
#' \eqn{|x|^{-\nu r - 1}} in the lower tail and like
#' \eqn{x^{-\nu (n - r + 1) - 1}} in the upper tail, so
#' \eqn{E|X_{(r)}|^k < \infty} if and only if \eqn{\nu_L r > k} and
#' \eqn{\nu_U (n - r + 1) > k}. \code{os_moments()} applies this condition with
#' \eqn{k = 2} to every requested rank and refuses a moment that does not exist.
#' The package cannot check the declared indices. The default
#' \code{tail = c(Inf, Inf)} declares that every moment exists, so the indices
#' of a heavy-tailed distribution must be given.
#'
#' \strong{Accuracy requirement.} The log-tails \code{p(x, log.p = TRUE)} and
#' \code{p(x, lower.tail = FALSE, log.p = TRUE)} must keep relative precision
#' deep in both tails, as R's base distributions do. A function that computes
#' \code{log(p(x))} from an underflowing probability loses the tail, and the
#' integrals with it. \code{os_dist()} checks that the three functions honor
#' their arguments at the median and are vectorized, and that \code{d}
#' integrates to \eqn{1/2} between the quartiles, where \code{p} rises by
#' \eqn{1/2}. \code{os_moments()} checks that the log-tails and the log-density
#' are finite at the median and at the cut points of the smallest and the
#' largest requested rank, and that the density of every order statistic it
#' integrates has mass 1. None of these certifies accuracy.
#'
#' \strong{A distribution far from the origin.} The integrals are taken about
#' the median of the supplied distribution, in units of its spread (its
#' interquartile range divided by that of the standard normal), so the scale on
#' which the distribution is supplied does not matter: the package's tests
#' check a normal supplied with standard deviations from \eqn{10^{-3}} to
#' \eqn{10^6}. Its position does: a distribution whose median is large against
#' its spread, say \eqn{10^8} with unit spread, loses the digits of
#' \eqn{x - Q(1/2)} in double precision and the computation stops with an
#' error. Such a distribution is supplied near the origin, and the shift is
#' passed to \code{os_moments()} as \code{os_shift}, which is applied exactly.
#' This also keeps \code{os_blue()} exact.
#'
#' \strong{A density that is infinite at an endpoint.} When the density is
#' infinite at a finite lower endpoint of the support, as for a power or a
#' gamma-type density with a shape below 1, \code{os_moments()} integrates next
#' to that endpoint in the distance from it, raised to a power that it
#' estimates from the distribution function. The distance is represented
#' without rounding when the endpoint is 0, so such a distribution is
#' supplied with its endpoint at 0. An infinite density at an upper endpoint
#' has no such treatment for a user's distribution, and the order statistics
#' next to it are computed less accurately or refused. Supplying the reflected
#' variable, whose singular endpoint is the lower one, brings that endpoint
#' under the treatment above.
#'
#' \strong{A density that is not smooth inside the support.} The quadrature
#' rule assumes an integrand that is smooth inside each piece of the
#' integration range. A point where the density has a kink or a jump, such as
#' the median of the Laplace distribution, is therefore declared in
#' \code{breaks}, and \code{os_moments()} then cuts every single integral and
#' both levels of every double integral there, so that each piece is smooth.
#' Without the declaration the covariances lose accuracy: for the Laplace
#' example below, at \eqn{n = 5}, the sum of the covariance matrix misses its
#' exact value by about \eqn{10^{-8}}, and with \code{breaks = 0} by less than
#' \eqn{10^{-12}}. The means and the variances are affected less. The error
#' estimate \code{quad_err} does not show the loss, and the residual
#' \code{sum_B} of the run-time checks, which is available for a complete rank
#' set, does. The package cannot find such points itself. The built-in families
#' have none.
#'
#' @examples
#' # A Pareto distribution with shape 3 on [1, Inf): the upper tail index is 3,
#' # so the variance of the maximum exists (3 * 1 > 2)
#' pareto <- function(shape) os_dist(
#'   p = function(q, lower.tail = TRUE, log.p = FALSE) {
#'     ls <- -shape * log(pmax(q, 1))                  # log S(x) = -shape * log(x)
#'     if (lower.tail) { if (log.p) log(-expm1(ls)) else -expm1(ls) }
#'     else if (log.p) ls else exp(ls)
#'   },
#'   d = function(x, log = FALSE) {
#'     ld <- ifelse(x >= 1, log(shape) - (shape + 1) * log(x), -Inf)
#'     if (log) ld else exp(ld)
#'   },
#'   q = function(p) (1 - p)^(-1 / shape),
#'   support = c(1, Inf), tail = c(Inf, shape), name = "pareto")
#' os_moments(5, dist = pareto(3), cov = FALSE)
#'
#' # With shape 1.5 the maximum has no variance, and the request is refused
#' try(os_moments(5, ranks = 5, dist = pareto(1.5)))
#' os_moments(5, ranks = 1:4, dist = pareto(1.5), cov = FALSE)$var
#'
#' # The Laplace (double exponential) distribution, symmetric about 0, where its
#' # density is not smooth
#' laplace <- os_dist(
#'   p = function(q, lower.tail = TRUE, log.p = FALSE) {
#'     s <- if (lower.tail) q else -q                  # F(x) = S(-x) by symmetry
#'     lp <- ifelse(s <= 0, log(0.5) + s, log1p(-0.5 * exp(-abs(s))))
#'     if (log.p) lp else exp(lp)
#'   },
#'   d = function(x, log = FALSE) if (log) log(0.5) - abs(x) else 0.5 * exp(-abs(x)),
#'   q = function(p) ifelse(p <= 0.5, log(2 * p), -log(2 * (1 - p))),
#'   symmetric = TRUE, name = "laplace",
#'   breaks = 0)                                       # the density has a kink at 0
#' os_moments(5, dist = laplace, os_shift = 10, os_multiplier = 2)
#' @seealso \code{\link{os_moments}}
#' @export
os_dist <- function(p, d, q, support = c(-Inf, Inf), symmetric = FALSE,
                    tail = c(Inf, Inf), name = "user", breaks = NULL) {
  obj <- structure(list(name = name, p = p, d = d, q = q, support = support,
                        symmetric = symmetric, tail = tail, center = NA_real_,
                        closed = list(mean = NULL, var = NULL, cov = NULL),
                        breaks = breaks),
                   class = "os_dist")
  .os_validate_dist(obj)
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' The built-in families: R's parameters and what each one does
#'
#' One entry per family: `par`, the parameters under R's own names with R's
#' defaults (NA for a parameter that R gives no default), and `label`, the name
#' the messages print. The families are the continuous distributions of the
#' stats package and the Gumbel. .os_family() turns a family and its
#' parameters into the standard member, whose moments the engine computes, and
#' the shift a0 and the multiplier b0 > 0 that carry the standard member Z to
#' the distribution asked for, X = a0 + b0 Z.
#' @keywords internal
#' @noRd
.os_family_spec <- list(
  norm    = list(par = c(mean = 0, sd = 1),            label = "normal"),
  unif    = list(par = c(min = 0, max = 1),            label = "uniform"),
  exp     = list(par = c(rate = 1),                    label = "exponential"),
  logis   = list(par = c(location = 0, scale = 1),     label = "logistic"),
  gumbel  = list(par = c(location = 0, scale = 1),     label = "Gumbel"),
  t       = list(par = c(df = NA),                     label = "Student t"),
  cauchy  = list(par = c(location = 0, scale = 1),     label = "Cauchy"),
  gamma   = list(par = c(shape = NA, rate = 1, scale = NA), label = "gamma"),
  chisq   = list(par = c(df = NA),                     label = "chi-square"),
  beta    = list(par = c(shape1 = NA, shape2 = NA),    label = "beta"),
  lnorm   = list(par = c(meanlog = 0, sdlog = 1),      label = "lognormal"),
  weibull = list(par = c(shape = NA, scale = 1),       label = "Weibull"),
  f       = list(par = c(df1 = NA, df2 = NA),          label = "F"))

#' A built-in family with its parameters, as the standard member and an affine map
#'
#' `args` is the list of the family's parameters that the user passed to
#' os_moments() by R's own names. Returns a list with `dist`, the os_dist object
#' of the standard member; `a0` and `b0`, the shift and the multiplier with
#' X = a0 + b0 Z for the standard member Z; and `par`, the parameters as
#' resolved, for printing.
#'
#' The standard member carries the shape parameters only. A location or a
#' scale parameter of R (the mean and sd of the normal, the rate of the
#' exponential, the scale of the gamma or the Weibull, exp(meanlog) of the
#' lognormal) is applied after the integration, exactly, so that a
#' distribution with a large mean or a small rate loses no digits:
#'   norm     N(0, 1);                           a0 = mean, b0 = sd
#'   unif     uniform on [0, 1];                 a0 = min,  b0 = max - min
#'   exp      unit rate;                         b0 = 1/rate
#'   logis    standard logistic;                 a0 = location, b0 = scale
#'   gumbel   F(x) = exp(-exp(-x));              a0 = location, b0 = scale
#'   t        Student t with df;                 tail indices (df, df)
#'   cauchy   standard Cauchy;                   a0 = location, b0 = scale; tail (1, 1)
#'   gamma    gamma(shape), rate 1;              b0 = scale = 1/rate
#'   chisq    chi-square with df
#'   beta     beta(shape1, shape2) on [0, 1]
#'   lnorm    lognormal(0, sdlog);               b0 = exp(meanlog)
#'   weibull  Weibull(shape), scale 1;           b0 = scale
#'   f        F(df1, df2);                       tail indices (Inf, df2/2)
#' The F density is proportional to x^(df1/2 - 1) (1 + df1 x/df2)^(-(df1 + df2)/2),
#' which decays like x^(-df2/2 - 1), hence its upper tail index df2/2.
#'
#' The five families of the first release are built exactly as before, with
#' their closed forms, which return the mean, the variance and the covariance
#' in the units of the standard member, or NA where the family has no closed
#' form for that rank (the Gumbel below the maximum). The beta carries in
#' `upper` the functions of the reflected variable 1 - X, which is the beta
#' with its two shapes exchanged, so that an infinite density at the upper
#' endpoint (shape2 < 1) is integrated in the distance from it (.os_end_view()).
#' A family whose density can be infinite at an endpoint carries in `end` the
#' index a of that endpoint, with the density behaving like distance^(a - 1):
#' the shape of the gamma and of the Weibull, df/2 of the chi-square, df1/2 of
#' the F, and the two shapes of the beta. The F takes its quantile function
#' from .qf_beta(), because R's qf() loses the digits of a small quantile. Its
#' df1 must be finite here, although R allows an infinite one.
#'
#' The noncentral versions (R's `ncp`) are refused. A parameter that is not
#' the family's own stops with a message that lists the family's parameters,
#' and for `location` or `scale` the message points at os_shift and
#' os_multiplier.
#' @keywords internal
#' @noRd
.os_family <- function(name, args = list()) {
  spec <- .os_family_spec[[name]]
  fail <- function(...) stop("os_moments(): ", ..., call. = FALSE)
  own <- names(spec$par)
  shown <- if (name == "gamma") "shape and rate (or scale)" else
    paste(own, collapse = if (length(own) == 2L) " and " else ", ")
  if (length(args)) {
    nms <- names(args)
    if (is.null(nms) || any(!nzchar(nms)))
      fail("the parameters of the \"", name, "\" family must be named: ", shown,
           ". Every argument of os_moments() after 'dist' is given by its full name")
    if ("ncp" %in% nms)
      fail("the noncentral ", spec$label, " distribution (R's 'ncp') is not built in; ",
           "it can be supplied through os_dist()")
    bad <- setdiff(nms, own)
    if (length(bad))
      fail("'", bad[1], "' is not a parameter of the \"", name, "\" family, whose parameter",
           if (length(own) > 1L) "s are " else " is ", shown,
           if (bad[1] %in% c("location", "scale"))
             ". The shift a and the multiplier b of a + bX are 'os_shift' and 'os_multiplier'"
           else ". Every argument of os_moments() after 'dist' is given by its full name")
    if (anyDuplicated(nms)) fail("the parameter '", nms[duplicated(nms)][1], "' is given twice")
    for (nm in nms) {
      v <- args[[nm]]
      if (!is.numeric(v) || length(v) != 1L || is.na(v))
        fail("'", nm, "' must be one number")
    }
  }
  par <- spec$par
  if (name == "gamma") {                      # rate or scale, not both, as in R
    if (!is.null(args$rate) && !is.null(args$scale))
      fail("give 'rate' or 'scale' for the \"gamma\" family, not both")
    if (!is.null(args$scale)) { args$rate <- 1 / args$scale; args$scale <- NULL }
    par <- par[c("shape", "rate")]
  }
  par[names(args)] <- unlist(args)
  if (anyNA(par))
    fail("the \"", name, "\" family needs '", names(par)[is.na(par)][1], "', which has no default")
  pos <- function(nm, finite = TRUE) if (!(par[[nm]] > 0) || (finite && !is.finite(par[[nm]])))
    fail("'", nm, "' must be ", if (finite) "one positive finite number" else "positive")
  fin <- function(nm) if (!is.finite(par[[nm]])) fail("'", nm, "' must be finite")
  # the p, d and q functions of R with the shape parameters fixed
  wrap <- function(P, D, Q, ...) list(
    p = function(q, lower.tail = TRUE, log.p = FALSE) P(q, ..., lower.tail = lower.tail, log.p = log.p),
    d = function(x, log = FALSE) D(x, ..., log = log),
    q = function(p) Q(p, ...))
  build <- function(w, ...) .os_validate_dist(structure(
    list(name = name, p = w$p, d = w$d, q = w$q, ..., center = NA_real_,
         closed = list(mean = NULL, var = NULL, cov = NULL)), class = "os_dist"), mass = FALSE)
  a0 <- 0; b0 <- 1
  obj <- switch(name,
    norm = { fin("mean"); pos("sd"); a0 <- par[["mean"]]; b0 <- par[["sd"]]
      os_dist(pnorm, dnorm, qnorm, support = c(-Inf, Inf), symmetric = TRUE, name = "norm") },
    unif = { fin("min"); fin("max")
      if (!(par[["max"]] > par[["min"]])) fail("'max' must be above 'min'")
      a0 <- par[["min"]]; b0 <- par[["max"]] - par[["min"]]
      os_dist(punif, dunif, qunif, support = c(0, 1), symmetric = TRUE, name = "unif") },
    exp = { pos("rate"); b0 <- 1 / par[["rate"]]
      os_dist(pexp, dexp, qexp, support = c(0, Inf), symmetric = FALSE, name = "exp") },
    logis = { fin("location"); pos("scale"); a0 <- par[["location"]]; b0 <- par[["scale"]]
      os_dist(plogis, dlogis, qlogis, support = c(-Inf, Inf), symmetric = TRUE, name = "logis") },
    gumbel = { fin("location"); pos("scale"); a0 <- par[["location"]]; b0 <- par[["scale"]]
      os_dist(.pgumbel, .dgumbel, .qgumbel, support = c(-Inf, Inf), symmetric = FALSE,
              name = "gumbel") },
    t = { pos("df", finite = FALSE)
      build(wrap(stats::pt, stats::dt, stats::qt, df = par[["df"]]), support = c(-Inf, Inf),
            symmetric = TRUE, tail = c(par[["df"]], par[["df"]])) },
    cauchy = { fin("location"); pos("scale"); a0 <- par[["location"]]; b0 <- par[["scale"]]
      build(list(p = stats::pcauchy, d = stats::dcauchy, q = stats::qcauchy),
            support = c(-Inf, Inf), symmetric = TRUE, tail = c(1, 1)) },
    gamma = { pos("shape"); pos("rate"); b0 <- 1 / par[["rate"]]
      build(wrap(stats::pgamma, stats::dgamma, stats::qgamma, shape = par[["shape"]]),
            support = c(0, Inf), symmetric = FALSE, tail = c(Inf, Inf),
            end = c(lo = par[["shape"]], hi = NA_real_)) },
    chisq = { pos("df")
      build(wrap(stats::pchisq, stats::dchisq, stats::qchisq, df = par[["df"]]),
            support = c(0, Inf), symmetric = FALSE, tail = c(Inf, Inf),
            end = c(lo = par[["df"]] / 2, hi = NA_real_)) },
    beta = { pos("shape1"); pos("shape2")
      o <- build(wrap(stats::pbeta, stats::dbeta, stats::qbeta, shape1 = par[["shape1"]],
                      shape2 = par[["shape2"]]), support = c(0, 1),
                 symmetric = par[["shape1"]] == par[["shape2"]], tail = c(Inf, Inf),
                 end = c(lo = par[["shape1"]], hi = par[["shape2"]]))
      o$upper <- wrap(stats::pbeta, stats::dbeta, stats::qbeta, shape1 = par[["shape2"]],
                      shape2 = par[["shape1"]])   # 1 - X is the beta with the shapes exchanged
      o },
    lnorm = { fin("meanlog"); pos("sdlog"); b0 <- exp(par[["meanlog"]])
      if (!is.finite(b0) || b0 <= 0) fail("exp(meanlog) is not a positive finite number")
      build(wrap(stats::plnorm, stats::dlnorm, stats::qlnorm, sdlog = par[["sdlog"]]),
            support = c(0, Inf), symmetric = FALSE, tail = c(Inf, Inf)) },
    weibull = { pos("shape"); pos("scale"); b0 <- par[["scale"]]
      build(wrap(stats::pweibull, stats::dweibull, stats::qweibull, shape = par[["shape"]]),
            support = c(0, Inf), symmetric = FALSE, tail = c(Inf, Inf),
            end = c(lo = par[["shape"]], hi = NA_real_)) },
    f = { pos("df1"); pos("df2", finite = FALSE)
      w <- wrap(stats::pf, stats::df, stats::qf, df1 = par[["df1"]], df2 = par[["df2"]])
      if (is.finite(par[["df2"]])) w$q <- .qf_beta(par[["df1"]], par[["df2"]])
      build(w, support = c(0, Inf), symmetric = FALSE, tail = c(Inf, par[["df2"]] / 2),
            end = c(lo = par[["df1"]] / 2, hi = NA_real_)) })
  obj$closed <- switch(name,
    unif = list(mean = .os_closed_unif_mean, var = .os_closed_unif_var,
                cov = .os_closed_unif_cov),
    exp = list(mean = .os_closed_exp_mean, var = .os_closed_exp_var,
               cov = .os_closed_exp_cov),
    logis = list(mean = .os_closed_logis_mean, var = .os_closed_logis_var, cov = NULL),
    gumbel = list(mean = .os_closed_gumbel_mean, var = .os_closed_gumbel_var, cov = NULL),
    list(mean = NULL, var = NULL, cov = NULL))
  list(dist = obj, a0 = a0, b0 = b0, par = par)
}

#' A family name with its parameters, or an os_dist object, resolved
#'
#' Returns the list of .os_family(): the os_dist object whose moments the
#' engine computes, the shift a0 and the multiplier b0 that the family's
#' parameters give (0 and 1 for a user's object), and the parameters. With an
#' os_dist object, `args` must be empty.
#' @keywords internal
#' @noRd
.os_resolve_dist <- function(dist, args = list()) {
  families <- names(.os_family_spec)
  if (inherits(dist, "os_dist")) {
    if (length(args))
      stop("os_moments(): the argument", if (length(args) > 1L) "s " else " ",
           paste0("'", if (is.null(names(args))) "..." else names(args), "'", collapse = ", "),
           " cannot be used with a distribution made by os_dist(), which has its parameters ",
           "inside its functions. The shift a and the multiplier b of a + bX are 'os_shift' ",
           "and 'os_multiplier', and every argument after 'dist' is given by its full name",
           call. = FALSE)
    return(list(dist = dist, a0 = 0, b0 = 1, par = NULL))
  }
  if (is.character(dist) && length(dist) == 1L && dist %in% families)
    return(.os_family(dist, args))
  stop("os_moments(): 'dist' must be one of ", paste0("\"", families, "\"", collapse = ", "),
       ", or an object made by os_dist()", call. = FALSE)
}

#' Validate a distribution object and compute its center
#'
#' Checks the types, the support, the tail indices, the points where the density
#' is not smooth, and that the three functions
#' honor their arguments at the median and are vectorized. Returns the object
#' with its center Q(1/2) filled in and its breaks sorted. `mass = FALSE` skips the integral of the
#' density between the quartiles, for the built-in families.
#' @keywords internal
#' @noRd
.os_validate_dist <- function(obj, mass = TRUE) {
  fail <- function(...) stop("os_dist(): ", ..., call. = FALSE)
  if (!is.function(obj$p) || !is.function(obj$d) || !is.function(obj$q))
    fail("'p', 'd' and 'q' must be functions")
  if (!is.numeric(obj$support) || length(obj$support) != 2L || anyNA(obj$support) ||
      obj$support[1] >= obj$support[2])
    fail("'support' must be two ordered endpoints, finite or infinite")
  if (!is.logical(obj$symmetric) || length(obj$symmetric) != 1L || is.na(obj$symmetric))
    fail("'symmetric' must be TRUE or FALSE")
  if (!is.numeric(obj$tail) || length(obj$tail) != 2L || anyNA(obj$tail) || any(obj$tail <= 0))
    fail("'tail' must be two positive tail indices (Inf allowed)")
  if (!is.character(obj$name) || length(obj$name) != 1L)
    fail("'name' must be a single string")
  # Each call names the signature it needs when it fails, and each value is checked for
  # being one finite number before it is compared, so that a missing argument or an NA
  # gives an error that says what to fix.
  sig_p <- "p(q, lower.tail = TRUE, log.p = FALSE)"
  sig_d <- "d(x, log = FALSE)"
  sig_q <- "q(p)"
  call_fn <- function(what, sig, expr) tryCatch(expr, error = function(e)
    fail(what, " failed (", conditionMessage(e), "); the function must have the form ", sig))
  near <- function(v, target) is.numeric(v) && length(v) == 1L && is.finite(v) &&
    abs(v - target) <= tol
  show <- function(v) if (is.numeric(v) && length(v) == 1L) signif(v, 8) else deparse(v)
  tol <- 1e-8
  # the median, and the three functions evaluated there
  cen <- call_fn("q(0.5)", sig_q, obj$q(0.5))
  if (!is.numeric(cen) || length(cen) != 1L || !is.finite(cen))
    fail("q(0.5) must return one finite number")
  if (cen < obj$support[1] || cen > obj$support[2])
    fail("the median q(0.5) = ", signif(cen, 6), " lies outside 'support'")
  pm <- call_fn("p(q(0.5))", sig_p, obj$p(cen))
  if (!near(pm, 0.5))
    fail("p(q(0.5)) must be 0.5, got ", show(pm), ": 'p' and 'q' disagree")
  # the two tails differ at the lower quartile, where the median cannot tell them apart
  x1 <- call_fn("q(0.25)", sig_q, obj$q(0.25))
  if (!is.numeric(x1) || length(x1) != 1L || !is.finite(x1))
    fail("q(0.25) must return one finite number")
  if (!near(call_fn("p(x)", sig_p, obj$p(x1)), 0.25))
    fail("p(q(0.25)) must be 0.25: 'p' and 'q' disagree")
  if (!near(call_fn("p(x, lower.tail = FALSE)", sig_p, obj$p(x1, lower.tail = FALSE)), 0.75))
    fail("p(x, lower.tail = FALSE) does not return the upper tail 1 - F(x)")
  if (!near(call_fn("p(x, log.p = TRUE)", sig_p, obj$p(x1, log.p = TRUE)), log(0.25)))
    fail("p(x, log.p = TRUE) does not return the log probability log F(x)")
  if (!near(call_fn("p(x, lower.tail = FALSE, log.p = TRUE)", sig_p,
                    obj$p(x1, lower.tail = FALSE, log.p = TRUE)), log(0.75)))
    fail("p(x, lower.tail = FALSE, log.p = TRUE) does not return the log upper tail")
  dm <- call_fn("d(x)", sig_d, obj$d(cen))
  dl <- call_fn("d(x, log = TRUE)", sig_d, obj$d(cen, log = TRUE))
  if (!is.numeric(dm) || length(dm) != 1L || !is.finite(dm) || dm <= 0)
    fail("d() must be finite and positive at the median")
  if (!near(dl, log(dm)))
    fail("d(x, log = TRUE) does not return the log density at the median")
  # the density against the distribution function: d integrates to F(q3) - F(q1) = 1/2
  # between the quartiles (a density with the wrong constant or Jacobian, 1.001 dnorm,
  # passed every check at the median and gave a variance of 1.001)
  # (skipped for a built-in family, whose density and distribution function are R's own:
  # for a gamma with shape 0.05 the quartiles are fifteen orders of magnitude apart and
  # the plain integral between them does not converge)
  x3 <- call_fn("q(0.75)", sig_q, obj$q(0.75))
  # the points where the density is not smooth: finite, strictly inside the support,
  # with a density that is not infinite, and symmetric about the median when the
  # distribution is. The symmetry is judged in units of the interquartile range, so that it
  # does not depend on the scale on which the distribution is supplied, with an allowance
  # for the rounding of the points and of the median.
  br <- obj$breaks
  if (!is.null(br)) {
    if (!is.numeric(br) || !length(br) || anyNA(br) || !all(is.finite(br)))
      fail("'breaks' must be NULL or a vector of finite numbers")
    br <- sort(unique(as.numeric(br)))
    if (any(br <= obj$support[1]) || any(br >= obj$support[2]))
      fail("'breaks' must lie strictly inside 'support'")
    db <- call_fn("d(breaks, log = TRUE)", sig_d, obj$d(br, log = TRUE))
    if (!is.numeric(db) || length(db) != length(br) || anyNA(db) || any(db == Inf))
      fail("d() must be finite at 'breaks': the points are for a kink or a jump of the ",
           "density, not for a point where it is infinite")
    if (obj$symmetric) {
      round_off <- 8 * .Machine$double.eps * max(abs(cen), abs(br))
      # a median so large against the spread that the rounding of a point is no longer
      # small against the tolerance: the symmetry cannot be judged, and os_moments() would
      # lose the digits of x - median as well
      if (round_off > 1e-6 * (x3 - x1))
        fail("the median ", signif(cen, 6), " is too large against the spread to place ",
             "'breaks'; supply the distribution near the origin and pass the shift to ",
             "os_moments() as os_shift")
      if (any(abs(rev(2 * cen - br) - br) > 1e-8 * (x3 - x1) + round_off))
        fail("'breaks' must be symmetric about the median ", signif(cen, 6),
             " when symmetric = TRUE")
    }
    obj$breaks <- br
  }
  if (mass) {
    # the integral is cut at the declared points between the quartiles, where the density
    # is not smooth (a narrow component between two declared jumps was otherwise missed)
    # A piece whose two limits have one sign and lie more than a factor 10 apart is
    # integrated in log|x|, where a power-law density is smooth (the quartiles of a density
    # that is infinite at 0 can lie many orders of magnitude apart, and the plain integral
    # then returned 0.74 for a density that is right).
    cuts <- c(x1, br[br > x1 & br < x3], x3)
    piece <- function(a, b) {
      if (a > 0 && b > 10 * a)
        return(integrate(function(s) obj$d(exp(s)) * exp(s), log(a), log(b), rel.tol = 1e-8)$value)
      if (b < 0 && a < 10 * b)
        return(integrate(function(s) obj$d(-exp(s)) * exp(s), log(-b), log(-a), rel.tol = 1e-8)$value)
      integrate(function(x) obj$d(x), a, b, rel.tol = 1e-8)$value
    }
    mass <- tryCatch(sum(vapply(seq_len(length(cuts) - 1L), function(i)
                       piece(cuts[i], cuts[i + 1L]), numeric(1))),
                     error = function(e) fail("d() could not be integrated between the quartiles: ",
                                              conditionMessage(e)))
    if (!is.finite(mass) || abs(mass - 0.5) > 1e-6)
      fail("d() integrates to ", signif(mass, 7), " between the quartiles, where F rises by 0.5: ",
           "'d' and 'p' disagree (a normalizing constant or a Jacobian is off)")
  }
  # vectorization
  xx <- c(cen, cen, cen)
  if (length(call_fn("p(x)", sig_p, obj$p(xx))) != 3L ||
      length(call_fn("d(x)", sig_d, obj$d(xx))) != 3L)
    fail("'p' and 'd' must be vectorized in their first argument")
  qq <- call_fn("q(p)", sig_q, obj$q(c(0.25, 0.75)))
  if (length(qq) != 2L || !all(is.finite(qq)) || qq[1] >= qq[2])
    fail("'q' must be vectorized and increasing in p")
  obj$center <- cen
  obj
}

#' The quantile function of the F distribution through the beta quantile of each tail
#'
#' With B = (df1 X/df2)/(1 + df1 X/df2) a beta(df1/2, df2/2) variable,
#' X = (df2/df1) B/(1 - B). R's qf() takes the beta quantile from the upper
#' tail for every p and forms 1/q - 1, which loses the digits of a small
#' quantile: for df1 = 0.5 and df2 = 10 it returns 0 at p = 1e-4, and
#' pf(qf(p)) differs from p by 2e-5 at p = 1e-3. Here the quantile below the
#' median is formed from the lower-tail beta quantile b as (df2/df1) b/(1 - b),
#' and the quantile above the median from c, the quantile of 1 - B, a
#' beta(df2/2, df1/2) variable, as (df2/df1) (1 - c)/c, so that each is a ratio
#' of numbers that keep their relative precision.
#' @keywords internal
#' @noRd
.qf_beta <- function(df1, df2) {
  ratio <- df2 / df1
  function(p) {
    out <- numeric(length(p))
    lo <- p <= 0.5
    if (any(lo)) { b <- stats::qbeta(p[lo], df1 / 2, df2 / 2); out[lo] <- ratio * b / (1 - b) }
    if (any(!lo)) {
      cc <- stats::qbeta(p[!lo], df2 / 2, df1 / 2, lower.tail = FALSE)
      out[!lo] <- ratio * (1 - cc) / cc
    }
    out
  }
}

#' The standard Gumbel distribution in the maximum form, F(x) = exp(-exp(-x))
#'
#' log F = -exp(-x) directly. log S = log(1 - F) as log(-expm1(log F)) up to
#' x = 30, and as the series -x + log(1 - exp(-x)/2) beyond, where exp(-x) is
#' near the double epsilon and -expm1(log F) loses its relative precision. These
#' are the two tail repairs of the feasibility check of the literature report.
#' @keywords internal
#' @noRd
.pgumbel <- function(q, lower.tail = TRUE, log.p = FALSE) {
  lF <- -exp(-q)
  if (lower.tail) {
    v <- lF
  } else {
    v <- log(-expm1(lF))
    big <- q > 30
    v[big] <- -q[big] + log1p(-exp(-q[big]) / 2)
  }
  if (log.p) v else exp(v)
}

#' @keywords internal
#' @noRd
.dgumbel <- function(x, log = FALSE) {
  v <- -x - exp(-x)
  if (log) v else exp(v)
}

#' @keywords internal
#' @noRd
.qgumbel <- function(p) -log(-log(p))
