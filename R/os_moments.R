# ==============================================================================
# File: os_moments.R
# Purpose: The main function of orderStatMoments: os_moments() returns the exact mean
#          vector and covariance matrix of the order statistics of a continuous
#          distribution at the requested ranks, from the closed forms where the
#          family has them and from the integrators of os_engine.R otherwise,
#          for a built-in family given by R's own parameters or a user's
#          distribution, with the shift and the multiplier, the symmetry halving and the
#          run-time checks, and keeps the moments of the standard member in the
#          object for the linear layer. Also the print method of the returned object.
#
# EXPORTED FUNCTIONS:
#   - os_moments()        Exact alpha, variances and B at the requested ranks
#   - print.orderStatMoments()   Print method
#
# INTERNAL HELPERS:
#   - .os_validate_args() Argument checks of os_moments()
#   - .os_mass_check()    The density of X_(r) integrates to 1, or the quadrature missed a mode
#   - .os_is_normal()     Whether the distribution's functions are R's normal functions
#   - .os_closed_value()  A closed-form hook's value, or NA when absent or not wanted

#################### EXPORTED MAIN FUNCTIONS ####################

#' Exact moments of order statistics
#'
#' Computes the means \eqn{\alpha_i = E[X_{(r_i)}]}, the variances and the
#' covariance matrix \eqn{B_{ij} = \mathrm{Cov}(X_{(r_i)}, X_{(r_j)})} of the
#' order statistics \eqn{X_{(1)} \le \dots \le X_{(n)}} of an i.i.d. sample of
#' size \code{n} from a continuous distribution, at the ranks \code{ranks}, by
#' numerical evaluation of the exact integrals, or from a closed form where the
#' family has one. The distribution is one of thirteen built-in families, given
#' as in R by its name and its parameters, or a distribution that the user
#' defines with \code{\link{os_dist}}.
#'
#' @param n The sample size, an integer \eqn{\ge 1}.
#' @param ranks The ranks, strictly increasing integers in \eqn{1, \dots, n}.
#'   The default is every rank, which gives the full covariance matrix.
#' @param dist The distribution: the name of a built-in family, \code{"norm"},
#'   \code{"unif"}, \code{"exp"}, \code{"logis"}, \code{"gumbel"}, \code{"t"},
#'   \code{"cauchy"}, \code{"gamma"}, \code{"chisq"}, \code{"beta"},
#'   \code{"lnorm"}, \code{"weibull"} or \code{"f"}, or an object made by
#'   \code{\link{os_dist}}.
#' @param ... The parameters of a built-in family, under the names and with the
#'   defaults of R's own distribution functions (the table in Details), for
#'   example \code{mean} and \code{sd} for \code{"norm"}, \code{shape} and
#'   \code{rate} for \code{"gamma"}, \code{df} for \code{"t"}. A parameter
#'   without a default in R must be given. Nothing is passed here when
#'   \code{dist} is an object made by \code{os_dist()}.
#' @param os_shift,os_multiplier The constants \eqn{a} and \eqn{b > 0} of
#'   \eqn{Y = a + bX}. The function returns the moments of the order statistics
#'   of \eqn{Y}, where \eqn{X} has the distribution that \code{dist} and its
#'   parameters give: the means become \eqn{a + b\,\alpha} and the covariances
#'   \eqn{b^2 B}. They serve where the distribution has no parameter of its
#'   own for a shift or a multiplier: a \eqn{t} distribution with a location
#'   and a scale, a shifted exponential, a three-parameter gamma, Weibull or
#'   lognormal, or a user's distribution.
#' @param cov \code{TRUE} for the full covariance matrix, \code{FALSE} for the
#'   means and the variances only, which need single integrals and no double
#'   integral.
#' @param method \code{"auto"} takes the closed form where the family has one
#'   and the quadrature elsewhere, \code{"integrate"} forces the quadrature,
#'   \code{"closed"} requires a closed form for every requested quantity and
#'   stops otherwise.
#' @param rel.tol The relative tolerance of the quadrature. The single integrals
#'   and the inner integral of a product moment run at \code{rel.tol}, the
#'   outer integral at \code{100 * rel.tol}.
#' @param check Whether to run the run-time checks of the Details and warn when
#'   a residual exceeds \eqn{10^{-6}}.
#'
#' @return An object of class \code{"orderStatMoments"}, a list with the components
#'   \code{n}, \code{ranks}, \code{alpha} (the means, named by rank),
#'   \code{var} (the variances), \code{B} (the covariance matrix, \code{NULL}
#'   when \code{cov = FALSE}), \code{quad_err} (the summed error estimates of
#'   the quadrature, see Details), \code{standard} (a list with the
#'   \code{alpha}, \code{var} and \code{B} of the standard member of the family,
#'   which \code{\link{os_blue}} reads), \code{dist} (the \code{os_dist} object
#'   of the standard member), \code{par} (the parameters of a built-in family,
#'   \code{NULL} for a user's distribution), \code{os_shift},
#'   \code{os_multiplier}, \code{affine} (the shift and the multiplier that
#'   carry the standard member to \eqn{Y}, which combine the family's own
#'   parameters with \code{os_shift} and \code{os_multiplier}), \code{cov},
#'   \code{method}, \code{route} (\code{"closed"}, \code{"integrate"} or
#'   \code{"mixed"}), \code{checks} (the residuals of the run-time checks,
#'   \code{NA} where a check did not apply) and \code{call}.
#'
#' @details
#' \strong{The built-in families.} A family is given as in R: \code{dist} is
#' the name that R's functions carry after \code{p}, \code{d} and \code{q}
#' (\code{"gamma"} for \code{pgamma}), and its parameters are passed by R's
#' own names. Every argument of \code{os_moments()} after \code{dist} is given
#' by its full name.
#' \tabular{lll}{
#'   \code{dist} \tab Parameters, with R's defaults \tab Closed forms returned by \code{method = "auto"} \cr
#'   \code{"norm"} \tab \code{mean = 0}, \code{sd = 1} \tab none \cr
#'   \code{"unif"} \tab \code{min = 0}, \code{max = 1} \tab every mean, variance and covariance \cr
#'   \code{"exp"} \tab \code{rate = 1} \tab every mean, variance and covariance \cr
#'   \code{"logis"} \tab \code{location = 0}, \code{scale = 1} \tab the means and the variances \cr
#'   \code{"gumbel"} \tab \code{location = 0}, \code{scale = 1} \tab the mean and the variance of the maximum \cr
#'   \code{"t"} \tab \code{df} \tab none \cr
#'   \code{"cauchy"} \tab \code{location = 0}, \code{scale = 1} \tab none \cr
#'   \code{"gamma"} \tab \code{shape}, and \code{rate = 1} or \code{scale} \tab none \cr
#'   \code{"chisq"} \tab \code{df} \tab none \cr
#'   \code{"beta"} \tab \code{shape1}, \code{shape2} \tab none \cr
#'   \code{"lnorm"} \tab \code{meanlog = 0}, \code{sdlog = 1} \tab none \cr
#'   \code{"weibull"} \tab \code{shape}, \code{scale = 1} \tab none \cr
#'   \code{"f"} \tab \code{df1}, \code{df2} \tab none
#' }
#' The Gumbel is in the maximum form, \eqn{F(x) = \exp(-e^{-(x - \mathrm{location})/\mathrm{scale}})},
#' and is the one family that R does not provide. The noncentral versions of
#' the \eqn{t}, the chi-square, the beta and the \eqn{F} (R's \code{ncp}) are not
#' built in, and \code{df1} of the \eqn{F} must be finite. The package computes the moments of the standard member of each
#' family, the member whose shape parameters are those given and whose other
#' parameters are at R's defaults, and applies the remaining parameters
#' afterwards, exactly: the mean and the standard deviation of the normal, the
#' ends of the uniform, the rate of the exponential, the gamma's rate or
#' scale, the Weibull's scale, and \eqn{\exp(\mathrm{meanlog})} of the
#' lognormal are a shift and a multiplier of the standard member. A
#' distribution with a large mean or a small rate therefore loses no digits.
#'
#' \strong{The integrals.} With cdf \eqn{F}, density \eqn{f} and
#' \eqn{y < z},
#' \deqn{f_{r:n}(x) = \frac{n!}{(r-1)!\,(n-r)!} F(x)^{r-1} \{1 - F(x)\}^{n-r} f(x),}
#' \deqn{f_{r,s:n}(y, z) = \frac{n!}{(r-1)!\,(s-r-1)!\,(n-s)!} F(y)^{r-1}
#'   \{F(z) - F(y)\}^{s-r-1} \{1 - F(z)\}^{n-s} f(y) f(z),}
#' and the moments are
#' \deqn{\mu_{r:n} = \int x f_{r:n}(x)\,dx, \qquad
#'   \mu^{(2)}_{r:n} = \int x^2 f_{r:n}(x)\,dx, \qquad
#'   \mu_{r,s:n} = \int\!\!\int_{y < z} y z f_{r,s:n}(y, z)\,dy\,dz,}
#' with \eqn{B_{ij} = \mu_{r_i, r_j:n} - \alpha_i \alpha_j}. The integrals are
#' evaluated by adaptive quadrature (\code{\link[stats]{integrate}}). The
#' numerical design has three parts.
#'
#' \emph{What is computed.} The integration is in standardized units of the
#' standard member, \eqn{t = \{x - Q(1/2)\}/\tilde\sigma}, where the spread
#' \eqn{\tilde\sigma} is the interquartile range divided by that of the
#' standard normal (so \eqn{\tilde\sigma = 1} for the normal). The absolute
#' tolerances and the cut points then mean the same for a distribution on
#' any scale. Each moment is taken about the center of its own order
#' statistic, \eqn{c_r = Q(p_r)} with \eqn{p_r = r/(n+1)}, so that a variance
#' is not the small difference of two large moments: taken about the median,
#' the variance of the uniform minimum kept seven digits at \eqn{n = 1000}
#' and none at \eqn{n = 10^5}.
#'
#' \emph{How it is integrated.} The range is cut at \eqn{c_r \pm h_r}, with
#' \eqn{h_r = d_r + 4/\sqrt{n+1}}, \eqn{d_r = 10 s_r / f(c_r)} and
#' \eqn{s_r = \sqrt{p_r(1-p_r)/(n+2)}}, so that the narrow peak of the order
#' statistic's density at large \eqn{n} lies inside one piece whose ends the
#' adaptive rule sees. The second term is a minimum half-width. The first
#' moment, whose integrand changes sign, is cut at \eqn{c_r} as well, and the
#' outer range of a product moment is cut at the points of both ranks. The
#' absolute tolerances are scaled to the width \eqn{d_r} of the order
#' statistic. A piece that the integrator flags is accepted only when its
#' error estimate is negligible, and the computation stops otherwise (see
#' "What stops with an error").
#'
#' \emph{How the integrand is evaluated.} The integrand is assembled on the
#' logarithmic scale, from \code{lbeta()}, the log-tails \code{p(x, log.p = TRUE)}
#' and \code{p(x, lower.tail = FALSE, log.p = TRUE)} and the log-density, so
#' that the powers of the tails never underflow. The difference
#' \eqn{F(z) - F(y)} is formed on the side where the two numbers are small, so
#' that nothing cancels.
#'
#' \strong{A density that is infinite at an endpoint.} A gamma, a chi-square,
#' a Weibull, a beta or an \eqn{F} distribution with a shape below 1 has a
#' density that grows like \eqn{e^{a - 1}} in the distance \eqn{e} from a
#' finite endpoint of its support. Next to such an endpoint the package
#' integrates in the variable \eqn{v = e^a}, in which the integrand is
#' bounded. The distance from the endpoint is then the integration variable
#' itself, so an order statistic that sits far closer to the endpoint than
#' the centered variable \eqn{t} can resolve is still placed. The range is
#' also cut around the peak of the order statistic in that variable. The
#' package's tests assert the result against exact values: the minimum of a
#' million observations from a beta(a, 1) at \eqn{a = 0.5}, 0.1 and 0.05,
#' within \eqn{10^{-10}} for the mean and \eqn{10^{-9}} for the variance in
#' relative terms, and the covariances next to the endpoint at
#' \eqn{n = 1001} within \eqn{10^{-7}}. For a built-in family
#' \eqn{a} is the family's own shape. For a user's distribution whose density
#' is infinite at a lower endpoint, \eqn{a} is estimated from the distribution
#' function, and the endpoint should be at 0, where the distance is
#' represented without rounding. The upper endpoint of the beta is treated in
#' the same way through the reflected variable \eqn{1 - X}. For a user's distribution with an infinite density at
#' an upper endpoint, the reflected variable is not available: the order
#' statistics next to that endpoint are computed less accurately (about
#' \eqn{10^{-10}} for the maximum of a beta(1/2, 1/2) supplied through
#' \code{os_dist()} at \eqn{n = 100}), and the computation stops once they
#' concentrate against the endpoint. Such a distribution is better supplied
#' reflected, with the endpoint at 0.
#'
#' \strong{The moment condition.} A moment is computed only if it exists:
#' with the tail indices \eqn{(\nu_L, \nu_U)} of the distribution,
#' \eqn{E|X_{(r)}|^2 < \infty} if and only if \eqn{\nu_L r > 2} and
#' \eqn{\nu_U (n - r + 1) > 2}, and a rank that fails stops the computation
#' with an error. The indices of the built-in families are set by the package:
#' \eqn{(\mathrm{df}, \mathrm{df})} for \code{"t"}, \eqn{(1, 1)} for
#' \code{"cauchy"}, \eqn{(\infty, \mathrm{df2}/2)} for \code{"f"}, and
#' infinite for the others, which always pass. The median of five Cauchy
#' observations has a variance, and the minimum and the maximum have no mean.
#'
#' \strong{The error estimate.} \code{quad_err} is the sum of the error
#' estimates that \code{integrate()} reports for the single integrals and for
#' the outer integrals of the product moments, in the units of the standard
#' member. It is an estimate, not a bound: the inner integral's estimate is
#' not carried into the outer one, and the errors of the means are not
#' propagated into the covariances. The accuracy evidence is the agreement
#' with the closed forms, the identities and the printed tables, documented in
#' the package's tests.
#'
#' \strong{Accuracy at large n.} With \code{method = "integrate"}, the means
#' and the variances of the uniform and the exponential at the minimum, the
#' quartiles, the median and the maximum agree with their closed forms to
#' within \eqn{10^{-9}} in relative terms at \eqn{n = 10^4}, \eqn{10^5} and
#' \eqn{10^6}, and so does the correlation of every pair of them, in absolute
#' terms. Up to the relative errors of the two variances, the relative error
#' of a covariance is the error of its correlation divided by the
#' correlation, so the covariance of two nearly independent
#' order statistics keeps fewer digits: the two extremes, whose correlation is
#' of order \eqn{1/n}, keep five to ten significant digits of their covariance
#' at those sizes, depending on the distribution and on \eqn{n}. The package's
#' tests assert the variances and the covariances of adjacent ranks at
#' \eqn{10^{-9}} and the extremes' covariance at \eqn{10^{-3}} beyond
#' \eqn{n = 10^4}. A distribution as heavy-tailed as the lognormal with
#' \code{sdlog = 2} is accurate to about \eqn{10^{-9}} in relative terms rather
#' than \eqn{10^{-11}}, because the far tail of the inner integral of a
#' covariance is resolved at a fixed tolerance. These are the cases that were
#' measured, not a bound for every distribution.
#'
#' \strong{What has been tested.} The tests cover three groups of cases. (1)
#' The normal, uniform,
#' exponential, logistic and Gumbel distributions: the closed forms of the
#' last four at sample sizes up to 1001, the identities on full covariance
#' matrices up to \eqn{n = 100}, the printed tables of the literature, and
#' the large samples of the paragraph above. (2) The other eight families at
#' ordinary parameter values: the minimum, the quartiles, the median and the
#' maximum at sample sizes up to 1001, at the ranks whose variance exists,
#' against an integral in the probability scale that uses none of the
#' package's integrands, and the identities on full matrices at \eqn{n = 10}
#' and 25 for the distributions that have a variance. (3) Cases that were
#' chosen because they are hard. The first kind is a density that is
#' infinite at an endpoint. Gamma, chi-square, Weibull, beta and \eqn{F}
#' distributions with \eqn{a} from 0.8 down to 0.05 were run at sample sizes
#' up to \eqn{10^6 + 1}. The beta(\eqn{1/K}, 1) and beta(1, \eqn{1/K})
#' distributions, whose moments are exact finite products, were run for
#' \eqn{K} from 2 to 40 (\eqn{a} down to 0.025), at sample sizes from
#' \eqn{10^3} to \eqn{10^7} and at ranks from the extreme to \eqn{0.999 n}.
#' The beta(\eqn{1/K}, 0.8) distribution, whose density is infinite at both
#' endpoints, was run for \eqn{K} from 20 to 100 next to its lower endpoint,
#' and its reflection, the beta(0.8, \eqn{1/K}) distribution, for \eqn{K} of
#' 20 and 30.
#' The second kind is a variance that barely exists: a \eqn{t} distribution
#' with 1.005 to 3 degrees of freedom and an \eqn{F} distribution with
#' \code{df2} from 4.02 to 6, at the ranks where the variance exists and at
#' sample sizes up to \eqn{10^6 + 1}. In group (3) the package was asked for
#' the mean and the variance of 25898 order statistics. It returned 25595 of
#' them, each within \eqn{10^{-9}} of its reference (a variance in relative
#' terms, a mean relative to the larger of the mean and the standard
#' deviation of its order statistic), and it stopped with an error for the
#' other 303. Of those, 206 have a variance that rounds to 0 in double
#' precision, and 67 are beta distributions that are refused when they are
#' built (see "What stops with an error"). A variance below the smallest
#' normal double, \eqn{2.2 \times 10^{-308}}, is returned to the precision of
#' a subnormal number only, about three digits at \eqn{10^{-320}}. Of 594
#' covariances next to an endpoint with an infinite density, at sample sizes
#' from 1000 to \eqn{10^6}, the package returned 590, each within
#' \eqn{10^{-9}} of its exact value in units of the correlation (the error
#' divided by the two standard deviations) and in relative terms, and it
#' stopped for four. The largest sample size in these tests is
#' \eqn{10^6 + 1}, apart from the beta(\eqn{1/K}, 1) and beta(1, \eqn{1/K})
#' distributions, which were run up to \eqn{n = 10^7}.
#'
#' \strong{What stops with an error.} The quadrature follows the peak of each
#' order statistic's density and the tails of the distribution. It does not
#' follow a second mode: for a mixture of two narrow normals, the mode that
#' holds a tenth of the mass fell between the nodes, and the moments came back
#' as those of the other mode. The package therefore integrates the density of
#' every order statistic it computes by quadrature and stops when the mass is
#' not 1 within \eqn{10^{-7}}, naming the rank and the mass it found. It also
#' stops when a piece of the range cannot be resolved. This happens for a
#' moment that does not exist when the tail indices of a user's distribution
#' were not declared, for a tail so heavy that the mass of the integrand lies
#' far out (the covariances of a lognormal with \code{sdlog = 3} or of a
#' Weibull with \code{shape = 0.2}), for the covariance of an extreme order
#' statistic whose variance barely exists (the minimum of 20 observations from
#' a \eqn{t} distribution with 2.05 to 2.3 degrees of freedom), and for a
#' variance that comes out zero or negative. A beta distribution whose
#' second shape is below about 0.03 can be refused when it is built, with a
#' message that \code{p} and \code{q} disagree: its median lies so close to
#' 1 that the rounding of the median moves the distribution function, and
#' the two cannot be checked against each other there. Its order statistics
#' are one minus those of the beta distribution with the two shapes
#' exchanged, taken at the ranks \eqn{n + 1 - r}, which the package
#' computes. A wrong number is never
#' returned knowingly. Extreme cases cannot be checked exhaustively, and
#' further development can address the most difficult of them.
#'
#' \strong{The run-time checks.} For a symmetric distribution and a rank set
#' closed under \eqn{r \mapsto n + 1 - r}, the means are antisymmetric and the
#' second moments symmetric about the median, which the integrators do not
#' impose; the product moments are then computed for the pairs with
#' \eqn{r + s \le n + 1} only and reflected. For a complete rank set, the order
#' statistics are a permutation of the sample, so \eqn{\sum_r \alpha_r = n E[X]},
#' \eqn{\sum_r \mu^{(2)}_{r:n} = n E[X^2]} and \eqn{\sum_{r,s} B_{rs} = n\,\mathrm{Var}(X)}
#' for any distribution with a finite second moment, and for the normal
#' (recognized by its functions \code{pnorm}, \code{dnorm} and \code{qnorm}, not
#' by its name) every row of \eqn{B} sums to \eqn{\sigma^2} as well. The
#' residuals are stored in \code{checks}, in the standardized units of the
#' integration, so that the threshold means the same on every scale, and a
#' residual above \eqn{10^{-6}} raises a warning, as does a covariance matrix
#' with a negative eigenvalue beyond rounding. The checks detect some failures
#' of the quadrature and do not certify its accuracy: a failure that every
#' integral shares, such as a peak missed in the same way at every rank, can
#' satisfy them.
#'
#' \strong{Cost.} A single integral takes a fraction of a millisecond and a
#' double integral a few hundredths of a second. When the covariances are not
#' closed forms, a full matrix needs one double integral per pair,
#' \eqn{n(n-1)/2} of them, halved for a symmetric family.
#'
#' @examples
#' m <- os_moments(5)            # all ranks of n = 5, standard normal
#' m$alpha
#' m$B
#' os_moments(20, ranks = c(1, 5, 10, 20), dist = "exp", rate = 0.5)
#' os_moments(100, ranks = c(1, 100), dist = "gumbel", cov = FALSE)
#'
#' # a gamma as in rgamma(), and a t distribution with a location and a scale
#' os_moments(10, ranks = c(1, 5, 10), dist = "gamma", shape = 2, rate = 3)
#' os_moments(7, dist = "t", df = 5, os_shift = 3, os_multiplier = 2, cov = FALSE)
#'
#' # the Cauchy: the median of five has a variance, the minimum has none
#' os_moments(5, ranks = 3, dist = "cauchy")
#' try(os_moments(5, ranks = 1, dist = "cauchy"))
#' @references
#' Wu, Z. and Tong, T. (2026). orderStatMoments: Exact means, variances, and covariances
#' of order statistics in R. Manuscript. The computation, its verification and
#' the comparison with other software.
#'
#' Royston, J. P. (1982). Algorithm AS 177: Expected normal order statistics
#' (exact and approximate). \emph{Applied Statistics}, 31(2), 161-165.
#' @seealso \code{\link{os_dist}}
#' @export
os_moments <- function(n, ranks = seq_len(n), dist = "norm", ...,
                       os_shift = 0, os_multiplier = 1,
                       cov = TRUE, method = c("auto", "integrate", "closed"),
                       rel.tol = 1e-11, check = TRUE) {
  method <- match.arg(method)
  .os_validate_args(n, ranks, os_shift, os_multiplier, cov, rel.tol, check)
  n <- as.integer(n); ranks <- as.integer(ranks); k <- length(ranks)
  fam <- .os_resolve_dist(dist, list(...))   # the standard member, and the shift a0 and the
  dist <- fam$dist                           # multiplier b0 that the family's parameters give
  ctrl <- .os_ctrl(rel.tol)
  .os_moment_condition(ranks, n, dist$tail, k = 2)
  cd <- .os_centered(dist)
  .os_tail_probe(ranks, n, cd, ctrl, dist)
  c0 <- dist$center
  sp <- cd$spread     # the integrals are in units of the spread; sp and sp^2 bring them back
  pr <- ranks / (n + 1)
  cr <- cd$q(pr)                # each moment is taken about its own location, in those units
  m1 <- numeric(k)              # E[X_(r)] - c_r in those units, the first centered moment
  cl <- dist$closed
  err <- 0
  used_closed <- logical(0)                   # one entry per quantity: closed form or not
  need_quad <- function(what, r, s = NULL) {  # stop under method = "closed"
    if (method == "closed")
      stop("os_moments(): no closed form for the ", what, " of rank ",
           if (is.null(s)) r else paste(r, "and", s), " of the \"", dist$name,
           "\" family; use method = \"auto\" or \"integrate\"", call. = FALSE)
  }

  # ---- means and variances, in centered units of the standard member ----
  a_c <- numeric(k); v <- numeric(k)
  # With an endpoint view (a density that is infinite at a finite endpoint, .os_end_view()),
  # an order statistic can sit so close to that endpoint that its distance from it is far
  # below the rounding of the centered coordinate. Its mean is then assembled from the
  # endpoint, a_end, and its first centered moment is the integral's own value, so that
  # neither passes through the centered coordinate. With a view at both endpoints the
  # mean is assembled from the endpoint the order statistic is nearer to (.os_side()), not
  # from the lower one whenever its probability is below one half: assembled from the upper
  # endpoint, the mean of the rank 510000 of a million from a beta(0.015, 0.8), 4.7e-20,
  # came back as 1 - 1 = 0. Without a view both are formed as before, to the last bit.
  views <- !is.null(cd$lo) || !is.null(cd$hi)
  a_end <- rep(NA_real_, k)
  for (i in seq_len(k)) {
    r <- ranks[i]
    m <- .os_closed_value(cl$mean, method, r, n)
    if (is.na(m)) {
      need_quad("mean", r)
      .os_mass_check(r, n, cd, ctrl, cr[i], pr[i])
      o <- .os_single(r, n, 1, cd, ctrl, cr[i], pr[i]); a_c[i] <- sp * (cr[i] + o$value)
      err <- err + sp * o$abs.error
      used_closed <- c(used_closed, FALSE)
      if (views)
        a_end[i] <- if (.os_side(cd, pr[i], nearest = FALSE) == "lo")
          dist$support[1] + sp * (cd$lo$dq(pr[i]) + o$value)
        else
          dist$support[2] - sp * (cd$hi$dq(pr[i]) - o$value)
    } else {
      a_c[i] <- m - c0; used_closed <- c(used_closed, TRUE)
    }
    m1[i] <- if (views && is.na(m)) o$value else a_c[i] / sp - cr[i]
    vv <- .os_closed_value(cl$var, method, r, n)
    if (is.na(vv)) {
      need_quad("variance", r)
      o <- .os_single(r, n, 2, cd, ctrl, cr[i], pr[i]); v[i] <- sp * (sp * (o$value - m1[i]^2))
      err <- err + sp * (sp * o$abs.error)
      used_closed <- c(used_closed, FALSE)
    } else {
      v[i] <- vv; used_closed <- c(used_closed, TRUE)
    }
  }

  # ---- a variance is positive; one that is not was lost by the quadrature ----
  bad <- which(!is.finite(v) | v <= 0)
  if (length(bad))
    stop("os_moments(): the variance of X_(", ranks[bad[1]], ") came out as ",
         signif(v[bad[1]], 3), ", so the quadrature could not resolve it. This happens when ",
         "the order statistic is concentrated on a scale far below the spread of the ",
         "distribution: an extreme rank of a bounded or semi-bounded distribution at a ",
         "large n on the integration route (method = \"auto\" uses the closed form where one ",
         "exists)", call. = FALSE)

  # ---- covariances: the pairs with r + s <= n + 1 when the parent is symmetric ----
  # and the rank set is closed under reflection, the mirror pairs filled afterwards;
  # the rank sums are formed in doubles, so that they cannot overflow the integer range
  reflect <- dist$symmetric && all(as.numeric(ranks) + rev(ranks) == n + 1)
  B_c <- NULL
  if (cov && k >= 2L) {
    B_c <- matrix(NA_real_, k, k); diag(B_c) <- v
    for (i in seq_len(k - 1L)) for (j in seq.int(i + 1L, k)) {
      r <- ranks[i]; s <- ranks[j]
      if (reflect && as.numeric(r) + s > n + 1) next
      cv <- .os_closed_value(cl$cov, method, r, s, n)
      if (is.na(cv)) {
        need_quad("covariance", r, s)
        o <- .os_prod(r, s, n, cd, ctrl, cr[i], cr[j], pr[i], pr[j])
        cv <- sp * (sp * (o$value - m1[i] * m1[j]))
        err <- err + sp * (sp * o$abs.error)
        used_closed <- c(used_closed, FALSE)
      } else {
        used_closed <- c(used_closed, TRUE)
      }
      B_c[i, j] <- B_c[j, i] <- cv
    }
    if (reflect) for (i in seq_len(k - 1L)) for (j in seq.int(i + 1L, k))
      if (as.numeric(ranks[i]) + ranks[j] > n + 1)  # Cov(X_(r), X_(s)) = Cov(X_(n+1-s), X_(n+1-r))
        B_c[i, j] <- B_c[j, i] <- B_c[k + 1L - j, k + 1L - i]
  } else if (cov) {
    B_c <- matrix(v, 1L, 1L)
  }
  route <- if (method == "integrate" || !any(used_closed)) "integrate"
           else if (all(used_closed)) "closed" else "mixed"

  # ---- run-time checks, in the standardized units of the engine, so that the warning
  # threshold means the same on every scale (a normal supplied with standard deviation
  # 1e6 has residuals of 1e-3 in its own units, 1e-9 in standardized ones) ----
  checks <- list(reflection = NA_real_, identities = NA)
  if (check) {
    a_s <- a_c / sp; v_s <- v / (sp * sp)
    B_s <- if (is.null(B_c)) NULL else B_c / (sp * sp)
    if (reflect) checks$reflection <- .os_reflection_check(a_s, v_s + a_s^2)
    if (k == n) {
      pm <- .os_parent_moments(dist, cd, ctrl)
      pm_s <- list(m1 = pm$m1 / sp, m2 = pm$m2 / (sp * sp))
      checks$identities <- .os_identity_check(a_s, v_s, B_s, n, pm_s, .os_is_normal(dist))
    }
    resid <- c(checks$reflection, checks$identities)
    resid <- resid[!is.na(resid)]
    if (length(resid) && max(resid) > 1e-6)
      warning("os_moments(): a run-time check has residual ", signif(max(resid), 3),
              " at n = ", n, " for \"", dist$name, "\"; the quadrature is suspect",
              call. = FALSE)
    # a covariance matrix is positive semidefinite; rounding may leave a negative
    # eigenvalue of relative size near the double epsilon, a failed quadrature more
    if (!is.null(B_c) && k >= 2L) {
      ev <- eigen(B_c, symmetric = TRUE, only.values = TRUE)$values
      if (min(ev) < -sqrt(.Machine$double.eps) * max(ev))
        warning("os_moments(): the covariance matrix has the negative eigenvalue ",
                signif(min(ev), 3), " (largest ", signif(max(ev), 3), ") at n = ", n,
                " for \"", dist$name, "\"; the quadrature is suspect", call. = FALSE)
    }
  }

  # ---- add the center back, apply the shift and the multiplier, name by rank ----
  # Y = os_shift + os_multiplier X and X = a0 + b0 Z for the standard member Z, so
  # Y = a + b Z with the a and b below. The multiplier is applied to the variances and
  # covariances one factor at a time, so that b^2 cannot overflow or underflow where the
  # products themselves are representable.
  nm <- paste0("r", ranks)
  alpha0 <- ifelse(is.na(a_end), a_c + c0, a_end)     # the standard member's moments
  names(alpha0) <- nm
  var0 <- v; names(var0) <- nm
  B0 <- if (is.null(B_c)) NULL else { Bm <- B_c; dimnames(Bm) <- list(nm, nm); Bm }
  a <- os_shift + os_multiplier * fam$a0
  b <- os_multiplier * fam$b0
  if (!is.finite(a) || !is.finite(b) || b <= 0)
    stop("os_moments(): the shift ", format(a), " and the multiplier ", format(b),
         " of the standard member, from the family's parameters and 'os_shift' and ",
         "'os_multiplier', are not a finite number and a positive finite number", call. = FALSE)
  alpha <- a + b * alpha0
  var <- b * (b * var0)
  B <- if (is.null(B0)) NULL else b * (b * B0)
  if (!all(is.finite(alpha)) || !all(is.finite(var)) || (!is.null(B) && !all(is.finite(B))))
    stop("os_moments(): the moments overflow double precision at the shift ", format(a),
         " and the multiplier ", format(b), " of the standard member", call. = FALSE)
  structure(list(n = n, ranks = ranks, alpha = alpha, var = var, B = B, quad_err = err,
                 standard = list(alpha = alpha0, var = var0, B = B0),
                 dist = dist, par = fam$par, os_shift = os_shift, os_multiplier = os_multiplier,
                 affine = c(a = a, b = b), cov = cov,
                 method = method, route = route, checks = checks, call = match.call()),
            class = "orderStatMoments")
}

#' Print an orderStatMoments object
#'
#' @param x An object of class \code{"orderStatMoments"}.
#' @param digits The number of significant digits.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.orderStatMoments <- function(x, digits = 7L, ...) {
  par <- if (length(x$par)) paste0(" (", paste0(names(x$par), " = ", format(x$par, digits = digits),
                                               collapse = ", "), ")") else ""
  ab <- if (x$os_shift != 0 || x$os_multiplier != 1)
    paste0(", os_shift = ", format(x$os_shift), ", os_multiplier = ", format(x$os_multiplier)) else ""
  cat("Exact moments of order statistics: ", x$dist$name, par, ab, ", n = ", x$n,
      ", ranks ", paste(x$ranks, collapse = " "), "\n", sep = "")
  cat("alpha (means):\n"); print(signif(x$alpha, digits))
  if (is.null(x$B)) {
    cat("var (variances; cov = FALSE):\n"); print(signif(x$var, digits))
  } else {
    cat("B (covariance matrix):\n"); print(signif(x$B, digits))
  }
  cat("route: ", x$route, "; quadrature error estimate: ", format(signif(x$quad_err, 3)),
      "\n", sep = "")
  ck <- x$checks
  cat("checks: reflection ",
      if (is.na(ck$reflection)) "n/a" else format(signif(ck$reflection, 3)),
      "; identities ",
      if (length(ck$identities) == 1L && is.na(ck$identities)) "n/a"
      else paste(names(ck$identities), format(signif(ck$identities, 3)), collapse = ", "),
      "\n", sep = "")
  invisible(x)
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' Argument checks of os_moments()
#' @keywords internal
#' @noRd
.os_validate_args <- function(n, ranks, os_shift, os_multiplier, cov, rel.tol, check) {
  fail <- function(...) stop("os_moments(): ", ..., call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n != round(n) || n < 1)
    fail("'n' must be one integer >= 1")
  if (n > .Machine$integer.max)
    fail("'n' must be at most ", .Machine$integer.max, ", the largest integer R stores")
  if (!is.numeric(ranks) || length(ranks) < 1L || anyNA(ranks) || any(ranks != round(ranks)))
    fail("'ranks' must be integers")
  if (any(ranks < 1) || any(ranks > n))
    fail("'ranks' must lie in 1, ..., n = ", n)
  if (length(ranks) > 1L && any(diff(ranks) <= 0))
    fail("'ranks' must be strictly increasing")
  if (!is.numeric(os_shift) || length(os_shift) != 1L || !is.finite(os_shift))
    fail("'os_shift' must be one finite number")
  if (!is.numeric(os_multiplier) || length(os_multiplier) != 1L || !is.finite(os_multiplier) ||
      os_multiplier <= 0)
    fail("'os_multiplier' must be one positive number")
  if (!is.logical(cov) || length(cov) != 1L || is.na(cov))
    fail("'cov' must be TRUE or FALSE")
  if (!is.numeric(rel.tol) || length(rel.tol) != 1L || !is.finite(rel.tol) || rel.tol <= 0)
    fail("'rel.tol' must be one positive number")
  if (!is.logical(check) || length(check) != 1L || is.na(check))
    fail("'check' must be TRUE or FALSE")
  invisible(TRUE)
}

#' The mass of the density of X_(r), against 1
#'
#' The range is cut around the peak of each order statistic and follows the
#' parent's tails, and a second mode of the parent can fall between the nodes of
#' a piece: for the mixture 0.9 N(0, 0.01) + 0.1 N(10, 0.01) the mass at 10 was
#' missed entirely at n = 1, and the mean came back as 0.0014 instead of 1 with
#' every identity satisfied, because the parent moments are computed by the same
#' quadrature. The zeroth moment, which costs one single integral per rank, is
#' the one quantity with a known value, so it is checked before the moments are
#' trusted.
#' @keywords internal
#' @noRd
.os_mass_check <- function(r, n, cd, ctrl, shift, u = NULL) {
  mass <- .os_single(r, n, 0, cd, ctrl, shift, u)$value
  if (!is.finite(mass) || abs(mass - 1) > 1e-7)
    stop("os_moments(): the density of X_(", r, ") integrates to ", signif(mass, 6),
         " instead of 1, so the quadrature missed part of its mass. The integration range ",
         "is split around the peak of each order statistic's density, which does not follow ",
         "a parent with more than one mode, or with mass far from its median on the scale ",
         "of its spread", call. = FALSE)
  invisible(mass)
}

#' Whether a distribution object is the normal, by its functions
#'
#' The unit row sums of B hold for the normal parent only, so the identity check
#' applies them when the three functions are R's own normal functions, and never
#' on the strength of the label, which is the user's to choose.
#' @keywords internal
#' @noRd
.os_is_normal <- function(dist) {
  identical(dist$p, stats::pnorm) && identical(dist$d, stats::dnorm) &&
    identical(dist$q, stats::qnorm)
}

#' A closed-form hook's value, or NA when the hook is absent, not wanted, or
#' returns NA for this rank
#' @keywords internal
#' @noRd
.os_closed_value <- function(fn, method, ...) {
  if (method == "integrate" || is.null(fn)) return(NA_real_)
  v <- fn(...)
  if (is.null(v) || length(v) != 1L || is.na(v)) NA_real_ else as.numeric(v)
}
