# =============================================================================
# utils.R
# Mathematical primitives for sequential forecast comparison
#
# Implements the low-level boundary functions and helper utilities called by
# confidence_sequences.R, etests.R, and all downstream modules.
#
# All formulae are referenced to their source:
#   CR24  Choe & Ramdas (2024), Operations Research 72(4), 1368-1387
#   H21   Howard, Ramdas, McAuliffe & Sekhon (2021), Ann. Statist. 49(2)
#
# Dependencies: lamW (Lambert W function)
# Install via: install.packages("lamW")
# =============================================================================


# =============================================================================
# 1. TUNING PARAMETER: rho <-> v_opt mapping
# =============================================================================
#
# The GE mixture boundary and the CM boundary both take a tuning parameter rho.
# The user specifies v_opt (the intrinsic time at which the boundary is
# tightest), and rho is derived from it.
#
# Source: Howard et al. (2021) Proposition 3.
# Formula: rho = v_opt / (-W_{-1}(-alpha^2 / e) - 1)
# where W_{-1} is the lower branch of the Lambert W function.
#
# The Python comparecast package uses a different approximation:
#   rho_py = v_opt / (2*log(1/alpha_os) + log(1 + 2*log(1/alpha_os)))
#   where alpha_os = 2 * alpha_opt
#
# Both are implemented below. rho_from_vopt() uses the Lambert W formula
# (paper-exact). rho_from_vopt_python() mirrors the Python convention.
# Numerical agreement is checked in tests/test_utils.R.
# -----------------------------------------------------------------------------

#' Convert optimal intrinsic time to rho
#'
#' Maps the user-specified intrinsic time `v_opt` (the point at which a
#' boundary is tightest) to the `rho` tuning parameter, via the Lambert W
#' formula of Howard et al. (2021), Proposition 3 (paper-exact).
#'
#' @param v_opt  Numeric > 0. Intrinsic time at which the boundary is tightest.
#'               Recommended default from CR24: 10.
#' @param alpha  Numeric in (0,1). Significance level (one-sided).
#'               For a two-sided boundary at level alpha, pass alpha/2 here.
#'
#' @return Numeric > 0. The rho tuning parameter.
#'
#' @details
#' \deqn{\rho = \frac{v_{opt}}{-W_{-1}(-\alpha^2 / e) - 1}}
#' The lower branch \eqn{W_{-1}} is defined for `x` in `[-1/e, 0)` and returns
#' values `<= -1`. For `alpha` in `(0, 1)`, `-alpha^2/e` is always in
#' `(-1/e, 0)`, so the branch is well-defined.
#'
#' @examples
#' rho_from_vopt(v_opt = 10, alpha = 0.025)
#'
#' @export
rho_from_vopt <- function(v_opt = 10, alpha = 0.025) {
  stopifnot(v_opt > 0, alpha > 0, alpha < 1)
  arg    <- -(alpha^2) / exp(1)           # in (-1/e, 0) for alpha in (0,1)
  w_val  <- lamW::lambertWm1(arg)         # lower branch, returns value <= -1
  denom  <- -w_val - 1                    # > 0
  rho    <- v_opt / denom
  return(rho)
}

#' Convert v_opt to rho via Python comparecast approximation
#'
#' @param v_opt     Numeric > 0.
#' @param alpha_opt Numeric in (0,1). One-sided alpha used inside the e-process
#'                  (comparecast default: alpha/2 for two-sided tests).
#'
#' @return Numeric > 0.
#'
#' @details
#'   Python formula: rho = v_opt / (2*log(1/a) + log(1 + 2*log(1/a)))
#'   where a = 2 * alpha_opt.
#'   Used only for cross-validation against Python output; not used in
#'   production code.
#'
#' @keywords internal
#' @noRd
rho_from_vopt_python <- function(v_opt = 10, alpha_opt = 0.025) {
  stopifnot(v_opt > 0, alpha_opt > 0, alpha_opt < 1)
  a     <- 2 * alpha_opt
  denom <- 2 * log(1 / a) + log(1 + 2 * log(1 / a))
  rho   <- v_opt / denom
  return(rho)
}


# =============================================================================
# 2. BOUNDARY 1: Normal Mixture (Conjugate-Mixture, CM)
# =============================================================================
#
# Used for: Theorem 1 (Hoeffding CS), evaluated at v = t.
#
# Source: H21 Equation (3.7) / Proposition S1; CR24 Section 4.3.3.
#
# Formula:
#   u_CM(v; rho, alpha) = sqrt((v + rho) * log((v + rho) / (alpha^2 * rho)))
#
# Notes:
#   - alpha here is the ONE-SIDED alpha (pass alpha/2 for a two-sided CS).
#   - The scale c is implicit: the formula assumes the process has been scaled
#     so that deviations are 1-sub-Gaussian. In practice for CR24 Theorem 1
#     with |delta_i| <= c, the intrinsic time is v = c^2 * t. See Section 4.3.
#   - The formula is valid at v = 0: (0 + rho)*log(rho/(alpha^2*rho)) =
#     rho*log(1/alpha^2) > 0. No guard needed.
# -----------------------------------------------------------------------------

#' Normal mixture (CM) boundary
#'
#' @param v      Numeric vector >= 0. Intrinsic time values (V_t or t).
#' @param alpha  Numeric in (0,1). ONE-SIDED significance level.
#' @param rho    Numeric > 0. Tuning parameter. Obtain via rho_from_vopt().
#'
#' @return Numeric vector of boundary values (cumulative-sum scale, before /t).
#'
#' @examples
#'   rho <- rho_from_vopt(v_opt = 10, alpha = 0.025)
#'   u   <- cm_boundary(v = 1:500, alpha = 0.025, rho = rho)
#'
#' @seealso [rho_from_vopt()] to compute `rho` from a target `v_opt`.
#'
#' @export
cm_boundary <- function(v, alpha, rho) {
  stopifnot(all(v >= 0), alpha > 0, alpha < 1, rho > 0)
  sqrt((v + rho) * log((v + rho) / (alpha^2 * rho)))
}


# =============================================================================
# 3. BOUNDARY 2: Gamma-Exponential Mixture (GE)
# =============================================================================
#
# Used for: Theorem 2 (EB CS) boundary; Theorem 3 e-process value.
#
# Source: H21 Proposition S5; CR24 Proposition EC.3 / EC.2.1
#
# The function m(s, v) serves dual purpose:
#   (a) E-process value:    E_t^mix = m(S_t, V_hat_t)
#   (b) CS boundary:        u_GE(v; alpha) = sup{s : m(s,v) < 1/alpha}
#       solved numerically via ge_boundary() using uniroot().
#
# Formula (m(s, v)):
#   Let a = rho / c^2   (shape parameter)
#   Let x = (c*s + v + rho) / c^2
#
#   m(s, v) = [a^a / (Gamma(a) * gamma_reg(a, a))]
#             * Gamma(v/c^2 + a)                        [note: see (*) below]
#             * gamma_reg(v/c^2 + a,  x)
#             * x^(-(v/c^2 + a))
#             * exp((c*s + v) / c^2)
#
# (*) IMPORTANT: The second Gamma() call takes argument (v/c^2 + a),
#     not (v + rho/c^2). These are the same when a = rho/c^2, since
#     v/c^2 + rho/c^2 = (v + rho)/c^2. We use the fully expanded form
#     for clarity.
#
# gamma_reg(a, x) = pgamma(x, shape=a, rate=1) in R  [regularised lower IGF]
#
# Note on the denominator normalisation constant:
#   Gamma(a) * gamma_reg(a, a) = the UNREGULARISED lower incomplete gamma
#   evaluated at (a, a). This equals pgamma(a, shape=a, rate=1) * gamma(a).
#   In R: gamma(a) * pgamma(a, shape=a, rate=1).
# -----------------------------------------------------------------------------

#' Gamma-Exponential mixture function m(s, v)
#'
#' Core function for both the GE boundary (via uniroot) and the e-process
#' value (direct evaluation).
#'
#' @param s    Numeric. The "wealth" argument (cumulative sum S_t for e-process;
#'             boundary search variable for CS).
#' @param v    Numeric >= 0. Intrinsic time (V_hat_t).
#' @param rho  Numeric > 0. Tuning parameter from rho_from_vopt().
#' @param c    Numeric > 0. Sub-exponential scale. For Theorems 2&3 with
#'             |delta_i| <= c/2, use c = hi - lo (e.g. c=2 for scores in
#'             `[-1,1]`). For Proposition EC.4 Winkler scores, use c = 2.
#'
#' @return Numeric. Value of m(s, v). Used directly as E_t^mix when
#'         s = S_t and v = V_hat_t. Used inside ge_boundary() for root finding.
#'
#' @details
#'   Computed in log space for numerical stability, then exponentiated.
#'   The formula is valid for s < v/c + rho/c (otherwise x <= 0, log undefined).
#'   For large s (well past rejection threshold), the function returns Inf
#'   which the e-process wrapper clips to 1e7.
#'
#' @keywords internal
#' @noRd
ge_mixture <- function(s, v, rho, c) {
  stopifnot(rho > 0, c > 0, v >= 0)

  a  <- rho / c^2          # shape parameter
  x  <- (c * s + v + rho) / c^2   # upper limit for incomplete gamma

  # Guard: if x <= 0, m(s,v) is not defined (s is too negative)
  if (x <= 0) return(.Machine$double.eps)

  # Log-space computation for numerical stability
  # log m(s,v) = a*log(a) - log(Gamma(a)) - log(gamma_reg(a,a))
  #            + log(Gamma(v/c^2 + a))
  #            + log(gamma_reg(v/c^2 + a, x))
  #            - (v/c^2 + a) * log(x)
  #            + (c*s + v) / c^2

  shape2 <- v / c^2 + a   # = (v + rho) / c^2

  log_norm_const <- a * log(a) -
    lgamma(a) -
    log(pgamma(a, shape = a, rate = 1))   # log(gamma_reg(a,a))

  log_m <- log_norm_const +
    lgamma(shape2) +
    log(pgamma(x, shape = shape2, rate = 1)) -
    shape2 * log(x) +
    (c * s + v) / c^2

  return(exp(log_m))
}

#' Vectorised GE mixture
#'
#' Applies \code{ge_mixture()} element-wise over parallel vectors \code{s}
#' and \code{v} via \code{mapply}.
#'
#' @param s Numeric vector. Cumulative sum values.
#' @param v Numeric vector. Intrinsic time values, same length as \code{s}.
#' @param rho Numeric > 0. GE mixture tuning parameter.
#' @param c Numeric > 0. Sub-exponential scale.
#'
#' @return Numeric vector of GE mixture values.
#'
#' @keywords internal
#' @noRd
ge_mixture_vec <- function(s, v, rho, c) {
  mapply(ge_mixture, s = s, v = v, MoreArgs = list(rho = rho, c = c))
}


#' Gamma-exponential mixture boundary
#'
#' @param v      Numeric vector >= 0. Intrinsic time values.
#' @param alpha  Numeric in (0,1). ONE-SIDED significance level.
#' @param rho    Numeric > 0. From rho_from_vopt().
#' @param c      Numeric > 0. Sub-exponential scale.
#' @param s_lo   Numeric. Lower search bound for uniroot. Default: -10.
#' @param s_hi   Numeric. Upper search bound for uniroot. Default: 500.
#'
#' @return Numeric vector of boundary values (cumulative-sum scale).
#'
#' @importFrom stats pgamma uniroot
#'
#' @details
#' Computes \eqn{u_{GE}(v; \alpha, \rho, c) = \sup\{s : m(s,v) < 1/\alpha\}} by
#' solving `m(s, v) = 1/alpha` numerically for `s` via [uniroot()], separately
#' for each `v_i`.
#'
#' Root-finding fallback: the search starts in `[s_lo, s_hi]`; if
#' `m(s_hi, v_i)` has not yet crossed the target, `s_hi` is doubled once and
#' retried. If it still fails, a warning is issued and `s_hi` is returned as
#' a conservative fallback value. Increase `s_hi` directly if this warning
#' appears often (e.g. at large `v` or small `alpha`); increase `abs(s_lo)`
#' if no root is found at small `v`.
#'
#' Computed elementwise; can be slow for long vectors — consider caching
#' boundary values when the same `(alpha, rho, c)` are reused.
#'
#' @examples
#' rho <- rho_from_vopt(v_opt = 10, alpha = 0.025)
#' ge_boundary(v = 1:3, alpha = 0.025, rho = rho, c = 2)
#'
#' @export
ge_boundary <- function(v, alpha, rho, c, s_lo = -10, s_hi = 500) {
  stopifnot(alpha > 0, alpha < 1, rho > 0, c > 0)

  target <- 1 / alpha
  c_sq <- c^2
  a <- rho / c_sq
  log_xmax <- log(.Machine$double.xmax)

  log_norm_const <- a * log(a) -
    lgamma(a) -
    log(stats::pgamma(a, shape = a, rate = 1))

  vapply(v, function(vi) {
    shape2 <- vi / c_sq + a

    ge_value <- function(s) {
      x <- (c * s + vi + rho) / c_sq

      if (x <= 0) {
        return(.Machine$double.eps)
      }

      log_m <- log_norm_const +
        lgamma(shape2) +
        log(stats::pgamma(x, shape = shape2, rate = 1)) -
        shape2 * log(x) +
        (c * s + vi) / c_sq

      if (is.nan(log_m)) {
        return(NaN)
      }

      if (log_m > log_xmax) {
        return(.Machine$double.xmax)
      }

      exp(log_m)
    }

    s_hi_i <- s_hi
    m_hi <- ge_value(s_hi_i)

    if (is.nan(m_hi) || m_hi < target) {
      s_hi_try <- s_hi_i * 2
      m_hi_try <- ge_value(s_hi_try)

      if (!is.nan(m_hi_try) && m_hi_try >= target) {
        s_hi_i <- s_hi_try
      } else {
        warning(
          "ge_boundary: upper search bound may be too small at v = ", vi,
          ". Returning s_hi = ", s_hi_i, "."
        )

        return(s_hi_i)
      }
    }

    tryCatch(
      uniroot(
        function(s) ge_value(s) - target,
        lower = s_lo,
        upper = s_hi_i,
        tol = 1e-8
      )$root,
      error = function(e) {
        warning("ge_boundary: uniroot failed at v = ", vi, ": ", e$message)
        NA_real_
      }
    )
  }, FUN.VALUE = numeric(1))
}



# =============================================================================
# 4. BOUNDARY 3: Polynomial Stitched (PS)
# =============================================================================
#
# Used for: Alternative to CM/GE for both Theorems 1 and 2.
#           NOT the recommended primary boundary (CM/GE are tighter in CR24).
#           Included for completeness and for cross-checking CR24 Table results.
#
# Source: H21 Theorem 1 / Equation (3.3); CR24 Section 4.3.2.
#
# Formula (scalar, l0 = 1):
#   S_alpha(v) = k1 * sqrt(v * L(v)) + c * k2 * L(v)
#   where L(v) = s * log(log(eta*v/m)) + log(zeta(s) / (alpha * log(eta)^s))
#   k1 = (eta^(1/4) + eta^(-1/4)) / sqrt(2)
#   k2 = (sqrt(eta) + 1) / 2
#   zeta(s) = Riemann zeta function
#
# Evaluated as S_alpha(max(v, m)) to ensure log(log(.)) is well-defined.
# The inner argument of the outer log is eta*(v vee m)/m >= eta > 1,
# so log(eta*(v vee m)/m) >= log(eta) > 0, and log(log(.)) is finite.
#
# Recommended defaults CR24: s = 1.4, eta = 2, v_opt (= m) = 10.
#
# The Riemann zeta function at non-integer arguments requires an external
# package. We use the 'VGAM' package, or fall back to a precomputed value
# for the recommended s = 1.4.
# -----------------------------------------------------------------------------

# Precomputed zeta values for recommended fixed s values.
# Avoids VGAM dependency for the standard use case.
.ZETA_PRECOMPUTED <- c(
  "1.4" = 3.6028,  # zeta(1.4), accurate to 4 decimal places
  "1.5" = 2.6124,  # zeta(1.5)
  "2.0" = 1.6449   # zeta(2.0) = pi^2/6
)

#' Riemann zeta function
#'
#' Uses VGAM if available, otherwise uses precomputed values for common s.
#'
#' @param s Numeric > 1.
#' @return Numeric. zeta(s).
#'
#' @keywords internal
#' @noRd
.zeta <- function(s) {
  key <- as.character(round(s, 1))
  if (!is.null(.ZETA_PRECOMPUTED[[key]])) {
    return(.ZETA_PRECOMPUTED[[key]])
  }
  if (requireNamespace("VGAM", quietly = TRUE)) {
    return(VGAM::zeta(s))
  }
  stop(
    "Riemann zeta(", s, ") not precomputed and 'VGAM' package not available.\n",
    "Install VGAM with: install.packages('VGAM'), or use s = 1.4 (default)."
  )
}

#' Polynomial stitched (PS) boundary
#'
#' Alternative boundary for both Theorem 1 and Theorem 2 constructions,
#' included for completeness and for cross-checking CR24 Table results.
#'
#' @param v      Numeric vector >= 0. Intrinsic time values.
#' @param alpha  Numeric in (0,1). ONE-SIDED significance level.
#' @param v_opt  Numeric > 0. Optimal intrinsic time (= m in H21).
#'               Default: 10.
#' @param c      Numeric > 0. Sub-exponential scale.
#' @param s      Numeric > 1. Stitching parameter. Default: 1.4.
#' @param eta    Numeric > 1. Geometric spacing. Default: 2.
#'
#' @return Numeric vector of boundary values (cumulative-sum scale).
#'
#' @details
#' This is **not** the recommended primary boundary: the CM/GE mixture
#' boundaries (`cm_boundary()`, `ge_boundary()`) are tighter in CR24 and are
#' used by default throughout `seqcomp`. Use `ps_boundary()` only when you
#' specifically need the polynomial-stitched construction.
#'
#' @examples
#' ps_boundary(v = 1:5, alpha = 0.025, v_opt = 10, c = 1)
#'
#' @export
ps_boundary <- function(v, alpha, v_opt = 10, c = 1, s = 1.4, eta = 2) {
  stopifnot(all(v >= 0), alpha > 0, alpha < 1, v_opt > 0, c > 0, s > 1, eta > 1)

  m   <- v_opt
  v_  <- pmax(v, m)   # v vee m: ensures log(log(.)) well-defined

  k1  <- (eta^(1/4) + eta^(-1/4)) / sqrt(2)
  k2  <- (sqrt(eta) + 1) / 2
  z_s <- .zeta(s)

  # Inner quantity L(v)
  L   <- s * log(log(eta * v_ / m)) + log(z_s / (alpha * log(eta)^s))

  boundary <- k1 * sqrt(v_ * L) + c * k2 * L
  return(boundary)
}


#' Hardcoded 95% Empirical Bernstein boundary from CR24
#'
#' Equation provided directly in Choe & Ramdas (2024), EC.2.2 as an example.
#' Valid only for alpha = 0.05, c = 1. For other alpha, use ps_boundary()
#' or ge_boundary().
#'
#' @param v Numeric vector >= 0. Intrinsic time (V_hat_t).
#' @return Numeric vector of boundary values (cumulative-sum scale).
#'
#' @keywords internal
#' @noRd
cs_boundary_cr23_hardcoded <- function(v) {
  v_ <- pmax(v, 1)
  1.7 * sqrt(v_) * (log(log(2 * v_)) + 3.8) +
    3.4 * log(log(2 * v_)) + 13
}

# =============================================================================
# 5. INTRINSIC TIME HELPERS
# =============================================================================

#' Predictable centering sequence (default gamma_i)
#'
#' gamma_1 = 0, gamma_t = hat_Delta_{t-1} for t >= 2.
#' This is the lagged running mean of the score differences.
#'
#' @param xs Numeric vector of score differences (hat_delta_i).
#' @param lag Integer >= 1. Lag for the centering sequence. Default 1.
#'            For lag > 1 (stream-splitting setting), the first `lag` entries
#'            are set to 0, and the running mean is lagged by `lag` steps.
#'
#' @return Numeric vector of same length as xs.
#'
#' @details
#'   For lag = 1 (standard case):
#'     gamma_t = (1/(t-1)) * sum_{i=1}^{t-1} delta_i,  t >= 2
#'     gamma_1 = 0
#'   This matches both the paper default CR24 and the comparecast convention.
#'
#' @keywords internal
#' @noRd
make_gammas <- function(xs, lag = 1) {
  T_   <- length(xs)
  mus  <- cumsum(xs) / seq_along(xs)   # running mean
  gam  <- mus                           # copy
  if (lag >= T_) {
    gam[] <- 0
  } else {
    gam[(lag + 1):T_] <- mus[1:(T_ - lag)]
    gam[1:lag]        <- 0
  }
  return(gam)
}

#' Empirical intrinsic time V_hat_t
#'
#' V_hat_t = sum_{i=1}^t (delta_i - gamma_i)^2
#'
#' @param xs     Numeric vector of score differences.
#' @param gammas Numeric vector of predictable centres (same length as xs).
#'               If NULL, constructed via make_gammas(xs).
#' @param floor  Logical. If TRUE (default), apply pmax(1, V_hat_t) to prevent
#'               log(0) in boundary formulas. Set FALSE for e-process computation
#'               (the GE mixture does not need the floor).
#'
#' @return Numeric vector of cumulative intrinsic time values.
#'
#' @keywords internal
#' @noRd
intrinsic_time <- function(xs, gammas = NULL, floor = TRUE) {
  if (is.null(gammas)) gammas <- make_gammas(xs)
  vt <- cumsum((xs - gammas)^2)
  if (floor) vt <- pmax(1, vt)
  return(vt)
}


# =============================================================================
# 6. E-PROCESS LOG-SPACE COMPUTATION AND CLIPPING
# =============================================================================

#' Sub-exponential CGF-like function psi_{E,c}(lambda)
#'
#' psi_{E,c}(lambda) = (-log(1 - c*lambda) - c*lambda) / c^2
#'
#' Source: CR24 Theorem 3 formula.
#'
#' @param lambda Numeric. Must satisfy 0 <= lambda < 1/c.
#' @param c      Numeric > 0. Sub-exponential scale.
#'
#' @return Numeric.
#'
#' @keywords internal
#' @noRd
psi_e <- function(lambda, c) {
  stopifnot(lambda >= 0, lambda < 1 / c, c > 0)
  (-log(1 - c * lambda) - c * lambda) / c^2
}

#' Log fixed-lambda e-process from a score-difference stream
#'
#' Internal helper for computing the cumulative log e-process under a fixed
#' betting parameter \eqn{\lambda}, allowing a time-varying scale sequence
#' \eqn{c_1, c_2, \ldots}.
#'
#' At each step \eqn{i}, the log increment is
#' \deqn{\lambda x_i - \psi_E(\lambda, c_i)(x_i - \gamma_i)^2.}
#'
#' The key difference from \code{log_eprocess_fixed()} is that \code{c_seq}
#' is a vector. Passing \code{c_seq = rep(c0, length(xs))} recovers the
#' constant-scale case.
#'
#' @param xs Numeric vector. Score-difference stream.
#' @param c_seq Numeric vector, same length as \code{xs}. Time-varying
#'   sub-exponential scale parameters.
#' @param lambda Numeric scalar. Fixed betting parameter satisfying
#'   \eqn{0 \le \lambda < 1 / \max(c_seq)}.
#' @param gammas Numeric vector, same length as \code{xs}, or \code{NULL}.
#'   Predictable centering sequence. If \code{NULL}, constructed using
#'   \code{make_gammas(xs)}.
#'
#' @return Numeric vector of cumulative log e-process values.
#'
#' @keywords internal
#' @noRd
log_eprocess_fixed_predictable <- function(xs, c_seq, lambda, gammas = NULL) {
  stopifnot(length(xs) == length(c_seq), all(c_seq > 0))
  if (is.null(gammas)) gammas <- make_gammas(xs)
  stopifnot(length(gammas) == length(xs))

  c_0 <- max(c_seq)
  stopifnot(lambda >= 0, lambda < 1 / c_0)

  log_increments <- lambda * xs - psi_e(lambda, c_seq) * (xs - gammas)^2
  cumsum(log_increments)
}

#' Fixed-lambda log e-process
#'
#' log E_t(lambda) = lambda * S_t - psi_{E,c}(lambda) * V_hat_t
#'
#' @param xs     Numeric vector of score differences.
#' @param lambda Numeric in `[0, 1/c)`. Betting parameter.
#' @param c      Numeric > 0. Sub-exponential scale.
#' @param gammas Numeric vector of predictable centres, or NULL for default.
#'
#' @return Numeric vector of log e-process values.
#'
#' @keywords internal
#' @noRd
log_eprocess_fixed <- function(xs, lambda, c, gammas = NULL) {
  if (is.null(gammas)) gammas <- make_gammas(xs)
  S_t  <- cumsum(xs)
  V_t  <- intrinsic_time(xs, gammas, floor = FALSE)   # no floor for e-process
  lambda * S_t - psi_e(lambda, c) * V_t
}

#' Log GE mixture from cumulative sum and intrinsic time
#'
#' Internal helper for computing
#' `log(m(S_t, V_t))` from precomputed cumulative sums and intrinsic times.
#'
#' This is used when the caller has already constructed the signed cumulative
#' sum process and the matching intrinsic time process. It centralises the
#' numerical guard used throughout the package: non-positive or `NaN` GE mixture
#' values are replaced by `log(.Machine$double.eps)`.
#'
#' @param S_t Numeric vector. Cumulative sum process.
#' @param V_t Numeric vector. Intrinsic time process, same length as `S_t`.
#' @param rho Numeric > 0. GE mixture tuning parameter.
#' @param c Numeric > 0. Sub-exponential scale.
#'
#' @return Numeric vector of log GE mixture values.
#'
#' @keywords internal
#' @noRd
log_ge_mixture_from_sv <- function(S_t, V_t, rho, c) {
  stopifnot(length(S_t) == length(V_t), rho > 0, c > 0)

  mapply(
    function(s, v) {
      val <- ge_mixture(s, v, rho, c)
      if (val <= 0 || is.nan(val) || is.infinite(val)) return(log(.Machine$double.eps))
      log(val)
    },
    s = S_t,
    v = V_t
  )
}

#' Mixture log e-process via GE mixture
#'
#' log E_t^mix = log m(S_t, V_hat_t)
#'
#' @param xs      Numeric vector of score differences.
#' @param rho     Numeric > 0. Tuning parameter.
#' @param c       Numeric > 0. Sub-exponential scale.
#' @param gammas  Numeric vector or NULL.
#'
#' @return Numeric vector of log e-process values.
#'
#' @keywords internal
#' @noRd
log_eprocess_mixture <- function(xs, rho, c, gammas = NULL) {
  if (is.null(gammas)) gammas <- make_gammas(xs)
  S_t <- cumsum(xs)
  V_t <- intrinsic_time(xs, gammas, floor = FALSE)
  log_ge_mixture_from_sv(S_t, V_t, rho, c)
}

#' Clip and exponentiate log e-process
#'
#' Clips log_e to `[log(clip_min), log(clip_max)]` before exponentiating.
#' Matches Python comparecast default (clip_max = 1e7).
#'
#' @param log_e    Numeric vector of log e-process values.
#' @param clip_max Numeric > 0. Maximum e-process value. Default: 1e7.
#' @param clip_min Numeric or NULL. Minimum (optional; NULL = no lower clip).
#'
#' @return Numeric vector of e-process values.
#'
#' @keywords internal
#' @noRd
clip_eprocess <- function(log_e, clip_max = 1e7, clip_min = NULL) {
  log_e <- pmin(log_e, log(clip_max))
  if (!is.null(clip_min) && clip_min > 0) {
    log_e <- pmax(log_e, log(clip_min))
  }
  exp(log_e)
}


# =============================================================================
# 7. UTILITY HELPERS
# =============================================================================

#' Running mean
#' @param xs Numeric vector.
#' @return Numeric vector of cumulative means.
#'
#' @keywords internal
#' @noRd
running_mean <- function(xs) cumsum(xs) / seq_along(xs)

#' Score difference bounds -> sub-Gaussian / sub-exponential scale
#'
#' Given score-difference bounds `[lo, hi]`, computes the two scale constants
#' used elsewhere in `seqcomp`:
#' * `c_thm1`  = `(hi - lo) / 2` — sub-Gaussian scale for Theorem 1, where
#'   `|delta_i| <= c_thm1`.
#' * `c_thm23` = `hi - lo` — sub-exponential scale for Theorems 2 & 3, where
#'   `|delta_i| <= c_thm23 / 2`.
#'
#' Both conventions bound the same quantity: after centering, `delta_i` lies
#' in `[-(hi-lo)/2, (hi-lo)/2]`, so `max|delta_i| = (hi-lo)/2`, which equals
#' both `c_thm1` and `c_thm23 / 2`.
#'
#' @param lo Numeric. Lower bound of score difference (usually a - b).
#' @param hi Numeric. Upper bound of score difference (usually b - a).
#'
#' @return Named list with elements c_thm1 and c_thm23.
#'
#' @examples
#' score_diff_scales(lo = -1, hi = 1)
#'
#' @export
score_diff_scales <- function(lo, hi) {
  stopifnot(hi > lo)
  c_val <- hi - lo
  list(
    c_thm1  = c_val / 2,  # Theorem 1: |delta_i| <= c, so c = (hi-lo)/2
    c_thm23 = c_val       # Theorems 2 & 3: |delta_i| <= c/2, so c = hi-lo
  )
}

#' Safe log (avoids log(0) = -Inf warnings in calibration)
#' @param x Numeric.
#' @param eps Numeric. Small positive guard. Default: 1e-16.
#'
#' @keywords internal
#' @noRd
.safe_log <- function(x, eps = 1e-16) log(pmax(x, eps))


#' Slice global or pair-specific parameters for SMCS loops
#'
#' @keywords internal
#' @noRd
.slice_param <- function(param, Tt, i, j) {
  if (is.null(param)) return(NULL)

  # 3D Array: T x m x m (e.g., dynamic betting fractions)
  if (is.array(param) && length(dim(param)) == 3) {
    return(param[, i, j])
  }
  # 2D Matrix: m x m (e.g., fixed conservative bounds per pair)
  if (is.matrix(param)) {
    return(param[i, j])
  }
  # Scalar or T-length vector (e.g., global c = 2)
  return(param)
}

# =============================================================================
# 8. NUMERICAL VERIFICATION: rho formula cross-check
# =============================================================================
#
# Run this once to check if Lambert W and Python approximation agree.
# Called explicitly in tests/testthat/test_utils.R, not on package load.
# -----------------------------------------------------------------------------

#' Verify that the two rho formulas agree numerically
#'
#' @param v_opt     Numeric. Default: 10.
#' @param alpha_vec Numeric vector of alpha values to check. Default: standard grid.
#' @param tol       Relative tolerance. Default: 0.01 (1%).
#'
#' @return Invisibly returns a data.frame of results; prints a summary.
#'
#' @keywords internal
#' @noRd
verify_rho_formulas <- function(
    v_opt     = 10,
    alpha_vec = c(0.005, 0.010, 0.025, 0.050, 0.100),
    tol       = 0.01
) {
  results <- lapply(alpha_vec, function(a) {
    rho_lw  <- rho_from_vopt(v_opt = v_opt, alpha = a)
    rho_py  <- rho_from_vopt_python(v_opt = v_opt, alpha_opt = a)
    rel_diff <- abs(rho_lw - rho_py) / rho_lw
    data.frame(
      alpha    = a,
      rho_lambertW = rho_lw,
      rho_python   = rho_py,
      rel_diff     = rel_diff,
      agree        = rel_diff < tol
    )
  })
  out <- do.call(rbind, results)
  cat("=== rho formula cross-check (v_opt =", v_opt, ") ===\n")
  print(out, row.names = FALSE, digits = 6)
  if (all(out$agree)) {
    cat("OK: All values agree within", tol * 100, "% relative tolerance.\n")
  } else {
    cat("WARNING: Some values differ by more than", tol * 100, "%.\n")
    cat("Using Lambert W formula (paper-exact) in production code.\n")
  }
  invisible(out)
}





