# =============================================================================
# confidence_sequences.R
# Anytime-valid confidence sequences for comparing two forecasters
#
# Implements:
#   cs_hoeffding()  — Theorem 1 CR24: sub-Gaussian (Hoeffding-style) CS
#   cs_bernstein()  — Theorem 2 CR24: empirical Bernstein (variance-adaptive) CS
#   cs_asymptotic() — EC.3, Eq. EC.29 CR24: asymptotic CS for unbounded differences
#
# All three return a data.frame with columns:
#   t        — time index
#   estimate — running mean hat_Delta_t
#   lower    — lower confidence bound
#   upper    — upper confidence bound
#
# Coverage guarantee for cs_hoeffding and cs_bernstein:
#   P(forall t >= 1: Delta_t in [lower_t, upper_t]) >= 1 - alpha
#
# By design, cs_asymptotic does not have a finite-sample coverage guarantee.
# It is valid in the limit as m -> infinity under the assumption of finite
# variance. Waudby-Smith et al. (2021) define:
# lim inf_{m -> infinity} P(forall t >= m: mu_t in C^(m)_t) >= 1 - alpha
# m being the "peeking time".
#
# References:
#   CR24  Choe & Ramdas (2024), Operations Research 72(4), 1368-1387
#   H21   Howard et al. (2021), Annals of Statistics 49(2), 1055-1080
#   WS21  Waudby-Smith et al. (2021), Annals of Statistics 52(6), 2613-2640
# =============================================================================

#' Hoeffding-style confidence sequence (Theorem 1, Choe & Ramdas 2024)
#'
#' Constructs a time-uniform confidence sequence for the mean score difference
#' \eqn{\Delta_t = \frac{1}{t} \sum_{i=1}^t E[\hat{\delta_i} \mid \mathcal{F}_{i-1}]}.
#'
#' @param scores1 Numeric vector. Scores S(p_t, y_t) for forecaster 1.
#' @param scores2 Numeric vector. Scores S(q_t, y_t) for forecaster 2.
#' @param alpha   Numeric in (0,1). Significance level. The CS has coverage
#'                1 - alpha uniformly over all t. Default: 0.05.
#' @param c       Numeric > 0. Sub-Gaussian scale. The process must satisfy
#'                |hat_delta_i| <= c for all i. For scores in `[a,b]`, the
#'                difference is in `[a-b, b-a]`, so c = b - a.
#'                Default: 1 (appropriate for Brier score differences in `[-1,1]`).
#' @param v_opt   Numeric > 0. Intrinsic time at which the CS is tightest.
#'                Default: 10 (recommended by CR24).
#' @param boundary Character. "mixture" (default, recommended) or "stitching".
#'
#' @return data.frame with columns t, estimate, lower, upper.
#'
#' @section Assumption:
#' Requires `hat_delta_i` to be c-sub-Gaussian given \eqn{\mathcal{F}_{i-1}},
#' i.e. `|hat_delta_i| <= c` for all `i`.
#'
#' @section Boundary:
#' \deqn{C_t^H = \hat\Delta_t \pm u^{CM}_{\alpha/2}(c^2 t; \rho) / t}
#' where \eqn{u^{CM}} is the normal mixture boundary and \eqn{c^2 t} is the
#' intrinsic time for a c-sub-Gaussian process with deterministic variance
#' proxy. The intrinsic time for Theorem 1 is `v_t = c^2 * t`, not `v_t = t`:
#' the CM boundary implicitly assumes 1-sub-Gaussian inputs, so the `c^2`
#' scaling must be applied explicitly. This matches the H21 convention,
#' where the boundary absorbs the sub-Gaussian parameter via the variance
#' process definition.
#'
#' Relation to Python comparecast: Python uses `v_t = sigma * t` where
#' `sigma = (hi - lo)/2 = c`. This is equivalent to our `c^2 * t` only when
#' `c = 1`. For `c != 1` the parametrisations differ; we follow the paper.
#'
#' @section Output:
#' Returns a `data.frame` with one row per `t` and columns `t`, `estimate`
#' (the running mean `hat_Delta_t`), `lower`, and `upper`, with coverage
#' guarantee \eqn{P(\forall t \geq 1 : \Delta_t \in [\text{lower}_t,
#' \text{upper}_t]) \geq 1 - \alpha}.
#'
#' @examples
#' scores1 <- c(-0.04, -0.09, -0.01, -0.16)
#' scores2 <- c(-0.09, -0.16, -0.04, -0.25)
#' cs_hoeffding(scores1, scores2, alpha = 0.05)
#'
#' @export
cs_hoeffding <- function(scores1, scores2,
                         alpha    = 0.05,
                         c        = 1,
                         v_opt    = 10,
                         boundary = "mixture") {

  stopifnot(
    length(scores1) == length(scores2),
    length(scores1) >= 1,
    alpha > 0, alpha < 1,
    c > 0,
    v_opt > 0,
    boundary %in% c("mixture", "stitching")
  )

  xs  <- scores1 - scores2          # hat_delta_t: score differences
  T_  <- length(xs)
  ts  <- seq_len(T_)                # time indices 1, 2, ..., T

  # Running mean: hat_Delta_t = (1/t) * sum_{i=1}^t hat_delta_i
  mus <- cumsum(xs) / ts

  # Intrinsic time for Theorem 1: deterministic, v_t = c^2 * t
  # For a c-sub-Gaussian process, the variance proxy is c^2 per observation.
  vs  <- c^2 * ts

  # One-sided alpha for two-sided CS: each side has crossing prob alpha/2
  alpha_os <- alpha / 2

  # Tuning parameter rho from v_opt
  rho <- rho_from_vopt(v_opt = v_opt, alpha = alpha_os)

  # Boundary values (cumulative-sum scale)
  if (boundary == "mixture") {
    u_t <- cm_boundary(v = vs, alpha = alpha_os, rho = rho)
  } else {
    # Polynomial stitching: c parameter enters boundary explicitly
    u_t <- ps_boundary(v = vs, alpha = alpha_os, v_opt = v_opt * c^2,
                       c = c, s = 1.4, eta = 2)
  }

  # CS radius: divide boundary by t to convert from cumulative-sum to mean scale
  radii <- u_t / ts

  data.frame(
    t        = ts,
    estimate = mus,
    lower    = mus - radii,
    upper    = mus + radii
  )
}

#' Empirical Bernstein confidence sequence (Theorem 2, Choe & Ramdas 2024)
#'
#' Constructs a variance-adaptive time-uniform CS using empirical intrinsic
#' time \eqn{\hat{V}_t = \sum_{i=1}^t (\hat{\delta}_i - \gamma_i)^2}.
#' Tighter than the Hoeffding CS when score differences have low variance.
#'
#' The CS is:
#' \deqn{C_t^{EB} = \hat{\Delta}_t \pm u_{\alpha/2}^{GE}(\hat{V}_t;\, \rho, c) \;/\; t}
#'
#' @param scores1  Numeric vector. Scores for forecaster 1.
#' @param scores2  Numeric vector. Scores for forecaster 2.
#' @param alpha    Numeric in (0,1). Significance level. Default: 0.05.
#' @param c        Numeric > 0. Sub-exponential scale. The process must satisfy
#'                 |hat_delta_i| <= c/2. For score differences in `[a-b, b-a]`,
#'                 c = b - a (e.g. c = 2 for Brier score differences in `[-1,1]`).
#'                 Default: 2.
#' @param v_opt    Numeric > 0. Optimal intrinsic time. Default: 10.
#' @param boundary Character. "mixture" (default, GE mixture) or "stitching"
#'                 (polynomial stitched) or "hardcoded" (CR24 example formula,
#'                 only valid for alpha=0.05, c=1).
#' @param gammas   Numeric vector or NULL. Predictable centering sequence.
#'                 If NULL, constructed as lagged running mean (default).
#' @param lcb_only Logical. If TRUE, return lower CS only: `[lower, +Inf)`.
#'                 Requires finite lower bound on hat_delta_i; provide c.
#' @param ucb_only Logical. If TRUE, return upper CS only: `(-Inf, upper]`.
#'
#' @return data.frame with columns t, estimate, lower, upper.
#'         lower = -Inf if ucb_only = TRUE; upper = Inf if lcb_only = TRUE.
#'
#' @examples
#' scores1 <- c(-0.04, -0.09, -0.01, -0.16)
#' scores2 <- c(-0.09, -0.16, -0.04, -0.25)
#' cs_bernstein(scores1, scores2, alpha = 0.05)
#'
#' @export
cs_bernstein <- function(scores1, scores2,
                         alpha    = 0.05,
                         c        = 2,
                         v_opt    = 10,
                         boundary = "mixture",
                         gammas   = NULL,
                         lcb_only = FALSE,
                         ucb_only = FALSE) {

  stopifnot(
    length(scores1) == length(scores2),
    length(scores1) >= 1,
    alpha > 0, alpha < 1,
    c > 0,
    v_opt > 0,
    boundary %in% c("mixture", "stitching", "hardcoded"),
    !(lcb_only && ucb_only)   # cannot request both simultaneously
  )

  if (boundary == "hardcoded" && (abs(alpha - 0.05) > 1e-10 || abs(c - 1) > 1e-10)) {
    warning(
      "boundary='hardcoded' is only valid for alpha=0.05 and c=1. ",
      "Supplied alpha=", alpha, ", c=", c, ". ",
      "Switching to boundary='mixture'."
    )
    boundary <- "mixture"
  }

  xs  <- scores1 - scores2
  T_  <- length(xs)
  ts  <- seq_len(T_)

  # Running mean
  mus <- cumsum(xs) / ts

  # Predictable centering sequence gamma_i
  if (is.null(gammas)) {
    gammas <- make_gammas(xs, lag = 1)
  } else {
    stopifnot(length(gammas) == T_)
  }

  # One-sided gamma adjustment for one-sided CS (mirrors Python comparecast)
  # lcb_only: center is clipped from above to prevent upward drift of gammas
  # ucb_only: center is clipped from below
  # These adjustments ensure the variance process remains valid for one-sided CSs.
  if (lcb_only) {
    # For lower CS: gammas clipped to [-c/2, c/2] from above
    gammas <- pmin(c / 2, gammas)
  }
  if (ucb_only) {
    # For upper CS: gammas clipped from below
    gammas <- pmax(-c / 2, gammas)
  }

  # Empirical intrinsic time: V_hat_t = sum_{i=1}^t (hat_delta_i - gamma_i)^2
  # The GE mixture is analytically valid at v=0.
  # However, polynomial stitching and hardcoded boundaries require v >= 1
  # to prevent NaN values in their log-based formulas.
  needs_floor <- (boundary != "mixture")
  vs <- intrinsic_time(xs, gammas, floor = needs_floor)

  # Alpha for one-sided vs two-sided
  alpha_use <- if (lcb_only || ucb_only) alpha else alpha / 2

  # Compute boundary values
  if (boundary == "mixture") {
    rho <- rho_from_vopt(v_opt = v_opt, alpha = alpha_use)
    u_t <- ge_boundary(v = vs, alpha = alpha_use, rho = rho, c = c)
  } else if (boundary == "stitching") {
    u_t <- ps_boundary(v = vs, alpha = alpha_use, v_opt = v_opt,
                       c = c, s = 1.4, eta = 2)
  } else {
    # hardcoded: CR24 example formula (alpha=0.05, c=1 only)
    u_t <- cs_boundary_cr23_hardcoded(vs)
  }

  # CS radius (mean scale)
  radii <- u_t / ts

  # Build output: replace one endpoint with +/-Inf for one-sided CSs
  lower <- if (ucb_only) rep(-Inf, T_) else mus - radii
  upper <- if (lcb_only) rep( Inf, T_) else mus + radii

  data.frame(
    t        = ts,
    estimate = mus,
    lower    = lower,
    upper    = upper
  )
}

#' Asymptotic confidence sequence (EC.3, Eq. EC.29, Choe & Ramdas 2024)
#'
#' **Asymptotic, not finite-sample**: coverage `>= 1 - alpha` holds only as
#' `t -> infinity`. Valid without requiring bounded score differences —
#' requires only that `hat_delta_t` has finite variance — so it's appropriate
#' for tick loss and other scoring rules where hard bounds depend on
#' unbounded realised values. Suitable for large evaluation windows.
#'
#' @param scores1 Numeric vector. Scores for forecaster 1.
#' @param scores2 Numeric vector. Scores for forecaster 2.
#' @param alpha   Numeric in (0,1). Significance level. Default: 0.05.
#' @param t_star  Numeric > 0. Sample size at which CS is tightest.
#'                Default: length(scores1) (tightest at end of sample).
#'
#' @return data.frame with columns t, estimate, lower, upper.
#'
#' @details
#' \deqn{C_t^A = \hat\Delta_t \pm \sqrt{
#'   \frac{2(t \sigma^2_t \rho^2 + 1)}{t^2 \rho^2}
#'   \log\frac{\sqrt{t \sigma^2_t \rho^2 + 1}}{\alpha}}}
#' where \eqn{\sigma^2_t = \frac{1}{t}\sum_{i=1}^t (\hat\delta_i - \hat\Delta_{i-1})^2}
#' and `rho` is tuned to be tightest at `t_star`:
#' \deqn{\rho(t_{star}) = \sqrt{\frac{2\log(1/\alpha) + \log(1 + 2\log(1/\alpha))}{t_{star}}}}
#'
#' The running variance estimator uses the predictable mean `hat_Delta_{t-1}`
#' (not the current mean `hat_Delta_t`) to maintain predictability, with
#' `hat_Delta_0 := 0`.
#'
#' @examples
#' scores1 <- c(-0.4, -0.2, -0.3, -0.1, -0.2)
#' scores2 <- c(-0.5, -0.3, -0.4, -0.2, -0.3)
#' cs_asymptotic(scores1, scores2, alpha = 0.05)
#'
#' @export
cs_asymptotic <- function(scores1, scores2,
                          alpha  = 0.05,
                          t_star = NULL) {

  stopifnot(
    length(scores1) == length(scores2),
    length(scores1) >= 1,
    alpha > 0, alpha < 1
  )

  xs  <- scores1 - scores2
  T_  <- length(xs)
  ts  <- seq_len(T_)

  # Default: tightest at end of sample
  if (is.null(t_star)) t_star <- T_

  # Tuning parameter rho
  rho <- sqrt((2 * log(1 / alpha) + log(1 + 2 * log(1 / alpha))) / t_star)

  # Running mean hat_Delta_t
  mus <- cumsum(xs) / ts

  # Predictable running mean hat_Delta_{t-1}: lag by 1, initialise at 0
  mus_lag <- c(0, mus[-T_])

  # Running variance estimator using predictable mean
  # sigma2_t = (1/t) * sum_{i=1}^t (hat_delta_i - hat_Delta_{i-1})^2
  sigma2_t <- cumsum((xs - mus_lag)^2) / ts

  # CS half-width (radius)
  # inner = t * sigma2_t * rho^2 + 1
  inner  <- ts * sigma2_t * rho^2 + 1

  # radius = sqrt( (2 * inner / (t^2 * rho^2)) * log(sqrt(inner) / alpha) )
  radius <- sqrt(
    (2 * inner / (ts^2 * rho^2)) *
      log(sqrt(inner) / alpha)
  )

  data.frame(
    t        = ts,
    estimate = mus,
    lower    = mus - radius,
    upper    = mus + radius
  )
}
