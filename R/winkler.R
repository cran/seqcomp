# =============================================================================
# winkler.R
# Proposition 4 pipeline: Winkler score extension for unbounded scoring rules
#
# Choe & Ramdas (2023), Proposition 4:
#   For a proper scoring rule S on binary outcomes, the normalised Winkler
#   score w(p, q, y) is bounded above by 1 and satisfies the sub-exponential
#   condition required for Theorems 2 and 3 with scale c = 2.
#
# This file provides:
#   winkler_cs()   — one-sided EB confidence sequence for Winkler scores
#   winkler_etest() — e-process for Winkler scores
#   winkler_compare() — combined pipeline returning both
#
# IMPORTANT RESTRICTIONS (from Proposition 4 and Section G of CR23):
#   1. Binary outcomes only: y in {0, 1}, p, q in (0, 1).
#   2. The CS is ONE-SIDED (upper bound only) because the Winkler score is
#      unbounded below. Output: (-Inf, U_t].
#   3. Two-sided CS requires a finite analytical lower bound (Corollary 2),
#      which must be derived per scoring rule and forecaster pair.
#   4. NOT applicable to QLIKE or other continuous-outcome scoring rules.
#
# References:
#   CR23  Choe & Ramdas (2023), Operations Research 72(4)
# =============================================================================

#' One-sided empirical Bernstein CS for Winkler scores (Proposition 4)
#'
#' Applies the Winkler normalisation and constructs a one-sided upper
#' confidence sequence for the mean Winkler score W_t = (1/t)*sum w_i.
#' The CS takes the form `(-Inf, U_t]`, valid uniformly over all t >= 1.
#'
#' @param p          Numeric vector in (0,1). Forecasts from model 1.
#' @param q          Numeric vector in (0,1). Forecasts from model 2.
#' @param y          Numeric vector containing only 0 and 1. Binary outcomes.
#' @param alpha      Numeric in (0,1). Significance level. Default: 0.05.
#' @param base_score Function. Underlying scoring rule. Default: log_score.
#' @param v_opt      Numeric > 0. Optimal intrinsic time. Default: 10.
#' @param lower_bound Numeric or NULL. Analytical lower bound on w_i for
#'                   two-sided CS via Corollary 2. If NULL (default), returns
#'                   one-sided CS only. If supplied, must satisfy w_i >= lower_bound
#'                   for all i almost surely.
#'
#' @return data.frame with columns t, estimate, lower, upper.
#'         lower = -Inf always (one-sided) unless lower_bound is supplied.
#'
#' @section Interpretation:
#' If `U_t < 0` for some `t`, this is time-uniform evidence that forecaster 1
#' (`p`) is worse than forecaster 2 (`q`) on average — i.e. a rejection is
#' evidence against `p`, not for it. More generally, `W_t > 0` suggests `p`
#' outperforms `q`; `W_t < 0` suggests `q` outperforms `p`.
#'
#' @details
#' Scale convention: Winkler score bounded above by 1, so c/2 = 1, c = 2.
#' This is hardcoded — do not change c without re-deriving the bound.
#'
#' @examples
#' p <- c(0.7, 0.6, 0.8, 0.65)
#' q <- c(0.5, 0.7, 0.6, 0.55)
#' y <- c(1, 1, 0, 1)
#' winkler_cs(p, q, y, alpha = 0.05)
#'
#' @export
winkler_cs <- function(p, q, y,
                       alpha       = 0.05,
                       base_score  = log_score,
                       v_opt       = 10,
                       lower_bound = NULL) {

  stopifnot(
    length(p) == length(q), length(p) == length(y),
    all(p > 0 & p < 1), all(q > 0 & q < 1),
    all(y %in% c(0, 1)),
    alpha > 0, alpha < 1
  )

  # Step 1: compute Winkler scores
  ws <- winkler_score(p, q, y, base_score = base_score)

  # Step 2: scale parameter — Winkler bounded above by 1, so c = 2
  c_wink <- 2

  if (is.null(lower_bound)) {
    # One-sided upper CS: (-Inf, U_t]
    # ucb_only = TRUE in cs_bernstein returns upper bound only
    cs <- cs_bernstein(
      scores1  = ws,
      scores2  = rep(0, length(ws)),
      alpha    = alpha,
      c        = c_wink,
      v_opt    = v_opt,
      ucb_only = TRUE
    )
  } else {
    # Two-sided CS via a finite lower bound.
    # For raw Winkler scores in [lower_bound, 1], our cs_bernstein() wrapper
    # expects a symmetric bound of the form |x_i| <= c/2.
    # A conservative compatible choice is therefore:
    #   c/2 >= max(1, -lower_bound)
    # i.e. c = 2 * max(1, -lower_bound).
    stopifnot(is.numeric(lower_bound), length(lower_bound) == 1, lower_bound < 1)

    if (any(ws < lower_bound - 1e-8)) {
      warning(
        "Some Winkler scores fall below the supplied lower_bound (",
        lower_bound, "). The two-sided CS validity may be compromised."
      )
    }

    c_twosided <- 2 * max(1, -lower_bound)

    cs <- cs_bernstein(
      scores1 = ws,
      scores2 = rep(0, length(ws)),
      alpha   = alpha,
      c       = c_twosided,
      v_opt   = v_opt
    )
  }

  cs
}

#' E-process for Winkler scores (Proposition 4 + Theorem 3)
#'
#' Tests whether the mean Winkler score W_t >= 0 for all t.
#' A rejection provides time-uniform evidence that forecaster 1 (p)
#' is worse than forecaster 2 (q) under the base scoring rule.
#'
#' @param p         Numeric vector in (0,1). Forecasts from model 1.
#' @param q         Numeric vector in (0,1). Forecasts from model 2.
#' @param y         Numeric vector containing only 0 and 1. Binary outcomes.
#' @param alpha     Numeric in (0,1). Significance level. Default: 0.05.
#' @param base_score Function. Underlying scoring rule. Default: log_score.
#' @param v_opt     Numeric > 0. Default: 10.
#' @param clip_max  Numeric. Maximum e-process value before clipping.
#'                  Default: 1e7.
#'
#' @return data.frame with columns t, e, log_e.
#'
#' @section Rejection rule:
#' Reject at level `alpha` when `e >= 1 / alpha`; this provides time-uniform
#' evidence that `p` is worse than `q`.
#'
#' @examples
#' p <- c(0.7, 0.6, 0.8, 0.65)
#' q <- c(0.5, 0.7, 0.6, 0.55)
#' y <- c(1, 1, 0, 1)
#' winkler_etest(p, q, y, alpha = 0.05)
#'
#' @export
winkler_etest <- function(p, q, y,
                          alpha      = 0.05,
                          base_score = log_score,
                          v_opt      = 10,
                          clip_max   = 1e7) {

  stopifnot(
    length(p) == length(q), length(p) == length(y),
    all(p > 0 & p < 1), all(q > 0 & q < 1),
    all(y %in% c(0, 1)),
    alpha > 0, alpha < 1,
    v_opt > 0,
    clip_max > 0
  )

  ws <- winkler_score(p, q, y, base_score = base_score)

  # Proposition 4 applies because the Winkler score is bounded above by 1.
  # To use the lower-bounded sub-exponential e-process machinery, negate it:
  #   x_i = -w_i >= -1
  xs <- -ws

  # Predictable centering sequence
  gammas <- make_gammas(xs, lag = 1)

  # Shared intrinsic time (no floor for e-processes)
  V_hat <- intrinsic_time(xs, gammas, floor = FALSE)
  S_t   <- cumsum(xs)

  # One-sided tuning parameter
  rho <- rho_from_vopt(v_opt = v_opt, alpha = alpha)

  log_e <- log_ge_mixture_from_sv(S_t, V_hat, rho, c = 2)

  e <- clip_eprocess(log_e, clip_max = clip_max)

  data.frame(
    t     = seq_along(xs),
    e     = e,
    log_e = pmin(log_e, log(clip_max))
  )
}

#' Summarise one-sided Winkler e-process rejections
#'
#' Internal helper used by `winkler_compare()`.
#'
#' @param ep Data frame returned by `winkler_etest()`, with an `e` column.
#' @param alpha Numeric in `(0, 1)`. Significance level.
#'
#' @return A list with elements `threshold`, `tau`, and `reject`.
#'
#' @keywords internal
#' @noRd
winkler_etest_rejections <- function(ep, alpha = 0.05) {
  threshold <- 1 / alpha
  tau <- which(ep$e >= threshold)

  list(
    threshold = threshold,
    tau       = if (length(tau) > 0) tau[1] else NA_integer_,
    reject    = length(tau) > 0
  )
}

#' Full Winkler comparison pipeline (Proposition 4)
#'
#' Convenience wrapper that computes Winkler scores, one-sided CS, and
#' e-process in a single call.
#'
#' @param p          Numeric vector in (0,1).
#' @param q          Numeric vector in (0,1).
#' @param y          Numeric vector containing only 0 and 1. Binary outcomes.
#' @param alpha      Numeric in (0,1). Default: 0.05.
#' @param base_score Function. Default: log_score.
#' @param v_opt      Numeric > 0. Default: 10.
#' @param lower_bound Numeric or NULL. See winkler_cs().
#'
#' @return Named list with elements:
#' * `winkler_scores` — raw Winkler score vector.
#' * `cs` — `data.frame` from `winkler_cs()`.
#' * `etest_p_worse` — one-sided e-process testing whether `p` is worse than `q`.
#' * `etest_q_worse` — one-sided e-process testing whether `q` is worse than `p`.
#' * `rejections` — list of one-sided rejection summaries.
#'
#' @examples
#' p <- c(0.7, 0.6, 0.8, 0.65)
#' q <- c(0.5, 0.7, 0.6, 0.55)
#' y <- c(1, 1, 0, 1)
#' winkler_compare(p, q, y, alpha = 0.05)
#'
#' @export
winkler_compare <- function(p, q, y,
                            alpha       = 0.05,
                            base_score  = log_score,
                            v_opt       = 10,
                            lower_bound = NULL) {

  ws_pq <- winkler_score(p, q, y, base_score = base_score)
  cs    <- winkler_cs(p, q, y, alpha = alpha, base_score = base_score,
                      v_opt = v_opt, lower_bound = lower_bound)

  ep_p_worse <- winkler_etest(
    p = p, q = q, y = y,
    alpha = alpha, base_score = base_score, v_opt = v_opt
  )

  ep_q_worse <- winkler_etest(
    p = q, q = p, y = y,
    alpha = alpha, base_score = base_score, v_opt = v_opt
  )

  rej_p_worse <- winkler_etest_rejections(ep_p_worse, alpha = alpha)
  rej_q_worse <- winkler_etest_rejections(ep_q_worse, alpha = alpha)

  list(
    winkler_scores = ws_pq,
    cs             = cs,
    etest_p_worse  = ep_p_worse,
    etest_q_worse  = ep_q_worse,
    rejections     = list(
      p_worse_than_q = rej_p_worse,
      q_worse_than_p = rej_q_worse
    )
  )
}
