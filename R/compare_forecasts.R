# =============================================================================
# compare_forecasts.R
# High-level wrapper for sequential comparison of two probabilistic forecasters
#
# This file provides a user-facing pipeline around the lower-level seqcomp
# building blocks:
#   1. compute pointwise scores,
#   2. compute score differences,
#   3. optionally construct a confidence sequence,
#   4. optionally construct two one-sided e-processes.
#
# The wrapper is intentionally restricted to the Choe & Ramdas framework
# implemented elsewhere in the package.
# =============================================================================

#' Compare Two Sequential Forecasters
#'
#' Computes pointwise scores for two probabilistic forecasters and compares
#' them sequentially using confidence sequences and, when valid finite-sample
#' bounds are available, e-processes.
#'
#' This is a convenience wrapper around [brier_score()], [spherical_score()],
#' [log_score()], [cs_hoeffding()], [cs_bernstein()], [cs_asymptotic()], and
#' [eprocess()]. It is designed for the common workflow where the user has two
#' forecast streams `p` and `q`, an outcome stream `y`, and wants a single tidy
#' output object.
#'
#' @param p Forecasts from forecaster 1. For binary outcomes, a numeric vector
#'   of probabilities for event `y = 1`. For categorical outcomes, a numeric
#'   matrix whose rows are probability vectors.
#' @param q Forecasts from forecaster 2, in the same format as `p`.
#' @param y Outcomes. For binary vector forecasts, a numeric vector in `{0, 1}`.
#'   For categorical matrix forecasts, integer class labels in `{1, ..., K}`.
#' @param scoring_rule Character. Scoring rule used to compare forecasts.
#'   Currently supports `"brier"`, `"spherical"`, and `"log"`.
#' @param alpha Numeric in `(0, 1)`. Significance level. Default is `0.05`.
#' @param cs_type Character or `NULL`. Confidence sequence type:
#'   `"bernstein"`, `"hoeffding"`, `"asymptotic"`, or `"none"`.
#'   If `NULL`, the wrapper uses `"bernstein"` for bounded scoring rules
#'   (`"brier"` and `"spherical"`) and `"asymptotic"` for `"log"`.
#' @param compute_cs Logical. If `TRUE`, compute a confidence sequence.
#'   Default is `TRUE`.
#' @param compute_e Logical. If `TRUE`, compute two one-sided e-processes.
#'   Default is `TRUE`. This is only allowed for bounded score differences
#'   under the current wrapper, namely `"brier"` and `"spherical"`.
#' @param v_opt Numeric > 0. Intrinsic time at which the mixture boundary or
#'   e-process is tuned to be tightest. Default is `10`.
#' @param boundary Character. Boundary type passed to [cs_hoeffding()] or
#'   [cs_bernstein()]. Default is `"mixture"`.
#' @param lcb_only Logical. If `TRUE`, compute a lower one-sided empirical
#'   Bernstein CS. Only used when `cs_type = "bernstein"`.
#' @param ucb_only Logical. If `TRUE`, compute an upper one-sided empirical
#'   Bernstein CS. Only used when `cs_type = "bernstein"`.
#' @param eps Numeric. Probability floor passed to [log_score()] when
#'   `scoring_rule = "log"`. Default is `1e-15`.
#' @param clip_max Numeric. Maximum e-process value before clipping. Passed to
#'   [eprocess()]. Default is `1e7`.
#'
#' @return A `data.frame` with one row per time point and columns:
#' \describe{
#'   \item{`t`}{Time index.}
#'   \item{`score_p`}{Pointwise score of forecaster `p`.}
#'   \item{`score_q`}{Pointwise score of forecaster `q`.}
#'   \item{`delta`}{Pointwise score difference, `score_p - score_q`.}
#'   \item{`estimate`}{Running mean score difference. Positive values favour
#'     forecaster `p`; negative values favour forecaster `q`.}
#'   \item{`lower`, `upper`}{Confidence sequence bounds. These are `NA` if
#'     `compute_cs = FALSE` or `cs_type = "none"`.}
#'   \item{`e_pq`, `e_qp`}{One-sided e-processes. `e_pq` tests whether
#'     forecaster `p` outperforms `q`; `e_qp` tests the reverse direction.
#'     These are `NA` if `compute_e = FALSE`.}
#' }
#'
#' @details
#' All scoring rules in `seqcomp` are positively oriented: higher scores are
#' better. Therefore
#' \deqn{\hat{\delta}_t = S(p_t, y_t) - S(q_t, y_t)}
#' is positive when forecaster `p` performs better than forecaster `q` at time
#' `t`.
#'
#' For `"brier"` and `"spherical"`, score differences are bounded in `[-1, 1]`.
#' The wrapper therefore uses `c = 1` for Hoeffding-style confidence sequences
#' and `c = 2` for empirical Bernstein confidence sequences and e-processes.
#'
#' For `"log"`, score differences are unbounded. The wrapper therefore defaults
#' to [cs_asymptotic()] and refuses to compute finite-sample e-processes. For
#' binary log-score comparisons where the Winkler construction is appropriate,
#' use [winkler_compare()] instead.
#'
#' @section Interpretation:
#' The confidence sequence estimates the running average score advantage of
#' `p` over `q`. If the whole interval lies above zero, the data favour `p`;
#' if the whole interval lies below zero, the data favour `q`.
#'
#' The e-processes are evidence processes for one-sided null hypotheses. At
#' level `alpha`, the two-sided rejection threshold used by [eprocess()] is
#' `2 / alpha`.
#'
#' @examples
#' set.seed(1)
#' y <- rbinom(200, 1, 0.5)
#' p <- rep(0.5, 200)
#' q <- runif(200)
#'
#' out <- compare_forecasts(p, q, y, scoring_rule = "brier")
#' tail(out)
#'
#' @export
compare_forecasts <- function(p, q, y,
                              scoring_rule = c("brier", "spherical", "log"),
                              alpha      = 0.05,
                              cs_type    = NULL,
                              compute_cs = TRUE,
                              compute_e  = TRUE,
                              v_opt      = 10,
                              boundary   = "mixture",
                              lcb_only   = FALSE,
                              ucb_only   = FALSE,
                              eps        = 1e-15,
                              clip_max   = 1e7) {

  scoring_rule <- match.arg(scoring_rule)

  stopifnot(
    alpha > 0, alpha < 1,
    is.logical(compute_cs), length(compute_cs) == 1,
    is.logical(compute_e),  length(compute_e)  == 1,
    v_opt > 0,
    eps >= 0,
    clip_max > 0
  )

  forecast_length <- function(x) {
    if (is.matrix(x)) nrow(x) else length(x)
  }

  T_ <- forecast_length(p)

  if (forecast_length(q) != T_ || length(y) != T_) {
    stop("p, q, and y must have the same number of forecast/outcome rows.")
  }

  if (is.matrix(p) != is.matrix(q)) {
    stop("p and q must have the same shape: both vectors or both matrices.")
  }

  if (is.matrix(p) && ncol(p) != ncol(q)) {
    stop("p and q must have the same number of columns.")
  }

  score_p <- switch(
    scoring_rule,
    "brier"     = brier_score(p, y),
    "spherical" = spherical_score(p, y),
    "log"       = log_score(p, y, eps = eps)
  )

  score_q <- switch(
    scoring_rule,
    "brier"     = brier_score(q, y),
    "spherical" = spherical_score(q, y),
    "log"       = log_score(q, y, eps = eps)
  )

  delta <- score_p - score_q
  estimate <- cumsum(delta) / seq_along(delta)

  bounded_rule <- scoring_rule %in% c("brier", "spherical")

  if (is.null(cs_type)) {
    cs_type <- if (bounded_rule) "bernstein" else "asymptotic"
  } else {
    cs_type <- match.arg(cs_type, c("bernstein", "hoeffding", "asymptotic", "none"))
  }

  lower <- rep(NA_real_, T_)
  upper <- rep(NA_real_, T_)

  if (compute_cs && cs_type != "none") {
    if (cs_type %in% c("bernstein", "hoeffding") && !bounded_rule) {
      stop(
        "Finite-sample Hoeffding/Bernstein confidence sequences require ",
        "bounded score differences. Use cs_type = 'asymptotic' for scoring_rule = 'log'."
      )
    }

    cs <- switch(
      cs_type,

      "bernstein" = cs_bernstein(
        scores1  = score_p,
        scores2  = score_q,
        alpha    = alpha,
        c        = 2,
        v_opt    = v_opt,
        boundary = boundary,
        lcb_only = lcb_only,
        ucb_only = ucb_only
      ),

      "hoeffding" = cs_hoeffding(
        scores1  = score_p,
        scores2  = score_q,
        alpha    = alpha,
        c        = 1,
        v_opt    = v_opt,
        boundary = boundary
      ),

      "asymptotic" = cs_asymptotic(
        scores1 = score_p,
        scores2 = score_q,
        alpha   = alpha
      )
    )

    lower <- cs$lower
    upper <- cs$upper
  }

  e_pq <- rep(NA_real_, T_)
  e_qp <- rep(NA_real_, T_)

  if (compute_e) {
    if (!bounded_rule) {
      stop(
        "e-processes in compare_forecasts() require bounded score differences. ",
        "Use compute_e = FALSE for scoring_rule = 'log', or use winkler_compare() ",
        "for binary log-score comparisons when the Winkler construction applies."
      )
    }

    ep <- eprocess(
      scores1  = score_p,
      scores2  = score_q,
      alpha    = alpha,
      c        = 2,
      v_opt    = v_opt,
      clip_max = clip_max
    )

    e_pq <- ep$e_pq
    e_qp <- ep$e_qp
  }

  out <- data.frame(
    t        = seq_len(T_),
    score_p  = score_p,
    score_q  = score_q,
    delta    = delta,
    estimate = estimate,
    lower    = lower,
    upper    = upper,
    e_pq     = e_pq,
    e_qp     = e_qp
  )

  attr(out, "scoring_rule") <- scoring_rule
  attr(out, "cs_type")      <- cs_type
  attr(out, "alpha")        <- alpha

  out
}
