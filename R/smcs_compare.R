# =============================================================================
# smcs_compare.R
# High-level wrapper for sequential comparison of multiple probabilistic forecasters
#
# This file provides a user-facing pipeline around the multi-model SMCS building
# blocks. It abstracts away the matrix/array wrangling and directly returns
# inclusion matrices for the evaluated models over time.
#
# 1. Computes pointwise scores for all models,
# 2. Automatically constructs bounding arrays for conditionally bounded rules (e.g. tick),
# 3. Constructs the Sequential Model Confidence Set (SMCS) under the strong null,
# 4. Constructs the SMCS under the weak null (for uniformly bounded rules).
# =============================================================================

#' Compare Multiple Sequential Forecasters (SMCS)
#'
#' Evaluates an arbitrary number of candidate forecasters simultaneously,
#' constructing Sequential Model Confidence Sets (SMCS) that maintain
#' family-wise error rate control over time.
#'
#' This is a high-level wrapper that automates pointwise score calculation,
#' boundary generation, and multiplicity corrections via [smcs_strong()] and
#' [smcs_weak()]. For uniformly bounded scoring rules, it returns SMCSs under
#' the strong, uniformly weak, and time-varying weak null hypotheses. For
#' conditionally bounded rules like `"tick"` loss, it automatically builds the
#' dynamic 3D arrays required for strong-null adaptive betting.
#'
#' @param forecasts A \eqn{T \times m} matrix of forecasts. For binary/categorical
#'   probability forecasts, these should be matrices of probabilities. For quantile
#'   forecasts, these should be raw predicted quantiles.
#'   *Note:* to replicate the log-scale bounds of the Arnold et al. (2026)
#'   Covid-19 study, pass log-transformed forecasts and outcomes.
#' @param outcomes A numeric vector of \eqn{T} realised outcomes.
#' @param scoring_rule Character. Scoring rule used to compare forecasts.
#'   Currently supports `"brier"`, `"spherical"`, and `"tick"`.
#' @param cs_method Character. Confidence sequence method for the weak null:
#'   `"bernstein"` or `"hoeffding"`. Default is `"bernstein"`.
#' @param tau Numeric in `(0, 1)`. The quantile level. Required only if
#'   `scoring_rule = "tick"`.
#' @param alpha Numeric in `(0, 1)`. Family-wise significance level. Default is `0.05`.
#' @param v_opt Numeric > 0. Intrinsic time at which the weak-null confidence sequence
#'   is tuned to be tightest. Default is `10`.
#' @param clip_max Numeric. Maximum e-process value before clipping in the strong-null
#'   test. Default is `1e7`.
#' @param betting_rule Character. Strong-null betting-fraction rule. `"default"`
#'   preserves the existing behavior: `"naive"` for `"brier"` and `"spherical"`,
#'   and `"arnold"` for `"tick"`. `"naive"`, `"agrapa"`, and `"ons"` are
#'   available for all scoring rules; `"arnold"` is available only for `"tick"`.
#' @param period Positive integer. Number of interleaved periodic sub-streams
#'   used by `"agrapa"` or `"ons"`. Must be `1` for `"naive"` and `"arnold"`.
#' @param ... Additional arguments passed to [build_agrapa_betting_array()] or
#'   [build_ons_betting_array()] for the selected adaptive rule, such as
#'   `kappa`, `prior_mean`, `prior_variance`, `fake_obs`, or `eta`.
#'
#' @return A list containing:
#' \describe{
#'   \item{`scores`}{A \eqn{T \times m} matrix of evaluated pointwise scores.}
#'   \item{`smcs_strong`}{A \eqn{T \times m} logical matrix tracking inclusion in the
#'     strong-null SMCS over time (permanent exclusions).}
#'   \item{`smcs_uniform_weak`}{A \eqn{T \times m} logical matrix tracking inclusion
#'     in the uniformly weak SMCS over time (permanent exclusions). Currently `NULL`
#'     for `"tick"` loss.}
#'   \item{`smcs_weak`}{A \eqn{T \times m} logical matrix tracking inclusion in the
#'     time-varying weak-null SMCS over time (models can exit and re-enter).
#'     Currently `NULL` for `"tick"` loss.}
#'   \item{`betting_rule`}{The resolved strong-null betting rule actually used.}
#'   \item{`period`}{The period supplied for an adaptive betting rule.}
#'   \item{`betting_args`}{The additional adaptive-rule settings supplied through
#'     `...`.}
#' }
#'
#' @examples
#' set.seed(42)
#' T_sim <- 100
#' y <- rbinom(T_sim, 1, 0.5)
#'
#' # Create 3 forecasters:
#' # M1 is a perfect oracle (always predicts the true y)
#' # M2 is slightly noisy (adds small uniform noise to y)
#' # M3 is an anti-oracle (predicts the exact opposite of y)
#' fcsts <- matrix(NA, nrow = T_sim, ncol = 3)
#' fcsts[, 1] <- y
#' fcsts[, 2] <- abs(y - runif(T_sim, 0, 0.1))
#' fcsts[, 3] <- 1 - y
#' colnames(fcsts) <- c("M1", "M2", "M3")
#'
#' out <- smcs_compare(fcsts, y, scoring_rule = "brier")
#'
#' # Print the object to see exclusions (M3 will be dropped rapidly)
#' out
#'
#' # View how the set sizes shrink over time
#' summary(out)
#'
#' # Additionally, use an adaptive strong-null betting rule
#' out_agrapa <- smcs_compare(
#'   fcsts, y, scoring_rule = "brier", betting_rule = "agrapa"
#' )
#'
#' out_agrapa
#'
#' @export
smcs_compare <- function(forecasts, outcomes,
                         scoring_rule = c("brier", "spherical", "tick"),
                         cs_method = c("bernstein", "hoeffding"),
                         tau = NULL, alpha = 0.05,
                         v_opt = 10, clip_max = 1e7,
                         betting_rule = c("default", "naive", "agrapa", "ons", "arnold"),
                         period = 1, ...) {
  scoring_rule <- match.arg(scoring_rule)
  cs_method <- match.arg(cs_method)
  betting_rule <- match.arg(betting_rule)
  betting_args <- list(...)

  if (!is.numeric(period) || length(period) != 1L || !is.finite(period) ||
      period < 1 || period != floor(period) ||
      period > .Machine$integer.max) {
    stop("period must be a positive integer.", call. = FALSE)
  }
  period <- as.integer(period)

  if (betting_rule == "default") {
    betting_rule <- if (scoring_rule == "tick") "arnold" else "naive"
  }

  if (betting_rule == "arnold" && scoring_rule != "tick") {
    stop(
      "betting_rule = 'arnold' is available only when scoring_rule = 'tick'.",
      call. = FALSE
    )
  }

  if (betting_rule %in% c("naive", "arnold") && period != 1L) {
    stop(
      "period is available only when betting_rule is 'agrapa' or 'ons'.",
      call. = FALSE
    )
  }

  if (betting_rule %in% c("naive", "arnold") && length(betting_args) > 0L) {
    stop(
      "Additional betting arguments are available only when betting_rule is 'agrapa' or 'ons'.",
      call. = FALSE
    )
  }

  Tt <- nrow(forecasts)
  m <- ncol(forecasts)
  stopifnot(length(outcomes) == Tt)

  scores_mat <- matrix(0, nrow = Tt, ncol = m)
  colnames(scores_mat) <- colnames(forecasts)
  if (is.null(colnames(scores_mat))) {
    colnames(scores_mat) <- paste0("Model_", seq_len(m))
  }

  for (i in seq_len(m)) {
    if (scoring_rule == "brier") {
      scores_mat[, i] <- brier_score(forecasts[, i], outcomes)
    } else if (scoring_rule == "spherical") {
      scores_mat[, i] <- spherical_score(forecasts[, i], outcomes)
    } else if (scoring_rule == "tick") {
      if (is.null(tau)) stop("tau must be provided for tick_loss.")
      scores_mat[, i] <- tick_loss(forecasts[, i], outcomes, tau)
    }
  }

  if (scoring_rule %in% c("brier", "spherical")) {
    c_param_strong <- 2
    c_param_weak <- if (cs_method == "hoeffding") 1 else 2

    lambda_param <- if (betting_rule == "naive") {
      NULL
    } else if (betting_rule == "agrapa") {
      build_agrapa_betting_array(
        scores_mat, c_mat = c_param_strong, period = period, ...
      )
    } else if (betting_rule == "ons") {
      build_ons_betting_array(
        scores_mat, c_mat = c_param_strong, period = period, ...
      )
    } else {
      stop("Internal error: unsupported betting rule.", call. = FALSE)
    }

    smcs_s <- smcs_strong(
      scores_mat, alpha = alpha, method = "betting",
      c_param = c_param_strong, lambda_param = lambda_param,
      clip_max = clip_max
    )
    smcs_uw <- smcs_strong(
      scores_mat, alpha = alpha, method = "mixture",
      c_param = c_param_strong, v_opt = v_opt, clip_max = clip_max
    )
    smcs_w <- smcs_weak(
      scores_mat, alpha = alpha, cs_method = cs_method,
      c_param = c_param_weak, v_opt = v_opt
    )

  } else if (scoring_rule == "tick") {
    bnds <- build_quantile_betting_arrays(forecasts, scores_mat, tau)

    lambda_param <- if (betting_rule == "naive") {
      NULL
    } else if (betting_rule == "arnold") {
      bnds$lambda_array
    } else if (betting_rule == "agrapa") {
      build_agrapa_betting_array(
        scores_mat, c_mat = bnds$c_array, period = period, ...
      )
    } else if (betting_rule == "ons") {
      build_ons_betting_array(
        scores_mat, c_mat = bnds$c_array, period = period, ...
      )
    } else {
      stop("Internal error: unsupported betting rule.", call. = FALSE)
    }

    smcs_s <- smcs_strong(
      scores_mat, alpha = alpha, method = "betting",
      c_param = bnds$c_array, lambda_param = lambda_param,
      clip_max = clip_max
    )

    warning(
      "smcs_uniform_weak and smcs_weak are currently omitted for unbounded ",
      "tick loss in this wrapper (require transformation)."
    )
    smcs_uw <- NULL
    smcs_w <- NULL
  }

  result <- list(
    scores            = scores_mat,
    smcs_strong       = smcs_s$smcs,
    smcs_uniform_weak = if (!is.null(smcs_uw)) smcs_uw$smcs else NULL,
    smcs_weak         = if (!is.null(smcs_w)) smcs_w$smcs else NULL,
    alpha             = alpha,
    scoring_rule      = scoring_rule,
    betting_rule      = betting_rule,
    period            = period,
    betting_args      = betting_args
  )
  class(result) <- "seqcomp_multi"
  return(result)
}

#' Print method for seqcomp_multi objects
#'
#' @param x A `seqcomp_multi` object.
#' @param ... Additional arguments passed to print.
#' @export
#' @noRd
print.seqcomp_multi <- function(x, ...) {
  Tt <- nrow(x$scores)
  m  <- ncol(x$scores)
  model_names <- colnames(x$scores)
  if (is.null(model_names)) model_names <- paste0("Model_", seq_len(m))

  betting_rule <- if (is.null(x$betting_rule)) "not recorded" else x$betting_rule

  cat(sprintf("<seqcomp multi-model comparison>\n"))
  cat(sprintf(
    "  %d models, %d time steps, alpha = %.3g, scoring rule = '%s', strong betting rule = '%s'\n\n",
    m, Tt, x$alpha, x$scoring_rule, betting_rule
  ))

  first_excluded <- function(mat) {
    apply(mat, 2, function(col) {
      idx <- which(!col)
      if (length(idx) == 0) NA_integer_ else idx[1]
    })
  }

  summarise_one <- function(mat, label, is_weak = FALSE) {
    if (is.null(mat)) return(invisible(NULL))
    fe <- first_excluded(mat)

    status <- ifelse(mat[Tt, ], "included", "excluded")

    # Highlight re-entry for the weak null
    if (is_weak) {
      reentered <- mat[Tt, ] & !is.na(fe)
      status[reentered] <- "included (re-entered)"
    }

    df <- data.frame(
      model         = model_names,
      status_at_T   = status,
      first_dropped = ifelse(is.na(fe), "-", fe)
    )

    cat(sprintf("-- %s --\n", label))
    print(df, row.names = FALSE)
    cat(sprintf("   final set size: %d -> %d\n\n", m, sum(mat[Tt, ])))
  }

  summarise_one(x$smcs_strong, "Strong-null SMCS (Permanent Exclusion)", is_weak = FALSE)
  summarise_one(
    x$smcs_uniform_weak,
    "Uniformly-weak SMCS (Permanent Exclusion)",
    is_weak = FALSE
  )
  summarise_one(x$smcs_weak,   "Weak-null SMCS (Dynamic Re-entry)", is_weak = TRUE)

  cat(
    "Use `x$smcs_strong`, `x$smcs_uniform_weak`, or `x$smcs_weak` for full ",
    "inclusion matrices, or `x$scores` for pointwise scores.\n"
  )
  invisible(x)
}

#' Summary method for seqcomp_multi objects
#'
#' @param object A `seqcomp_multi` object.
#' @param checkpoints Numeric vector of quantiles for the timeline.
#' @param ... Additional arguments.
#' @export
#' @noRd
summary.seqcomp_multi <- function(object, checkpoints = c(0, 0.25, 0.5, 0.75, 1), ...) {
  Tt  <- nrow(object$scores)
  idx <- unique(pmax(1, round(checkpoints * Tt)))

  set_size <- function(mat) {
    if (is.null(mat)) rep(NA_integer_, length(idx)) else rowSums(mat[idx, , drop = FALSE])
  }

  df <- data.frame(
    t                     = idx,
    strong_set_size       = set_size(object$smcs_strong),
    uniform_weak_set_size = set_size(object$smcs_uniform_weak),
    weak_set_size         = set_size(object$smcs_weak)
  )
  return(df)
}
