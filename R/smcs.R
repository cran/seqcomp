# =============================================================================
# smcs.R
# Sequential Model Confidence Sets (SMCS)
#
# Implements multi-model extensions for evaluating an arbitrary number of
# forecasters sequentially, using closure principles and joint confidence sequences.
#
# Implements:
#   vovk_wang_merge()               — Accelerated O(m log m) closed-testing e-value merge
#   smcs_strong()                   — SMCS under the strong and uniformly weak nulls via
#                                     closure principle
#   smcs_weak()                     — SMCS under the weak null via joint confidence sequences
#   build_quantile_betting_arrays() — 3D array builder for quantile betting
#   build_agrapa_betting_array()    — 3D array builder for aGRAPA Betting
#   build_ons_betting_array()       — 3D array builder for ONS-m Betting
#
# By convention, the SMCS maintains a running intersection over time: once a
# model is excluded at time t, it remains excluded for all T > t.
#
# References:
#   A26   Arnold et al. (2026), J R Stat Soc B 88, qkag066
#   VW21  Vovk & Wang (2021), Annals of Statistics 49(3), 1736-1754
#   WSR24  Waudby-Smith & Ramdas (2024), J R Stat Soc B 86(1), 1-27
# =============================================================================

#' Arithmetic-mean closed-testing e-value merge
#'
#' Given a vector of e-values \eqn{e_1, \ldots, e_m} (where each \eqn{e_i} represents
#' the evidence against the intersection null hypothesis for model \eqn{i}), computes
#' the closed-testing adjusted e-values using the arithmetic mean as the e-merging
#' function.
#'
#' This implements the Accelerated E-Value Calibration algorithm (Tim Stephan,
#' ETH Zurich; Section H of the supplementary material to Arnold et al., 2026),
#' which solves the Vovk & Wang (2021) closed-testing minimization in \eqn{O(m \log m)}
#' time rather than the naive \eqn{O(m^2)}.
#'
#' @param e_values Numeric vector of non-negative e-values, one per model.
#'
#' @return A numeric vector of the same length, containing the closed-testing
#'   adjusted e-values \eqn{e^\star_i}, in the original (unsorted) order.
#'
#' @examples
#' raw_evalues <- c(10, 2, 1)
#' # Model 1 has strong evidence against it, Model 3 has none.
#' vovk_wang_merge(raw_evalues)
#'
#' @export
vovk_wang_merge <- function(e_values) {
  m <- length(e_values)
  if (m == 1L) return(e_values)
  if (any(e_values < 0)) stop("e_values must be nonnegative.")

  ord <- order(e_values)
  e_sorted <- e_values[ord]
  S <- cumsum(e_sorted)

  e_star_sorted <- numeric(m)
  e_star_sorted[1] <- e_sorted[1]
  i <- 1L

  for (k in 2:m) {
    cand <- (e_sorted[k] + S[i]) / (i + 1)
    e_star_sorted[k] <- cand
    while (i < k - 1L) {
      i <- i + 1L
      e_temp <- (e_sorted[k] + S[i]) / (i + 1)
      if (e_temp <= e_star_sorted[k]) {
        e_star_sorted[k] <- e_temp
      } else {
        i <- i - 1L
        break
      }
    }
  }

  e_star <- numeric(m)
  e_star[ord] <- e_star_sorted
  e_star
}


# At m = 2, vovk_wang_merge is NOT a no-op.


#' Sequential Model Confidence Set (Strong & Uniformly Weak Null)
#'
#' Constructs a Sequential Model Confidence Set (SMCS) evaluating a family-wise
#' intersection null hypothesis. By maintaining a running intersection over time,
#' any model excluded from the set is permanently eliminated.
#'
#' Depending on the `method` chosen, this function tests different hypotheses:
#' \itemize{
#'   \item \strong{`method = "betting"`}: Tests the \strong{Strong Null} hypothesis
#'     (conditional step-by-step superiority). Uses a product-form betting martingale.
#'   \item \strong{`method = "mixture"`}: Tests the \strong{Uniformly Weak Null}
#'     hypothesis (average superiority over time). Uses an exponential-mixture martingale.
#'     Because strong superiority implies uniform weak superiority, feeding `"mixture"`
#'     into this closed-testing machinery yields a valid (though strictly testing the
#'     uniformly weak null) SMCS.
#' }
#'
#' The function computes pairwise e-processes between all models, constructs an
#' intersection e-process for each model, applies a closed-testing multiplicity
#' adjustment via [vovk_wang_merge()], and permanently excludes models when their
#' adjusted e-value exceeds \eqn{1/\alpha}.
#'
#' @param scores A \eqn{T \times m} matrix of positively-oriented scores.
#' @param alpha Numeric in `(0, 1)`. Family-wise significance level. Default is `0.05`.
#' @param method Character. `"betting"` (uses [eprocess_betting()], true strong null)
#'   or `"mixture"` (uses [eprocess()], tests uniformly weak null). Default is `"betting"`.
#' @param c_param Numeric scalar, \eqn{m \times m} matrix, or \eqn{T \times m \times m} array.
#'   The predictable bound parameter. Required. Time-varying arrays are only allowed
#'   if `method = "betting"`.
#' @param lambda_param Optional parameter for betting fractions, matching the shape
#'   allowed for `c_param`. Only used if `method = "betting"`.
#' @param ... Additional arguments passed to the underlying pairwise e-process function
#'   (e.g., `v_opt` and `clip_max` for `"mixture"`; `clip_max` for `"betting"`).
#'
#' @return A list containing:
#' \describe{
#'   \item{`E_i_dot`}{A \eqn{T \times m} matrix of unadjusted intersection e-processes.}
#'   \item{`E_star`}{A \eqn{T \times m} matrix of closed-testing adjusted e-processes.}
#'   \item{`smcs`}{A \eqn{T \times m} logical matrix. `TRUE` indicates the model remains
#'     in the SMCS at time `t`.}
#' }
#'
#' @examples
#' set.seed(1)
#' # 3 models, 100 time steps. Model 3 is artificially much worse.
#' scores <- matrix(runif(300, -0.5, 0), nrow = 100, ncol = 3)
#' scores[, 3] <- scores[, 3] - 0.5
#' colnames(scores) <- c("M1", "M2", "M3")
#'
#' # Using the betting method with a global bound c = 2
#' res <- smcs_strong(scores, alpha = 0.05, method = "betting", c_param = 2)
#' tail(res$smcs)
#'
#' @export
smcs_strong <- function(scores, alpha = 0.05,
                        method = c("betting", "mixture"),
                        c_param = NULL, lambda_param = NULL, ...) {
  method <- match.arg(method)
  Tt <- nrow(scores)
  m <- ncol(scores)
  if (m < 2) stop("Need at least 2 models.")
  if (is.null(c_param)) stop("c_param must be provided.")

  is_scalar <- length(c_param) == 1
  is_mxm_matrix <- is.matrix(c_param) && nrow(c_param) == m && ncol(c_param) == m

  if (method == "mixture" && !(is_scalar || is_mxm_matrix)) {
    stop(
      "method = 'mixture' requires c_param to be a scalar or an m x m matrix. ",
      "eprocess() only accepts a single scalar sub-exponential scale c; a ",
      "time-varying (length-T vector or T x m x m array) c_param is only ",
      "supported for method = 'betting', via eprocess_betting()'s vectorised c_t."
    )
  }

  E_i_dot <- matrix(0, Tt, m)
  for (i in seq_len(m)) {
    e_sum <- rep(0, Tt)
    for (j in seq_len(m)) {
      if (i == j) next

      c_ij <- .slice_param(c_param, Tt, i, j)
      lam_ij <- .slice_param(lambda_param, Tt, i, j)

      # We need evidence that j beats i (i.e. e_qp where p=i, q=j)
      e_ij <- if (method == "mixture") {
        eprocess(scores[, i], scores[, j], alpha = alpha, c = c_ij, ...)$e_qp
      } else {
        eprocess_betting(scores[, i], scores[, j], c_t = c_ij, lambda_t = lam_ij, ...)$e_qp
      }
      e_sum <- e_sum + e_ij
    }
    E_i_dot[, i] <- e_sum / (m - 1)
  }

  E_star <- matrix(0, Tt, m)
  for (t in seq_len(Tt)) {
    E_star[t, ] <- vovk_wang_merge(E_i_dot[t, ])
  }

  in_set <- E_star < (1 / alpha)
  smcs <- apply(in_set, 2, function(x) as.logical(cummin(x)))
  colnames(smcs) <- colnames(scores)

  list(E_i_dot = E_i_dot, E_star = E_star, smcs = smcs)
}

#' Sequential Model Confidence Set (Weak Null)
#'
#' Constructs a Sequential Model Confidence Set (SMCS) evaluating the weak null
#' hypothesis that a model outperforms all other candidate models on average
#' over time.
#'
#' Uses a joint confidence sequence decoupling result: a model \eqn{i} remains in
#' the SMCS at time \eqn{t} if and only if, for every competitor \eqn{j}, the
#' pairwise \eqn{(1 - \alpha/(m(m-1)))}-confidence sequence for the average score
#' difference does not strictly rule out that \eqn{i} is better than \eqn{j}.
#'
#' Unlike the strong null, this SMCS does not maintain a strict running
#' intersection; a model's average score can recover over time, allowing it to
#' dynamically exit and re-enter the confidence set.
#'
#' @param scores A \eqn{T \times m} matrix of positively-oriented scores.
#' @param alpha Numeric in `(0, 1)`. Family-wise significance level. Default is `0.05`.
#' @param cs_method Character. `"bernstein"` (uses [cs_bernstein()]) or `"hoeffding"`
#'   (uses [cs_hoeffding()]).
#' @param c_param Numeric scalar or \eqn{m \times m} matrix. The uniform bound
#'   parameter. Must be constant over time.
#' @param ... Additional arguments passed to the chosen CS function (e.g., `v_opt`).
#'
#' @return A list containing:
#' \describe{
#'   \item{`smcs`}{A \eqn{T \times m} logical matrix. `TRUE` indicates the model is
#'     in the weakly superior set at time `t`.}
#'   \item{`alpha_adjusted`}{The Bonferroni-adjusted significance level applied to
#'     each pairwise sequence.}
#' }
#'
#' @examples
#' set.seed(2)
#' scores <- matrix(runif(300, -0.5, 0), nrow = 100, ncol = 3)
#' scores[, 3] <- scores[, 3] - 0.5
#' colnames(scores) <- c("M1", "M2", "M3")
#'
#' res <- smcs_weak(scores, alpha = 0.05, cs_method = "bernstein", c_param = 2)
#' tail(res$smcs)
#'
#' @export
smcs_weak <- function(scores, alpha = 0.05,
                      cs_method = c("bernstein", "hoeffding"),
                      c_param = NULL, ...) {
  cs_method <- match.arg(cs_method)
  Tt <- nrow(scores)
  m <- ncol(scores)
  if (m < 2) stop("Need at least 2 models.")
  if (is.null(c_param)) stop("c_param must be provided.")

  is_scalar <- length(c_param) == 1
  is_mxm_matrix <- is.matrix(c_param) && nrow(c_param) == m && ncol(c_param) == m
  if (!(is_scalar || is_mxm_matrix)) {
    stop(
      "smcs_weak() requires c_param to be a scalar or an m x m matrix (constant ",
      "over time); cs_bernstein()/cs_hoeffding() only accept a single scalar c."
    )
  }

  alpha_adj <- alpha / (m * (m - 1))
  cs_fun <- if (cs_method == "bernstein") cs_bernstein else cs_hoeffding

  excluded <- matrix(FALSE, Tt, m)
  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      if (i == j) next

      c_ij <- .slice_param(c_param, Tt, i, j)
      cs <- cs_fun(scores[, i], scores[, j], alpha = alpha_adj, c = c_ij, ...)

      excluded[, i] <- excluded[, i] | (cs$upper < 0)
    }
  }

  smcs <- !excluded
  colnames(smcs) <- colnames(scores)

  list(smcs = smcs, alpha_adjusted = alpha_adj)
}



#' Build 3D Parameter Arrays for Quantile Betting
#'
#' Pre-computes the dynamic, pair-specific bounding arrays (`c_array`) and betting
#' fractions (`lambda_array`) required to evaluate the Sequential Model Confidence
#' Set for quantile forecasts under the strong null hypothesis.
#'
#' @note
#' **Scale Translation:** This function assumes `forecasts` and `scores` are
#' evaluated on the raw, linear scale. To replicate the exact log-scale bounds
#' used in the Arnold et al. (2026) Covid-19 case study, the forecast matrix and
#' outcomes must be log-transformed prior to passing them to this pipeline.
#'
#' @param forecasts A \eqn{T \times m} matrix of raw quantile forecasts.
#' @param scores A \eqn{T \times m} matrix of positively-oriented tick-loss scores.
#' @param tau Numeric scalar in `(0, 1)`. The quantile level.
#' @param eps Positive numeric scalar. Used both as the denominator safeguard
#'   in [lambda_betting_quantile()] and as the minimum predictable bound
#'   returned in `c_array`. Default: `1e-8`.
#'
#' @return A list containing two \eqn{T \times m \times m} numeric arrays: `c_array`
#'   and `lambda_array`.
#'
#' @examples
#' set.seed(3)
#' fcsts <- matrix(runif(150), nrow = 50, ncol = 3)
#' y <- rbinom(50, 1, 0.5)
#' scores <- matrix(0, nrow = 50, ncol = 3)
#' for(i in 1:3) scores[, i] <- tick_loss(fcsts[, i], y, tau = 0.5)
#'
#' arrays <- build_quantile_betting_arrays(fcsts, scores, tau = 0.5)
#' dim(arrays$c_array) # 50 x 3 x 3
#'
#' @export
build_quantile_betting_arrays <- function(forecasts, scores, tau, eps = 1e-8) {
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a finite positive scalar.", call. = FALSE)
  }
  Tt <- nrow(forecasts)
  m <- ncol(forecasts)

  c_array <- array(0, dim = c(Tt, m, m))
  lambda_array <- array(0, dim = c(Tt, m, m))

  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      if (i == j) next

      # Extract score difference (p = i, q = j)
      xs <- scores[, i] - scores[, j]
      delta_hat_lag1 <- c(0, utils::head(xs, -1))

      bnds <- lambda_betting_quantile(
        p_t = forecasts[, i],
        q_t = forecasts[, j],
        tau = tau,
        delta_hat_lag1 = delta_hat_lag1,
        eps = eps
      )

      # `eps` is both the lambda denominator safeguard and the bound floor.
      c_safe <- pmax(bnds$c_t, eps)

      c_array[, i, j] <- c_safe
      lambda_array[, i, j] <- pmin(bnds$lambda_t, 1 / c_safe)
    }
  }
  list(c_array = c_array, lambda_array = lambda_array)
}

#' Build 3D Parameter Array for aGRAPA Betting
#'
#' Pre-computes the dynamic, pair-specific betting fractions required to
#' evaluate the Sequential Model Confidence Set under the strong null
#' hypothesis, using the aGRAPA algorithm.
#'
#' @param scores A \eqn{T \times m} matrix of positively-oriented scores.
#' @param c_mat A numeric scalar, an \eqn{m \times m} matrix, or a
#'   \eqn{T \times m \times m} array of predictable bounds (matching the
#'   `c_param` argument of [smcs_strong()]).
#' @param period Integer. If greater than 1, uses the periodic betting variant
#'   of aGRAPA, which updates the betting fraction only every `period` time steps.
#' @param ... Additional arguments passed to [lambda_betting_agrapa()] (e.g.,
#'   `kappa`, `prior_mean`, `prior_variance`, `fake_obs`).
#'
#' @return A \eqn{T \times m \times m} numeric array of predictable betting
#'   fractions, suitable for passing to the `lambda_param` argument of
#'   [smcs_strong()].
#'
#' @seealso [lambda_betting_agrapa()], [eprocess_betting()], [smcs_strong()]
#'
#' @examples
#' set.seed(42)
#' y <- rbinom(50, 1, 0.5)
#' forecasts <- matrix(runif(150), nrow = 50, ncol = 3)
#' scores <- matrix(0, nrow = 50, ncol = 3)
#' for (i in 1:3) scores[, i] <- brier_score(forecasts[, i], y)
#'
#' # For Brier scores, differences are bounded in [-1, 1], so c = 2
#' lam_array <- build_agrapa_betting_array(scores, c_mat = 2)
#' dim(lam_array) # 50 x 3 x 3
#'
#' @export
build_agrapa_betting_array <- function(scores, c_mat, period = 1, ...) {
  Tt <- nrow(scores)
  m <- ncol(scores)
  lambda_array <- array(0, dim = c(Tt, m, m))

  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      if (i == j) next

      xs_j_beats_i <- scores[, j] - scores[, i]
      c_ij <- .slice_param(c_mat, Tt, i, j)

      if (period > 1) {
        lambda_array[, i, j] <- lambda_betting_periodic(
          xs = xs_j_beats_i, c = c_ij, period = period,
          base_fun = lambda_betting_agrapa, ...
        )
      } else {
        lambda_array[, i, j] <- lambda_betting_agrapa(
          xs = xs_j_beats_i, c = c_ij, ...
        )
      }
    }
  }
  lambda_array
}

#' Build 3D Parameter Array for ONS-m Betting
#'
#' Pre-computes the dynamic, pair-specific betting fractions required to
#' evaluate the Sequential Model Confidence Set under the strong null
#' hypothesis, using the Online Newton Step (ONS-m) algorithm.
#'
#' @param scores A \eqn{T \times m} matrix of positively-oriented scores.
#' @param c_mat A numeric scalar, an \eqn{m \times m} matrix, or a
#'   \eqn{T \times m \times m} array of predictable bounds (matching the
#'   `c_param` argument of [smcs_strong()]).
#' @param period Integer. If greater than 1, uses the periodic betting variant
#'   of ONS-m, which updates the betting fraction only every `period` time steps.
#' @param ... Additional arguments passed to [lambda_betting_ons()] (e.g.,
#'   `eta`).
#'
#' @return A \eqn{T \times m \times m} numeric array of predictable betting
#'   fractions, suitable for passing to the `lambda_param` argument of
#'   [smcs_strong()].
#'
#' @seealso [lambda_betting_ons()], [eprocess_betting()], [smcs_strong()]
#'
#' @examples
#' set.seed(43)
#' y <- rbinom(50, 1, 0.5)
#' forecasts <- matrix(runif(150), nrow = 50, ncol = 3)
#' scores <- matrix(0, nrow = 50, ncol = 3)
#' for (i in 1:3) scores[, i] <- brier_score(forecasts[, i], y)
#'
#' # For Brier scores, differences are bounded in [-1, 1], so c = 2
#' lam_array <- build_ons_betting_array(scores, c_mat = 2)
#' dim(lam_array) # 50 x 3 x 3
#'
#' @export
build_ons_betting_array <- function(scores, c_mat, period = 1, ...) {
  Tt <- nrow(scores)
  m <- ncol(scores)

  lambda_array <- array(0, dim = c(Tt, m, m))

  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      if (i == j) next

      # Optimize bet for j beating i to align with e_qp in smcs_strong
      xs_j_beats_i <- scores[, j] - scores[, i]
      c_ij <- .slice_param(c_mat, Tt, i, j)

      if (period > 1) {
        lambda_array[, i, j] <- lambda_betting_periodic(
          xs = xs_j_beats_i, c = c_ij, period = period,
          base_fun = lambda_betting_ons, ...
        )
      } else {
        lambda_array[, i, j] <- lambda_betting_ons(
          xs = xs_j_beats_i,
          c  = c_ij,
          ...
        )
      }
    }
  }
  lambda_array
}
