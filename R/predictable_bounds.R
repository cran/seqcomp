# =============================================================================
# predictable_bounds.R
# Proposition EC.7: e-process with time-varying predictable bounds
#
# When score differences have bounds that vary over time but are known
# one step ahead (F_{t-1}-measurable), the fixed-c assumption of Theorems
# 2 & 3 can be replaced by a predictable sequence (c_i).
#
# Proposition EC.7 CR24:
#   If |hat_delta_i| <= c_i/2 and (c_i) is strictly positive and
#   F_{i-1}-measurable, then for fixed lambda in [0, 1/c_0):
#
#   E_t(lambda) = prod_{i=1}^t exp{lambda*hat_delta_i - psi_{E,c_i}(lambda)
#                                  * (hat_delta_i - gamma_i)^2}
#
#   is a valid e-process for H_0^w(p, q).
#
# IMPORTANT LIMITATIONS:
#   - No closed-form mixture exists when c_i varies over time.
#   - Lambda must be fixed (no mixture over lambda).
#   - c_0 = sup_i(c_i) governs the domain: lambda in [0, 1/c_0).
#   - If c_0 = Inf, lambda must be 0 (trivial test). In practice,
#     c_i must be bounded from above for the test to have power.
#
# References:
#   CR24  Choe & Ramdas (2024), Operations Research 72(4), Proposition EC.7
# =============================================================================

#' Fixed-lambda e-process with predictable bounds (Proposition EC.7)
#'
#' Constructs a valid e-process when score difference bounds vary over time
#' but are predictable (known at time i-1 before observing hat_delta_i).
#'
#' @param scores1  Numeric vector. Scores for forecaster 1.
#' @param scores2  Numeric vector. Scores for forecaster 2.
#' @param c_seq    Numeric vector. Predictable bound sequence (c_i), same
#'                 length as scores1. Must satisfy `|scores1[i]-scores2[i]| <= c_i/2`
#'                 and c_i > 0 for all i.
#' @param lambda   Numeric in `[0, 1/c_0)`. Betting parameter. Must be
#'                 strictly less than 1/c_0 where c_0 = max(c_seq).
#'                 If NULL, uses the recommended default lambda = 0.5/c_0.
#' @param alpha    Numeric in (0,1). Significance level for rejection rule.
#'                 Default: 0.05. Not used in computation, only for API
#'                 consistency. Pass the same value to [predictable_rejections()]
#'                 when evaluating rejection.
#' @param gammas   Numeric vector or NULL. Predictable centering sequence.
#'                 If NULL, constructed as lagged running mean.
#' @param clip_max Numeric. Maximum e-process value. Default: 1e7.
#' @param strict   Logical. If TRUE, any violation of the bound condition at
#'                 any time point will stop execution with an error. If FALSE
#'                 (default), a warning is issued but the e-process is still computed.
#'                 Note that violations invalidate the e-process guarantee,
#'                 so strict = TRUE is recommended for formal inference.
#'
#' @return data.frame with columns:
#'   t, e_pq, e_qp, log_e_pq, log_e_qp, c_seq, lambda_used
#'
#' @section Predictability:
#' The bound sequence `c_seq` (and the centering sequence `gammas`) must be
#' predictable: `c_i` is fixed and known at time `i - 1`, before
#' `scores1[i]`/`scores2[i]` (and hence `hat_delta_i`) are observed —
#' formally, \eqn{c_i} is \eqn{\mathcal{F}_{i-1}}-measurable. A bound chosen
#' after seeing `hat_delta_i` (e.g. derived from the realised data range)
#' invalidates the e-process guarantee, even if it numerically satisfies
#' `|hat_delta_i| <= c_i/2`.
#'
#' @details
#'   The e-process is computed as:
#'   \deqn{\log E_t(\lambda) = \sum_{i=1}^t \Bigl[\lambda\,\hat{\delta}_i
#'         - \psi_{E,c_i}(\lambda)\,(\hat{\delta}_i - \gamma_i)^2\Bigr]}
#'
#'   where
#'   \deqn{\psi_{E,c}(\lambda) = \frac{-\log(1 - c\lambda) - c\lambda}{c^2}}
#'   is evaluated at each step with the current \eqn{c_i}.
#'
#'   LAMBDA CHOICE: lambda = 0.5/c_0 is a conservative default that stays
#'   well within the valid domain `[0, 1/c_0)`. For better power, lambda can
#'   be tuned to the expected signal size, but must never reach 1/c_0.
#'
#'   VALIDITY CHECK: The function verifies |hat_delta_i| <= c_i/2 at each
#'   step and warns if violated. Violations invalidate the e-process guarantee.
#'
#' @examples
#' scores1 <- c(0.10, 0.20, 0.15, 0.25)
#' scores2 <- c(0.05, 0.10, 0.10, 0.20)
#' c_seq <- rep(1, length(scores1))
#' ep <- eprocess_predictable(scores1, scores2, c_seq = c_seq)
#' head(ep)
#'
#' @export
eprocess_predictable <- function(scores1, scores2,
                                 c_seq,
                                 lambda    = NULL,
                                 alpha     = 0.05,
                                 gammas    = NULL,
                                 clip_max  = 1e7,
                                 strict    = FALSE) {

  T_ <- length(scores1)
  stopifnot(
    length(scores2) == T_,
    length(c_seq)   == T_,
    all(c_seq > 0),
    alpha > 0, alpha < 1,
    clip_max > 0
  )

  xs  <- scores1 - scores2
  c_0 <- max(c_seq)

  # Validity check: |hat_delta_i| <= c_i/2 for all i
  violations <- which(abs(xs) > c_seq / 2 + 1e-8)
  valid_bound <- length(violations) == 0

  if (!valid_bound) {
    msg <- paste0(
      length(violations), " observation(s) violate |hat_delta_i| <= c_i/2. ",
      "First violation at t = ", violations[1], ". ",
      "The e-process validity guarantee is compromised. ",
      "Consider increasing c_seq at those time points."
    )
    if (strict) stop(msg) else warning(msg)
  }

  if (is.null(lambda)) {
    if (is.infinite(c_0))
      stop("c_0 = sup(c_seq) is infinite. Cannot construct a non-trivial e-process.",
           " Ensure c_seq is bounded from above.")
    lambda <- 0.5 / c_0
  }

  stopifnot(lambda >= 0, lambda < 1 / c_0)

  if (is.null(gammas)) {
    gammas <- make_gammas(xs, lag = 1)
  } else {
    stopifnot(length(gammas) == T_)
  }

  log_e_pq <- log_eprocess_fixed_predictable(xs, c_seq, lambda, gammas)
  log_e_qp <- log_eprocess_fixed_predictable(-xs, c_seq, lambda, -gammas)

  log_e_pq_clipped <- pmin(log_e_pq, log(clip_max))
  log_e_qp_clipped <- pmin(log_e_qp, log(clip_max))

  e_pq <- pmax(exp(log_e_pq_clipped), .Machine$double.eps)
  e_qp <- pmax(exp(log_e_qp_clipped), .Machine$double.eps)

  data.frame(
    t           = seq_len(T_),
    e_pq        = e_pq,
    e_qp        = e_qp,
    log_e_pq    = log(e_pq),
    log_e_qp    = log(e_qp),
    c_seq       = c_seq,
    lambda_used = lambda
  )
}

#' Summarise predictable bounds e-process
#'
#' @param ep    data.frame. Output of eprocess_predictable().
#' @param alpha Numeric. Significance level.
#'
#' @return Named list matching the `eprocess_rejections()` format
#' (`threshold`, `tau_pq`, `tau_qp`, `reject_pq`, `reject_qp`), plus:
#' * `c_range` — range of `c_seq` used.
#' * `lambda` — `lambda` value used.
#'
#' @examples
#' scores1 <- c(0.10, 0.20, 0.15, 0.25)
#' scores2 <- c(0.05, 0.10, 0.10, 0.20)
#' c_seq <- rep(1, length(scores1))
#' ep <- eprocess_predictable(scores1, scores2, c_seq = c_seq)
#' predictable_rejections(ep, alpha = 0.05)
#'
#' @export
predictable_rejections <- function(ep, alpha = 0.05) {
  threshold <- 2 / alpha
  tau_pq    <- which(ep$e_pq >= threshold)
  tau_qp    <- which(ep$e_qp >= threshold)

  list(
    threshold  = threshold,
    tau_pq     = if (length(tau_pq) > 0) tau_pq[1] else NA_integer_,
    tau_qp     = if (length(tau_qp) > 0) tau_qp[1] else NA_integer_,
    reject_pq  = length(tau_pq) > 0,
    reject_qp  = length(tau_qp) > 0,
    c_range    = range(ep$c_seq),
    lambda     = ep$lambda_used[1]
  )
}
