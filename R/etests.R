# =============================================================================
# etests.R
# Sequential e-processes for testing the weak null hypothesis
#
# Implements:
#   eprocess()  — Theorem 3 CR23: sub-exponential mixture e-process
#
# The weak null being tested:
#   H_0^w(p, q): Delta_t = (1/t) * sum_{i=1}^t E[hat_delta_i | F_{i-1}] <= 0
#   for all t = 1, 2, ...
#   i.e. forecaster 1 is no better than forecaster 2 on average.
#
# Two one-sided e-processes are run simultaneously:
#   e_pq: tests H_0^w(p, q) — is p better than q?
#   e_qp: tests H_0^w(q, p) — is q better than p?
#
# Rejection rule (two-sided, level alpha):
#   Reject H_0^w(p, q) if e_pq >= 2/alpha  [p outperforms q]
#   Reject H_0^w(q, p) if e_qp >= 2/alpha  [q outperforms p]
#
# Returns a data.frame with columns:
#   t      — time index
#   e_pq   — e-process value for H_0^w(p, q)
#   e_qp   — e-process value for H_0^w(q, p)
#   log_e_pq, log_e_qp — log-scale values (for plotting)
#
# References:
#   CR23  Choe & Ramdas (2023), Operations Research 72(4)
#   H21   Howard et al. (2021), Annals of Statistics 49(2)
# =============================================================================

#' Sub-exponential mixture e-process (Theorem 3, Choe & Ramdas 2023)
#'
#' Constructs two simultaneous one-sided e-processes for sequentially testing
#' whether forecaster 1 (p) outperforms forecaster 2 (q) or vice versa.
#'
#' The mixture e-process at time t is:
#' \deqn{E_t^{\mathrm{mix}} = m(S_t, \hat{V}_t)}
#' where \eqn{S_t = \sum_{i=1}^t \hat{\delta}_i},
#' \eqn{\hat{V}_t = \sum_{i=1}^t (\hat{\delta}_i - \gamma_i)^2},
#' and \eqn{m(s, v)} is the Gamma-Exponential mixture function (Proposition 3, CR23).
#'
#' @param scores1   Numeric vector. Scores S(p_t, y_t) for forecaster 1.
#' @param scores2   Numeric vector. Scores S(q_t, y_t) for forecaster 2.
#' @param alpha     Numeric in (0,1). Significance level. Rejection threshold
#'                  is 2/alpha for the two-sided test. Default: 0.05.
#' @param c         Numeric > 0. Sub-exponential scale. Must satisfy
#'                  |hat_delta_i| <= c/2 for all i.
#'                  For score differences in `[-(b-a), b-a]`: c = 2*(b-a).
#'                  Default: 2 (for Brier score differences in `[-1,1]`).
#' @param v_opt     Numeric > 0. Intrinsic time at which e-process grows
#'                  fastest. Default: 10 (recommended by CR23).
#' @param alpha_opt Numeric in (0,1). One-sided alpha used to compute rho.
#'                  Default: alpha/2 (matches comparecast two-sided convention).
#' @param gammas    Numeric vector or NULL. Predictable centering sequence.
#'                  If NULL, constructed as lagged running mean.
#' @param clip_max  Numeric. Maximum e-process value before clipping.
#'                  Default: 1e7 (matches Python comparecast).
#'
#' @return data.frame with columns t, e_pq, e_qp, log_e_pq, log_e_qp.
#'
#' @section Rejection rule:
#' At level `alpha`: reject \eqn{H_0^w(p, q)} (conclude `p` outperforms `q`)
#' when `e_pq >= 2/alpha`; reject \eqn{H_0^w(q, p)} (conclude `q` outperforms
#' `p`) when `e_qp >= 2/alpha`. Use `eprocess_rejections()` to extract the
#' first crossing time for each.
#'
#' @details
#'   VARIANCE PROCESS: The intrinsic time V_hat_t uses NO floor (unlike the
#'   EB CS). The GE mixture m(s, v) is well-defined at v=0 (returns 1 when
#'   s=0), so no floor is needed. Adding a floor would distort e-values.
#'
#'   SCALE CONVENTION: c is the sub-exponential scale parameter such that
#'   |hat_delta_i| <= c/2. This is the Theorems 2 & 3 convention from CR23.
#'   For Brier score differences in `[-1,1]`: c = 2.
#'   For Winkler scores (bounded above by 1): c = 2.
#'
#'   LOG-SPACE: E-process values are computed in log-space and clipped before
#'   exponentiating to avoid numerical overflow.
#'
#' @examples
#' scores1 <- c(-0.04, -0.09, -0.01, -0.16)
#' scores2 <- c(-0.09, -0.16, -0.04, -0.25)
#' ep <- eprocess(scores1, scores2, alpha = 0.05)
#' head(ep)
#'
#' @export
eprocess <- function(scores1, scores2,
                     alpha     = 0.05,
                     c         = 2,
                     v_opt     = 10,
                     alpha_opt = NULL,
                     gammas    = NULL,
                     clip_max  = 1e7) {

  stopifnot(
    length(scores1) == length(scores2),
    length(scores1) >= 1,
    alpha > 0, alpha < 1,
    c > 0,
    v_opt > 0,
    clip_max > 0
  )

  # Default alpha_opt: half of alpha for two-sided test
  if (is.null(alpha_opt)) alpha_opt <- alpha / 2

  xs <- scores1 - scores2
  T_ <- length(xs)

  # Predictable centering sequence
  if (is.null(gammas)) {
    gammas <- make_gammas(xs, lag = 1)
  } else {
    stopifnot(length(gammas) == T_)
  }

  # Tuning parameter rho from v_opt and alpha_opt
  rho <- rho_from_vopt(v_opt = v_opt, alpha = alpha_opt)

  # Compute shared variance process (sign-invariant)
  V_shared <- intrinsic_time(xs, gammas, floor = FALSE)
  S_pq     <- cumsum(xs)
  S_qp     <- -S_pq

  log_e_pq <- log_ge_mixture_from_sv(S_pq, V_shared, rho, c)
  log_e_qp <- log_ge_mixture_from_sv(S_qp, V_shared, rho, c)

  # Clip and exponentiate
  e_pq <- clip_eprocess(log_e_pq, clip_max = clip_max)
  e_qp <- clip_eprocess(log_e_qp, clip_max = clip_max)

  data.frame(
    t        = seq_len(T_),
    e_pq     = e_pq,
    e_qp     = e_qp,
    log_e_pq = pmin(log_e_pq, log(clip_max)),
    log_e_qp = pmin(log_e_qp, log(clip_max))
  )
}

#' Determine rejection times for an e-process output
#'
#' @param ep      data.frame. Output of eprocess().
#' @param alpha   Numeric. Significance level. Threshold is 2/alpha.
#'
#' @return Named list with elements:
#' * `threshold` — rejection threshold (`2 / alpha`).
#' * `tau_pq` — first `t` where `e_pq >= threshold` (`NA` if never crossed).
#' * `tau_qp` — first `t` where `e_qp >= threshold` (`NA` if never crossed).
#' * `reject_pq` — logical: was \eqn{H_0^w(p,q)} ever rejected?
#' * `reject_qp` — logical: was \eqn{H_0^w(q,p)} ever rejected?
#'
#' @examples
#' scores1 <- c(-0.04, -0.09, -0.01, -0.16)
#' scores2 <- c(-0.09, -0.16, -0.04, -0.25)
#' ep <- eprocess(scores1, scores2, alpha = 0.05)
#' eprocess_rejections(ep, alpha = 0.05)
#'
#' @export
eprocess_rejections <- function(ep, alpha = 0.05) {
  threshold <- 2 / alpha

  tau_pq <- which(ep$e_pq >= threshold)
  tau_qp <- which(ep$e_qp >= threshold)

  list(
    threshold  = threshold,
    tau_pq     = if (length(tau_pq) > 0) tau_pq[1] else NA_integer_,
    tau_qp     = if (length(tau_qp) > 0) tau_qp[1] else NA_integer_,
    reject_pq  = length(tau_pq) > 0,
    reject_qp  = length(tau_qp) > 0
  )
}
