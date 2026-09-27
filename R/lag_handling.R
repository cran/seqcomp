# =============================================================================
# lag_handling.R
# Propositions 5 & 6: sequential forecast comparison for lag-h forecasts
#
# When forecasts are made h steps ahead (h >= 2), the score difference
# process is no longer a martingale difference sequence with respect to
# the natural filtration. Stream splitting partitions the data into h
# interleaved subsequences, each of which IS a valid martingale sequence,
# then combines the resulting e-processes.
#
#   This module is implemented for theoretical completeness and to cover
#   Propositions 5 & 6 in the thesis.
#
# Two combination rules:
#   "pw"  — period-wise weak null: average e-processes across streams
#   "w"   — standard weak null: minimum e-process across streams
#
# Both use a p-to-e calibration step.
#
# References:
#   CR24  Choe & Ramdas (2024), Operations Research 72(4), Section 4.4
# =============================================================================

#' Split a sequence into h interleaved lag streams
#'
#' For lag `h`, the k-th stream (`k = 1, ..., h`) contains indices
#' \eqn{\{k,\, k+h,\, k+2h,\, \ldots\}}, following the CR24 convention.
#'
#' @param xs  Numeric vector. Score differences hat_delta_t.
#' @param h   Integer >= 1. Lag (number of steps ahead).
#'
#' @return List of length h. Each element is a numeric vector
#'         containing the score differences for that stream.
#'
#' @examples
#'   split_streams(1:10, h = 3)
#'   # stream 1: indices 1, 4, 7, 10
#'   # stream 2: indices 2, 5, 8
#'   # stream 3: indices 3, 6, 9
#'
#' @export
split_streams <- function(xs, h) {
  stopifnot(h >= 1, length(xs) >= h)
  T_ <- length(xs)
  lapply(seq_len(h), function(k) xs[seq(k, T_, by = h)])
}

#' Unroll a stream-wise quantity back to the original time scale
#'
#' After computing a per-stream cumulative quantity (e.g. e-process values),
#' restores them to the original length `T` by zero-padding the first
#' `k - 1` positions, repeating each stream value `h` times, then truncating
#' to length `T`.
#'
#' @param stream_vals Numeric vector. Values for stream k (length ~ T/h).
#' @param k           Integer. Stream index (1-based).
#' @param h           Integer. Lag.
#' @param T_          Integer. Total original sequence length.
#'
#' @return Numeric vector of length T_.
#'
#' @section Alignment only, not theoretical updating:
#' For lagged forecasts (`h >= 2`), the returned series is aligned to the
#' evaluated score-difference index after stream splitting. It should
#' **not** be interpreted as a process that updates at the original
#' forecast-issuance time. The unrolled process is for visualization and
#' alignment only; the theoretical validity argument relies strictly on the
#' streamwise sub-filtrations, not on this unrolled representation.
#'
#' @examples
#' unroll_stream(c(1, 2, 3), k = 2, h = 2, T_ = 6)
#'
#' @export
unroll_stream <- function(stream_vals, k, h, T_) {
  padded    <- c(rep(0, k - 1), rep(stream_vals, each = h))
  length(padded) <- T_   # truncate or extend with NA then replace
  padded[is.na(padded)] <- stream_vals[length(stream_vals)]
  padded[seq_len(T_)]
}

#' P-to-e calibrator
#'
#' Converts anytime-valid p-values to e-values using the mixture or simple
#' calibrator, as used in CR24 Section 4.4.
#'
#' @param p        Numeric vector of p-values in (0, 1].
#' @param strategy Character. `"mixture"` (default, from Vovk & Wang 2021) or
#'                 `"simple"`. See Details for the formulas.
#' @param eps      Numeric. Numerical guard for log(0). Default: 1e-16.
#'
#' @return Numeric vector of e-values >= 0.
#'
#' @details
#' Mixture calibrator (default, matches Python comparecast behaviour):
#' \deqn{f(p) = \frac{1 - p + p\log(p)}{p\,(\log p)^2}}
#'
#' Simple calibrator (`strategy = "simple"`):
#' \deqn{f(p) = \frac{1}{2\sqrt{p}}}
#'
#' @examples
#' p <- c(0.5, 0.1, 0.01)
#' calibrate_p_to_e(p)
#' calibrate_p_to_e(p, strategy = "simple")
#'
#' @export
calibrate_p_to_e <- function(p, strategy = "mixture", eps = 1e-16) {
  stopifnot(all(p >= 0), all(p <= 1 + 1e-10))
  p <- pmin(pmax(p, eps), 1)

  if (strategy == "mixture") {
    num   <- 1 - p + p * log(p)
    denom <- p * log(p)^2 + eps
    pmax(num / denom, 0)
  } else if (strategy == "simple") {
    1 / (2 * sqrt(p) + eps)
  } else {
    stop("Unknown strategy: '", strategy, "'. Use 'mixture' or 'simple'.")
  }
}

#' Lag-h e-process for sequential forecast comparison (Propositions 5 & 6)
#'
#' For h-step-ahead forecasts, constructs an anytime-valid e-process by
#' stream splitting and combining h individual e-processes.
#'
#' @param scores1     Numeric vector. Scores for forecaster 1.
#' @param scores2     Numeric vector. Scores for forecaster 2.
#' @param h           Integer >= 1. Forecast lag.
#'                    For h=1, reduces to standard eprocess() — no splitting.
#' @param alpha       Numeric in (0,1). Significance level. Default: 0.05.
#' @param c           Numeric > 0. Sub-exponential scale. Default: 2.
#' @param v_opt       Numeric > 0. Default: 10.
#' @param null        Character. Null hypothesis type:
#' * `"pw"` — period-wise weak null (average combination).
#' * `"w"` — standard weak null (minimum combination).
#' @param calibrate   Logical. Apply p-to-e calibration. Default: TRUE.
#' @param cal_strategy Character. "mixture" (default) or "simple".
#'
#' @return data.frame with columns t, e_pq, e_qp, log_e_pq, log_e_qp.
#'
#' @details
#'   For h = 1: calls eprocess() directly and returns its output unchanged.
#'
#'   For h >= 2:
#'     1. Split xs into h streams
#'     2. Compute e-process on each stream independently
#'     3. Combine using the appropriate null rule
#'     4. Convert to p-process, combine, calibrate back to e-process
#'     5. Unroll to original time scale
#'
#'   The period-wise ("pw") null is less conservative than the standard ("w")
#'   null but tests a different (weaker) hypothesis. See CR24 Section 4.4.
#'
#' @examples
#' scores1 <- c(-0.04, -0.09, -0.01, -0.16, -0.04, -0.09)
#' scores2 <- c(-0.09, -0.16, -0.04, -0.25, -0.09, -0.16)
#' ep <- eprocess_lag(scores1, scores2, h = 2, alpha = 0.05)
#' head(ep)
#'
#' @export
eprocess_lag <- function(scores1, scores2,
                         h            = 1,
                         alpha        = 0.05,
                         c            = 2,
                         v_opt        = 10,
                         null         = "pw",
                         calibrate    = TRUE,
                         cal_strategy = "mixture") {

  stopifnot(
    length(scores1) == length(scores2),
    h >= 1,
    null %in% c("pw", "w")
  )

  # h = 1: trivial case, full martingale property holds
  if (h == 1L) {
    return(eprocess(scores1, scores2, alpha = alpha, c = c, v_opt = v_opt))
  }

  xs <- scores1 - scores2
  T_ <- length(xs)

  alpha_opt <- alpha / 2   # two-sided
  rho       <- rho_from_vopt(v_opt = v_opt, alpha = alpha_opt)

  # Split into h streams
  streams <- split_streams(xs, h)

  # Compute e-process on each stream
  # Each stream is a valid martingale difference sequence
  stream_epq <- vector("list", h)
  stream_eqp <- vector("list", h)

  for (k in seq_len(h)) {
    xs_k   <- streams[[k]]
    gam_k  <- make_gammas(xs_k, lag = 1)
    V_k    <- intrinsic_time(xs_k, gam_k, floor = FALSE)
    S_pq_k <- cumsum(xs_k)
    S_qp_k <- -S_pq_k

    stream_epq[[k]] <- clip_eprocess(
      log_ge_mixture_from_sv(S_pq_k, V_k, rho, c)
    )

    stream_eqp[[k]] <- clip_eprocess(
      log_ge_mixture_from_sv(S_qp_k, V_k, rho, c)
    ) }

  # Combine streams
  combine_streams <- function(stream_evals) {
    # Running maximum within each stream (anytime-valid p-value construction)
    run_max_list <- lapply(stream_evals, cummax)

    if (null == "pw") {
      # Period-wise: average running maxima, then convert to p
      run_max_mat <- do.call(cbind, lapply(seq_len(h), function(k) {
        unroll_stream(run_max_list[[k]], k = k, h = h, T_ = T_)
      }))
      mean_max_e <- pmax(rowMeans(run_max_mat), .Machine$double.eps)
      combined_p    <- pmin(1, exp(1) * log(max(h, 2)) / pmax(mean_max_e, .Machine$double.eps))

    } else {
      # Standard weak null: max p-value across streams (most conservative)
      p_per_stream <- lapply(seq_len(h), function(k) {
        pmin(1, 1 / pmax(
          unroll_stream(run_max_list[[k]], k = k, h = h, T_ = T_),
          .Machine$double.eps
        ))
      })
      p_mat         <- do.call(cbind, p_per_stream)
      combined_p    <- apply(p_mat, 1, max)
    }

    if (calibrate) {
      pmax(calibrate_p_to_e(combined_p, strategy = cal_strategy), .Machine$double.eps)
    } else {
      combined_p
    }
  }

  e_pq_combined <- combine_streams(stream_epq)
  e_qp_combined <- combine_streams(stream_eqp)

  # Clip final values
  e_pq_combined <- pmin(e_pq_combined, 1e7)
  e_qp_combined <- pmin(e_qp_combined, 1e7)

  data.frame(
    t        = seq_len(T_),
    e_pq     = e_pq_combined,
    e_qp     = e_qp_combined,
    log_e_pq = log(pmax(e_pq_combined, .Machine$double.eps)),
    log_e_qp = log(pmax(e_qp_combined, .Machine$double.eps))
  )
}
