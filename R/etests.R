# =============================================================================
# etests.R
# Sequential e-processes for testing the weak and strong null hypotheses
#
# Implements:
#   eprocess()                — Theorem 3 CR24: sub-exponential mixture e-process
#                               (weak null)
#   eprocess_betting()        — Proposition 3.2 A26: product-form betting
#                               e-process (strong null)
#   lambda_betting_quantile() — Section 5.1 A26: adaptive betting fraction for
#                               quantile forecasts (strong null)
#   lambda_betting_agrapa()   — adaptive betting fraction for the strong null,
#                               adapting WSR24's aGRAPA plug-in estimator
#   lambda_betting_ons()      — adaptive betting fraction for the strong null,
#                               adapting WSR24's ONS-m algorithm
#
#   lambda_betting_agrapa() and lambda_betting_ons() are original adaptations
#   to the bounded strong-null betting setting; the derivation is not given
#   in WSR24 itself.
#
# The weak null being tested:
#   H_0^w(p, q): Delta_t = (1/t) * sum_{i=1}^t E[hat_delta_i | F_{i-1}] <= 0
#   for all t = 1, 2, ...
#   i.e. forecaster 1 is no better than forecaster 2 on average.
#   Tested via exponential-mixture martingales.
#
# The strong null being tested:
#   H_0^s(p, q): mu_t = E[hat_delta_t | F_{t-1}] <= 0 for all t = 1, 2, ...
#   i.e. forecaster 1 is no better than forecaster 2 at every time step.
#   Tested via product-form predictable betting martingales.
#
# References:
#   CR24   Choe & Ramdas (2024), Operations Research 72(4)
#   H21    Howard et al. (2021), Annals of Statistics 49(2)
#   A26    Arnold et al. (2026), J R Stat Soc B 88, qkag066
#   WSR24  Waudby-Smith & Ramdas (2024), J R Stat Soc B 86(1), 1-27
# =============================================================================

#' Sub-exponential mixture e-process (Theorem 3, Choe & Ramdas 2024)
#'
#' Constructs two simultaneous one-sided e-processes for sequentially testing
#' whether forecaster 1 (p) outperforms forecaster 2 (q) or vice versa.
#'
#' The mixture e-process at time t is:
#' \deqn{E_t^{\mathrm{mix}} = m(S_t, \hat{V}_t)}
#' where \eqn{S_t = \sum_{i=1}^t \hat{\delta}_i},
#' \eqn{\hat{V}_t = \sum_{i=1}^t (\hat{\delta}_i - \gamma_i)^2},
#' and \eqn{m(s, v)} is the Gamma-Exponential mixture function (Proposition EC.3,
#' CR24).
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
#'                  fastest. Default: 10 (recommended by CR24).
#' @param alpha_opt Numeric in (0,1). One-sided alpha used to compute rho.
#'                  Default: alpha/2 (matches comparecast two-sided convention).
#' @param gammas    Numeric vector or NULL. Predictable centering sequence.
#'                  If NULL, constructed as lagged running mean.
#' @param clip_max  Numeric. Maximum e-process value before clipping.
#'                  Default: 1e7 (matches Python comparecast).
#'
#' @return data.frame with the following columns:
#' \describe{
#'  \item{`t`}{Time index.}
#'  \item{`e_pq`, `e_qp`}{One-sided e-processes. `e_pq` tests H_0^w(p, q):
#'  whether forecaster `p` outperforms `q`; `e_qp` tests H_0^w(q, p): whether
#'  forecaster `q` outperforms `p`.}
#'  \item{`log_e_pq`, `log_e_qp`}{Log-scale values of the e-processes, clipped
#'  at log(clip_max).}
#'  }
#'
#'
#' @section Rejection rule:
#' At level `alpha`: reject \eqn{H_0^w(p, q)} (conclude `p` outperforms `q`)
#' when `e_pq >= 2/alpha`; reject \eqn{H_0^w(q, p)} (conclude `q` outperforms
#' `p`) when `e_qp >= 2/alpha`. Use `eprocess_rejections()` to extract the
#' first crossing time for each.
#'
#' @details
#'   VARIANCE PROCESS: The intrinsic time V_hat_t uses NO floor. The GE mixture
#'   m(s, v) is well-defined at v=0 (returns 1 when s=0), so no floor is needed.
#'   Adding a floor would yield less power in the e-process.
#'
#'   SCALE CONVENTION: c is the sub-exponential scale parameter such that
#'   |hat_delta_i| <= c/2. This is the Theorems 2 & 3 convention from CR24.
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

#' Betting-style e-process for the strong null hypothesis
#'
#' Implements the product-form test martingale
#' \deqn{E_t = \prod_{r=1}^{t} (1 + \lambda_r \hat\delta_r)}
#' for testing the strong null \eqn{H_0^s(p, q): \delta_t \le 0} for all t.
#'
#' @param scores1 Numeric vector. Scores S(p_t, y_t) for forecaster 1.
#' @param scores2 Numeric vector. Scores S(q_t, y_t) for forecaster 2.
#' @param c_t Numeric scalar or vector (same length as scores) of predictable bounds
#'   such that |scores1 - scores2| <= c_t / 2 almost surely at every step.
#' @param lambda_t Optional numeric vector of predictable betting fractions in
#'   `[0, 1/c_t]`. If `NULL` (default), uses the fixed fraction `lambda_t = 1 / (2 * c_t)`.
#' @param clip_max Numeric. Maximum e-process value before clipping. Default: 1e7.
#'
#' @return data.frame with columns t, e_pq, e_qp, log_e_pq, log_e_qp.
#'
#' @examples
#' set.seed(123)
#' T_sim <- 100
#' y <- rbinom(T_sim, size = 1, prob = 0.7)
#'
#' # Forecaster 1 has a genuine edge (closer to true probability 0.7)
#' p1 <- runif(T_sim, 0.5, 0.9)
#' # Forecaster 2 is just guessing uniformly
#' p2 <- runif(T_sim, 0.1, 0.9)
#'
#' # 1. Compute positively oriented Brier scores
#' s1 <- brier_score(p1, y)
#' s2 <- brier_score(p2, y)
#'
#' # 2. Establish the bound.
#' # Brier scores lie in [-1, 0], so the maximum absolute difference is 1.
#' # eprocess_betting() requires |s1 - s2| <= c_t / 2, so c_t = 2 is globally valid.
#' res <- eprocess_betting(scores1 = s1, scores2 = s2, c_t = 2)
#'
#' # The e-process e_pq tests the null that Forecaster 1 is NOT better than Forecaster 2.
#' # Because Forecaster 1 is genuinely better, e_pq accumulates massive evidence.
#' head(res)
#' tail(res, 3)
#'
#' @export
eprocess_betting <- function(scores1, scores2, c_t, lambda_t = NULL, clip_max = 1e7) {
  stopifnot(
    length(scores1) == length(scores2),
    length(scores1) >= 1,
    clip_max > 0
  )

  t_len <- length(scores1)
  xs <- scores1 - scores2

  if (length(c_t) == 1L) {
    c_t <- rep(c_t, t_len)
  }

  if (length(c_t) != t_len) {
    stop("c_t must be a scalar or a vector of the same length as scores1/scores2.")
  }
  if (any(c_t <= 0)) {
    stop("c_t must be strictly positive at every step.")
  }

  if (is.null(lambda_t)) {
    lambda_t <- 1 / (2 * c_t)
  } else {
    if (length(lambda_t) == 1L) {
      lambda_t <- rep(lambda_t, t_len)
    }
  }
  stopifnot(length(lambda_t) == t_len)

  if (any(lambda_t < 0) || any(lambda_t > 1 / c_t + 1e-8)) {
    stop("lambda_t must lie in [0, 1/c_t] at every step.")
  }

  terms_pq <- 1 + lambda_t * xs
  terms_qp <- 1 - lambda_t * xs

  if (any(terms_pq <= 0) || any(terms_qp <= 0)) {
    stop("A betting factor was non-positive; check that |scores1 - scores2| <= c_t / 2 at every step.")
  }

  log_e_pq <- cumsum(log(terms_pq))
  log_e_qp <- cumsum(log(terms_qp))

  data.frame(
    t        = seq_len(t_len),
    e_pq     = clip_eprocess(log_e_pq, clip_max = clip_max),
    e_qp     = clip_eprocess(log_e_qp, clip_max = clip_max),
    log_e_pq = pmin(log_e_pq, log(clip_max)),
    log_e_qp = pmin(log_e_qp, log(clip_max))
  )
}

#' Adaptive betting fraction for quantile-forecast strong-null tests
#'
#' Implements the quantile-specific adaptive betting scheme of Arnold et al. (2026),
#' translating their log-scale loss convention to seqcomp's positively-oriented scores.
#'
#' @note
#' **Scale Translation:** This function assumes `p_t`, `q_t`, and `delta_hat_lag1`
#' are calculated on the raw, linear scale. To replicate the exact log-scale bounds
#' used in the Arnold et al. (2026) Covid-19 case study, the forecast vectors
#' passed to this function must be log-transformed prior to evaluation.
#'
#' @details
#' If two forecasts are identical at a time point, their tick-loss difference
#' and its analytic bound are both zero. To accommodate the strictly positive
#' bound required by [eprocess_betting()] and the adaptive betting rules, this
#' function replaces such bounds (and any smaller bounds) by `eps`. This is a
#' conservative predictable enlargement of the bound; the corresponding
#' Arnold betting fraction is capped at `1 / c_safe`.
#'
#' @param p_t Numeric vector. Quantile forecasts of the first forecaster.
#' @param q_t Numeric vector. Quantile forecasts of the second forecaster.
#' @param tau Numeric scalar in (0, 1). The quantile level.
#' @param delta_hat_lag1 Numeric vector. The *previous* step's score difference.
#'   Must have `delta_hat_lag1[1] = 0`.
#' @param eps Positive numeric scalar. Used as a denominator safeguard against
#'   zero. Default: `1e-8`.
#'
#' @return A list with `c_t` and `lambda_t` vectors for `eprocess_betting()`.
#'
#' @examples
#' set.seed(456)
#' T_sim <- 100
#' y <- rnorm(T_sim)
#' tau <- 0.90
#'
#' # Forecaster 1 correctly predicts the true 90th percentile (~1.28)
#' p_t <- rep(qnorm(tau), T_sim)
#' # Forecaster 2 is biased and incorrectly predicts the median (0.0)
#' q_t <- rep(0.0, T_sim)
#'
#' # 1. Compute pointwise tick loss (positively oriented)
#' s_p <- tick_loss(p_t, y, tau)
#' s_q <- tick_loss(q_t, y, tau)
#'
#' # 2. Compute the lagged score difference required by Arnold's heuristic
#' xs <- s_p - s_q
#' delta_lag1 <- c(0, head(xs, -1))
#'
#' # 3. Generate the predictable bounds (c_t) and betting fractions (lambda_t)
#' bnds <- lambda_betting_quantile(p_t, q_t, tau, delta_hat_lag1 = delta_lag1)
#'
#' # 4. Plug these directly into the betting e-process
#' res <- eprocess_betting(
#'   scores1 = s_p,
#'   scores2 = s_q,
#'   c_t = bnds$c_t,
#'   lambda_t = bnds$lambda_t
#' )
#'
#' # View the final evidence accumulation
#' tail(res[, c("t", "e_pq", "e_qp")], 3)
#'
#' @export
lambda_betting_quantile <- function(p_t, q_t, tau, delta_hat_lag1 = NULL, eps = 1e-8) {
  t_len <- length(p_t)
  if (length(q_t) != t_len) stop("p_t and q_t must have the same length.")
  if (tau <= 0 || tau >= 1) stop("tau must lie strictly between 0 and 1.")

  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a finite positive scalar.", call. = FALSE)
  }

  if (is.null(delta_hat_lag1)) {
    delta_hat_lag1 <- c(0, rep(NA_real_, t_len - 1))
  }
  if (length(delta_hat_lag1) != t_len) stop("delta_hat_lag1 must have the same length as p_t.")
  if (t_len > 1 && anyNA(delta_hat_lag1[-1])) {
    stop("delta_hat_lag1 must supply the previous step's score difference for t >= 2.")
  }

  c_t <- 2 * max(tau, 1 - tau) * abs(p_t - q_t)
  c_safe <- pmax(c_t, eps)

  u <- abs(tau - 0.5)
  K_t <- ((2 - u) / (1 + u)) * ((3 * pi / 2 + atan(delta_hat_lag1)) / pi)

  lambda_t <- pmin(1 / (K_t * c_t + eps), 1 / c_safe)

  list(c_t = c_safe, lambda_t = lambda_t)
}

#' Approximate GRAPA (aGRAPA) betting fractions for a constant bound
#'
#' Constructs a predictable betting-fraction sequence for the product-form
#' strong-null e-process ([eprocess_betting()]), by mapping the
#' constant-bound score-difference stream onto \eqn{[0,1]} and applying the
#' aGRAPA plug-in estimator of Waudby-Smith & Ramdas (2024), Online
#' Supplementary Material, Section B.3, evaluated at the fixed null
#' \eqn{m = 1/2}.
#'
#' @param xs Numeric vector. Score-difference stream
#'   \eqn{\hat\delta_t = S(p_t,y_t) - S(q_t,y_t)}.
#' @param c  Numeric > 0. Either a scalar (constant bound) or a
#'   length-`length(xs)` predictable vector `c_t` (must satisfy `|xs_t| <= c_t/2`
#'   pointwise, with `c_t` known before `xs_t` is observed). Scalars are
#'   recycled to a vector of the same length as `xs`.
#' @param kappa Numeric in (0, 1]. WSR's internal truncation safeguard on the
#'   `Y`-scale (`[-2*kappa, 2*kappa]`). Default `0.5`. **Note on saturation:**
#'   Because `seqcomp` applies a strict `[0, 1]` explicit projection *after*
#'   WSR's native truncation to ensure one-sided null validity, any `kappa >= 0.5`
#'   is mathematically equivalent. The lower bound `-2*kappa` is always superseded
#'   by the `0` floor, and for `kappa >= 0.5`, the upper bound `2*kappa >= 1` is
#'   superseded by the `1` cap. `kappa` only actively restricts the bet size when
#'   `kappa < 0.5` (e.g., `kappa = 0.1`).
#' @param prior_mean Numeric. Regularization prior for the running mean
#'   estimator on the `Y`-scale. Default `0.5` (WSR's own default; also
#'   happens to equal the fixed null `m = 1/2` here, which is why the
#'   very first returned value is exactly `0`).
#' @param prior_variance Numeric in (0, 0.25]. Regularization prior for the
#'   running variance estimator. Default `0.25` (WSR's own default).
#' @param fake_obs Numeric > 0. Number of "fake observations" for
#'   regularization. Default `1` (WSR's own default).
#'
#' @return Numeric vector of length `length(xs)`: predictable
#'   \eqn{\lambda_t} values in `[0, 1/c]`, for use as `eprocess_betting()`'s
#'   `lambda_t` argument together with `c_t = c`.
#'
#' @details
#' Derivation (constant `c` only): map
#' \eqn{Y_t = \hat\delta_t / c + 1/2 \in [0,1]}, so the null
#' \eqn{\mu_t \le 0} becomes \eqn{\mathbb{E}[Y_t \mid \mathcal{F}_{t-1}] \le 1/2}.
#' WSR's aGRAPA plug-in is run exactly as published, fixed at `m = 1/2`,
#' on the `Y`-scale running mean/variance. Two clips are then applied in
#' sequence: WSR's own native truncation `[-2*kappa, 2*kappa]`, followed by
#' an explicit floor/cap onto `[0, 1]`. The floor is important for the one-sided
#' composite null `mu <= 0` to remain valid, the played `lambda_t` must be
#' nonnegative on *every* round. WSR's raw aGRAPA value routinely goes negative
#' whenever the running mean currently sits below the null, which happens
#' routinely under the null itself. The cap enforces Arnold's own Prop 3.2 bound
#' `lambda_d <= 1/c`, which is strictly tighter than WSR's native
#' bankruptcy-avoidance bound.
#'
#' The floor/cap projection and the `m_0 = 1/2` null hold identically for a
#' predictable, time-varying `c_t`. Because `c_t` is \eqn{F_{t-1}}-measurable,
#' \eqn{Y_t = x_t / c_t + 1/2} still satisfies
#' \eqn{\mathbb{E}[Y_t \mid \mathcal{F}_{t-1}] \le 1/2} under the strong null,
#' and the same aGRAPA plug-in can be run unmodified at the fixed null `m = 1/2`;
#' only the final \eqn{Y \mapsto d} conversion (`lam_Y / c`) becomes pointwise.
#' This time-varying extension is original to this package.
#'
#'
#' @references
#' Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded
#' random variables by betting. Journal of the Royal Statistical Society Series
#' B: Statistical Methodology, 86(1), 1–27.
#'
#' @examples
#' xs <- c(0.6, -0.2, 0.4)
#' lambda_betting_agrapa(xs, c = 2)
#'
#' @export
lambda_betting_agrapa <- function(xs, c,
                                  kappa          = 0.5,
                                  prior_mean     = 0.5,
                                  prior_variance = 0.25,
                                  fake_obs       = 1) {
  T_ <- length(xs)
  if (length(c) == 1L) {
    c <- rep(c, T_)
  }
  stopifnot(
    length(c) == T_, all(c > 0),
    kappa > 0, kappa <= 1,
    prior_variance > 0, prior_variance <= 0.25,
    fake_obs > 0
  )
  if (any(abs(xs) > c / 2 + 1e-8)) {
    stop("xs contains values with |xs_t| > c_t/2; the Y_t = xs_t/c_t + 1/2 mapping ",
         "requires this bound to hold pointwise for a (possibly time-varying, ",
         "predictable) c_t.")
  }

  ts <- seq_len(T_)

  # Map to Y-scale: Y_t = xs_t / c + 1/2 in [0,1]; fixed null m0 = 1/2.
  ys <- xs / c + 0.5
  m0 <- 0.5

  # Predictable running mean/variance on the Y-scale (WSR aGRAPA plug-in,
  # Sec. B.3 / Eq. 15), lagged by one step, exactly as in WSR's own code.
  mu_hat_t    <- (fake_obs * prior_mean + cumsum(ys)) / (ts + fake_obs)
  mu_hat_lag1 <- c(prior_mean, mu_hat_t[-T_])

  sigma2_t    <- (fake_obs * prior_variance + cumsum((ys - mu_hat_t)^2)) / (ts + fake_obs)
  sigma2_lag1 <- c(prior_variance, sigma2_t[-T_])

  # aGRAPA plug-in at m = m0, WSR's own truncation [-2*kappa, 2*kappa]
  raw_num <- mu_hat_lag1 - m0
  raw     <- raw_num / (sigma2_lag1 + raw_num^2)
  lam_Y   <- pmin(pmax(raw, -kappa / (1 - m0)), kappa / m0)

  # Explicit projection onto Arnold's Prop 3.2 feasible region (Y-scale [0,1]).
  lam_Y <- pmin(pmax(lam_Y, 0), 1)

  # Map back to the d-scale.
  lam_Y / c
}


#' Online Newton Step (ONS-m) betting fractions for a constant bound
#'
#' Constructs a predictable betting-fraction sequence for the product-form
#' strong-null e-process ([eprocess_betting()]), by mapping the
#' constant-bound score-difference stream onto \eqn{[0,1]} and running
#' WSR's ONS-m algorithm (Waudby-Smith & Ramdas 2024, Online Supplementary
#' Material, Algorithm 1, Section B.5), projected each round onto `[0,1]`
#' instead of their native symmetric box, and evaluated at the fixed null
#' \eqn{m = 1/2}.
#'
#' @param xs Numeric vector. Score-difference stream
#'   \eqn{\hat\delta_t = S(p_t,y_t) - S(q_t,y_t)}.
#' @param c  Numeric > 0. Either a scalar (constant bound) or a
#'   length-`length(xs)` predictable vector `c_t` (must satisfy `|xs_t| <= c_t/2`
#'   pointwise, with `c_t` known before `xs_t` is observed). Scalars are
#'   recycled to a vector of the same length as `xs`.
#' @param eta Numeric > 0. WSR's ONS step-size constant. Default
#'   `2 / (2 - log(3))`, WSR's own stated value (natural log).
#'
#' @return Numeric vector of length `length(xs)`: predictable
#'   \eqn{\lambda_t} values in `[0, 1/c]`, for use as `eprocess_betting()`'s
#'   `lambda_t` argument together with `c_t = c`.
#'
#' @details
#' Same `Y_t = xs_t/c + 1/2` mapping as [lambda_betting_agrapa()]. WSR's
#' Algorithm 1 is run exactly as stated, but using the standard OCO gradient
#' for minimizing the negative log-wealth \eqn{-\log(1 + \lambda y_t)}, and
#' replacing WSR's own projection target `[-c/(1-m), c/m]` with `[0, 1]`.
#' Following generic online-convex-optimization theory, \eqn{\lambda_t^O} is
#' projected onto `[0,1]` before being used to compute the next round's
#' gradient. Validity of the resulting e-process only requires the played
#' `lambda_d,t` to lie in `[0, 1/c]` regardless of derivation.
#'
#' This modification is original to this package.
#'
#' @references
#' Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded
#' random variables by betting. Journal of the Royal Statistical Society Series
#' B: Statistical Methodology, 86(1), 1–27.
#'
#' @examples
#' xs <- c(0.6, -0.2, 0.4)
#' lambda_betting_ons(xs, c = 2)
#'
#' @export
lambda_betting_ons <- function(xs, c, eta = 2 / (2 - log(3))) {
  T_ <- length(xs)
  if (length(c) == 1L) {
    c <- rep(c, T_)
  }
  stopifnot(length(c) == T_, all(c > 0), eta > 0)
  if (any(abs(xs) > c / 2 + 1e-8)) {
    stop("xs contains values with |xs_t| > c_t/2; the Y_t = xs_t/c_t + 1/2 mapping ",
         "requires this bound to hold pointwise for a (possibly time-varying, ",
         "predictable) c_t.")
  }

  ys  <- xs / c   # y_t = Y_t - 1/2, elementwise if c is a vector

  lam_Y     <- numeric(T_)
  lam_state <- 0
  A         <- 1

  for (t in seq_len(T_)) {
    lam_Y[t] <- lam_state

    y_t <- ys[t]

    # Standard OCO gradient for minimizing -log(1 + lambda * y_t)
    grad <- -y_t / (1 + lam_state * y_t)
    A    <- A + grad^2

    raw       <- lam_state - eta * grad / A
    lam_state <- min(1, max(0, raw))
  }

  lam_Y / c
}

#' Period-conditional predictable betting fractions
#'
#' Wraps a base betting-fraction rule (e.g. lambda_betting_agrapa,
#' lambda_betting_ons) by running it independently on each of `period`
#' interleaved sub-streams, matching split_streams()'s indexing convention.
#' Each sub-stream's lambda_t is computed using only prior values from the
#' SAME sub-stream, so predictability (F_{t-1}-measurability) is preserved.
#' Useful for known periodic anomalies (e.g., day-of-week seasonality).
#'
#' @keywords internal
#' @noRd
lambda_betting_periodic <- function(xs, c, period, base_fun = lambda_betting_agrapa, ...) {
  T_ <- length(xs)
  if (length(c) == 1L) c <- rep(c, T_)
  stopifnot(length(c) == T_, period >= 1, T_ >= period)

  lambda_out <- numeric(T_)
  for (k in seq_len(period)) {
    idx <- seq(k, T_, by = period)
    lambda_out[idx] <- base_fun(xs = xs[idx], c = c[idx], ...)
  }
  lambda_out
}
