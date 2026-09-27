# =============================================================================
# scores.R
# Scoring rules for probabilistic forecast evaluation
#
# All scoring rules are POSITIVELY ORIENTED: higher values = better forecast.
# This matches the CR24 convention:
#   hat_delta_t = S(p_t, y_t) - S(q_t, y_t)
# is positive when forecaster p outperforms forecaster q.
#
# Binary/categorical convention:
#   - Vector probability input is treated as binary and uses y in {0, 1}.
#   - Matrix probability input is treated as categorical and uses y in {1, ..., K}
#     with R-native one-based class labels.
#
# Scoring rules implemented:
#   1. brier_score()      — Brier/quadratic score, binary or categorical,
#                           bounded in [-1, 0]
#   2. log_score()        — Logarithmic score, binary or categorical,
#                           unbounded below
#   3. spherical_score()  — Spherical score, binary or categorical,
#                           bounded in [0, 1]
#   4. tick_loss()        — Quantile/tick loss, negated
#   5. qlike_score()      — QLIKE for variance forecasts, negated
#   6. winkler_score()    — Winkler normalised score for binary outcomes
#   7. CRPS wrappers      — CRPS wrappers over scoringRules, negated
#
# References:
#   GR07  Gneiting & Raftery (2007), JASA 102(477)
#   CR24  Choe & Ramdas (2024), Operations Research 72(4)
#   JKL19 Jordan, Krueger & Lerch (2019), Journal of Statistical Software, 90(12)
# =============================================================================


#' Preprocess categorical probability forecasts
#'
#' Internal helper used by matrix-input branches of `brier_score()`,
#' `spherical_score()`, and `log_score()`.
#'
#' Vector input is converted to a two-column probability matrix using the
#' convention that the vector gives `P(Y = class 2)`. Matrix input is treated
#' as categorical probability forecasts with rows on the probability simplex.
#'
#' @param p Numeric vector or matrix of forecast probabilities.
#' @param y Integer/numeric vector of one-based class labels.
#'
#' @return A list with elements `p_mat`, `y_int`, `T`, and `K`.
#'
#' @keywords internal
#' @noRd
.preprocess_cat <- function(p, y) {

  # ── 1. Shape detection ──────────────────────────────────────────────────
  if (is.vector(p) && !is.matrix(p)) {
    # Binary shorthand: p is a length-T vector of P(Y = class 2).
    # NOTE: when called from brier_score / spherical_score / log_score,
    # the vector path is handled by the binary code branch in those
    # functions and .preprocess_cat() is NOT called.  This branch exists
    # for standalone / pipeline use where matrix form is always desired.
    if (!is.numeric(p))
      stop("p must be numeric.")
    T_ <- length(p)
    K  <- 2L
    p_mat <- cbind(1 - p, p)
  } else if (is.matrix(p)) {
    if (!is.numeric(p))
      stop("p must be a numeric matrix.")
    T_ <- nrow(p)
    K  <- ncol(p)
    p_mat <- p
  } else {
    stop("p must be a numeric vector (binary shorthand) or a numeric matrix (categorical).")
  }

  # ── 2. Dimension sanity ─────────────────────────────────────────────────
  if (as.integer(K) < 2L)
    stop("p must have at least 2 columns (K >= 2).")

  # ── 3. Probability validity ─────────────────────────────────────────────
  if (any(!is.finite(p_mat)))
    stop("p contains non-finite values (NA, NaN, Inf).")
  if (any(p_mat < -1e-9 | p_mat > 1 + 1e-9))
    stop("All entries of p must be in [0, 1].")
  p_mat <- pmin(pmax(p_mat, 0), 1)   # clip negligible floating-point deviations

  row_sums <- rowSums(p_mat)
  if (any(abs(row_sums - 1) > 1e-6)) {
    stop(sprintf(
      "Rows of p must sum to 1; max deviation is %.2e.",
      max(abs(row_sums - 1))
    ))
  }

  # ── 4. Strict integer validation for y ─────────────────────────────────
  if (!is.numeric(y) && !is.integer(y))
    stop("y must be an integer/numeric vector of class labels.")
  if (any(!is.finite(y)))
    stop("y contains non-finite values.")
  if (any(abs(y - round(y)) > 0))
    stop("y must contain integer-valued class labels (e.g. y = 1.9 is not accepted).")
  y <- as.integer(y)

  if (length(y) != T_)
    stop(sprintf(
      "length(y) = %d but nrow(p) = %d; they must agree.",
      length(y), T_
    ))
  if (any(y < 1L | y > as.integer(K)))
    stop(sprintf(
      "y values must be integers in 1:%d (R 1-based labels).",
      K
    ))

  list(p_mat = p_mat, y_int = y, T = T_, K = as.integer(K))
}

#' Brier score for binary and categorical forecasts
#'
#' Computes the positively oriented Brier/quadratic score. Vector probability
#' input is treated as binary; matrix probability input is treated as
#' categorical.
#'
#' @param p Numeric vector in `[0, 1]` for binary forecasts, or a numeric
#'   matrix whose rows are probability vectors for categorical forecasts.
#' @param y For binary vector input, numeric vector in `{0, 1}`. For
#'   categorical matrix input, integer vector in `{1, ..., K}`, where
#'   `K = ncol(p)`.
#'
#' @return Numeric vector of scores in `[-1, 0]`. Higher is better.
#'
#' @details
#' For binary forecasts, this computes
#' \deqn{S(p, y) = -(p-y)^2.}
#'
#' For categorical forecasts, this computes
#' \deqn{S(\mathbf{p}, y) = -\frac{1}{2}\|\mathbf{p} - e_y\|_2^2,}
#' where `e_y` is the one-hot vector of the realised category.
#'
#' With the convention that category 2 corresponds to the binary event
#' `y = 1`, the categorical formula recovers the binary formula exactly
#' when `K = 2`.
#'
#' @section Bounds:
#' Score differences lie in `[-1, 1]`, so use `c = 1` for Theorem 1 and
#' `c = 2` for Theorems 2 and 3.
#'
#' @examples
#' p <- c(0.2, 0.7, 0.9)
#' y <- c(0, 1, 1)
#' brier_score(p, y)
#'
#' @export
brier_score <- function(p, y) {
  if (is.matrix(p)) {
    prep <- .preprocess_cat(p, y)
    idx  <- cbind(seq_len(prep$T), prep$y_int)
    p_y  <- prep$p_mat[idx]

    return(p_y - 0.5 * rowSums(prep$p_mat^2) - 0.5)
  }

  stopifnot(
    length(p) == length(y),
    all(p >= 0 & p <= 1),
    all(y %in% c(0, 1))
  )

  -(p - y)^2
}

#' Logarithmic score for binary and categorical forecasts
#'
#' Computes the positively oriented logarithmic score. Vector probability input
#' is treated as binary; matrix probability input is treated as categorical.
#'
#' @param p Numeric vector in `[0, 1]` for binary forecasts, or a numeric
#'   matrix whose rows are probability vectors for categorical forecasts.
#' @param y For binary vector input, numeric vector in `{0, 1}`. For
#'   categorical matrix input, integer vector in `{1, ..., K}`, where
#'   `K = ncol(p)`.
#' @param eps Numeric. Probability floor used before taking logarithms.
#'   Default is `1e-15`. Set to `0` to disable clipping.
#'
#' @return Numeric vector of scores in `(-Inf, 0]`. Higher is better.
#'
#' @details
#' For binary forecasts, this computes
#' \deqn{S(p, y) = y\log(p) + (1-y)\log(1-p).}
#'
#' For categorical forecasts, this computes
#' \deqn{S(\mathbf{p}, y) = \log(p_y),}
#' where `p_y` is the forecast probability assigned to the realised category.
#'
#' @section Use with seqcomp:
#' The logarithmic score is unbounded below. It should not be used directly
#' with the finite-sample bounded-difference confidence sequences or
#' e-processes. For binary outcomes, use `winkler_score()` and `winkler_cs()`
#' when the Winkler construction is appropriate. For unbounded score
#' differences, use `cs_asymptotic()` or supply genuine predictable bounds
#' to `eprocess_predictable()`.
#'
#' @examples
#' p <- c(0.2, 0.7, 0.9)
#' y <- c(0, 1, 1)
#' log_score(p, y)
#'
#' @export
log_score <- function(p, y, eps = 1e-15) {
  if (is.matrix(p)) {
    prep <- .preprocess_cat(p, y)
    idx  <- cbind(seq_len(prep$T), prep$y_int)
    p_y  <- prep$p_mat[idx]

    return(log(pmax(p_y, eps)))
  }

  stopifnot(
    length(p) == length(y),
    all(is.finite(p)),
    all(y %in% c(0, 1)),
    eps >= 0
  )

  p <- pmax(pmin(p, 1 - eps), eps)
  y * log(p) + (1 - y) * log(1 - p)
}

#' Negated tick loss for quantile forecasts
#'
#' Computes the positively oriented (negated) tick/quantile loss
#' (Koenker & Bassett, 1978).
#'
#' @param q     Numeric vector. Quantile forecasts at level tau.
#' @param y     Numeric vector. Realised outcomes.
#' @param tau Numeric in (0,1). Quantile level.
#'
#' @return Numeric vector of negated tick loss scores. Higher = better.
#'
#' @details
#' The standard tick loss is
#' \deqn{\rho_\tau(u) = u \left(\tau - \mathbb{1}(u < 0)\right),}
#' where \eqn{u = y - q_\tau} is the forecast error. This is loss-oriented
#' (lower = better), so the function negates it:
#' \deqn{S_T(q, y; \tau) = -(y - q)\left(\tau - \mathbb{1}(y < q)\right).}
#'
#' Tick loss is unbounded on general real-valued outcomes. Bounds derived from
#' an empirical data range are ex-post and do not provide theorem-valid
#' constants for finite-sample Hoeffding/Bernstein confidence sequences or
#' e-processes.
#'
#' Sign convention: the negation means `hat_delta_t > 0` when forecaster `p`
#' has smaller tick loss, hence a better quantile forecast, than forecaster `q`.
#'
#' @examples
#' q <- c(1.0, 1.5, 2.0)
#' y <- c(1.2, 1.4, 2.3)
#' tick_loss(q, y, tau = 0.5)
#'
#' @export
tick_loss <- function(q, y, tau) {
  stopifnot(tau > 0, tau < 1, length(q) == length(y))
  u <- y - q
  -(u * (tau - as.numeric(u < 0)))
}

#' Negated QLIKE score for variance forecasts
#'
#' Computes the positively oriented (negated) QLIKE quasi-likelihood loss for
#' variance forecasts.
#'
#' @param sigma2_hat Numeric vector. Forecast variance (strictly positive).
#' @param sigma2     Numeric vector. Realised variance (strictly positive).
#'
#' @return Numeric vector of negated QLIKE scores. Higher is better.
#'   Maximum value is 0, achieved at a perfect forecast `sigma2_hat = sigma2`.
#'   Unbounded below.
#'
#' @details
#' Standard QLIKE loss is
#' \deqn{L_{QL}(\hat\sigma^2, \sigma^2) = \frac{\sigma^2}{\hat\sigma^2} -
#'   \log\frac{\sigma^2}{\hat\sigma^2} - 1.}
#' This is loss-oriented (lower = better, minimum 0 at a perfect forecast), so
#' the function negates it: \eqn{S_{QL} = -L_{QL}}.
#'
#' Literature note: some sources define QLIKE as
#' `log(sigma2_hat) + sigma2 / sigma2_hat`, which differs by constants from
#' the form above. Here the loss is normalised to have minimum 0 and is then
#' negated for positive orientation.
#'
#' @section Unbounded below:
#' QLIKE is unbounded below. It should not be used directly with the
#' finite-sample bounded-difference confidence sequences or e-processes.
#' Use `cs_asymptotic()` for QLIKE-based confidence sequences, or use
#' `eprocess_predictable()` only when genuine ex ante predictable bounds are
#' available. QLIKE is not compatible with the Winkler construction because
#' Winkler scores are restricted to binary outcomes and probability forecasts.
#'
#' @examples
#' sigma2_hat <- c(1.0, 1.5, 2.0)
#' sigma2 <- c(1.1, 1.4, 2.2)
#' qlike_score(sigma2_hat, sigma2)
#'
#' @export
qlike_score <- function(sigma2_hat, sigma2) {
  stopifnot(all(sigma2_hat > 0), all(sigma2 > 0))
  ratio <- sigma2 / sigma2_hat
  -(ratio - log(ratio) - 1)
}

#' Winkler-normalized binary score
#'
#' Normalises the score difference S(p,y) - S(q,y) by the maximum possible
#' score difference given the forecaster ordering, mapping the result to
#' `(-Inf, 1]` (Proposition EC.4, Choe & Ramdas 2024). Used to apply Theorems 2 & 3
#' to unbounded scoring rules on binary outcomes.
#'
#' @param p         Numeric vector in (0,1). Forecasts from model 1.
#' @param q         Numeric vector in (0,1). Forecasts from model 2.
#' @param y         Numeric vector containing only 0 and 1. Binary outcomes.
#' @param base_score Function. The underlying scoring rule S(p, y).
#'                   Must accept two arguments: forecast probability and outcome.
#'                   Default: log_score (with eps clipping).
#' @param eps       Numeric. Zero-protection for the normaliser denominator.
#'                  Default: 1e-8 (matches Python comparecast convention).
#'
#' @return Numeric vector. Winkler scores in `(-Inf, 1]`.
#'         Upper bound of 1 is tight: w = 1 when y = 1(p > q).
#'
#' @details
#' \deqn{w(p, q, y) = \frac{S(p,y) - S(q,y)}{S(p, \mathbb{1}(p>q)) - S(q, \mathbb{1}(p>q))}}
#' with the convention 0/0 := 0.
#'
#' The lower bound is problem-dependent (depends on how extreme p and q can
#' be). For a two-sided CS via Corollary 2, the user must establish a finite
#' lower bound analytically. If no finite lower bound can be guaranteed, use
#' the one-sided (upper) CS only, as in the CR24 MLB experiments.
#'
#' @section When to use:
#' Strictly limited to binary outcomes `y` in `{0, 1}` and probability
#' forecasts `p`, `q` in `(0, 1)`. Not applicable to QLIKE or other
#' continuous-outcome scoring rules. See CR24 EC.7 for discussion.
#'
#' For use in Theorems 2 & 3: upper bound = 1 implies c/2 = 1, so use `c = 2`
#' in all GE boundary and e-process calls.
#'
#' @examples
#' p <- c(0.7, 0.6, 0.8, 0.65)
#' q <- c(0.5, 0.7, 0.6, 0.55)
#' y <- c(1, 1, 0, 1)
#' winkler_score(p, q, y)
#'
#' @export
winkler_score <- function(p, q, y,
                          base_score = log_score,
                          eps = 1e-8) {
  stopifnot(
    all(p > 0 & p < 1),
    all(q > 0 & q < 1),
    all(y %in% c(0, 1)),
    length(p) == length(q),
    length(p) == length(y)
  )

  # Numerator: actual score difference
  num <- base_score(p, y) - base_score(q, y)

  # Denominator: score difference at the "favourable" outcome
  # If p > q: favourable outcome is y=1 (both forecasters assigned more
  #           probability to 1, but p assigned more, so p wins when y=1)
  # If p < q: favourable outcome is y=0
  y_fav <- as.numeric(p > q)          # 1 where p > q, 0 where p <= q

  denom <- base_score(p, y_fav) - base_score(q, y_fav)

  # Zero-protection: replace |denom| < eps with eps (preserving sign)
  # When p == q exactly, both scores are identical and the ratio is 0/0 := 0.
  # The eps guard handles near-ties numerically.
  denom_safe <- ifelse(abs(denom) < eps, eps, denom)

  result <- num / denom_safe

  # Enforce 0/0 := 0 convention from CR24
  result[abs(num) < eps & abs(denom) < eps] <- 0

  return(result)
}


# -----------------------------------------------------------------------------
# CRPS wrappers
#
# scoringRules convention: lower score = better forecast (loss-oriented).
# seqcomp convention:      higher score = better forecast (reward-oriented).
#
# Each wrapper checks inputs and negates the scoringRules output.
# -----------------------------------------------------------------------------


#' Negated CRPS for normal predictive distributions
#'
#' Computes the Continuous Ranked Probability Score for a normal predictive
#' distribution using `scoringRules::crps_norm()` and negates it so that higher
#' values are better.
#'
#' @param mu    Numeric vector. Location parameters (conditional means).
#' @param sigma Numeric vector. Scale parameters (conditional SDs, > 0).
#' @param x     Numeric vector. Realised observations.
#'
#' @return Numeric vector of CRPS values in `(-Inf, 0]` (negated loss).
#'
#' @details
#'   Calls \code{scoringRules::crps_norm(y = x, mean = mu, sd = sigma)} and
#'   negates. Use for GARCH(1,1)-norm forecasts where mu is the conditional
#'   mean and sigma is the conditional standard deviation.
#'
#' @examples
#' if (requireNamespace("scoringRules", quietly = TRUE)) {
#'   crps_normal(mu = c(0, 1), sigma = c(1, 2), x = c(0.2, 1.3))
#' }
#'
#' @export
crps_normal <- function(mu, sigma, x) {
  if (!requireNamespace("scoringRules", quietly = TRUE))
    stop("Package 'scoringRules' required. Install with: ",
         "install.packages('scoringRules')")
  stopifnot(
    length(mu)    == length(x),
    length(sigma) == length(x),
    all(sigma > 0)
  )
  -scoringRules::crps_norm(y = x, mean = mu, sd = sigma)
}


#' Negated CRPS for empirical predictive distributions
#'
#' Wrapper over \code{scoringRules::crps_sample} using method = "edf"
#' (empirical distribution function, O(n log n) via quantile decomposition
#' of Laio & Tamea, 2007). Positively oriented: higher = better.
#'
#' @param ensemble Matrix. T x n matrix of forecast draws. Each row
#'                 corresponds to one observation in y and comprises
#'                 n simulation draws from the predictive distribution.
#'                 For Historical Simulation: each row is the past
#'                 WINDOW returns.
#' @param y        Numeric vector of length T. Realised observations.
#'
#' @return Numeric vector of length T of CRPS values in `(-Inf, 0]`
#'         (negated loss).
#'
#' @details
#'   Requires \code{nrow(ensemble) == length(y)}. Passes
#'   \code{dat = ensemble} directly to \code{crps_sample} which handles
#'   vectorisation over rows natively. \code{show_messages} is suppressed
#'   as the "edf" method requires no bandwidth selection messages.
#'
#' @examples
#' if (requireNamespace("scoringRules", quietly = TRUE)) {
#'   ensemble <- matrix(c(0.1, 0.2, 0.3, 1.0, 1.1, 1.2), nrow = 2, byrow = TRUE)
#'   crps_empirical(ensemble, y = c(0.25, 1.05))
#' }
#'
#' @export
crps_empirical <- function(ensemble, y) {
  if (!requireNamespace("scoringRules", quietly = TRUE))
    stop("Package 'scoringRules' required. Install with: ",
         "install.packages('scoringRules')")
  stopifnot(
    is.matrix(ensemble),
    nrow(ensemble) == length(y)
  )
  -scoringRules::crps_sample(y   = y,
                             dat = ensemble,
                             method        = "edf",
                             show_messages = FALSE)
}


#' Negated CRPS for Student-t predictive distributions
#'
#' Wrapper over \code{scoringRules::crps_t}. Positively oriented:
#' higher = better. The dof > 2 constraint ensures finite variance,
#' which is required for the CRPS to be well-defined for the t-distribution.
#'
#' @param mu    Numeric vector. Location parameters (conditional means).
#' @param sigma Numeric vector. Scale parameters (conditional SDs, > 0).
#' @param dof   Numeric vector or scalar. Degrees of freedom (> 2).
#'              May be scalar if constant across all observations (e.g.
#'              estimated once per rolling window).
#' @param x     Numeric vector. Realised observations.
#'
#' @return Numeric vector of CRPS values in `(-Inf, 0]` (negated loss).
#'
#' @details
#'   Calls \code{scoringRules::crps_t(y = x, df = dof, location = mu,
#'   scale = sigma)} and negates. Use for GARCH(1,1)-std forecasts where
#'   dof is the estimated degrees-of-freedom parameter from ugarchroll.
#'
#' @examples
#' if (requireNamespace("scoringRules", quietly = TRUE)) {
#'   crps_std(mu = c(0, 1), sigma = c(1, 2), dof = 5, x = c(0.2, 1.3))
#' }
#'
#' @export
crps_std <- function(mu, sigma, dof, x) {
  if (!requireNamespace("scoringRules", quietly = TRUE))
    stop("Package 'scoringRules' required. Install with: ",
         "install.packages('scoringRules')")
  dof <- rep_len(dof, length(x))
  stopifnot(
    length(mu)    == length(x),
    length(sigma) == length(x),
    all(sigma > 0),
    all(dof   > 2)
  )
  -scoringRules::crps_t(y        = x,
                        df       = dof,
                        location = mu,
                        scale    = sigma)
}

#' Spherical score for binary and categorical forecasts
#'
#' Computes the positively oriented spherical score. Vector probability input
#' is treated as binary; matrix probability input is treated as categorical.
#'
#' @param p Numeric vector in `[0, 1]` for binary forecasts, or a numeric
#'   matrix whose rows are probability vectors for categorical forecasts.
#' @param y For binary vector input, numeric vector in `{0, 1}`. For
#'   categorical matrix input, integer vector in `{1, ..., K}`, where
#'   `K = ncol(p)`.
#'
#' @return Numeric vector of scores in `[0, 1]`. Higher is better.
#'
#' @details For binary forecasts, this computes \deqn{S(p, y) =
#' \frac{py + (1-p)(1-y)}{\sqrt{p^2 + (1-p)^2}}.}
#' For categorical forecasts, this computes \deqn{S(\mathbf{p}, y) =
#' \frac{p_y}{\|\mathbf{p}\|_2},}
#' where `p_y` is the forecast probability assigned to the realised category.
#'
#' Score differences lie in `[-1, 1]`, so use `c = 1` for Theorem 1 and
#' `c = 2` for Theorems 2 and 3.
#'
#' @examples
#' p <- c(0.2, 0.7, 0.9)
#' y <- c(0, 1, 1)
#' spherical_score(p, y)
#'
#' @export
spherical_score <- function(p, y) {
  if (is.matrix(p)) {
    prep <- .preprocess_cat(p, y)
    idx  <- cbind(seq_len(prep$T), prep$y_int)
    p_y  <- prep$p_mat[idx]

    return(p_y / sqrt(rowSums(prep$p_mat^2)))
  }

  stopifnot(
    length(p) == length(y),
    all(p >= 0 & p <= 1),
    all(y %in% c(0, 1))
  )

  (p * y + (1 - p) * (1 - y)) / sqrt(p^2 + (1 - p)^2)
}

#' Score difference bounds for a named scoring rule
#'
#' Returns lo, hi and the derived scale parameters c_thm1, c_thm23 for the
#' score difference process
#'   hat_delta_t = S(p, y) - S(q, y),
#' in those cases where a genuine, theorem-valid bound is available.
#'
#' Convention (utils.R::score_diff_scales):
#'   c_thm1  = (hi - lo) / 2   # Theorem 1:      |delta_i| <= c
#'   c_thm23 = hi - lo        # Theorems 2 & 3: |delta_i| <= c/2
#'
#' @param scoring_rule Character. One of:
#' * `"brier"`, `"spherical"` — bounded, exact finite-sample `c`.
#' * `"winkler"` — descriptive helper for the one-sided CS on the log score.
#' * `"tick"` — unbounded; returns `NULL` with guidance.
#' * `"crps"`, `"crps_normal"`, `"crps_empirical"`, `"crps_std"` — unbounded;
#'   returns `NULL` with guidance.
#' * `"log"`, `"qlike"` — unbounded; returns `NULL` with guidance.
#'
#' @return Named list with elements lo, hi, c_thm1, c_thm23 for bounded
#'         rules, or NULL for unbounded rules (with an informative message).
#'
#' @section Per-rule notes:
#' * **Brier / Spherical** — individual scores lie in `[-1, 0]` (Brier) or
#'   `[0, 1]` (Spherical), so score differences lie in `[-1, 1]` either way.
#'   This bound is exact and yields finite-sample anytime-valid CS via
#'   Hoeffding/Bernstein.
#' * **Winkler** — bounded above by 1; the lower bound is problem-dependent,
#'   so `lo = -Inf` and only `hi = 1` is used, as a descriptive helper for
#'   the one-sided CS wrapper `winkler_cs()`. Not intended for generic
#'   Hoeffding/Bernstein use (Theorem 1 requires a finite symmetric interval).
#' * **Tick loss** — unbounded on general financial returns. Any bound
#'   derived from an empirical data range is ex-post and not
#'   filtration-respecting, so it cannot justify finite-sample anytime
#'   validity. Use `cs_asymptotic()` for tick comparisons.
#' * **CRPS** (normal, t, empirical) — unbounded, since both the predictive
#'   distributions and the realised outcomes are unbounded. A historical
#'   data range is again an ex-post surrogate and does not provide a
#'   theorem-valid `c` for Hoeffding/Bernstein. Use `cs_asymptotic()`, or
#'   supply genuine ex ante bounds in problem-specific code if available.
#' * **Log / QLIKE** — both unbounded. For binary log-score comparisons, use
#'   `winkler_score()` + `winkler_cs()` when the Winkler construction is
#'   appropriate. For categorical log-score, QLIKE, and other unbounded
#'   score differences, use `cs_asymptotic()`, or `eprocess_predictable()`
#'   only with genuine ex ante predictable bounds.
#'
#' @examples
#' score_bounds("brier")
#' score_bounds("winkler")
#'
#' @export
score_bounds <- function(scoring_rule) {

  # Helper: add lo, hi to the scales from score_diff_scales()
  full_bounds <- function(lo, hi) {
    c(list(lo = lo, hi = hi), score_diff_scales(lo, hi))
  }

  switch(
    scoring_rule,

    "spherical" = ,
    "brier" = {
      # Individual Brier score in [-1, 0], Spherical in [0, 1].
      # Both have score differences bounded in [-1, 1].
      full_bounds(lo = -1, hi = 1)
    },

    "winkler" = {
      # Descriptive only: used by winkler_cs() for one-sided log-score CS.
      # lo = -Inf, hi = 1; c_thm1 is not meaningful here.
      list(lo = -Inf, hi = 1, c_thm1 = NA, c_thm23 = 2)
    },

    "tick" = {
      message(
        "Tick loss is unbounded on general outcomes. ",
        "There is no theorem-valid fixed c for Hoeffding/Bernstein. ",
        "Use cs_asymptotic() for tick-based confidence sequences."
      )
      invisible(NULL)
    },

    "crps"           = ,
    "crps_normal"    = ,
    "crps_empirical" = ,
    "crps_std"       = {
      message(
        "CRPS differences are unbounded when both forecasts and outcomes ",
        "have unbounded support. ",
        "Empirical data ranges cannot be used to justify theorem-valid ",
        "Hoeffding/Bernstein bounds. ",
        "Use cs_asymptotic() for CRPS-based confidence sequences."
      )
      invisible(NULL)
    },

    "log" = {
      message(
        "Logarithmic score is unbounded. ",
        "For binary outcomes use winkler_score() + winkler_cs() when appropriate. ",
        "For categorical log-score or other unbounded score differences, use ",
        "cs_asymptotic(), or eprocess_predictable() only with genuine predictable bounds."
      )
      invisible(NULL)
    },

    "qlike" = {
      message(
        "QLIKE is unbounded. ",
        "Use cs_asymptotic() for QLIKE-based confidence sequences."
      )
      invisible(NULL)
    },

    {
      stop(
        "Unknown scoring_rule: '", scoring_rule, "'. ",
        "Supported: 'brier', 'spherical', 'winkler', 'tick', 'crps', 'crps_normal', ",
        "'crps_empirical', 'crps_std', 'log', 'qlike'."
      )
    }
  )
}
