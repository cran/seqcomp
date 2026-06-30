# =============================================================================
# tests/testthat/test_predictable_bounds.R
# =============================================================================

set.seed(42)
T_ <- 300

# ── A. STRUCTURAL TESTS ──────────────────────────────────────────────────────

xs       <- sample(c(-1, 1), T_, replace = TRUE)
c_const  <- rep(2, T_)    # |xs| <= c/2 = 1
lambda_c <- 0.5 / 2       # = 0.25, default for c_0 = 2

ep_pred <- eprocess_predictable(
  scores1 = xs,
  scores2 = rep(0, T_),
  c_seq   = c_const,
  lambda  = lambda_c
)

test_that("eprocess_predictable returns a well-formed data.frame", {
  expect_true(is.data.frame(ep_pred))
  expect_true(all(c("t", "e_pq", "e_qp", "log_e_pq", "log_e_qp",
                    "c_seq", "lambda_used") %in% names(ep_pred)))
  expect_equal(nrow(ep_pred), T_)
})

test_that("eprocess_predictable e-values are strictly positive and respect the clip ceiling", {
  expect_true(all(ep_pred$e_pq > 0) && all(ep_pred$e_qp > 0))
  expect_true(all(ep_pred$e_pq <= 1e7) && all(ep_pred$e_qp <= 1e7))
  expect_true(all(is.finite(ep_pred$log_e_pq)) && all(is.finite(ep_pred$log_e_qp)))
})

test_that("eprocess_predictable records lambda_used and c_seq correctly", {
  expect_true(all(ep_pred$lambda_used == lambda_c))
  expect_true(all(ep_pred$c_seq == 2))
})

test_that("eprocess_predictable satisfies sign symmetry", {
  ep_rev <- eprocess_predictable(
    scores1 = -xs,
    scores2 = rep(0, T_),
    c_seq   = c_const,
    lambda  = lambda_c
  )

  expect_equal(ep_pred$e_pq, ep_rev$e_qp, tolerance = 1e-8)
  expect_equal(ep_pred$e_qp, ep_rev$e_pq, tolerance = 1e-8)
})

test_that("NULL lambda defaults to 0.5 / c_0", {
  ep_default <- eprocess_predictable(
    scores1 = xs,
    scores2 = rep(0, T_),
    c_seq   = c_const
  )
  expect_true(all(ep_default$lambda_used == 0.25))
})

test_that("eprocess_predictable warns on predictable-bound violations", {
  xs_bad <- c(1.5, rep(0.5, T_ - 1))   # first obs violates |xs| <= c/2 = 1
  c_bad  <- rep(2, T_)

  expect_warning(
    eprocess_predictable(
      scores1 = xs_bad,
      scores2 = rep(0, T_),
      c_seq   = c_bad,
      lambda  = 0.1
    )
  )
})

test_that("strict = TRUE turns a bound violation into an error", {
  xs_bad <- c(1.5, rep(0.5, T_ - 1))
  c_bad  <- rep(2, T_)

  expect_error(
    eprocess_predictable(
      scores1 = xs_bad,
      scores2 = rep(0, T_),
      c_seq   = c_bad,
      lambda  = 0.1,
      strict  = TRUE
    )
  )
})

test_that("lambda must be strictly less than 1 / max(c_seq)", {
  expect_error(
    eprocess_predictable(
      scores1 = xs,
      scores2 = rep(0, T_),
      c_seq   = c_const,
      lambda  = 0.5   # invalid: 1 / max(c_seq) = 0.5
    )
  )
})

test_that("negative lambda raises an error", {
  expect_error(
    eprocess_predictable(
      scores1 = xs,
      scores2 = rep(0, T_),
      c_seq   = c_const,
      lambda  = -0.1
    )
  )
})

test_that("different c_seq produces a different e-process path under the same lambda", {
  set.seed(10)
  xs_pos <- abs(rnorm(T_, 0.3, 0.3))
  xs_pos <- pmin(xs_pos, 0.9)
  c_tight <- rep(2,  T_)
  c_loose <- rep(10, T_)
  lam_fair <- 0.05   # must satisfy lambda < 1 / max(c_seq) for both

  ep_tight <- eprocess_predictable(
    scores1 = xs_pos, scores2 = rep(0, T_), c_seq = c_tight, lambda = lam_fair
  )
  ep_loose <- eprocess_predictable(
    scores1 = xs_pos, scores2 = rep(0, T_), c_seq = c_loose, lambda = lam_fair
  )

  expect_gt(max(abs(ep_tight$e_pq - ep_loose$e_pq)), 1e-12)
})

# ── B. REJECTION SUMMARY ─────────────────────────────────────────────────────

test_that("predictable_rejections returns a correctly structured two-sided summary", {
  rej <- predictable_rejections(ep_pred, alpha = 0.05)

  expect_true(is.list(rej))
  expect_true(all(c("threshold", "tau_pq", "tau_qp", "reject_pq", "reject_qp",
                    "c_range", "lambda") %in% names(rej)))
  expect_equal(rej$threshold, 40)
  expect_true(is.logical(rej$reject_pq) && length(rej$reject_pq) == 1)
  expect_true(is.logical(rej$reject_qp) && length(rej$reject_qp) == 1)
  expect_identical(as.numeric(rej$c_range), c(2, 2))
  expect_equal(rej$lambda, lambda_c)
})

# ── C. OPTIONAL LIGHT SANITY CHECK UNDER NULL ────────────────────────────────
#
# N=100, T=100 is two orders of magnitude cheaper than the other coverage
# simulations in this package, so this stays at normal test speed.

test_that("light null rejection rate is not obviously inflated under the predictable e-process", {
  set.seed(123)
  N_sim     <- 100
  T_sim     <- 100
  alpha_sim <- 0.05
  threshold <- 2 / alpha_sim
  n_rej     <- 0L

  for (i in seq_len(N_sim)) {
    xs_i <- sample(c(-1, 1), T_sim, replace = TRUE)
    c_i  <- rep(2, T_sim)
    ep_i <- eprocess_predictable(
      scores1 = xs_i, scores2 = rep(0, T_sim), c_seq = c_i, lambda = 0.2
    )
    if (any(pmax(ep_i$e_pq, ep_i$e_qp) >= threshold)) n_rej <- n_rej + 1L
  }
  type1_light <- n_rej / N_sim

  expect_lte(type1_light, 0.15)
})
