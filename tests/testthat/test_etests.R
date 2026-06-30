# =============================================================================
# tests/testthat/test_etests.R
#
# Tests:
#   A. Structural: correct output, initialisation at 1, monotonicity under alt
#   B. Type I error simulation: e-process under null never exceeds 2/alpha
#      in more than alpha fraction of paths (Markov inequality guarantee)
#   C. Power: e-process under alternative reaches threshold faster than chance
#
# B and C are gated behind SEQCOMP_RUN_SLOW_TESTS=true; run locally with:
#   SEQCOMP_RUN_SLOW_TESTS=true Rscript -e 'devtools::test()'
# =============================================================================

set.seed(42)
s1 <- brier_score(runif(300), rbinom(300, 1, 0.5))
s2 <- brier_score(runif(300), rbinom(300, 1, 0.5))

# ── A. STRUCTURAL TESTS ──────────────────────────────────────────────────────

test_that("eprocess returns a well-formed data.frame", {
  ep <- eprocess(s1, s2, alpha = 0.05, c = 2)

  expect_true(is.data.frame(ep))
  expect_true(all(c("t", "e_pq", "e_qp", "log_e_pq", "log_e_qp") %in% names(ep)))
  expect_equal(nrow(ep), 300)
  expect_equal(ep$t, 1:300)
})

test_that("eprocess e-values are positive and respect the clip ceiling", {
  ep <- eprocess(s1, s2, alpha = 0.05, c = 2)

  expect_true(all(ep$e_pq > 0) && all(ep$e_qp > 0))
  expect_true(all(ep$e_pq <= 1e7) && all(ep$e_qp <= 1e7))
})

test_that("log_e_pq matches log(e_pq)", {
  ep <- eprocess(s1, s2, alpha = 0.05, c = 2)
  expect_equal(ep$log_e_pq, log(ep$e_pq), tolerance = 1e-8)
})

test_that("e-process initialises near 1 for a near-zero first observation", {
  xs_small <- c(0.01, rep(0, 299))
  ep_small <- eprocess(xs_small, rep(0, 300), alpha = 0.05, c = 2)
  expect_equal(ep_small$e_pq[1], 1, tolerance = 0.1)
})

test_that("eprocess satisfies sign symmetry: e_pq(p,q) == e_qp(q,p)", {
  ep_fwd <- eprocess(s1, s2, alpha = 0.05, c = 2)
  ep_rev <- eprocess(s2, s1, alpha = 0.05, c = 2)

  expect_equal(ep_fwd$e_pq, ep_rev$e_qp, tolerance = 1e-8)
  expect_equal(ep_fwd$e_qp, ep_rev$e_pq, tolerance = 1e-8)
})

test_that("eprocess_rejections detects a clear alternative", {
  xs_alt <- rep(0.5, 300)   # constant positive differences: strong evidence for p
  ep_alt <- eprocess(xs_alt, rep(0, 300), alpha = 0.05, c = 2)
  rej    <- eprocess_rejections(ep_alt, alpha = 0.05)

  expect_equal(rej$threshold, 40)
  expect_true(rej$reject_pq)
  expect_false(is.na(rej$tau_pq))
  expect_false(rej$reject_qp)
})

test_that("eprocess_rejections does not reject on a null path", {
  xs_null_path <- rep(0, 300)
  ep_null_path <- eprocess(xs_null_path, rep(0, 300), alpha = 0.05, c = 2)
  rej_null     <- eprocess_rejections(ep_null_path, alpha = 0.05)

  expect_false(rej_null$reject_pq)
})

# ── B + C1. TYPE I ERROR AND POWER SIMULATION ────────────────────────────────
#
# Under H_0^w(p,q): Delta_t <= 0 for all t.
# By Markov's inequality: P(exists t: E_t >= 1/alpha_os) <= alpha_os = alpha/2
# Two-sided: P(max(e_pq, e_qp) >= 2/alpha) <= alpha.
#
# Power is then checked against the Type I error rate measured here, so both
# simulations live in one test (power_est > type1_pq mirrors the original
# script's direct comparison).

test_that("eprocess controls Type I error under the sharp null and has power exceeding it under the alternative", {
  skip_on_cran()
  skip_if(Sys.getenv("SEQCOMP_RUN_SLOW_TESTS") != "true",
          "slow Type I error / power simulation; set SEQCOMP_RUN_SLOW_TESTS=true to run")

  alpha_sim <- 0.05
  threshold <- 2 / alpha_sim   # = 40

  # --- Type I error under the sharp null (Rademacher, mean 0) ---
  N     <- 1000
  T_sim <- 500
  n_reject_pq <- 0L
  n_reject_qp <- 0L
  set.seed(42)
  for (i in seq_len(N)) {
    xs <- sample(c(-1, 1), T_sim, replace = TRUE)
    ep_i <- eprocess(scores1 = xs, scores2 = rep(0, T_sim),
                     alpha = alpha_sim, c = 2)
    if (any(ep_i$e_pq >= threshold)) n_reject_pq <- n_reject_pq + 1L
    if (any(ep_i$e_qp >= threshold)) n_reject_qp <- n_reject_qp + 1L
  }
  type1_pq <- n_reject_pq / N
  type1_qp <- n_reject_qp / N

  expect_lte(type1_pq, alpha_sim)
  expect_lte(type1_qp, alpha_sim)

  # --- Power under the alternative (delta = 0.3) ---
  N_pow <- 500
  T_pow <- 1000
  delta <- 0.3
  n_power <- 0L
  set.seed(99)
  for (i in seq_len(N_pow)) {
    xs <- runif(T_pow, -1 + delta, 1 + delta)
    xs <- pmax(pmin(xs, 1), -1)
    ep_i <- eprocess(scores1 = xs, scores2 = rep(0, T_pow),
                     alpha = alpha_sim, c = 2)
    if (any(ep_i$e_pq >= threshold)) n_power <- n_power + 1L
  }
  power_est <- n_power / N_pow

  expect_gt(power_est, 0.5)
  expect_gt(power_est, type1_pq)  # test is informative: power beats the null rejection rate
})

# ── C2. POWER CURVE ACROSS DELTA VALUES ──────────────────────────────────────

test_that("eprocess power increases with the effect size delta", {
  skip_on_cran()
  skip_if(Sys.getenv("SEQCOMP_RUN_SLOW_TESTS") != "true",
          "slow power curve simulation; set SEQCOMP_RUN_SLOW_TESTS=true to run")

  alpha_sim <- 0.05
  threshold <- 2 / alpha_sim
  deltas    <- c(0.1, 0.2, 0.3, 0.5)
  N_curve   <- 200
  T_curve   <- 1000

  power_curve <- numeric(length(deltas))
  set.seed(77)
  for (j in seq_along(deltas)) {
    d <- deltas[j]
    n_rej <- 0L
    for (i in seq_len(N_curve)) {
      xs   <- pmax(pmin(runif(T_curve, -1 + d, 1 + d), 1), -1)
      ep_i <- eprocess(xs, rep(0, T_curve), alpha = alpha_sim, c = 2)
      if (any(ep_i$e_pq >= threshold)) n_rej <- n_rej + 1L
    }
    power_curve[j] <- n_rej / N_curve
  }

  expect_true(all(power_curve >= 0 & power_curve <= 1))
  expect_gt(power_curve[length(power_curve)], power_curve[1])
})
