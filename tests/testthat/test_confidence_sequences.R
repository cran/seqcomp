# =============================================================================
# tests/testthat/test_confidence_sequences.R
#
# Tests fall into two categories:
#   A. Structural tests: correct output format, monotone width, etc.
#   B. Coverage simulation: the key anytime-valid guarantee.
#      Run 1000 paths under the null; verify that the CS never exits
#      in >= (1 - alpha) of paths, uniformly over all t.
#      These are gated behind SEQCOMP_RUN_SLOW_TESTS=true since they're
#      too heavy for ordinary devtools::check() runs.
# =============================================================================

set.seed(42)
s1 <- brier_score(runif(200), rbinom(200, 1, 0.5))
s2 <- brier_score(runif(200), rbinom(200, 1, 0.5))

# ── A. STRUCTURAL TESTS ──────────────────────────────────────────────────────

test_that("cs_hoeffding returns a well-formed data.frame", {
  cs_h <- cs_hoeffding(s1, s2, alpha = 0.05, c = 1)

  expect_true(is.data.frame(cs_h))
  expect_true(all(c("t", "estimate", "lower", "upper") %in% names(cs_h)))
  expect_equal(nrow(cs_h), 200)
  expect_equal(cs_h$t, 1:200)
})

test_that("cs_hoeffding bounds are ordered and have positive width", {
  cs_h <- cs_hoeffding(s1, s2, alpha = 0.05, c = 1)

  expect_true(all(cs_h$lower < cs_h$estimate & cs_h$estimate < cs_h$upper))
  expect_true(all(cs_h$upper - cs_h$lower > 0))
})

test_that("cs_hoeffding width shrinks over time", {
  cs_h <- cs_hoeffding(s1, s2, alpha = 0.05, c = 1)
  widths <- cs_h$upper - cs_h$lower
  # Allow for slight non-monotonicity early on (t=1 to t=2 can wobble due to rho);
  # check the broad trend instead of strict monotonicity.
  expect_lt(mean(widths[151:200]), mean(widths[1:50]))
})

test_that("cs_bernstein returns a well-formed data.frame", {
  cs_eb <- cs_bernstein(s1, s2, alpha = 0.05, c = 2)

  expect_true(is.data.frame(cs_eb))
  expect_true(all(c("t", "estimate", "lower", "upper") %in% names(cs_eb)))
  expect_equal(nrow(cs_eb), 200)
  expect_true(all(cs_eb$lower < cs_eb$estimate & cs_eb$estimate < cs_eb$upper))
  expect_true(all(cs_eb$upper - cs_eb$lower > 0))
})

test_that("EB CS is tighter than Hoeffding under low-variance score differences", {
  set.seed(7)
  p1 <- pmin(pmax(rnorm(500, 0.5, 0.05), 0.01), 0.99)
  p2 <- pmin(pmax(rnorm(500, 0.5, 0.05), 0.01), 0.99)
  y  <- rbinom(500, 1, 0.5)
  s1_lv <- brier_score(p1, y)
  s2_lv <- brier_score(p2, y)

  cs_h_lv  <- cs_hoeffding(s1_lv, s2_lv, alpha = 0.05, c = 1)
  cs_eb_lv <- cs_bernstein(s1_lv, s2_lv, alpha = 0.05, c = 2)

  width_h  <- cs_h_lv$upper  - cs_h_lv$lower
  width_eb <- cs_eb_lv$upper - cs_eb_lv$lower

  expect_lt(mean(width_eb), mean(width_h))
  expect_lt(width_eb[500], width_h[500])
})

test_that("lcb_only produces an upper bound of Inf and a finite lower bound", {
  cs_lcb <- cs_bernstein(s1, s2, alpha = 0.05, c = 2, lcb_only = TRUE)

  expect_true(all(is.infinite(cs_lcb$upper)))
  expect_true(all(is.finite(cs_lcb$lower)))
})

test_that("ucb_only produces a lower bound of -Inf and a finite upper bound", {
  cs_ucb <- cs_bernstein(s1, s2, alpha = 0.05, c = 2, ucb_only = TRUE)

  expect_true(all(is.infinite(cs_ucb$lower) & cs_ucb$lower < 0))
  expect_true(all(is.finite(cs_ucb$upper)))
})

test_that("hardcoded boundary (alpha=0.05, c=1) produces valid output", {
  cs_hc <- cs_bernstein(s1, s2, alpha = 0.05, c = 1, boundary = "hardcoded")

  expect_true(all(is.finite(cs_hc$lower)))
  expect_true(all(is.finite(cs_hc$upper)))
  expect_true(all(cs_hc$lower < cs_hc$upper))
})

test_that("hardcoded boundary warns and falls back to the mixture boundary for unsupported alpha/c", {
  expect_warning(
    cs_warn <- cs_bernstein(s1, s2, alpha = 0.10, c = 2, boundary = "hardcoded"),
    regexp = "hardcoded"
  )
  expect_true(all(is.finite(cs_warn$lower)))
})

test_that("cs_asymptotic returns a well-formed data.frame", {
  cs_asy <- cs_asymptotic(s1, s2, alpha = 0.05)

  expect_true(is.data.frame(cs_asy))
  expect_true(all(c("t", "estimate", "lower", "upper") %in% names(cs_asy)))
  expect_equal(nrow(cs_asy), 200)
  expect_true(all(cs_asy$lower < cs_asy$estimate & cs_asy$estimate < cs_asy$upper))
})

test_that("cs_asymptotic width shrinks over time and stays finite", {
  cs_asy <- cs_asymptotic(s1, s2, alpha = 0.05)

  expect_lt(
    mean(cs_asy$upper[151:200] - cs_asy$lower[151:200]),
    mean(cs_asy$upper[1:50]    - cs_asy$lower[1:50])
  )
  expect_true(all(is.finite(cs_asy$lower)))
  expect_true(all(is.finite(cs_asy$upper)))
})

test_that("cs_asymptotic produces a finite, positive-width CS under low variance", {
  set.seed(7)
  p1 <- pmin(pmax(rnorm(500, 0.5, 0.05), 0.01), 0.99)
  p2 <- pmin(pmax(rnorm(500, 0.5, 0.05), 0.01), 0.99)
  y  <- rbinom(500, 1, 0.5)
  s1_lv <- brier_score(p1, y)
  s2_lv <- brier_score(p2, y)

  cs_asy_lv <- cs_asymptotic(s1_lv, s2_lv, alpha = 0.05)
  width_asy <- cs_asy_lv$upper - cs_asy_lv$lower

  expect_true(all(width_asy > 0))
  expect_true(all(is.finite(width_asy)))
})

# ── B. COVERAGE SIMULATION ───────────────────────────────────────────────────
#
# Key property: for any stopping time tau, P(Delta_tau not in C_tau) <= alpha.
# Stronger: P(exists t: Delta_t not in C_t) <= alpha (uniform over all t).
#
# Simulation design:
#   - T = 500 observations per path
#   - N = 1000 paths
#   - Null: Delta_t = 0 (both forecasters identical, score diffs iid)
#   - Coverage = fraction of paths where CS contains 0 at ALL time points
#   - Expected: coverage >= 1 - alpha = 0.95
#
# Run locally with: SEQCOMP_RUN_SLOW_TESTS=true Rscript -e 'devtools::test()'

test_that("cs_hoeffding achieves nominal (1 - alpha) uniform coverage under the null", {
  skip_on_cran()
  skip_if(Sys.getenv("SEQCOMP_RUN_SLOW_TESTS") != "true",
          "slow coverage simulation; set SEQCOMP_RUN_SLOW_TESTS=true to run")

  N   <- 1000
  T_  <- 500
  alpha_sim <- 0.05
  c_sim     <- 1   # |delta_i| <= 1 for Brier differences

  n_covered_h <- 0L
  set.seed(100)
  for (i in seq_len(N)) {
    xs <- sample(c(-1, 1), T_, replace = TRUE)  # Rademacher: saturates c=1 bound
    cs_i <- cs_hoeffding(scores1 = xs, scores2 = rep(0, T_),
                         alpha = alpha_sim, c = c_sim)
    if (all(cs_i$lower <= 0 & 0 <= cs_i$upper)) n_covered_h <- n_covered_h + 1L
  }
  coverage_h <- n_covered_h / N

  expect_gte(coverage_h, 1 - alpha_sim)
  expect_lt(coverage_h, 1.00)  # not trivially conservative
})

test_that("cs_bernstein achieves nominal (1 - alpha) uniform coverage under the null", {
  skip_on_cran()
  skip_if(Sys.getenv("SEQCOMP_RUN_SLOW_TESTS") != "true",
          "slow coverage simulation; set SEQCOMP_RUN_SLOW_TESTS=true to run")

  N   <- 1000
  T_  <- 500
  alpha_sim <- 0.05

  n_covered_eb <- 0L
  set.seed(200)
  for (i in seq_len(N)) {
    xs <- sample(c(-1, 1), T_, replace = TRUE)
    cs_i <- cs_bernstein(scores1 = xs, scores2 = rep(0, T_),
                         alpha = alpha_sim, c = 2)  # c = hi - lo = 1 - (-1)
    if (all(cs_i$lower <= 0 & 0 <= cs_i$upper)) n_covered_eb <- n_covered_eb + 1L
  }
  coverage_eb <- n_covered_eb / N

  expect_gte(coverage_eb, 1 - alpha_sim)
  expect_lt(coverage_eb, 1.00)
})

test_that("EB CS is substantially tighter than Hoeffding under low-variance data", {
  T_ <- 500
  set.seed(300)
  xs_lv <- rnorm(T_, mean = 0, sd = 0.1)
  xs_lv <- pmax(pmin(xs_lv, 1), -1)

  cs_h_cmp  <- cs_hoeffding(xs_lv, rep(0, T_), alpha = 0.05, c = 1)
  cs_eb_cmp <- cs_bernstein(xs_lv, rep(0, T_), alpha = 0.05, c = 2)

  width_h_cmp  <- cs_h_cmp$upper  - cs_h_cmp$lower
  width_eb_cmp <- cs_eb_cmp$upper - cs_eb_cmp$lower

  expect_lt(mean(width_eb_cmp) / mean(width_h_cmp), 0.80)
})

test_that("cs_asymptotic achieves nominal (1 - alpha) uniform coverage under the null", {
  skip_on_cran()
  skip_if(Sys.getenv("SEQCOMP_RUN_SLOW_TESTS") != "true",
          "slow coverage simulation; set SEQCOMP_RUN_SLOW_TESTS=true to run")

  N   <- 1000
  T_  <- 500
  alpha_sim <- 0.05

  n_covered_asy <- 0L
  set.seed(400)
  for (i in seq_len(N)) {
    xs <- sample(c(-1, 1), T_, replace = TRUE)
    cs_i <- cs_asymptotic(xs, rep(0, T_), alpha = alpha_sim)
    if (all(cs_i$lower <= 0 & 0 <= cs_i$upper)) n_covered_asy <- n_covered_asy + 1L
  }
  coverage_asy <- n_covered_asy / N

  expect_gte(coverage_asy, 0.95)
})
