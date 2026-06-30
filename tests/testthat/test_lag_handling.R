# =============================================================================
# tests/testthat/test_lag_handling.R
# =============================================================================

# --- A. split_streams ---

test_that("split_streams partitions a sequence into h interleaved streams", {
  xs <- 1:12
  streams <- split_streams(xs, h = 3)

  expect_equal(length(streams), 3)
  expect_equal(streams[[1]], c(1, 4, 7, 10))
  expect_equal(streams[[2]], c(2, 5, 8, 11))
  expect_equal(streams[[3]], c(3, 6, 9, 12))
})

test_that("split_streams with h=1 returns a single stream equal to the input", {
  xs <- 1:12
  streams1 <- split_streams(xs, h = 1)

  expect_equal(length(streams1), 1)
  expect_equal(streams1[[1]], xs)
})

# --- B. calibrate_p_to_e ---

test_that("calibrate_p_to_e produces non-negative, decreasing e-values", {
  p_vals <- c(0.001, 0.01, 0.05, 0.1, 0.5, 1.0)
  e_vals <- calibrate_p_to_e(p_vals)

  expect_true(all(e_vals >= 0))
  expect_true(all(diff(e_vals) <= 0))
  expect_lt(e_vals[length(e_vals)], 1e-3)  # p=1 gives e near 0
})

test_that("calibrate_p_to_e with strategy='simple' produces non-negative, decreasing e-values", {
  p_vals <- c(0.001, 0.01, 0.05, 0.1, 0.5, 1.0)
  e_simple <- calibrate_p_to_e(p_vals, strategy = "simple")

  expect_true(all(e_simple >= 0))
  expect_true(all(diff(e_simple) <= 0))
})

# --- C. eprocess_lag h=1 matches eprocess ---

test_that("eprocess_lag with h=1 matches eprocess() exactly", {
  set.seed(1)
  T_ <- 200
  xs <- sample(c(-1, 1), T_, replace = TRUE)

  ep_standard <- eprocess(xs, rep(0, T_), alpha = 0.05, c = 2)
  ep_lag1     <- eprocess_lag(xs, rep(0, T_), h = 1, alpha = 0.05, c = 2)

  expect_equal(ep_lag1$e_pq, ep_standard$e_pq, tolerance = 1e-8)
  expect_equal(ep_lag1$e_qp, ep_standard$e_qp, tolerance = 1e-8)
})

# --- D. eprocess_lag h=3 structure ---

test_that("eprocess_lag with h=3 returns correctly shaped, positive e-values under both null types", {
  set.seed(1)
  T_ <- 200
  xs <- sample(c(-1, 1), T_, replace = TRUE)

  ep_lag3_pw <- eprocess_lag(xs, rep(0, T_), h = 3, alpha = 0.05, c = 2, null = "pw")
  ep_lag3_w  <- eprocess_lag(xs, rep(0, T_), h = 3, alpha = 0.05, c = 2, null = "w")

  expect_equal(nrow(ep_lag3_pw), T_)
  expect_equal(nrow(ep_lag3_w),  T_)
  expect_true(all(ep_lag3_pw$e_pq > 0))
  expect_true(all(ep_lag3_w$e_pq  > 0))
})

test_that("pw-null e-values are on average at least as large as w-null e-values (less conservative)", {
  set.seed(1)
  T_ <- 200
  xs <- sample(c(-1, 1), T_, replace = TRUE)

  ep_lag3_pw <- eprocess_lag(xs, rep(0, T_), h = 3, alpha = 0.05, c = 2, null = "pw")
  ep_lag3_w  <- eprocess_lag(xs, rep(0, T_), h = 3, alpha = 0.05, c = 2, null = "w")

  expect_gte(mean(ep_lag3_pw$e_pq), mean(ep_lag3_w$e_pq))
})

# --- E. Type I error h=3 simulation ---
#
# Gated: N=500 paths x two eprocess_lag() calls each is too slow for routine
# devtools::check(). Run locally with:
#   SEQCOMP_RUN_SLOW_TESTS=true Rscript -e 'devtools::test()'

test_that("eprocess_lag (h=3) controls Type I error under both null types", {
  skip_on_cran()
  skip_if(Sys.getenv("SEQCOMP_RUN_SLOW_TESTS") != "true",
          "slow Type I error simulation; set SEQCOMP_RUN_SLOW_TESTS=true to run")

  N_sim <- 500
  T_sim <- 300
  threshold <- 2 / 0.05
  n_rej_pw <- 0L
  n_rej_w  <- 0L

  set.seed(77)
  for (i in seq_len(N_sim)) {
    xs_i <- sample(c(-1, 1), T_sim, replace = TRUE)

    ep_pw <- eprocess_lag(xs_i, rep(0, T_sim), h = 3, alpha = 0.05, c = 2, null = "pw")
    ep_w  <- eprocess_lag(xs_i, rep(0, T_sim), h = 3, alpha = 0.05, c = 2, null = "w")

    if (any(ep_pw$e_pq >= threshold)) n_rej_pw <- n_rej_pw + 1L
    if (any(ep_w$e_pq  >= threshold)) n_rej_w  <- n_rej_w  + 1L
  }

  t1_pw <- n_rej_pw / N_sim
  t1_w  <- n_rej_w  / N_sim

  expect_lte(t1_pw, 0.05)
  expect_lte(t1_w,  0.05)
})

# --- F. pw combination grows under alternative (T=2000) ---

test_that("pw e_pq grows beyond 1 under a sustained alternative at T=2000", {
  set.seed(55)
  T_long <- 2000
  xs_alt_long <- pmax(pmin(runif(T_long, -0.7, 1.3), 1), -1)  # mean difference ~0.3

  ep_alt_long <- eprocess_lag(xs_alt_long, rep(0, T_long),
                              h = 3, alpha = 0.05, c = 2, null = "pw")

  expect_gt(max(ep_alt_long$e_pq), 1)
})
