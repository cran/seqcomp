# =============================================================================
# tests/testthat/test_scores.R
# =============================================================================

# --- 1. Brier score ---

test_that("brier_score: perfect forecasts score 0", {
  expect_equal(brier_score(1, 1), 0)
  expect_equal(brier_score(0, 0), 0)
})

test_that("brier_score: worst forecasts score -1", {
  expect_equal(brier_score(0, 1), -1)
  expect_equal(brier_score(1, 0), -1)
})

test_that("brier_score: p=0.5 gives -0.25", {
  expect_equal(brier_score(0.5, 1), -0.25)
})

test_that("brier_score is positively oriented", {
  expect_gt(brier_score(0.8, 1), brier_score(0.2, 1))
})

test_that("brier_score values are bounded in [-1, 0]", {
  set.seed(101)
  p <- seq(0, 1, by = 0.01)
  y <- rbinom(length(p), 1, 0.5)
  scores <- brier_score(p, y)
  expect_true(all(scores >= -1 & scores <= 0))
})

# --- 2. Log score ---

test_that("log_score: y=1 at p=1 is near 0", {
  expect_gt(log_score(1, 1), -1e-10)
})

test_that("log_score matches log(p) at p=0.5 for both outcomes", {
  expect_equal(log_score(0.5, 1), log(0.5), tolerance = 1e-10)
  expect_equal(log_score(0.5, 0), log(0.5), tolerance = 1e-10)
})

test_that("log_score is positively oriented", {
  expect_gt(log_score(0.8, 1), log_score(0.3, 1))
})

test_that("log_score eps clipping prevents -Inf", {
  expect_true(is.finite(log_score(0, 1, eps = 1e-15)))
})

# --- 3. Tick loss ---

test_that("tick_loss is finite for exceedance and non-exceedance", {
  q_test <- -0.02
  alpha_test <- 0.05
  expect_true(is.finite(tick_loss(q_test, -0.05, alpha_test)))
  expect_true(is.finite(tick_loss(q_test,  0.01, alpha_test)))
})

test_that("tick_loss is 0 when forecast equals outcome exactly", {
  expect_equal(tick_loss(0.0, 0.0, 0.05), 0, tolerance = 1e-12)
})

test_that("tick_loss propriety: correct quantile scores better than a biased one", {
  set.seed(1)
  n_prop <- 5000
  y_prop <- rnorm(n_prop)
  q_correct <- qnorm(0.05)
  q_wrong   <- qnorm(0.05) + 0.5

  mean_correct <- mean(tick_loss(rep(q_correct, n_prop), y_prop, 0.05))
  mean_wrong   <- mean(tick_loss(rep(q_wrong,   n_prop), y_prop, 0.05))

  expect_gt(mean_correct, mean_wrong)
})

# --- 4. QLIKE score ---

test_that("qlike_score: perfect forecast scores 0", {
  expect_equal(qlike_score(1.0, 1.0), 0, tolerance = 1e-12)
})

test_that("qlike_score is bounded above by 0", {
  expect_lte(qlike_score(2.0, 1.0), 0)
  expect_lte(qlike_score(0.5, 1.0), 0)
})

test_that("qlike_score is positively oriented", {
  sigma2_true <- 0.04
  expect_gt(qlike_score(sigma2_true, sigma2_true),
            qlike_score(sigma2_true * 2, sigma2_true))
})

test_that("qlike_score propriety: correct variance beats biased forecast on average", {
  set.seed(2)
  n_q <- 5000
  sigma2_true_val <- 0.04
  eps_q <- rnorm(n_q, sd = sqrt(sigma2_true_val))
  rv    <- eps_q^2

  mean_q_correct <- mean(qlike_score(rep(sigma2_true_val, n_q), rv))
  mean_q_wrong   <- mean(qlike_score(rep(sigma2_true_val * 1.5, n_q), rv))

  expect_gt(mean_q_correct, mean_q_wrong)
})

# --- 5. Winkler score ---

test_that("winkler_score equals 1 at the favourable outcome", {
  p_w <- c(0.8, 0.3, 0.6)
  q_w <- c(0.5, 0.6, 0.4)
  y_fav <- as.numeric(p_w > q_w)
  w_fav <- winkler_score(p_w, q_w, y_fav)
  expect_equal(w_fav, rep(1, 3), tolerance = 1e-8)
})

test_that("winkler_score is <= 0 at the unfavourable outcome", {
  p_w <- c(0.8, 0.3, 0.6)
  q_w <- c(0.5, 0.6, 0.4)
  y_fav   <- as.numeric(p_w > q_w)
  y_unfav <- 1 - y_fav
  w_unfav <- winkler_score(p_w, q_w, y_unfav)
  expect_true(all(w_unfav <= 1e-10))
})

test_that("winkler_score is bounded above by 1", {
  set.seed(3)
  p_rand <- runif(1000, 0.05, 0.95)
  q_rand <- runif(1000, 0.05, 0.95)
  y_rand <- rbinom(1000, 1, 0.5)
  w_rand <- winkler_score(p_rand, q_rand, y_rand)
  expect_true(all(w_rand <= 1 + 1e-10))
})

test_that("winkler_score follows the 0/0 convention when p == q", {
  p_eq <- rep(0.5, 5); q_eq <- rep(0.5, 5); y_eq <- c(0, 1, 0, 1, 0)
  w_eq <- winkler_score(p_eq, q_eq, y_eq)
  expect_equal(w_eq, rep(0, 5), tolerance = 1e-6)
})

# --- 6. crps_normal ---

test_that("crps_normal: perfect deterministic forecast is near 0", {
  skip_if_not_installed("scoringRules")
  expect_equal(crps_normal(mu = 0, sigma = 1e-6, x = 0), 0, tolerance = 1e-4)
})

test_that("crps_normal: all values are <= 0", {
  skip_if_not_installed("scoringRules")
  set.seed(5)
  mu_t  <- rnorm(100)
  sig_t <- abs(rnorm(100)) + 0.1
  x_t   <- rnorm(100)
  expect_true(all(crps_normal(mu_t, sig_t, x_t) <= 1e-10))
})

test_that("crps_normal: correct sigma scores better than overdispersed", {
  skip_if_not_installed("scoringRules")
  set.seed(6)
  y_t       <- rnorm(1000)
  crps_good <- mean(crps_normal(rep(0, 1000), rep(1, 1000), y_t))
  crps_bad  <- mean(crps_normal(rep(0, 1000), rep(3, 1000), y_t))
  expect_gt(crps_good, crps_bad)
})

test_that("crps_normal: unbiased forecast scores better than biased", {
  skip_if_not_installed("scoringRules")
  set.seed(6)
  y_t          <- rnorm(1000)
  crps_centred <- mean(crps_normal(y_t,     rep(1, 1000), y_t))
  crps_biased  <- mean(crps_normal(y_t + 1, rep(1, 1000), y_t))
  expect_gt(crps_centred, crps_biased)
})

test_that("crps_normal: sigma <= 0 raises an error", {
  skip_if_not_installed("scoringRules")
  expect_error(crps_normal(0, -1, 0))
})

test_that("crps_normal: length mismatch raises an error", {
  skip_if_not_installed("scoringRules")
  expect_error(crps_normal(c(0, 0), c(1, 1), c(0, 0, 0)))
})

# --- 7. crps_empirical ---

test_that("crps_empirical agrees with crps_normal for a large N(0,1) ensemble", {
  skip_if_not_installed("scoringRules")
  set.seed(7)
  n_ens   <- 10000
  y_obs   <- rnorm(20)
  ens_mat <- matrix(rnorm(20 * n_ens, mean = 0, sd = 1), nrow = 20, ncol = n_ens)

  crps_norm_vals <- crps_normal(rep(0, 20), rep(1, 20), y_obs)
  crps_emp_vals  <- crps_empirical(ens_mat, y_obs)

  expect_lt(max(abs(crps_emp_vals - crps_norm_vals)), 0.05)
})

test_that("crps_empirical: all values are <= 0", {
  skip_if_not_installed("scoringRules")
  set.seed(8)
  small_ens <- matrix(rnorm(50 * 250), nrow = 50, ncol = 250)
  y_small   <- rnorm(50)
  expect_true(all(crps_empirical(small_ens, y_small) <= 1e-10))
})

test_that("crps_empirical: centred ensemble scores better than shifted ensemble", {
  skip_if_not_installed("scoringRules")
  set.seed(9)
  y_ref     <- rnorm(100)
  ens_good  <- matrix(rnorm(100 * 2000, mean = 0, sd = 1), nrow = 100)
  ens_shift <- matrix(rnorm(100 * 2000, mean = 2, sd = 1), nrow = 100)

  expect_gt(mean(crps_empirical(ens_good,  y_ref)),
            mean(crps_empirical(ens_shift, y_ref)))
})

test_that("crps_empirical: non-matrix ensemble raises an error", {
  skip_if_not_installed("scoringRules")
  set.seed(9)
  y_ref    <- rnorm(100)
  ens_good <- matrix(rnorm(100 * 2000, mean = 0, sd = 1), nrow = 100)
  expect_error(crps_empirical(as.data.frame(ens_good), y_ref))
})

test_that("crps_empirical: nrow/length mismatch raises an error", {
  skip_if_not_installed("scoringRules")
  set.seed(9)
  y_ref    <- rnorm(100)
  ens_good <- matrix(rnorm(100 * 2000, mean = 0, sd = 1), nrow = 100)
  expect_error(crps_empirical(ens_good, y_ref[1:50]))
})

# --- 8. crps_std ---

test_that("crps_std: perfect deterministic forecast is near 0 at large df", {
  skip_if_not_installed("scoringRules")
  expect_equal(crps_std(mu = 0, sigma = 1e-6, dof = 30, x = 0), 0, tolerance = 1e-4)
})

test_that("crps_std: all values are <= 0", {
  skip_if_not_installed("scoringRules")
  set.seed(10)
  mu_s  <- rnorm(100)
  sig_s <- abs(rnorm(100)) + 0.1
  dof_s <- runif(100, min = 3, max = 30)
  x_s   <- rnorm(100)
  expect_true(all(crps_std(mu_s, sig_s, dof_s, x_s) <= 1e-10))
})

test_that("crps_std converges to crps_normal as dof grows large", {
  skip_if_not_installed("scoringRules")
  set.seed(11)
  y_conv   <- rnorm(200)
  crps_n   <- crps_normal(rep(0, 200), rep(1, 200), y_conv)
  crps_t30 <- crps_std(rep(0, 200), rep(1, 200), dof = 1000, y_conv)
  expect_lt(max(abs(crps_t30 - crps_n)), 1e-3)
})

test_that("crps_std: higher df scores better than lower df on normal data", {
  skip_if_not_installed("scoringRules")
  set.seed(12)
  y_norm  <- rnorm(2000)
  crps_hi <- mean(crps_std(rep(0, 2000), rep(1, 2000), dof = 100, y_norm))
  crps_lo <- mean(crps_std(rep(0, 2000), rep(1, 2000), dof = 3,   y_norm))
  expect_gt(crps_hi, crps_lo)
})

test_that("crps_std recycles a scalar dof to the input vector length", {
  skip_if_not_installed("scoringRules")
  set.seed(10)
  mu_s  <- rnorm(100)
  sig_s <- abs(rnorm(100)) + 0.1
  x_s   <- rnorm(100)
  expect_equal(length(crps_std(mu_s, sig_s, dof = 5, x_s)), 100)
})

test_that("crps_std: dof <= 2 raises an error", {
  skip_if_not_installed("scoringRules")
  expect_error(crps_std(0, 1, dof = 2, x = 0))
})

test_that("crps_std: sigma <= 0 raises an error", {
  skip_if_not_installed("scoringRules")
  expect_error(crps_std(0, -1, dof = 5, x = 0))
})

test_that("crps_std: length mismatch raises an error", {
  skip_if_not_installed("scoringRules")
  expect_error(crps_std(c(0, 0), c(1, 1), dof = 5, x = c(0, 0, 0)))
})

# --- 9. Spherical score ---

test_that("spherical_score: perfect forecasts score 1", {
  expect_equal(spherical_score(1, 1), 1)
  expect_equal(spherical_score(0, 0), 1)
})

test_that("spherical_score: worst forecasts score 0", {
  expect_equal(spherical_score(0, 1), 0)
  expect_equal(spherical_score(1, 0), 0)
})

test_that("spherical_score at p=0.5 equals sqrt(0.5)", {
  expect_equal(spherical_score(0.5, 1), sqrt(0.5), tolerance = 1e-10)
})

test_that("spherical_score is positively oriented", {
  expect_gt(spherical_score(0.8, 1), spherical_score(0.2, 1))
})

test_that("spherical_score values are bounded in [0, 1]", {
  set.seed(101)
  p_test <- seq(0, 1, by = 0.01)
  y_test <- rbinom(length(p_test), 1, 0.5)
  scores <- spherical_score(p_test, y_test)
  expect_true(all(scores >= 0 & scores <= 1))
})

# --- 10. Categorical checks ---

test_that("categorical brier_score matches hand-computed values", {
  P3 <- matrix(c(0.7, 0.1, 0.2,
                 0.2, 0.5, 0.3),
               nrow = 2, byrow = TRUE)
  y3 <- c(1L, 2L)
  expect_equal(brier_score(P3, y3), c(-0.07, -0.19), tolerance = 1e-12)
})

test_that("categorical spherical_score matches hand-computed values", {
  P3 <- matrix(c(0.7, 0.1, 0.2,
                 0.2, 0.5, 0.3),
               nrow = 2, byrow = TRUE)
  y3 <- c(1L, 2L)
  expect_equal(spherical_score(P3, y3),
               c(0.7 / sqrt(0.54), 0.5 / sqrt(0.38)), tolerance = 1e-12)
})

test_that("categorical log_score matches hand-computed values", {
  P3 <- matrix(c(0.7, 0.1, 0.2,
                 0.2, 0.5, 0.3),
               nrow = 2, byrow = TRUE)
  y3 <- c(1L, 2L)
  expect_equal(log_score(P3, y3), c(log(0.7), log(0.5)), tolerance = 1e-12)
})

test_that("binary vector path agrees with the K=2 categorical matrix path", {
  p_vec <- c(0.7, 0.2)
  y_01  <- c(1,   0)
  P2    <- cbind(1 - p_vec, p_vec)
  y_12  <- y_01 + 1L

  expect_equal(brier_score(p_vec, y_01),     brier_score(P2, y_12),     tolerance = 1e-12)
  expect_equal(spherical_score(p_vec, y_01), spherical_score(P2, y_12), tolerance = 1e-12)
  expect_equal(log_score(p_vec, y_01),       log_score(P2, y_12),       tolerance = 1e-12)
})

test_that("categorical preprocessing rejects non-integer y", {
  P3 <- matrix(c(0.7, 0.1, 0.2, 0.2, 0.5, 0.3), nrow = 2, byrow = TRUE)
  expect_error(brier_score(P3, c(1, 1.9)))
})

test_that("categorical preprocessing rejects NA y", {
  P3 <- matrix(c(0.7, 0.1, 0.2, 0.2, 0.5, 0.3), nrow = 2, byrow = TRUE)
  expect_error(brier_score(P3, c(1L, NA_integer_)))
})

test_that("categorical preprocessing rejects rows not summing to 1", {
  P3 <- matrix(c(0.7, 0.1, 0.2, 0.2, 0.5, 0.3), nrow = 2, byrow = TRUE)
  y3 <- c(1L, 2L)
  P_bad <- P3
  P_bad[1, 1] <- 0.8
  expect_error(brier_score(P_bad, y3))
})

test_that("categorical preprocessing rejects negative probabilities", {
  P3 <- matrix(c(0.7, 0.1, 0.2, 0.2, 0.5, 0.3), nrow = 2, byrow = TRUE)
  y3 <- c(1L, 2L)
  P_bad <- P3
  P_bad[1, 2] <- -0.01
  expect_error(brier_score(P_bad, y3))
})

test_that("categorical preprocessing rejects zero class labels", {
  P3 <- matrix(c(0.7, 0.1, 0.2, 0.2, 0.5, 0.3), nrow = 2, byrow = TRUE)
  expect_error(brier_score(P3, c(0L, 2L)))
})

test_that("categorical preprocessing rejects K=1 matrices", {
  P_bad <- matrix(1, nrow = 2, ncol = 1)
  expect_error(brier_score(P_bad, c(1L, 1L)))
})

# --- 11. score_bounds ---

test_that("score_bounds: Brier has exact bounds and correct constants", {
  sb_brier <- score_bounds("brier")
  expect_equal(sb_brier$lo, -1)
  expect_equal(sb_brier$hi,  1)
  expect_equal(sb_brier$c_thm1,  1)
  expect_equal(sb_brier$c_thm23, 2)
})

test_that("score_bounds: Brier constants are arithmetically consistent with lo/hi", {
  sb_brier <- score_bounds("brier")
  expect_equal(sb_brier$c_thm1,  (sb_brier$hi - sb_brier$lo) / 2, tolerance = 1e-12)
  expect_equal(sb_brier$c_thm23, (sb_brier$hi - sb_brier$lo),     tolerance = 1e-12)
})

test_that("score_bounds: spherical has the same score-difference interval as Brier", {
  sb_spherical <- score_bounds("spherical")
  expect_equal(sb_spherical$lo, -1)
  expect_equal(sb_spherical$hi,  1)
  expect_equal(sb_spherical$c_thm1,  1)
  expect_equal(sb_spherical$c_thm23, 2)
})

test_that("score_bounds: Winkler is a special case with lo = -Inf and c_thm1 = NA", {
  sb_wink <- score_bounds("winkler")
  expect_true(is.infinite(sb_wink$lo) && sb_wink$lo < 0)
  expect_equal(sb_wink$hi, 1)
  expect_true(is.na(sb_wink$c_thm1))
  expect_equal(sb_wink$c_thm23, 2)
})

unbounded_rules <- c("tick", "crps", "crps_normal", "crps_empirical", "crps_std",
                     "log", "qlike")
for (rule in unbounded_rules) {
  test_that(sprintf("score_bounds('%s') returns NULL with an 'unbounded' message", rule), {
    expect_message(out <- score_bounds(rule), "unbounded")
    expect_null(out)
  })
}

test_that("score_bounds raises an error for an unknown rule", {
  expect_error(score_bounds("rmse"))
})
