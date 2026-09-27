# =============================================================================
# tests/testthat/test_smcs.R
# Sequential Model Confidence Sets (Arnold et al., 2026)
# =============================================================================

# ── A. STRUCTURAL TESTS FOR vovk_wang_merge ──────────────────────────────────

test_that("vovk_wang_merge works for trivial cases", {
  expect_equal(vovk_wang_merge(5), 5)
  expect_error(vovk_wang_merge(c(1, -2, 3)), "nonnegative")
})

test_that("vovk_wang_merge accurately computes closed-testing minimums", {
  # Hand-calculated example: input c(10, 2, 1)
  # Subset containing 1: min is {1} -> 1
  # Subset containing 2: min of {2}, {1,2} -> (2+1)/2 = 1.5
  # Subset containing 10: min of {10}, {10,2}, {10,2,1} -> (10+2+1)/3 = 13/3 = 4.3333...

  raw_e <- c(10, 2, 1)
  adj_e <- vovk_wang_merge(raw_e)

  expect_equal(adj_e[1], 13/3, tolerance = 1e-6)
  expect_equal(adj_e[2], 1.5, tolerance = 1e-6)
  expect_equal(adj_e[3], 1.0, tolerance = 1e-6)
})

test_that("vovk_wang_merge output is always <= input (conservative adjustment)", {
  set.seed(99)
  raw_e <- rexp(50, rate = 0.1)
  adj_e <- vovk_wang_merge(raw_e)

  expect_true(all(adj_e <= raw_e + 1e-12))
})


# ── B. STRUCTURAL TESTS FOR smcs_strong & smcs_weak ──────────────────────────

set.seed(42)
T_sim <- 150
m_models <- 4
# Matrix of Brier scores (bounded in [-1, 0] so score diffs in [-1, 1], c=2)
scores_mat <- matrix(runif(T_sim * m_models, min = -0.5, max = 0), nrow = T_sim, ncol = m_models)
colnames(scores_mat) <- paste0("Model_", 1:m_models)

test_that("smcs_strong returns correct structure with mixture method", {
  res <- smcs_strong(scores_mat, alpha = 0.05, method = "mixture", c_param = 2, v_opt = 10)

  expect_true(is.list(res))
  expect_true(all(c("E_i_dot", "E_star", "smcs") %in% names(res)))
  expect_equal(dim(res$E_star), c(T_sim, m_models))
  expect_equal(dim(res$smcs), c(T_sim, m_models))
  expect_true(is.logical(res$smcs))
})

test_that("smcs_strong returns correct structure with betting method", {
  res <- smcs_strong(scores_mat, alpha = 0.05, method = "betting", c_param = 2)

  expect_equal(dim(res$E_star), c(T_sim, m_models))
  expect_true(is.logical(res$smcs))
})

test_that("smcs_weak returns correct structure", {
  res <- smcs_weak(scores_mat, alpha = 0.05, cs_method = "bernstein", c_param = 2, v_opt = 10)

  expect_true(is.list(res))
  expect_true(all(c("smcs", "alpha_adjusted") %in% names(res)))
  expect_equal(dim(res$smcs), c(T_sim, m_models))
  expect_equal(res$alpha_adjusted, 0.05 / (m_models * (m_models - 1)))
})

test_that("SMCS maintains running intersection for strong null, but not necessarily weak", {
  res_strong <- smcs_strong(scores_mat, alpha = 0.05, method = "betting", c_param = 2)
  res_weak   <- smcs_weak(scores_mat, alpha = 0.05, cs_method = "hoeffding", c_param = 1)

  # For any column, there should be no transitions from FALSE back to TRUE
  check_monotonicity <- function(mat) {
    all(apply(mat, 2, function(col) {
      false_idx <- which(!col)
      if (length(false_idx) == 0) return(TRUE)
      all(!col[false_idx[1]:length(col)])
    }))
  }

  # Strong null MUST be monotonic (permanent exclusion)
  expect_true(check_monotonicity(res_strong$smcs))

  # Weak null is point-in-time, so we do NOT enforce monotonicity on res_weak$smcs.
  # We just check that it cleanly outputs a logical matrix without crashing.
  expect_true(is.logical(res_weak$smcs))
})

# ── C. STATISTICAL PROPERTIES (GATED SIMULATIONS) ────────────────────────────

test_that("SMCS reliably eliminates a demonstrably inferior model", {
  # Model 4 is strictly worse than Models 1-3
  scores_alt <- scores_mat
  scores_alt[, 4] <- scores_alt[, 4] - 0.4 # Worse scores

  # Strong
  res_s <- smcs_strong(scores_alt, alpha = 0.05, method = "betting", c_param = 2)
  expect_false(res_s$smcs[T_sim, 4])          # Model 4 should be eliminated
  expect_true(all(res_s$smcs[T_sim, 1:3]))    # Others should survive

  # Weak
  res_w <- smcs_weak(scores_alt, alpha = 0.05, cs_method = "bernstein", c_param = 2)
  expect_false(res_w$smcs[T_sim, 4])
})

test_that("smcs_strong integrates perfectly with aGRAPA and ONS-m 3D arrays", {
  # Model 4 is strictly worse
  scores_alt <- scores_mat
  scores_alt[, 4] <- scores_alt[, 4] - 0.4

  # Build arrays (c = 2 for Brier score differences)
  lam_agrapa <- build_agrapa_betting_array(scores_alt, c_mat = 2)
  lam_ons    <- build_ons_betting_array(scores_alt, c_mat = 2)

  res_agrapa <- smcs_strong(scores_alt, alpha = 0.05, method = "betting",
                            c_param = 2, lambda_param = lam_agrapa)
  res_ons    <- smcs_strong(scores_alt, alpha = 0.05, method = "betting",
                            c_param = 2, lambda_param = lam_ons)

  # Check structure
  expect_equal(dim(res_agrapa$smcs), c(T_sim, m_models))
  expect_equal(dim(res_ons$smcs), c(T_sim, m_models))

  # Strong null MUST be monotonic (permanent exclusion)
  check_monotonicity <- function(mat) {
    all(apply(mat, 2, function(col) {
      false_idx <- which(!col)
      if (length(false_idx) == 0) return(TRUE)
      all(!col[false_idx[1]:length(col)])
    }))
  }
  expect_true(check_monotonicity(res_agrapa$smcs))
  expect_true(check_monotonicity(res_ons$smcs))

  # Model 4 should be eliminated
  expect_false(res_agrapa$smcs[T_sim, 4])
  expect_false(res_ons$smcs[T_sim, 4])
})

test_that("build_agrapa_betting_array routes the correct sign (prevents silent underpowering)", {
  # Trivial case: T=2, m=2
  # Model 2 always wins by the max margin (score diff = 1)
  scores <- matrix(c(0, 0,   # Model 1
                     1, 1),  # Model 2
                   nrow = 2, ncol = 2)
  c_mat <- 2

  # Compute the betting array
  lam_array <- build_agrapa_betting_array(scores, c_mat = c_mat)

  # We want lambda_array[, 1, 2], which is the bet for Model 2 beating Model 1
  lam_12 <- lam_array[, 1, 2]

  # Expected from hand-trace:
  # t=1: Y=(1,1). mu_lag=0.5. raw=0 -> lam_d=0
  # t=2: Y=(1,1). mu_lag=0.75, sig2_lag=0.15625. raw=0.25/0.21875 = 1.14 -> clipped to 1 -> lam_d=0.5
  expect_equal(lam_12[1], 0)
  expect_equal(lam_12[2], 0.5)

  # Verify downstream execution in smcs_strong
  # eprocess_betting evaluates: 1 - lam * (S_1 - S_2)
  # S_1 - S_2 = -1 everywhere.
  # t=1: 1 - 0 * (-1) = 1.0
  # t=2: 1 - 0.5 * (-1) = 1.5
  res <- smcs_strong(scores, alpha = 0.05, method = "betting", c_param = c_mat, lambda_param = lam_array)

  # For m=2, E_i_dot[, 1] is exactly the evidence that Model 2 beats Model 1
  expect_equal(res$E_i_dot[1, 1], 1.0)
  expect_equal(res$E_i_dot[2, 1], 1.5)
})

test_that("smcs_weak FWER is controlled under the global null", {
  skip_on_cran()
  skip_if(Sys.getenv("SEQCOMP_RUN_SLOW_TESTS") != "true",
          "slow FWER simulation; set SEQCOMP_RUN_SLOW_TESTS=true to run")

  N_sim <- 500
  T_fw <- 200
  m_fw <- 3
  alpha_sim <- 0.05
  n_fwer_violations <- 0L

  set.seed(2026)
  for (i in seq_len(N_sim)) {
    # Global null: all models output identically distributed noise
    null_scores <- matrix(rnorm(T_fw * m_fw, mean = 0, sd = 0.1), nrow = T_fw)
    # Clip to bounds [-0.5, 0.5] -> score differences in [-1, 1], c = 2
    null_scores <- pmax(pmin(null_scores, 0.5), -0.5)

    res <- smcs_weak(null_scores, alpha = alpha_sim, cs_method = "bernstein", c_param = 2)

    # A Type I error occurs if ANY model is falsely excluded under the global null
    if (any(!res$smcs[T_fw, ])) {
      n_fwer_violations <- n_fwer_violations + 1L
    }
  }

  fwer <- n_fwer_violations / N_sim
  expect_lte(fwer, alpha_sim)
})
