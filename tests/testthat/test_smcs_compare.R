# =============================================================================
# tests/testthat/test_multi_compare.R
# Tests for V2 SMCS parameter slicing, array builders, and the high-level wrapper
# =============================================================================

set.seed(42)
T_sim <- 50
m_models <- 3

outcomes <- rbinom(T_sim, 1, 0.5)
forecasts <- matrix(runif(T_sim * m_models), nrow = T_sim, ncol = m_models)
colnames(forecasts) <- paste0("M", 1:m_models)

# Create a dummy scores matrix to use for SMCS tests
scores_mat <- matrix(0, nrow = T_sim, ncol = m_models)
for(i in 1:m_models) scores_mat[, i] <- brier_score(forecasts[, i], outcomes)

# ── A. INTERNAL PARAMETER SLICING ────────────────────────────────────────────

test_that(".slice_param routes parameters correctly", {
  # 1. Scalar
  expect_equal(seqcomp:::.slice_param(2, T_sim, 1, 2), 2)

  # 2. Matrix (m x m)
  mat <- matrix(1:9, nrow = 3)
  expect_equal(seqcomp:::.slice_param(mat, T_sim, 1, 2), mat[1, 2])

  # 3. 3D Array (T x m x m)
  arr <- array(runif(T_sim * 3 * 3), dim = c(T_sim, 3, 3))
  expect_equal(seqcomp:::.slice_param(arr, T_sim, 1, 2), arr[, 1, 2])
  expect_length(seqcomp:::.slice_param(arr, T_sim, 1, 2), T_sim)

  # 4. NULL
  expect_null(seqcomp:::.slice_param(NULL, T_sim, 1, 2))
})

# ── B. DEFENSIVE ERROR HANDLERS ──────────────────────────────────────────────

test_that("smcs_strong and smcs_weak strictly enforce parameter shapes", {
  vec_c <- rep(2, T_sim)
  arr_c <- array(2, dim = c(T_sim, m_models, m_models))

  # smcs_weak should reject time-varying bounds
  expect_error(
    smcs_weak(scores_mat, c_param = vec_c),
    "requires c_param to be a scalar or an m x m matrix"
  )
  expect_error(
    smcs_weak(scores_mat, c_param = arr_c),
    "requires c_param to be a scalar or an m x m matrix"
  )

  # smcs_strong(method="mixture") should reject time-varying bounds
  expect_error(
    smcs_strong(scores_mat, method = "mixture", c_param = vec_c),
    "requires c_param to be a scalar or an m x m matrix"
  )
  expect_error(
    smcs_strong(scores_mat, method = "mixture", c_param = arr_c),
    "requires c_param to be a scalar or an m x m matrix"
  )

  # smcs_strong(method="betting") SHOULD accept time-varying bounds
  expect_error(smcs_strong(scores_mat, method = "betting", c_param = arr_c), NA)
})

# ── C. 3D ARRAY BUILDER ──────────────────────────────────────────────────────

test_that("build_quantile_betting_arrays outputs valid 3D arrays", {
  # Mock some tick loss scores (unbounded)
  tick_scores <- matrix(0, nrow = T_sim, ncol = m_models)
  for(i in 1:m_models) tick_scores[, i] <- tick_loss(forecasts[, i], outcomes, tau = 0.5)

  eps <- 1e-8
  bnds <- build_quantile_betting_arrays(
    forecasts, tick_scores, tau = 0.5, eps = eps
  )

  expect_true(is.list(bnds))
  expect_true(all(c("c_array", "lambda_array") %in% names(bnds)))
  expect_equal(dim(bnds$c_array), c(T_sim, m_models, m_models))
  expect_equal(dim(bnds$lambda_array), c(T_sim, m_models, m_models))

  # Check the returned arrays include the strictly positive bound floor.
  manual <- lambda_betting_quantile(
    p_t = forecasts[, 1],
    q_t = forecasts[, 2],
    tau = 0.5,
    delta_hat_lag1 = c(0, head(tick_scores[, 1] - tick_scores[, 2], -1)),
    eps = eps
  )
  c_safe <- pmax(manual$c_t, eps)
  lambda_safe <- pmin(manual$lambda_t, 1 / c_safe)

  expect_equal(bnds$c_array[, 1, 2], c_safe)
  expect_equal(bnds$lambda_array[, 1, 2], lambda_safe)
})

test_that("build_quantile_betting_arrays floors bounds for tied forecasts", {
  tied_forecasts <- forecasts
  tied_forecasts[, 2] <- tied_forecasts[, 1]

  tied_scores <- matrix(0, nrow = T_sim, ncol = m_models)
  for (i in seq_len(m_models)) {
    tied_scores[, i] <- tick_loss(tied_forecasts[, i], outcomes, tau = 0.5)
  }

  eps <- 1e-8
  bnds <- build_quantile_betting_arrays(
    tied_forecasts, tied_scores, tau = 0.5, eps = eps
  )

  expect_equal(bnds$c_array[, 1, 2], rep(eps, T_sim))
  expect_equal(bnds$c_array[, 2, 1], rep(eps, T_sim))
  expect_true(all(bnds$lambda_array[, 1, 2] <= 1 / bnds$c_array[, 1, 2]))
  expect_true(all(bnds$lambda_array[, 2, 1] <= 1 / bnds$c_array[, 2, 1]))
})

test_that("quantile betting eps must be a finite positive scalar", {
  expect_error(
    lambda_betting_quantile(
      p_t = c(0, 1),
      q_t = c(1, 0),
      tau = 0.5,
      delta_hat_lag1 = c(0, 0),
      eps = Inf
    ),
    "eps must be a finite positive scalar"
  )

  expect_error(
    build_quantile_betting_arrays(
      forecasts, scores_mat, tau = 0.5, eps = NA_real_
    ),
    "eps must be a finite positive scalar"
  )
})

# ── D. HIGH-LEVEL MULTI-MODEL WRAPPER ────────────────────────────────────────

test_that("smcs_compare handles bounded scoring rules (Brier)", {
  # We pass v_opt = 10 via '...' to prove it gets routed to the martingales properly
  # We also test the new cs_method argument
  res <- smcs_compare(
    forecasts, outcomes,
    scoring_rule = "brier",
    cs_method = "hoeffding",
    alpha = 0.05,
    v_opt = 10
  )

  expect_true(is.list(res))
  expect_true(all(c(
    "scores", "smcs_strong", "smcs_uniform_weak", "smcs_weak",
    "betting_rule", "period", "betting_args"
  ) %in% names(res)))
  expect_identical(res$betting_rule, "naive")
  expect_identical(res$period, 1L)
  expect_identical(res$betting_args, list())
  expect_equal(dim(res$scores), c(T_sim, m_models))
  expect_equal(dim(res$smcs_strong), c(T_sim, m_models))
  expect_equal(dim(res$smcs_uniform_weak), c(T_sim, m_models))
  expect_equal(dim(res$smcs_weak), c(T_sim, m_models))

  expect_true(is.logical(res$smcs_strong))
  expect_true(is.logical(res$smcs_uniform_weak))

  direct_uniform_weak <- smcs_strong(
    res$scores,
    alpha = 0.05,
    method = "mixture",
    c_param = 2,
    v_opt = 10,
    clip_max = 1e7
  )
  expect_identical(res$smcs_uniform_weak, direct_uniform_weak$smcs)
  expect_true(is.logical(res$smcs_weak))
})

test_that("smcs_compare handles conditionally bounded rules (Tick)", {
  # Tick loss triggers the 3D array builder and omits both weak-null SMCSs.
  expect_warning(
    res <- smcs_compare(forecasts, outcomes, scoring_rule = "tick", tau = 0.5),
    "smcs_uniform_weak and smcs_weak are currently omitted"
  )

  expect_equal(dim(res$scores), c(T_sim, m_models))
  expect_equal(dim(res$smcs_strong), c(T_sim, m_models))
  expect_null(res$smcs_weak)
  expect_null(res$smcs_uniform_weak)
  expect_identical(res$betting_rule, "arnold")
  expect_identical(res$period, 1L)
  expect_identical(res$betting_args, list())
})

test_that("smcs_compare routes adaptive betting rules", {
  res_agrapa <- smcs_compare(
    forecasts, outcomes,
    scoring_rule = "brier",
    betting_rule = "agrapa",
    period = 1,
    kappa = 0.25
  )

  expect_identical(res_agrapa$betting_rule, "agrapa")
  expect_identical(res_agrapa$period, 1L)
  expect_identical(res_agrapa$betting_args, list(kappa = 0.25))
  expect_equal(dim(res_agrapa$smcs_strong), c(T_sim, m_models))

  expect_warning(
    res_ons <- smcs_compare(
      forecasts, outcomes,
      scoring_rule = "tick",
      tau = 0.5,
      betting_rule = "ons",
      period = 5,
      eta = 1
    ),
    "smcs_uniform_weak and smcs_weak are currently omitted"
  )

  expect_identical(res_ons$betting_rule, "ons")
  expect_identical(res_ons$period, 5L)
  expect_identical(res_ons$betting_args, list(eta = 1))
  expect_equal(dim(res_ons$smcs_strong), c(T_sim, m_models))
  expect_null(res_ons$smcs_uniform_weak)
  expect_null(res_ons$smcs_weak)
})

test_that("smcs_compare rejects unsupported betting-rule combinations", {
  expect_error(
    smcs_compare(
      forecasts, outcomes,
      scoring_rule = "brier",
      betting_rule = "arnold"
    ),
    "available only when scoring_rule = 'tick'"
  )

  expect_error(
    smcs_compare(
      forecasts, outcomes,
      scoring_rule = "brier",
      betting_rule = "naive",
      period = 2
    ),
    "period is available only"
  )

  expect_error(
    smcs_compare(
      forecasts, outcomes,
      scoring_rule = "tick",
      tau = 0.5,
      betting_rule = "arnold",
      eta = 1
    ),
    "Additional betting arguments are available only"
  )

  expect_error(
    smcs_compare(
      forecasts, outcomes,
      scoring_rule = "brier",
      betting_rule = "agrpa"
    ),
    "should be one of"
  )
})

test_that("smcs_compare enforces inputs", {
  # Missing tau for tick loss
  expect_error(
    smcs_compare(forecasts, outcomes, scoring_rule = "tick"),
    "tau must be provided for tick_loss"
  )

  # Dimension mismatch between outcomes and forecasts
  expect_error(
    smcs_compare(forecasts, outcomes[1:10], scoring_rule = "brier"),
    "length\\(outcomes\\) == Tt is not TRUE"
  )
})

# ── E. S3 METHODS ────────────────────────────────────────────────────────────

test_that("S3 print and summary methods work for seqcomp_multi", {
  res <- smcs_compare(
    forecasts, outcomes, scoring_rule = "brier", alpha = 0.05
  )

  # 1. Test print (use capture.output so it doesn't spam the test logs)
  expect_error(out_str <- capture.output(print(res)), NA)
  expect_true(any(grepl("seqcomp multi-model comparison", out_str, fixed = TRUE)))
  expect_true(any(grepl("strong betting rule = 'naive'", out_str, fixed = TRUE)))
  expect_true(any(grepl("Strong-null SMCS", out_str, fixed = TRUE)))
  expect_true(any(grepl("Uniformly-weak SMCS", out_str, fixed = TRUE)))

  # 2. Test summary
  expect_error(summ <- summary(res, checkpoints = c(0, 0.5, 1)), NA)
  expect_true(is.data.frame(summ))
  expect_equal(nrow(summ), 3)
  expect_true(all(c(
    "t", "strong_set_size", "uniform_weak_set_size", "weak_set_size"
  ) %in% names(summ)))
})
