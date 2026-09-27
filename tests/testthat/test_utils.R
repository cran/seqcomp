# =============================================================================
# tests/testthat/test_utils.R
# Numerical verification for all primitives in R/utils.R
# =============================================================================

# ------------------------------------------------------------
# 1. rho_from_vopt (Lambert W formula)
# ------------------------------------------------------------

test_that("rho_from_vopt returns a positive, plausible value", {
  rho1 <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  expect_gt(rho1, 0)
  # known numerical range: -W_{-1}(-0.025^2/e) - 1 ~= 7.3, so rho ~= 1.37
  expect_gt(rho1, 0.5)
  expect_lt(rho1, 3.0)
})

test_that("rho_from_vopt increases with v_opt", {
  rho1      <- rho_from_vopt(v_opt = 10,  alpha = 0.025)
  rho_small <- rho_from_vopt(v_opt = 1,   alpha = 0.025)
  rho_large <- rho_from_vopt(v_opt = 100, alpha = 0.025)
  expect_lt(rho_small, rho1)
  expect_lt(rho1, rho_large)
})

test_that("rho_from_vopt changes with alpha", {
  rho_strict <- rho_from_vopt(v_opt = 10, alpha = 0.005)
  rho_loose  <- rho_from_vopt(v_opt = 10, alpha = 0.100)
  expect_false(isTRUE(all.equal(rho_strict, rho_loose)))
})

# ------------------------------------------------------------
# 2. rho formula cross-check (Lambert W vs Python approximation)
# ------------------------------------------------------------

test_that("Lambert W rho is always smaller than the Python approximation", {
  capture.output({
    result <- seqcomp:::verify_rho_formulas(v_opt = 10, tol = 0.02)
  })

  expect_true(all(result$rho_lambertW < result$rho_python))
})

# ------------------------------------------------------------
# 3. cm_boundary (Normal Mixture)
# ------------------------------------------------------------

test_that("cm_boundary returns positive, increasing values", {
  rho_cm <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  v_seq  <- c(1, 10, 100, 500)
  u_cm   <- cm_boundary(v = v_seq, alpha = 0.025, rho = rho_cm)

  expect_true(all(u_cm > 0))
  expect_true(all(diff(u_cm) > 0))
})

test_that("cm_boundary matches the closed-form formula at v=10", {
  rho_cm <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  u_cm10 <- cm_boundary(v = 10, alpha = 0.025, rho = rho_cm)
  manual <- sqrt((10 + rho_cm) * log((10 + rho_cm) / (0.025^2 * rho_cm)))

  expect_equal(u_cm10, manual, tolerance = 1e-10)
})

test_that("cm_boundary is finite at v=0", {
  rho_cm <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  u_cm0 <- cm_boundary(v = 0, alpha = 0.025, rho = rho_cm)
  expect_true(is.finite(u_cm0))
})

test_that("CM radius (u/t) shrinks as t grows", {
  rho_cm <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  t_seq  <- c(1, 10, 100, 500)
  radii  <- cm_boundary(v = t_seq, alpha = 0.025, rho = rho_cm) / t_seq
  expect_true(all(diff(radii) < 0))
})

# ------------------------------------------------------------
# 4. ge_mixture (GE mixture function m(s, v)) -- internal
# ------------------------------------------------------------

test_that("ge_mixture(0, 0) equals 1 (e-process starts at 1)", {
  rho_ge <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  m00 <- seqcomp:::ge_mixture(s = 0, v = 0, rho = rho_ge, c = 2)
  expect_equal(m00, 1, tolerance = 1e-6)
})

test_that("ge_mixture is positive and finite for positive s", {
  rho_ge <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  m_pos <- seqcomp:::ge_mixture(s = 5, v = 10, rho = rho_ge, c = 2)
  expect_gt(m_pos, 0)
  expect_true(is.finite(m_pos))
})

test_that("ge_mixture increases with s for fixed v", {
  rho_ge <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  m_s1 <- seqcomp:::ge_mixture(s = 1,  v = 10, rho = rho_ge, c = 2)
  m_s2 <- seqcomp:::ge_mixture(s = 5,  v = 10, rho = rho_ge, c = 2)
  m_s3 <- seqcomp:::ge_mixture(s = 10, v = 10, rho = rho_ge, c = 2)
  expect_lt(m_s1, m_s2)
  expect_lt(m_s2, m_s3)
})

test_that("ge_mixture(0, large v) stays finite and positive under the null", {
  rho_ge <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  m_null_large_v <- seqcomp:::ge_mixture(s = 0, v = 100, rho = rho_ge, c = 2)
  expect_true(is.finite(m_null_large_v))
  expect_gt(m_null_large_v, 0)
})

test_that("ge_mixture returns 0 for very negative s (x<=0 guard)", {
  rho_ge <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  m_neg <- seqcomp:::ge_mixture(s = -1000, v = 1, rho = rho_ge, c = 2)
  expect_lt(m_neg, 1e-12)
  expect_gt(m_neg, 0)
})

# ------------------------------------------------------------
# 5. ge_boundary (GE boundary via uniroot)
# ------------------------------------------------------------

test_that("ge_boundary is finite and increasing in v", {
  rho_ge2 <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  v_test  <- c(1, 5, 10, 50, 100)
  u_ge    <- ge_boundary(v = v_test, alpha = 0.025, rho = rho_ge2, c = 2)

  expect_true(all(is.finite(u_ge)))
  expect_true(all(diff(u_ge) > 0))
})

test_that("ge_boundary inverts ge_mixture to hit 1/alpha", {
  rho_ge2 <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  v_test  <- c(1, 5, 10, 50, 100)
  u_ge    <- ge_boundary(v = v_test, alpha = 0.025, rho = rho_ge2, c = 2)

  for (i in seq_along(v_test)) {
    m_check <- seqcomp:::ge_mixture(u_ge[i], v_test[i], rho = rho_ge2, c = 2)
    expect_equal(m_check, 40, tolerance = 0.1 / 40)
  }
})

test_that("ge_boundary avoids overflow warnings on small-c inputs", {
  alpha <- 0.025
  rho <- rho_from_vopt(v_opt = 1, alpha = alpha)
  v_test <- c(0, 1, 10, 1000)

  u_ge <- NULL

  expect_no_warning({
    u_ge <- ge_boundary(
      v = v_test,
      alpha = alpha,
      rho = rho,
      c = 0.01
    )
  })

  expect_true(all(is.finite(u_ge)))

  m_at_boundary <- vapply(
    seq_along(v_test),
    function(i) {
      seqcomp:::ge_mixture(
        s = u_ge[[i]],
        v = v_test[[i]],
        rho = rho,
        c = 0.01
      )
    },
    numeric(1)
  )

  expect_equal(
    m_at_boundary,
    rep(1 / alpha, length(v_test)),
    tolerance = 1e-5
  )
})

test_that("GE radius (u/t) shrinks as t grows", {
  rho_ge2  <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  t_test   <- c(1, 5, 10, 50, 100)
  u_ge2    <- ge_boundary(v = t_test, alpha = 0.025, rho = rho_ge2, c = 2)
  radii_ge <- u_ge2 / t_test
  expect_true(all(diff(radii_ge) < 0))
})

# ------------------------------------------------------------
# 6. ps_boundary (Polynomial Stitched)
# ------------------------------------------------------------

test_that("ps_boundary returns finite, positive, non-decreasing values", {
  u_ps <- ps_boundary(v = c(1, 10, 100, 500), alpha = 0.025, v_opt = 10,
                       c = 2, s = 1.4, eta = 2)

  expect_true(all(is.finite(u_ps) & u_ps > 0))
  # PS boundary is non-decreasing, not strictly increasing
  expect_true(all(diff(u_ps) >= 0))
})

test_that("hardcoded CR24 EC.27 formula matches the expected value at v=10", {
  hc_10 <- seqcomp:::cs_boundary_cr23_hardcoded(10)
  expect_lt(abs(hc_10 - 43.057), 0.01)
})

# ------------------------------------------------------------
# 7. intrinsic_time and make_gammas -- internal
# ------------------------------------------------------------

test_that("make_gammas produces the correctly lagged running mean", {
  set.seed(42)
  xs_test <- rnorm(100, mean = 0.1, sd = 0.5)
  gam <- seqcomp:::make_gammas(xs_test, lag = 1)

  expect_equal(length(gam), length(xs_test))
  expect_equal(gam[1], 0)
  expect_equal(gam[2], xs_test[1], tolerance = 1e-12)

  manual_gam <- c(0, cumsum(xs_test[-length(xs_test)]) /
                       seq_along(xs_test[-length(xs_test)]))
  expect_equal(gam, manual_gam, tolerance = 1e-12)
})

test_that("intrinsic_time behaves correctly with floor = TRUE", {
  set.seed(42)
  xs_test <- rnorm(100, mean = 0.1, sd = 0.5)
  vt_floor <- seqcomp:::intrinsic_time(xs_test, floor = TRUE)

  expect_equal(length(vt_floor), length(xs_test))
  expect_true(all(vt_floor >= 1))
  expect_true(all(diff(vt_floor) >= 0))
})

test_that("intrinsic_time runs without error when floor = FALSE", {
  set.seed(42)
  xs_test <- rnorm(100, mean = 0.1, sd = 0.5)
  expect_error(seqcomp:::intrinsic_time(xs_test, floor = FALSE), NA)
})

# ------------------------------------------------------------
# 8. psi_e (sub-exponential CGF) -- internal
# ------------------------------------------------------------

test_that("psi_e is positive for positive lambda and zero at lambda=0", {
  psi_val <- seqcomp:::psi_e(lambda = 0.1, c = 2)
  expect_gt(psi_val, 0)

  psi_zero <- seqcomp:::psi_e(lambda = 0, c = 2)
  expect_equal(psi_zero, 0, tolerance = 1e-12)
})

test_that("psi_e matches the manual formula at lambda=0.1, c=2", {
  psi_val <- seqcomp:::psi_e(lambda = 0.1, c = 2)
  manual_psi <- (-log(1 - 2 * 0.1) - 2 * 0.1) / 4
  expect_equal(psi_val, manual_psi, tolerance = 1e-12)
})

test_that("psi_e is increasing in lambda", {
  psi_seq <- vapply(
    seq(0.01, 0.4, by = 0.01),
    seqcomp:::psi_e,
    c = 2,
    FUN.VALUE = numeric(1)
  )
  expect_true(all(diff(psi_seq) > 0))
})

# ------------------------------------------------------------
# 9. log_eprocess_mixture and clip_eprocess -- internal
# ------------------------------------------------------------

test_that("log_eprocess_mixture is well-behaved under the null", {
  set.seed(99)
  xs_null <- rnorm(200, mean = 0, sd = 0.5)
  rho_ep <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  log_e_null <- seqcomp:::log_eprocess_mixture(xs_null, rho = rho_ep, c = 2)

  expect_equal(length(log_e_null), length(xs_null))
  expect_false(any(is.nan(log_e_null)))
})

test_that("log_eprocess_mixture grows more under the alternative than under the null", {
  set.seed(99)
  xs_null <- rnorm(200, mean = 0,   sd = 0.5)
  xs_alt  <- rnorm(200, mean = 0.3, sd = 0.5)
  rho_ep <- rho_from_vopt(v_opt = 10, alpha = 0.025)

  log_e_null <- seqcomp:::log_eprocess_mixture(xs_null, rho = rho_ep, c = 2)
  log_e_alt  <- seqcomp:::log_eprocess_mixture(xs_alt,  rho = rho_ep, c = 2)

  expect_gt(mean(log_e_alt), mean(log_e_null))
})

test_that("clip_eprocess caps values from above and keeps them positive", {
  set.seed(99)
  xs_null <- rnorm(200, mean = 0,   sd = 0.5)
  xs_alt  <- rnorm(200, mean = 0.3, sd = 0.5)
  rho_ep <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  log_e_alt <- seqcomp:::log_eprocess_mixture(xs_alt, rho = rho_ep, c = 2)

  e_clipped <- seqcomp:::clip_eprocess(log_e_alt, clip_max = 1e7)
  expect_true(all(e_clipped <= 1e7))
  expect_true(all(e_clipped > 0))
})

# ------------------------------------------------------------
# 10. score_diff_scales
# ------------------------------------------------------------

test_that("score_diff_scales computes correct constants for Brier score diffs", {
  # Brier score: bounds [-1,0], so score diff is in [-1,1]
  scales_brier <- score_diff_scales(lo = -1, hi = 1)
  expect_equal(scales_brier$c_thm1, 1)
  expect_equal(scales_brier$c_thm23, 2)
})

test_that("score_diff_scales computes correct constants for tick loss diffs", {
  # Tick loss at quantile 0.99: diff in [-0.99, 0.99]
  scales_tick <- score_diff_scales(lo = -0.99, hi = 0.99)
  expect_equal(scales_tick$c_thm23, 1.98, tolerance = 1e-10)
})
