# =============================================================================
# tests/testthat/test_internal_eprocess_helpers.R
# =============================================================================

test_that("log_ge_mixture_from_sv matches direct GE mixture loop", {
  set.seed(123)
  xs <- runif(100, -0.5, 0.5)
  gammas <- make_gammas(xs)
  V_t <- intrinsic_time(xs, gammas, floor = FALSE)
  S_t <- cumsum(xs)

  rho <- rho_from_vopt(v_opt = 10, alpha = 0.025)
  c <- 2

  direct <- mapply(
    function(s, v) {
      val <- ge_mixture(s, v, rho, c)
      if (val <= 0 || is.nan(val)) return(log(.Machine$double.eps))
      log(val)
    },
    s = S_t,
    v = V_t
  )

  helper <- log_ge_mixture_from_sv(S_t, V_t, rho, c)

  expect_equal(helper, direct, tolerance = 1e-12)
})

test_that("log_eprocess_fixed_predictable matches direct predictable loop", {
  set.seed(124)
  xs <- runif(100, -0.4, 0.4)
  c_seq <- rep(2, length(xs))
  lambda <- 0.2
  gammas <- make_gammas(xs)

  direct_increments <- numeric(length(xs))
  for (i in seq_along(xs)) {
    direct_increments[i] <-
      lambda * xs[i] - psi_e(lambda, c_seq[i]) * (xs[i] - gammas[i])^2
  }

  direct <- cumsum(direct_increments)
  helper <- log_eprocess_fixed_predictable(xs, c_seq, lambda, gammas)

  expect_equal(helper, direct, tolerance = 1e-12)
})
