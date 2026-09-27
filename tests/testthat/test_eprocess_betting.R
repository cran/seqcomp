# =============================================================================
# tests/testthat/test_eprocess_betting.R
# Numerical verification for betting-style e-processes and adaptive lambdas
# =============================================================================

set.seed(42)
T_ <- 200
p_forecasts <- runif(T_)
q_forecasts <- runif(T_)
y_outcomes  <- rbinom(T_, 1, 0.5)
s1 <- brier_score(p_forecasts, y_outcomes)
s2 <- brier_score(q_forecasts, y_outcomes)

# ── A. STRUCTURAL TESTS FOR eprocess_betting ─────────────────────────────────

test_that("eprocess_betting returns a well-formed data.frame", {
  ep_b <- eprocess_betting(s1, s2, c_t = 2)

  expect_true(is.data.frame(ep_b))
  expect_true(all(c("t", "e_pq", "e_qp", "log_e_pq", "log_e_qp") %in% names(ep_b)))
  expect_equal(nrow(ep_b), T_)
  expect_equal(ep_b$t, 1:T_)
})

test_that("eprocess_betting enforces lambda_t bounds and c_t positivity", {
  expect_error(
    eprocess_betting(s1, s2, c_t = -1),
    "c_t must be strictly positive"
  )

  expect_error(
    eprocess_betting(s1, s2, c_t = 2, lambda_t = rep(0.6, T_)), # 0.6 > 1/2
    "lambda_t must lie in"
  )

  expect_error(
    eprocess_betting(s1, s2, c_t = 2, lambda_t = rep(-0.1, T_)),
    "lambda_t must lie in"
  )
})

test_that("eprocess_betting satisfies sign symmetry", {
  ep_fwd <- eprocess_betting(s1, s2, c_t = 2)
  ep_rev <- eprocess_betting(s2, s1, c_t = 2)

  expect_equal(ep_fwd$e_pq, ep_rev$e_qp, tolerance = 1e-8)
  expect_equal(ep_fwd$e_qp, ep_rev$e_pq, tolerance = 1e-8)
})

test_that("eprocess_betting handles default lambda_t = 1/(2*c_t)", {
  ep_auto <- eprocess_betting(s1, s2, c_t = 2)
  ep_manual <- eprocess_betting(s1, s2, c_t = 2, lambda_t = rep(0.25, T_))

  expect_equal(ep_auto$e_pq, ep_manual$e_pq, tolerance = 1e-12)
})

# ── B. STRUCTURAL TESTS FOR lambda_betting_quantile ──────────────────────────

test_that("lambda_betting_quantile returns correctly shaped outputs", {
  # Mock the lagged score differences: start with 0, then the first T-1 differences
  mock_lag <- c(0, head(s1 - s2, -1))

  bnds <- lambda_betting_quantile(p_forecasts, q_forecasts, tau = 0.5, delta_hat_lag1 = mock_lag)

  expect_true(is.list(bnds))
  expect_true(all(c("c_t", "lambda_t") %in% names(bnds)))
  expect_equal(length(bnds$c_t), T_)
  expect_equal(length(bnds$lambda_t), T_)
})

test_that("lambda_betting_quantile enforces lambda_t <= 1/c_t everywhere", {
  mock_lag <- c(0, head(s1 - s2, -1))

  bnds <- lambda_betting_quantile(p_forecasts, q_forecasts, tau = 0.75, delta_hat_lag1 = mock_lag)

  # Because of the eps term and K_t >= 1, lambda_t should always be strictly < 1/c_t
  expect_true(all(bnds$lambda_t <= 1 / bnds$c_t))
})

test_that("lambda_betting_quantile floors tied forecasts for direct use", {
  eps <- 1e-8

  bnds <- lambda_betting_quantile(
    p_t = rep(0.5, 3),
    q_t = rep(0.5, 3),
    tau = 0.5,
    delta_hat_lag1 = rep(0, 3),
    eps = eps
  )

  expect_equal(bnds$c_t, rep(eps, 3))
  expect_true(all(bnds$lambda_t >= 0))
  expect_true(all(bnds$lambda_t <= 1 / bnds$c_t))

  expect_no_error(
    eprocess_betting(
      scores1 = rep(0, 3),
      scores2 = rep(0, 3),
      c_t = bnds$c_t,
      lambda_t = bnds$lambda_t
    )
  )
})

test_that("lambda_betting_quantile calculates correctly for a known input", {
  # Manual calculation:
  # p=0.8, q=0.5, tau=0.75, delta_lag1=0
  # c_t = 2 * max(0.75, 0.25) * |0.8 - 0.5| = 1.5 * 0.3 = 0.45
  # u = |0.75 - 0.5| = 0.25
  # K_t = ((2 - 0.25) / 1.25) * ((1.5*pi + 0) / pi) = (1.75 / 1.25) * 1.5 = 1.4 * 1.5 = 2.1
  # lambda_t = 1 / (2.1 * 0.45 + 1e-8)

  bnds <- lambda_betting_quantile(0.8, 0.5, tau = 0.75, delta_hat_lag1 = 0, eps = 1e-8)

  expect_equal(bnds$c_t, 0.45, tolerance = 1e-10)
  expect_equal(bnds$lambda_t, 1 / (2.1 * 0.45 + 1e-8), tolerance = 1e-10)
})

test_that("lambda_betting_quantile validates its inputs", {
  expect_error(lambda_betting_quantile(p_forecasts, q_forecasts, tau = -0.1), "tau must lie strictly")
  expect_error(lambda_betting_quantile(p_forecasts, q_forecasts[1:50], tau = 0.5), "same length")
  expect_error(lambda_betting_quantile(p_forecasts, q_forecasts, tau = 0.5, delta_hat_lag1 = c(0, NA)), "same length")
})

# ── C. STRUCTURAL TESTS FOR lambda_betting_agrapa AND lambda_betting_ons ─────

test_that("lambda_betting_ons matches hand-derived trace with constant-c", {
  xs <- c(0.6, -0.2, 0.4)
  lam <- lambda_betting_ons(xs, c = 2)
  # Derived with eta = 2 / (2 - log(3))
  # t=1: y=0.3, grad=-0.3, A=1.09 -> raw=0.6104781
  # t=2: y=-0.1, grad=0.1065016, A=1.101343 -> raw=0.3959897
  expected_Y <- c(0, 0.6104781, 0.3959897)
  expect_equal(lam, expected_Y / 2, tolerance = 1e-3)
})

test_that("lambda_betting_ons matches hand-derived trace with time-varying c_t", {
  xs <- c(0.6, -0.2, 0.4)
  c_t <- c(2, 4, 1)
  lam <- lambda_betting_ons(xs, c = c_t)

  # The raw Y-scale bets (lam_Y) before final division by c_t:
  # t=1: 0
  # t=2: 0 - 2.21807 * (-0.3 / 1.09) = 0.610478
  # t=3: 0.610478 - 2.21807 * (0.051574 / 1.092660) = 0.505784
  expected_Y <- c(0, 0.610478, 0.505784)

  expect_equal(lam, expected_Y / c_t, tolerance = 1e-3)
})

test_that("lambda_betting_periodic with base_fun=lambda_betting_ons matches hand-derived trace", {
  xs <- c(0.6, -0.2, 0.4, 0.5)
  c_t <- 2
  period <- 2
  lam <- lambda_betting_periodic(xs, c = c_t, period = period, base_fun = lambda_betting_ons)
  expected_Y <- c(0, 0, 0.6104, 0)

  expect_equal(lam, expected_Y / c_t, tolerance = 1e-3)
})

test_that("lambda_betting_agrapa and lambda_betting_ons stay within Prop 3.2 bound", {
  set.seed(1)
  xs <- runif(50, -1, 1)  # c = 2 => bound is [-1,1]
  expect_true(all(lambda_betting_agrapa(xs, c = 2) >= 0))
  expect_true(all(lambda_betting_agrapa(xs, c = 2) <= 0.5))
  expect_true(all(lambda_betting_ons(xs, c = 2) >= 0))
  expect_true(all(lambda_betting_ons(xs, c = 2) <= 0.5))
})

test_that("both lambda selectors feed eprocess_betting to a valid, finite e-process", {
  set.seed(2)
  n <- 100
  s1 <- runif(n, -0.5, 0.5)
  s2 <- runif(n, -0.5, 0.5)
  xs <- s1 - s2
  xs <- pmax(pmin(xs, 1), -1)  # enforce |xs| <= c/2 for c = 2
  lam <- lambda_betting_ons(xs, c = 2)
  ep <- eprocess_betting(s1, s2, c_t = 2, lambda_t = lam)
  expect_true(all(is.finite(ep$e_pq)))
  expect_true(all(is.finite(ep$e_qp)))
  expect_true(all(ep$e_pq >= 0))
  expect_true(all(ep$e_qp >= 0))
})
