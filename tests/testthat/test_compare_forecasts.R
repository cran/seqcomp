# =============================================================================
# tests/testthat/test_compare_forecasts.R
# =============================================================================

test_that("compare_forecasts returns a well-formed data.frame for Brier scores", {
  set.seed(101)
  T_ <- 100
  y <- rbinom(T_, 1, 0.5)
  p <- runif(T_)
  q <- runif(T_)

  out <- compare_forecasts(p, q, y, scoring_rule = "brier")

  expect_true(is.data.frame(out))
  expect_equal(nrow(out), T_)
  expect_equal(
    names(out),
    c("t", "score_p", "score_q", "delta", "estimate", "lower", "upper", "e_pq", "e_qp")
  )
  expect_equal(out$t, seq_len(T_))
  expect_true(all(is.finite(out$score_p)))
  expect_true(all(is.finite(out$score_q)))
  expect_true(all(is.finite(out$delta)))
  expect_true(all(is.finite(out$estimate)))
  expect_true(all(is.finite(out$lower)))
  expect_true(all(is.finite(out$upper)))
  expect_true(all(out$e_pq > 0))
  expect_true(all(out$e_qp > 0))
})

test_that("compare_forecasts computes the same Brier scores as brier_score", {
  set.seed(102)
  T_ <- 80
  y <- rbinom(T_, 1, 0.5)
  p <- runif(T_)
  q <- runif(T_)

  out <- compare_forecasts(p, q, y, scoring_rule = "brier")

  expect_equal(out$score_p, brier_score(p, y))
  expect_equal(out$score_q, brier_score(q, y))
  expect_equal(out$delta, brier_score(p, y) - brier_score(q, y))
  expect_equal(out$estimate, cumsum(out$delta) / seq_len(T_))
})

test_that("compare_forecasts agrees with direct cs_bernstein and eprocess calls", {
  set.seed(103)
  T_ <- 120
  y <- rbinom(T_, 1, 0.5)
  p <- runif(T_)
  q <- runif(T_)

  out <- compare_forecasts(
    p, q, y,
    scoring_rule = "brier",
    cs_type = "bernstein",
    compute_cs = TRUE,
    compute_e = TRUE,
    alpha = 0.05
  )

  sp <- brier_score(p, y)
  sq <- brier_score(q, y)

  cs_direct <- cs_bernstein(sp, sq, alpha = 0.05, c = 2)
  ep_direct <- eprocess(sp, sq, alpha = 0.05, c = 2)

  expect_equal(out$lower, cs_direct$lower, tolerance = 1e-10)
  expect_equal(out$upper, cs_direct$upper, tolerance = 1e-10)
  expect_equal(out$e_pq, ep_direct$e_pq, tolerance = 1e-10)
  expect_equal(out$e_qp, ep_direct$e_qp, tolerance = 1e-10)
})

test_that("compare_forecasts supports Hoeffding confidence sequences", {
  set.seed(104)
  T_ <- 100
  y <- rbinom(T_, 1, 0.5)
  p <- runif(T_)
  q <- runif(T_)

  out <- compare_forecasts(
    p, q, y,
    scoring_rule = "brier",
    cs_type = "hoeffding",
    compute_e = FALSE
  )

  sp <- brier_score(p, y)
  sq <- brier_score(q, y)
  cs_direct <- cs_hoeffding(sp, sq, alpha = 0.05, c = 1)

  expect_equal(out$lower, cs_direct$lower, tolerance = 1e-10)
  expect_equal(out$upper, cs_direct$upper, tolerance = 1e-10)
  expect_true(all(is.na(out$e_pq)))
  expect_true(all(is.na(out$e_qp)))
})

test_that("compare_forecasts supports spherical scores", {
  set.seed(105)
  T_ <- 100
  y <- rbinom(T_, 1, 0.5)
  p <- runif(T_)
  q <- runif(T_)

  out <- compare_forecasts(p, q, y, scoring_rule = "spherical")

  expect_equal(out$score_p, spherical_score(p, y))
  expect_equal(out$score_q, spherical_score(q, y))
  expect_true(all(is.finite(out$lower)))
  expect_true(all(is.finite(out$upper)))
  expect_true(all(out$e_pq > 0))
  expect_true(all(out$e_qp > 0))
})

test_that("compare_forecasts supports categorical matrix forecasts", {
  set.seed(106)
  T_ <- 60
  K <- 3

  raw_p <- matrix(runif(T_ * K), nrow = T_, ncol = K)
  raw_q <- matrix(runif(T_ * K), nrow = T_, ncol = K)

  p <- raw_p / rowSums(raw_p)
  q <- raw_q / rowSums(raw_q)
  y <- sample(seq_len(K), T_, replace = TRUE)

  out <- compare_forecasts(p, q, y, scoring_rule = "brier")

  expect_equal(out$score_p, brier_score(p, y))
  expect_equal(out$score_q, brier_score(q, y))
  expect_equal(nrow(out), T_)
  expect_true(all(is.finite(out$lower)))
  expect_true(all(is.finite(out$upper)))
})

test_that("compare_forecasts defaults to asymptotic CS for log scores", {
  set.seed(107)
  T_ <- 100
  y <- rbinom(T_, 1, 0.5)
  p <- runif(T_, 0.05, 0.95)
  q <- runif(T_, 0.05, 0.95)

  out <- compare_forecasts(
    p, q, y,
    scoring_rule = "log",
    compute_e = FALSE
  )

  sp <- log_score(p, y)
  sq <- log_score(q, y)
  cs_direct <- cs_asymptotic(sp, sq, alpha = 0.05)

  expect_equal(attr(out, "cs_type"), "asymptotic")
  expect_equal(out$score_p, sp)
  expect_equal(out$score_q, sq)
  expect_equal(out$lower, cs_direct$lower, tolerance = 1e-10)
  expect_equal(out$upper, cs_direct$upper, tolerance = 1e-10)
  expect_true(all(is.na(out$e_pq)))
  expect_true(all(is.na(out$e_qp)))
})

test_that("compare_forecasts refuses e-processes for unbounded log scores", {
  set.seed(108)
  T_ <- 50
  y <- rbinom(T_, 1, 0.5)
  p <- runif(T_, 0.05, 0.95)
  q <- runif(T_, 0.05, 0.95)

  expect_error(
    compare_forecasts(p, q, y, scoring_rule = "log", compute_e = TRUE),
    regexp = "bounded score differences"
  )
})

test_that("compare_forecasts refuses finite-sample CS for unbounded log scores", {
  set.seed(109)
  T_ <- 50
  y <- rbinom(T_, 1, 0.5)
  p <- runif(T_, 0.05, 0.95)
  q <- runif(T_, 0.05, 0.95)

  expect_error(
    compare_forecasts(
      p, q, y,
      scoring_rule = "log",
      cs_type = "bernstein",
      compute_e = FALSE
    ),
    regexp = "bounded score differences"
  )
})

test_that("compare_forecasts can suppress CS and e-process output", {
  set.seed(110)
  T_ <- 70
  y <- rbinom(T_, 1, 0.5)
  p <- runif(T_)
  q <- runif(T_)

  out <- compare_forecasts(
    p, q, y,
    scoring_rule = "brier",
    compute_cs = FALSE,
    compute_e = FALSE
  )

  expect_true(all(is.na(out$lower)))
  expect_true(all(is.na(out$upper)))
  expect_true(all(is.na(out$e_pq)))
  expect_true(all(is.na(out$e_qp)))
  expect_true(all(is.finite(out$estimate)))
})

test_that("compare_forecasts detects shape and length mismatches", {
  y <- c(0, 1, 0)
  p <- c(0.2, 0.7, 0.4)
  q <- c(0.3, 0.6)

  expect_error(
    compare_forecasts(p, q, y, scoring_rule = "brier"),
    regexp = "same number"
  )

  p_mat <- cbind(1 - p, p)

  expect_error(
    compare_forecasts(p_mat, p, y, scoring_rule = "brier"),
    regexp = "same shape"
  )
})
