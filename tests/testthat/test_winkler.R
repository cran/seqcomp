# =============================================================================
# tests/testthat/test_winkler.R
# =============================================================================

set.seed(42)
n <- 500
y <- rbinom(n, 1, 0.5)
# Two forecasters: p slightly better than q on average
p <- pmin(pmax(rbeta(n, 3, 2), 0.01), 0.99)   # skewed toward 1
q <- pmin(pmax(rbeta(n, 2, 2), 0.01), 0.99)   # symmetric

# --- A. winkler_cs: structure ---

test_that("winkler_cs returns a well-formed one-sided data.frame", {
  cs_w <- winkler_cs(p, q, y, alpha = 0.05)

  expect_true(is.data.frame(cs_w))
  expect_true(all(c("t", "estimate", "lower", "upper") %in% names(cs_w)))
  expect_equal(nrow(cs_w), n)
  expect_true(all(is.infinite(cs_w$lower) & cs_w$lower < 0))
  expect_true(all(is.finite(cs_w$upper)))
})

test_that("winkler_cs estimate is the running mean of the Winkler score", {
  cs_w <- winkler_cs(p, q, y, alpha = 0.05)
  ws <- winkler_score(p, q, y)
  expect_equal(cs_w$estimate, cumsum(ws) / seq_along(ws), tolerance = 1e-10)
})

test_that("winkler_cs upper bound always exceeds the estimate", {
  cs_w <- winkler_cs(p, q, y, alpha = 0.05)
  expect_true(all(cs_w$upper > cs_w$estimate))
})

# --- B. winkler_etest: one-sided structure ---

test_that("winkler_etest returns a well-formed one-sided e-process data.frame", {
  ep_w <- winkler_etest(p, q, y, alpha = 0.05)

  expect_true(is.data.frame(ep_w))
  expect_true(all(c("t", "e", "log_e") %in% names(ep_w)))
  expect_equal(nrow(ep_w), n)
  expect_true(all(ep_w$e > 0))
  expect_true(all(ep_w$e <= 1e7))
  expect_true(all(is.finite(ep_w$log_e)))
  expect_true(all(ep_w$log_e <= log(1e7) + 1e-12))
})

test_that("winkler_etest_rejections returns a correctly structured one-sided summary", {
  ep_w <- winkler_etest(p, q, y, alpha = 0.05)
  rej_w <- seqcomp:::winkler_etest_rejections(ep_w, alpha = 0.05)

  expect_equal(rej_w$threshold, 20)
  expect_true(all(c("threshold", "tau", "reject") %in% names(rej_w)))
  expect_true(is.logical(rej_w$reject) && length(rej_w$reject) == 1)
})

# --- C. winkler_compare: combined one-sided pipeline ---

test_that("winkler_compare returns the expected top-level structure", {
  wc <- winkler_compare(p, q, y, alpha = 0.05)
  expect_true(all(c("winkler_scores", "cs", "etest_p_worse", "etest_q_worse",
                    "rejections") %in% names(wc)))
})

test_that("winkler_compare winkler_scores are correctly sized and bounded above by 1", {
  wc <- winkler_compare(p, q, y, alpha = 0.05)
  expect_equal(length(wc$winkler_scores), n)
  expect_true(all(wc$winkler_scores <= 1 + 1e-10))
})

test_that("winkler_compare's cs matches a standalone winkler_cs() call", {
  wc <- winkler_compare(p, q, y, alpha = 0.05)
  cs_w <- winkler_cs(p, q, y, alpha = 0.05)
  expect_equal(wc$cs$upper, cs_w$upper, tolerance = 1e-10)
})

test_that("winkler_compare's etest components have the one-sided e-process columns", {
  wc <- winkler_compare(p, q, y, alpha = 0.05)

  expect_true(is.data.frame(wc$etest_p_worse))
  expect_true(all(c("t", "e", "log_e") %in% names(wc$etest_p_worse)))
  expect_true(is.data.frame(wc$etest_q_worse))
  expect_true(all(c("t", "e", "log_e") %in% names(wc$etest_q_worse)))
})

test_that("winkler_compare's rejections are correctly nested with threshold 20", {
  wc <- winkler_compare(p, q, y, alpha = 0.05)

  expect_true(all(c("p_worse_than_q", "q_worse_than_p") %in% names(wc$rejections)))
  expect_equal(wc$rejections$p_worse_than_q$threshold, 20)
  expect_equal(wc$rejections$q_worse_than_p$threshold, 20)
})

# --- D. one-sided guarantee: upper bound covers true mean ---
#
# Simulation: N paths under a symmetric null construction.
# The upper CS should contain 0 in at least about 1-alpha of paths.
# Gated: too slow for routine devtools::check(). Run locally with:
#   SEQCOMP_RUN_SLOW_TESTS=true Rscript -e 'devtools::test()'

test_that("winkler_cs one-sided upper bound achieves nominal coverage under a symmetric null", {
  skip_on_cran()
  skip_if(Sys.getenv("SEQCOMP_RUN_SLOW_TESTS") != "true",
          "slow coverage simulation; set SEQCOMP_RUN_SLOW_TESTS=true to run")

  N_sim <- 500
  T_sim <- 300
  n_covered <- 0L
  set.seed(123)
  for (i in seq_len(N_sim)) {
    p_i <- pmin(pmax(rbeta(T_sim, 2, 2), 0.01), 0.99)
    q_i <- pmin(pmax(rbeta(T_sim, 2, 2), 0.01), 0.99)
    y_i <- rbinom(T_sim, 1, 0.5)
    cs_i <- winkler_cs(p_i, q_i, y_i, alpha = 0.05)
    if (all(0 <= cs_i$upper)) n_covered <- n_covered + 1L
  }
  coverage_w <- n_covered / N_sim

  expect_gte(coverage_w, 0.95)
})

# --- E. two-sided CS with lower_bound ---
#
# A finite lower_bound activates the two-sided CS branch. brier_score keeps
# the underlying score bounded here. The supplied lower_bound is interpreted
# as a user-supplied analytical bound on the Winkler scores; the function
# checks observed violations but cannot verify the mathematical guarantee.
# For unbounded scores such as the log score, no finite lower bound can be
# guaranteed in general, so the one-sided CS is the safer/default procedure.

test_that("winkler_cs with a finite lower_bound produces a valid two-sided CS", {
  cs_2s <- winkler_cs(
    p = p, q = q, y = y, alpha = 0.05,
    base_score = brier_score,
    lower_bound = -10
  )

  expect_true(all(is.finite(cs_2s$lower)))
  expect_true(all(is.finite(cs_2s$upper)))
  expect_true(all(cs_2s$lower < cs_2s$upper))
})
