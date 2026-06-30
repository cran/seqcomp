## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----setup--------------------------------------------------------------------
library(seqcomp)

## -----------------------------------------------------------------------------
set.seed(1)

n <- 300
y <- rbinom(n, size = 1, prob = 0.55)

## -----------------------------------------------------------------------------
p <- ifelse(y == 1, 0.62, 0.38)
q <- rep(0.50, n)

head(data.frame(y = y, p = p, q = q))

## -----------------------------------------------------------------------------
cmp <- compare_forecasts(
  p = p,
  q = q,
  y = y,
  scoring_rule = "brier"
)

head(cmp)

## -----------------------------------------------------------------------------
plot(
  cmp$t, cmp$estimate,
  type = "l",
  ylim = range(c(cmp$lower, cmp$upper, 0), finite = TRUE),
  xlab = "Time",
  ylab = "Mean score difference",
  main = "Sequential comparison using the Brier score"
)
lines(cmp$t, cmp$lower, lty = 2)
lines(cmp$t, cmp$upper, lty = 2)
abline(h = 0, col = "gray50")

## -----------------------------------------------------------------------------
alpha <- 0.05
threshold <- 2 / alpha

threshold

## -----------------------------------------------------------------------------
plot(
  cmp$t, cmp$e_pq,
  type = "l",
  log = "y",
  xlab = "Time",
  ylab = "e-process value",
  main = "Evidence that p outperforms q"
)
abline(h = threshold, lty = 2, col = "gray50")

## -----------------------------------------------------------------------------
eprocess_rejections(cmp, alpha = alpha)

## -----------------------------------------------------------------------------
score_p <- brier_score(p, y)
score_q <- brier_score(q, y)

head(score_p)
head(score_q)

## -----------------------------------------------------------------------------
cs <- cs_bernstein(
  scores1 = score_p,
  scores2 = score_q,
  alpha = 0.05,
  c = 2
)

head(cs)

## -----------------------------------------------------------------------------
ep <- eprocess(
  scores1 = score_p,
  scores2 = score_q,
  alpha = 0.05,
  c = 2
)

head(ep)

## -----------------------------------------------------------------------------
cmp_brier <- compare_forecasts(p, q, y, scoring_rule = "brier")
tail(cmp_brier)

## -----------------------------------------------------------------------------
cmp_spherical <- compare_forecasts(p, q, y, scoring_rule = "spherical")
tail(cmp_spherical)

## -----------------------------------------------------------------------------
cmp_log <- compare_forecasts(
  p = p,
  q = q,
  y = y,
  scoring_rule = "log",
  compute_e = FALSE
)

tail(cmp_log)

## -----------------------------------------------------------------------------
wcmp <- winkler_compare(p, q, y)

names(wcmp)

