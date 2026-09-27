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
set.seed(2026)
T_sim <- 300
y <- rbinom(T_sim, size = 1, prob = 0.7)

# Create 3 forecasters
forecasts <- matrix(0, nrow = T_sim, ncol = 3)
colnames(forecasts) <- c("M1", "M2", "M3")

forecasts[, "M1"] <- ifelse(y == 1, 0.8, 0.2)   # Highly accurate
forecasts[, "M2"] <- runif(T_sim, 0.4, 0.6)     # Mediocre/random
forecasts[, "M3"] <- ifelse(y == 1, 0.2, 0.8)   # Consistently wrong

head(forecasts)


## -----------------------------------------------------------------------------
multi_cmp <- smcs_compare(
  forecasts = forecasts,
  outcomes = y,
  scoring_rule = "brier",
  alpha = 0.05
)

# Printing the object provides a clean summary of which models 
# were excluded and when they first dropped out of the set.
multi_cmp


## -----------------------------------------------------------------------------
tail(multi_cmp$smcs_strong)


## -----------------------------------------------------------------------------
par(mfrow = c(3, 1), mar = c(2, 4, 2, 1))
colors <- c("blue", "gray", "red")

for (i in 1:3) {
  plot(
    1:T_sim, multi_cmp$smcs_strong[, i], 
    type = "s", col = colors[i], lwd = 2,
    ylim = c(-0.1, 1.1), yaxt = "n", ylab = "In SMCS?",
    main = paste("Model:", colnames(forecasts)[i])
  )
  axis(2, at = c(0, 1), labels = c("Excluded", "Included"), las = 2)
}
par(mfrow = c(1, 1))


## -----------------------------------------------------------------------------
# Manually compute Brier scores
scores_mat <- matrix(0, nrow = T_sim, ncol = 3)
for(i in 1:3) scores_mat[, i] <- brier_score(forecasts[, i], y)

# Construct the Weak SMCS directly
# Brier score differences are bounded in [-1, 1], so we use c_param = 2
res_weak <- smcs_weak(
  scores = scores_mat,
  alpha = 0.05,
  cs_method = "bernstein",
  c_param = 2
)

tail(res_weak$smcs)


## -----------------------------------------------------------------------------
lam_agrapa <- build_agrapa_betting_array(scores_mat, c_mat = 2)

res_strong_agrapa <- smcs_strong(
  scores_mat,
  alpha        = 0.05,
  method       = "betting",
  c_param      = 2,
  lambda_param = lam_agrapa
)

tail(res_strong_agrapa$smcs)


## ----eval = FALSE-------------------------------------------------------------
# lam_agrapa_p7 <- build_agrapa_betting_array(scores_mat, c_mat = 2, period = 7)

## -----------------------------------------------------------------------------
set.seed(7)
tau <- 0.5

truth <- rnorm(T_sim)
q_forecasts <- cbind(
  M1 = truth + rnorm(T_sim, sd = 0.1),    # accurate
  M2 = rep(0, T_sim),                     # uninformative
  M3 = truth + rnorm(T_sim, sd = 0.1) + 1 # biased
)

tick_cmp <- smcs_compare(
  forecasts    = q_forecasts,
  outcomes     = truth,
  scoring_rule = "tick",
  tau          = tau,
  alpha        = 0.05
)

tail(tick_cmp$smcs_strong)

