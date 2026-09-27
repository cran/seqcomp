## ----include = FALSE----------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)


## ----setup--------------------------------------------------------------------
library(seqcomp)

## -----------------------------------------------------------------------------
set.seed(2027)
T_sim <- 1000
c_t <- rep(2, T_sim)
mu  <- ifelse(seq_len(T_sim) <= 600, 0.10, 0.40)
d_t <- pmin(pmax(mu + rnorm(T_sim), -c_t / 2), c_t / 2)

lam_agrapa <- lambda_betting_agrapa(d_t, c = c_t)
lam_ons    <- lambda_betting_ons(d_t, c = c_t)

plot(seq_len(T_sim), lam_agrapa, type = "l", col = "red",
     xlab = "t", ylab = expression(lambda[t]),
     main = "One simulated path: betting fraction around a break at t = 600")
lines(seq_len(T_sim), lam_ons, col = "orange")
abline(v = 600, col = "blue", lty = 3)
legend("topleft", legend = c("aGRAPA", "ONS-m"), col = c("red", "orange"),
       lwd = 1, bty = "n")


## ----eval = FALSE-------------------------------------------------------------
# simulate_dt <- function(T_, mu_fn, c_fn, rho = 0, sd_scale = 1) {
#   mu  <- vapply(seq_len(T_), mu_fn, numeric(1))
#   c_t <- vapply(seq_len(T_), c_fn, numeric(1))
#   eps <- numeric(T_)
#   eps[1] <- rnorm(1, sd = sd_scale)
#   for (t in 2:T_) {
#     eps[t] <- rho * eps[t - 1] + rnorm(1, sd = sd_scale * sqrt(1 - rho^2))
#   }
#   d_t <- pmin(pmax(mu + eps, -c_t / 2), c_t / 2)
#   list(d_t = d_t, c_t = c_t)
# }
# 
# mu_break   <- function(t) if (t <= 600) 0.10 else 0.40
# c_const_fn <- function(t) 2.0
# 
# N_runs <- 200
# T_sim  <- 1000
# lam_pre  <- matrix(NA, N_runs, 2, dimnames = list(NULL, c("aGRAPA", "ONS-m")))
# lam_post <- lam_pre
# 
# for (sim in seq_len(N_runs)) {
#   dat   <- simulate_dt(T_sim, mu_break, c_const_fn, rho = 0, sd_scale = 1)
#   lam_a <- lambda_betting_agrapa(dat$d_t, c = dat$c_t)
#   lam_o <- lambda_betting_ons(dat$d_t, c = dat$c_t)
#   lam_pre[sim, ]  <- c(mean(lam_a[400:600]),  mean(lam_o[400:600]))
#   lam_post[sim, ] <- c(mean(lam_a[800:1000]), mean(lam_o[800:1000]))
# }
# 
# colMeans(lam_pre)   # average bet just before the break
# colMeans(lam_post)  # average bet well after the break

## ----echo=FALSE---------------------------------------------------------------
tab <- data.frame(
  Method = c("Oracle", "ONS-m", "aGRAPA", "Naive"),
  Pre_break = c(0.132, 0.139, 0.132, 0.250),
  Post_break = c(0.494, 0.353, 0.255, 0.250),
  check.names = FALSE
)
knitr::kable(tab, format = "markdown")

## ----echo=FALSE---------------------------------------------------------------
tab2 <- data.frame(
  Break_t = c(100, 300, 600, 800),
  aGRAPA = c(1.975e37, 6.554e24, 1.567e11, 2.920e4),
  ONS_m = c(7.295e37, 8.168e27, 6.591e13, 3.626e5)
)
knitr::kable(tab2, format = "markdown", col.names = c(
  "Break t*", "aGRAPA post-break multiplier", "ONS-m post-break multiplier"))

## ----eval = FALSE-------------------------------------------------------------
# rho_vals <- c(0.0, 0.5, 0.9)
# N_runs <- 200
# T_sim  <- 1000
# get_mdd <- function(e) max(cummax(e) / e)
# 
# for (r in rho_vals) {
#   final_e <- matrix(NA, N_runs, 2, dimnames = list(NULL, c("aGRAPA", "ONS-m")))
#   mdd     <- final_e
#   for (sim in seq_len(N_runs)) {
#     dat   <- simulate_dt(T_sim, mu_fn = function(t) 0.05, c_fn = function(t) 2.0,
#                          rho = r, sd_scale = 1.0)
#     lam_a <- lambda_betting_agrapa(dat$d_t, c = dat$c_t)
#     lam_o <- lambda_betting_ons(dat$d_t, c = dat$c_t)
#     e_a   <- cumprod(1 + lam_a * dat$d_t)
#     e_o   <- cumprod(1 + lam_o * dat$d_t)
#     final_e[sim, ] <- c(e_a[T_sim], e_o[T_sim])
#     mdd[sim, ]     <- c(get_mdd(e_a), get_mdd(e_o))
#   }
#   # medians reported in the table below
# }

## ----echo=FALSE---------------------------------------------------------------
tab3 <- data.frame(
  rho = c(0.0, 0.5, 0.9),
  aGRAPA_final_wealth = c(0.44, 0.80, 1.08),
  ONS_m_final_wealth = c(0.36, 43.60, format(1158520000000, big.mark = ",")),
  aGRAPA_drawdown = c(11.52, 75.22, format(98600.77, big.mark = ",")),
  ONS_m_drawdown = c(26.29, 72.71, format(30731.16, big.mark = ","))
)
knitr::kable(tab3, format = "markdown", col.names = c(
  "$\\rho$", "aGRAPA final wealth", "ONS-m final wealth", "aGRAPA drawdown", 
  "ONS-m drawdown"), escape = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# simulate_garch <- function(T_, omega = 0.05, alpha = 0.15, beta = 0.80) {
#   y <- sigma2 <- numeric(T_)
#   sigma2[1] <- omega / (1 - alpha - beta)
#   y[1] <- rnorm(1, sd = sqrt(sigma2[1]))
#   for (t in 2:T_) {
#     sigma2[t] <- omega + alpha * y[t - 1]^2 + beta * sigma2[t - 1]
#     y[t] <- rnorm(1, sd = sqrt(sigma2[t]))
#   }
#   list(y = y, sigma = sqrt(sigma2))
# }
# 
# set.seed(2029)
# T_sim <- 1000
# dgp <- simulate_garch(T_sim)
# 
# q_oracle <- rep(0, T_sim)
# q_micro  <- q_oracle + 0.001 * rep(c(1, -1), T_sim / 2)
# 
# d_t <- tick_loss(q_oracle, dgp$y, tau = 0.5) - tick_loss(q_micro, dgp$y, tau = 0.5)
# 
# delta_lag1 <- c(0, head(d_t, -1))
# bnds <- lambda_betting_quantile(q_oracle, q_micro, tau = 0.5,
#                                 delta_hat_lag1 = delta_lag1)
# c_t <- pmax(bnds$c_t, 1e-8)
# 
# lam_naive  <- 1 / (2 * c_t)
# lam_arnold <- pmin(bnds$lambda_t, 1 / c_t)
# lam_agrapa <- lambda_betting_agrapa(d_t, c = c_t)
# lam_ons    <- lambda_betting_ons(d_t, c = c_t)

## ----echo=FALSE---------------------------------------------------------------
tab4 <- data.frame(
  Method = c("Naive", "Arnold", "aGRAPA", "ONS-m"),
  Avg_lambda = c(500.0, 333.3, 5.0, 70.8),
  Terminal_wealth = c(0.0, 0.0, 0.2, 0.0)
)
knitr::kable(tab4, format = "markdown", col.names = c(
  "Method", "Average $\\lambda_t$ played", "Terminal wealth"), escape = FALSE)

## ----echo=FALSE---------------------------------------------------------------
tab5 <- data.frame(
  Method = c("Naive", "Arnold", "aGRAPA", "ONS-m"),
  Terminal_wealth = c(1.207e16, 9.728e12, 1.465e15, 9.792e14),
  Maximum_drawdown = c(19.11, 5.51, 65.23, 69.90)
)
knitr::kable(tab5, format = "markdown", col.names = c(
  "Method", "Terminal wealth", "Maximum drawdown"), escape = FALSE)


## ----eval = FALSE-------------------------------------------------------------
# T_sim <- 1000
# alpha <- 0.10
# sundays <- seq(7, T_sim, by = 7)
# m_vals <- c(2, 5, 15, 30, 49)
# 
# for (m_curr in m_vals) {
#   scores_mat <- matrix(0, nrow = T_sim, ncol = m_curr)
#   scores_mat[, 1] <- 1.0
#   scores_mat[sundays, 1] <- 0.0   # fails every seventh round
#   scores_mat[, 2] <- 1.0          # never fails
#   if (m_curr > 2) scores_mat[, 3:m_curr] <- 0.0   # uninformative decoys
# 
#   C_mat <- matrix(2.0, nrow = m_curr, ncol = m_curr)
#   lam_agrapa <- build_agrapa_betting_array(scores_mat, c_mat = C_mat, period = 7)
#   lam_ons    <- build_ons_betting_array(scores_mat, c_mat = C_mat, period = 7)
# 
#   res_agrapa <- smcs_strong(scores_mat, alpha = alpha, method = "betting",
#                             c_param = C_mat, lambda_param = lam_agrapa)
#   res_ons    <- smcs_strong(scores_mat, alpha = alpha, method = "betting",
#                             c_param = C_mat, lambda_param = lam_ons)
#   # round of exclusion for model 1 reported below
# }

## ----echo=FALSE---------------------------------------------------------------
tab6 <- data.frame(
  m = c(2, 5, 15, 30, 49),
  aGRAPA_period7 = c(63, 84, 105, 119, 126),
  ONS_m_period7 = c(63, 84, 105, 119, 126),
  Naive = c(98, 140, 182, 203, 217)
)
knitr::kable(tab6, format = "markdown", col.names = c(
  "$m$", "aGRAPA (period 7)", "ONS-m (period 7)", "Naive"), escape = FALSE)

## ----echo=FALSE---------------------------------------------------------------
tab7 <- data.frame(
  m = c(2, 10, 25, 49),
  aGRAPA_median_t = c(217.0, 301.0, 346.5, 339.5),
  ONS_m_median_t = c(280.0, 420.0, 439.0, 472.5),
  Naive_median_t = c(497.0, 836.5, 868.5, ">1,000")
)
knitr::kable(tab7, format = "markdown", col.names = c(
  "$m$", "aGRAPA (median $t$)", "ONS-m (median $t$)", "Naive (median $t$)"), escape = FALSE)


## ----echo=FALSE---------------------------------------------------------------
tab8 <- data.frame(
  m = c(2, 10, 25, 49),
  Arnold_exclusion = c(">8,000", ">8,000", ">8,000", ">8,000"),
  aGRAPA_exclusion = c(format(1786, big.mark = ","), format(2516, big.mark = ","),
                       format(2540, big.mark = ","), format(2561, big.mark = ",")),
  aGRAPA_final_e = c(format(10000000, big.mark = ","),
                     format(1111111.38, big.mark = ","),
                     format(416666.94, big.mark = ","),
                     format(208333.60, big.mark = ","))
)
knitr::kable(tab8, format = "markdown", col.names = c(
  "$m$", "Arnold: round of exclusion", "aGRAPA: round of exclusion", 
  "aGRAPA final adjusted e-value"), escape = FALSE)


## ----echo=FALSE---------------------------------------------------------------
tab9 <- data.frame(
  Scenario = c("Autocorrelated or abrupt mean shift", "Nearly identical models (tiny $c_t$)",
               "Strong volatility clustering", "Smooth, slowly varying, mean-reverting quantile forecasts"),
  Recommended_rule = c("ONS-m", "aGRAPA", "Naive or Arnold's heuristic (if drawdown matters)", 
                       "Arnold's heuristic")
)
knitr::kable(tab9, format = "markdown", col.names = c(
  "Situation observed in these simulations", "Reasonable starting point"), escape = FALSE)


