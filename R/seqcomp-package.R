#' seqcomp: Sequential Comparison of Probabilistic Forecasts
#'
#' `seqcomp` provides tools for comparing probabilistic forecasters
#' sequentially, following the anytime-valid framework of Choe and Ramdas
#' (2024).
#'
#' The package is built around the score difference
#'
#' \deqn{\hat{\delta}_t = S(p_t, y_t) - S(q_t, y_t),}
#'
#' where scores are positively oriented, so larger values are better. Positive
#' score differences favour forecaster `p`; negative score differences favour
#' forecaster `q`.
#'
#' @section Main workflow:
#' For most applications, start with [compare_forecasts()]. It computes
#' pointwise scores, running mean score differences, confidence sequences, and
#' e-processes in one call.
#'
#' @section Scoring rules:
#' The package includes positively oriented scoring rules such as
#' [brier_score()], [log_score()], [spherical_score()], [tick_loss()],
#' [qlike_score()], [winkler_score()], [crps_normal()], [crps_empirical()],
#' and [crps_std()].
#'
#' @section Confidence sequences:
#' Use [cs_hoeffding()] for Hoeffding-style confidence sequences,
#' [cs_bernstein()] for empirical Bernstein confidence sequences, and
#' [cs_asymptotic()] for asymptotic confidence sequences when finite-sample
#' boundedness is not available.
#'
#' @section E-processes:
#' Use [eprocess()] for the main sub-exponential mixture e-process and
#' [eprocess_rejections()] to extract first rejection times. For multi-step
#' forecasts, see [eprocess_lag()]. For predictable time-varying bounds, see
#' [eprocess_predictable()].
#'
#' @section Winkler scores:
#' For binary probability forecasts with unbounded base scores, use
#' [winkler_score()], [winkler_cs()], [winkler_etest()], or
#' [winkler_compare()].
#'
#' @references
#' Choe, Y. J. and Ramdas, A. (2024). Comparing Sequential Forecasters.
#' Operations Research, 72(4), 1368-1387.
#'
#' Howard, S. R., Ramdas, A., McAuliffe, J. and Sekhon, J. (2021).
#' Time-uniform, nonparametric, nonasymptotic confidence sequences.
#' The Annals of Statistics, 49(2).
#'
#' @name seqcomp-package
#' @aliases seqcomp seqcomp-package
NULL
