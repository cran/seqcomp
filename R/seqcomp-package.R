#' seqcomp: Sequential Comparison of Probabilistic Forecasts
#'
#' `seqcomp` provides tools for comparing probabilistic forecasters
#' sequentially, following the anytime-valid framework of Choe and Ramdas
#' (2024). For three or more forecasters, the package additionally implements
#' Sequential Model Confidence Sets following Arnold, Gavrilopoulos, Schulz
#' and Ziegel (2026).
#'
#' The package is built around the score difference
#'
#' \deqn{\hat{\delta}_t = S(p_t, y_t) - S(q_t, y_t),}
#'
#' where scores are positively oriented, so larger values are better. Positive
#' score differences favor forecaster `p`; negative score differences favor
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
#' @section Multiple forecasters:
#' For three or more candidate forecasters, use [smcs_compare()]
#' as the main entry point. It extends the pairwise workflow above to a
#' Sequential Model Confidence Set (SMCS): the running set of models not yet
#' shown to underperform some competitor. Lower-level access is available via
#' [smcs_strong()] (strong or uniformly weak null, via closure testing) and
#' [smcs_weak()] (time-varying weak null, via joint confidence sequences),
#' with [vovk_wang_merge()] and [eprocess_betting()] as supporting building
#' blocks. To construct optimal adaptive betting fractions for the strong null,
#' use [build_agrapa_betting_array()] or [build_ons_betting_array()] for
#' constant-bound scores, and [build_quantile_betting_arrays()] for conditionally
#' bounded rules such as tick loss.
#'
#' @section Winkler scores:
#' For binary probability forecasts with unbounded base scores, use
#' [winkler_score()], [winkler_cs()], [winkler_etest()], or
#' [winkler_compare()].
#'
#' @references
#' Arnold, S., Gavrilopoulos, G., Schulz, B. and Ziegel, J. (2026). Sequential
#' model confidence sets. Journal of the Royal Statistical Society Series B:
#' Statistical Methodology, qkag066.
#'
#' Choe, Y. J. and Ramdas, A. (2024). Comparing Sequential Forecasters.
#' Operations Research, 72(4), 1368-1387.
#'
#' Howard, S. R., Ramdas, A., McAuliffe, J. and Sekhon, J. (2021).
#' Time-uniform, nonparametric, nonasymptotic confidence sequences.
#' The Annals of Statistics, 49(2), 1055-1080.
#'
#' Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded
#' random variables by betting. Journal of the Royal Statistical Society
#' Series B: Statistical Methodology, 86(1), 1-27.
#'
#' @name seqcomp-package
#' @aliases seqcomp seqcomp-package
NULL
