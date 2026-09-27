# seqcomp 0.3.0

## New Features
* Added `lambda_betting_agrapa()` and `lambda_betting_ons()`, adaptive
  betting-fraction rules for the strong-null e-process (`eprocess_betting()`),
  adapting the aGRAPA and ONS-m algorithms of Waudby-Smith and Ramdas (2024)
  to Arnold et al. (2026)'s bounded strong-null setting. These adaptations
  are original to `seqcomp`, not given in either source paper.
* Added `build_agrapa_betting_array()` and `build_ons_betting_array()` to
  construct the T x m x m `lambda_param` arrays these rules require for
  `smcs_strong()`, including a `period` argument for known periodic
  (e.g. day-of-week) structure.
* Added a new vignette, "Choosing an Adaptive Betting Strategy: Temporal
  Dependence and Pathologies", comparing naive, aGRAPA, ONS-m, and Arnold
  et al.'s quantile-specific betting rules across structural breaks,
  autocorrelated noise, near-degenerate bounds, volatility clustering, and
  closed-testing dilution.

## Renamed
* `compare_multiple_forecasts()` is renamed to `smcs_compare()`. This
  function was never released on CRAN (only available on GitHub as part of
  unreleased 0.2.0 development), so this is not a breaking change for any
  CRAN user.

## Minor Improvements & Fixes
* Removed the floor of 1 on the intrinsic time in `cs_bernstein()` after
  verifying it does not cause any numerical issues (as for validity,
  the floor was a conservative measure in the fist place, so nothing
  that affects validity).

# seqcomp 0.2.0

## New Features
* Introduced Sequential Model Confidence Sets (SMCS) for multi-model evaluations, based on Arnold et al. (2026).
* Added `smcs_compare()` as a high-level wrapper to sequentially compare 3 or more forecasters simultaneously.
* Added `smcs_strong()` and `smcs_weak()` to construct confidence sets under the strong, uniformly weak, and weak null hypotheses.
* Implemented predictable betting-style e-processes (`eprocess_betting()`) and adaptive betting fractions for quantile forecasts (`lambda_betting_quantile()`).
* Added `vovk_wang_merge()` for highly efficient $O(m \log m)$ closed-testing e-value merging.

## Minor Improvements & Fixes
* Updated documentation and citations to reflect the official publication of Choe & Ramdas (2024) in *Operations Research*, 72(4).

# seqcomp 0.1.0

* Initial CRAN submission.
