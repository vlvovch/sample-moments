# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] - 2025-12-09
- Added deterministic unit tests for `NumberStatistics` and `TwoNumberStatistics` (moments, cumulants, factorial moments/cumulants, ratios).
- Documented thread safety expectations and extended README with factorial moments/cumulants and build/test instructions.
- Added `CITATION.cff` for software citation metadata and included Kendall & Stuart reference.
- Expanded binomial examples to report factorial cumulants and ratios.
- Fixed buffer sizing in bulk observation loading, handled zero-division guards, and clarified ratio computations.
- History highlights leading up to 1.0.0:
  - 2021-01-10 to 2022-07-19: Built the header-only library with moments, central moments, cumulants, scaled cumulants, and factorial variants (single and joint); improved scaled cumulant error estimates; added covariance utilities and joint cumulants.
  - 2024-03-10 to 2024-08-08: Added ability to set moment sums externally, fixed joint moment ratio calculations, and expanded to factorial cumulants/moments; bug fixes for scaled cumulant ratios and external sums.
  - 2024-09-25 to 2024-09-26: CI/GitHub Actions updates, GoogleTest fixes, and matrix expansion.
  - 2025-12-09: Hardened robustness (off-by-one fixes, zero-event/zero-denominator/NaN guards), made mean-shift checks const, suppressed unused parameter warnings, made partition cache initialization thread-safe, added cumulant ratios to binomial example ahead of factorial extensions, added tests/docs/citation, and removed macOS-13 from CI.
