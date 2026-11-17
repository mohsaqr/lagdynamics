# Implementation notes removed from `?plot_forest` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/plot-forest.R); it is kept here for maintainers.

- **Estimate column**: `<metric>_observed`, or `observed` for counts;
  intervals `<metric>_ci_low` / `<metric>_ci_high`.
- **Significance**: `adj_res_stable` of the bootstrap, direction from the
  sign of `adj_res_observed`, for every metric.
- **Ordering**: spokes sorted alphabetically by `"from -> to"`; `n_top`
  keeps the largest |estimate| but preserves that order. Rows with a
  non-finite estimate or bound are dropped.
- A certainty result (`lsa_certainty` inherits `lsa_bootstrap`) is
  accepted and works with `metric = "prob"`.
