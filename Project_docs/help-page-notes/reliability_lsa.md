# Implementation notes removed from `?reliability_lsa` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/reliability_lsa.R); it is kept here for
maintainers.

- **Split**: half 1 is `sample.int(S, max(1, S %/% 2))` sequences drawn
  without replacement; half 2 is the rest (one more sequence when `S` is
  odd). Each half is refitted with `.refit_from_events()` and the fit's
  recipe.
- **Failed replicates**: a half made only of single-event sequences has
  zero transitions; the engine raises `lagdynamics_no_transitions`, which
  is caught for that replicate only (`tryCatch(...,
  lagdynamics_no_transitions = function(e) NULL)`) and the replicate is
  `NA`. Any other error propagates. A replicate is also `NA` when fewer
  than two cells are finite in both halves or either half has zero
  variance over those cells. `n_failed = sum(!is.finite(correlations))`
  counts both kinds.
- **Correlation**: over cells finite in both halves
  (`stats::cor(method = method)`), cells in `as.vector()` order.
- **Summary**: `mean` over finite correlations (`NA` if none); `sd` and
  the 2.5% / 97.5% quantiles need at least two finite correlations. The
  interval level is fixed at 95%; there is no argument for it.
- **Stored elements**: `correlations`, `mean`, `sd`, `ci_low`,
  `ci_high`, `R`, `weights`, `method`, `n_sequences`, `n_failed`, `fit`.
  Objects built before `n_failed` existed fall back to counting
  non-finite correlations (`%||%`) in `print()` and the summary table.
- **Grouped**: `reliability_lsa.lsa_group()` applies the method to each
  group with the same `...`; container class
  `c("lsa_reliability_group", "list")` with a `levels` attribute;
  `as.data.frame()` uses `.grouped_df()`, which forwards `...` (so
  `what = "summary"` reaches each group).
