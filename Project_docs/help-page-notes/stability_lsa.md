# Implementation notes removed from `?stability_lsa` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/stability_lsa.R); it is kept here for
maintainers.

- **Subset size**: `n_keep = max(1, round(S * proportion))` sequences,
  drawn with `sample.int(S, n_keep, replace = FALSE)`.
- **Single sequence**: a contiguous window of
  `m = max(2, round(T * proportion))` events starting at a uniform
  position in `1..(T - m + 1)`.
- **Decision per replicate**: `as.vector(fb$p < recipe$alpha)`, i.e.
  significance in either direction at the fit's alpha. `NA` p-values
  (non-estimable cells) give `NA`, later set to `FALSE`.
- **Failed replicates**: a subset of only single-event sequences has no
  transitions; `lagdynamics_no_transitions` is caught for that replicate
  only and it returns an all-`NA` row. `n_failed` counts all-`NA` rows;
  all `NA`s are then set to `FALSE`, so failed replicates lower every
  cell's stability. Other errors propagate.
- **observed_sig**: `is.finite(p) & p < alpha` on the original fit.
- **stability_matrix**: `R x K^2` logical, cells in `as.vector()` order.
- **Validation**: `proportion` in (0, 1); `min_stable` in (0, 1] (the
  old help text said (0, 1)). Errors from `stopifnot()` are unclassed.
- **print()**: `n_failed %||% 0` for objects created before the field
  existed.
