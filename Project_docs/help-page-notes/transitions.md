# Implementation notes removed from `?transitions` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/filter-transitions.R); it is kept here for maintainers.

- **Default order** (`sort = "none"`) is the column-major order of the
  `K x K` matrices (`to` varies slowest).
- **Sorting** puts non-finite keys last (`-Inf`) and orders decreasingly.
- **Inference methods** (R/transitions-inference.R) are thin filters over
  `as.data.frame()` of the result: `significant = TRUE` keeps rows where
  the class's flag column is `TRUE` (`adj_res_stable` for bootstrap,
  `stable` for certainty and stability, `significant` for permutation
  and comparisons); `sort = "strength"` orders by the absolute value of
  `adj_res_observed`, `adj_res_observed`, `stability`,
  `observed_adj_res`, `diff` respectively.
- **lsa_lags method** applies the filters per lag (each fit's own
  `alpha`) and row-binds.
- **Projection**: the integer ids and the duplicate `weight`/`edge`
  protocol columns are dropped; any extra engine column (beyond the
  standard schema) is kept.
