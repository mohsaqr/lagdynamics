# Implementation notes removed from `?compare_lsa` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/compare_lsa.R); it is kept here for
maintainers.

- **Permutation**: per pair, the per-sequence event vectors of both fits
  are pooled; each replicate draws `sample.int(S)` and assigns the first
  `n_a` sequences to pseudo-group a and the rest to b, so group sizes (in
  sequences) are preserved. Each pseudo-group is refitted with
  `.refit_from_events()` using the first fit's recipe.
- **Per-transition p-value**: two-sided,
  `p_perm = (1 + #{|d*| >= |d_obs|}) / (1 + n_finite)` where `n_finite`
  is the number of finite permuted differences for that cell (not a
  flat `R`); the old help text gave `(1 + R)` as the denominator, which
  is only true when every replicate is finite. Exceedances are counted
  with `na.rm = TRUE`. Cells with a non-finite observed difference, no
  finite replicate, or pooled count `obs_a + obs_b < min_count` get
  `p_perm = NA`.
- **Tested set**: `tested <- is.finite(p_perm)` is the single source of
  truth for both the family correction and the omnibus, so observed and
  permuted omnibus statistics use the same cells.
- **Omnibus**: `sum(diff[tested]^2)` (squared Frobenius norm of the
  difference matrix over tested cells); null from the same permutations;
  `p = (1 + #{null >= obs}) / (1 + R)`. `NA` when no cell is testable.
- **significant**: `is.finite(p_adj) & p_adj < alpha`, where `alpha` is
  the first fit's `params$alpha`.
- **Pairwise (> 2 groups)**: `utils::combn()` pairs; each pair is run by
  `.compare_two()` with `adjust = "none"`; the per-edge `p_perm` of all
  pairs are pooled and `p.adjust()`ed once (family-wise); the per-pair
  omnibus p-values are adjusted separately across pairs. The underlying
  two-group objects are kept in `comparisons`, named `"<a>_vs_<b>"`.
- **log_or** (`.lsa_log_or()`): for cell i -> j the 2 x 2 collapse is
  `a = obs[i, j]`, `b = row_i - a`, `c = col_j - a`,
  `d = N - row_i - col_j + a`; if any of the four is zero, 0.5 is added
  to all four (Haldane-Anscombe). Cells with `NA` adjusted residual in
  the engine output (structural zeros, zero-margin rows) are `NA`.
  Yule's Q saturates at +/-1 on zero cells, which is why `log_or` is the
  default.
- **Measure matrices** (`.lsa_compare_matrix()`): NA/Inf are kept, never
  coerced to 0; `lift = obs / exp`.
- **adj_res message**: emitted when `max(N) / min(N) > 1.5` over the
  fits' total transition counts.
- **Validation** (`.compare_validate_all()`): every fit must be
  event-level; labels, `params$engine`, `params$lag` and
  `params$structural_zeros` must be `identical()` to the first fit's.
- **Two bare fits**: group names are `c("a", "b")`; for an `lsa_group`
  without names they are `"g1"`, `"g2"`, ...
- **Interpretation caveat**: the odds ratio is non-collapsible; per-group
  log odds ratios are group-specific departures from independence and
  should not be pooled across groups with different marginal state
  distributions.
- **Two-group global** is stored as a list `(statistic, p_value, R)`;
  `as.data.frame(x, what = "global")` turns it into a one-row frame with
  `group_a`, `group_b` prepended (`.comparison_global_df()`).
