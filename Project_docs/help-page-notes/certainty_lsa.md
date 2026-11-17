# Implementation notes removed from `?certainty_lsa` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/certainty_lsa.R); it is kept here for
maintainers.

- **Origin**: ported from `Nestimate::certainty()` (Johnston & Jendoubi,
  2026) and adapted to the `obs` count matrix of an `lsa` fit.
- **Posterior**: for cell (i, j), `a = count[i, j] + prior` and
  `b = rowSums(count + prior)[i] - a`; the marginal posterior is
  `Beta(a, b)`. `prob_mean = a / (a + b)`,
  `prob_se = sqrt(a * b / ((a + b)^2 * (a + b + 1)))`, interval from
  `qbeta()` at `(1 - level_alpha) / 2` and its complement.
- **Decision, stability**: the band is
  `[min(W * r1, W * r2), max(W * r1, W * r2)]` with `W` the observed
  (MLE) probability `fit$prob`; `p_value = pbeta(lo) + 1 - pbeta(hi)`
  (posterior mass outside the band).
- **Decision, threshold**: `p_value = pbeta(edge_threshold)` (mass below
  the threshold). A `NULL` threshold is
  `quantile(W[W > 0 & finite], 0.10)`, or 0 when no such cell exists.
- **Certain**: `stable = valid & W > 0 & p_value < 1 - level_alpha`.
  `valid` excludes non-finite `fit$prob` (zero-margin rows) and declared
  structural zeros; for those cells `prob_mean`, `prob_se`, the interval
  and `p_value` are `NA`.
- **Class**: `c("lsa_certainty", "lsa_bootstrap", "list")` so
  `as.data.frame.lsa_bootstrap()`, `plot_forest()` and the bootstrap
  `transitions()` path apply. `adj_res_stable` is a copy of `stable`,
  carried only so `plot_forest()` can colour spokes. `R` is
  `NA_integer_` and `level` is `"analytic"`.
- **Stored matrices**: `mean`, `sd`, `ci_lower`, `ci_upper`, `p_values`
  (K x K), plus `prior`, `level_alpha`, `inference`,
  `consistency_range`, `edge_threshold` (the resolved value) and `fit`.
- **print()**: "certain edges: k of n" where n counts cells with a
  finite positive observed probability.
- **Grouped input**: `certainty_lsa()` is applied to each group fit; the
  container is `c("lsa_certainty_group", "list")`; `as.data.frame()`
  uses `.grouped_df()`.
- Runs in microseconds: no resampling.
