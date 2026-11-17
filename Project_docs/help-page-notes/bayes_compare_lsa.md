# Implementation notes removed from `?bayes_compare_lsa` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/bayes_compare_lsa.R); it is kept here for
maintainers.

- **Origin**: ported from `Nestimate::bayes_compare()` (Johnston &
  Jendoubi, 2026). Reuses `.compare_resolve()` and
  `.compare_validate_all()` from R/compare_lsa.R.
- **Posterior**: per group, `a = obs + prior`, row concentration
  `rowSums(a)`; `prob = a / rowSums(a)` (posterior mean, not the MLE).
  `diff = prob_a - prob_b` is exact.
- **Monte Carlo**: one vectorised `rbeta(K^2 * draws, ...)` per group,
  arranged `K^2 x draws`; `dd = da - db`. Interval from per-cell
  `quantile()` at `(1 - ci) / 2` and its complement.
  `pd = max(mean(dd > 0), 1 - mean(dd > 0))`; `p_value = 2 * (1 - pd)`;
  `effect_size = diff / sd(dd)` (population sd from the draws), set to 0
  when not finite.
- **Valid cells**: finite `fit$prob` in both groups and not a declared
  structural zero (first fit's `params$structural_zeros`). Invalid cells
  get `p_value = NA` and are never significant (two-group path).
- **significant (two groups)**: `valid & (ci_low > 0 | ci_high < 0) &
  |diff| > mean_threshold & min(|ci_low|, |ci_high|) > bound_threshold`.
- **significant (pairwise)**: same rule with `valid` replaced by
  `is.finite(p_value)`. A code comment says the flag is "recomputed under
  the family-wide p", but `p_adj` does not enter the rule: `adjust` only
  changes the `p_adj` column, in both paths. (Reported to the lead as a
  comment/behaviour mismatch; the help page documents the behaviour.)
- **Pairwise**: each pair through `.bayes_compare_two()` with
  `adjust = "none"`; `p_adj` computed once across all pairs' finite
  `p_value`s. `global` is a data frame `group_a, group_b, n_credible`;
  for two groups `global` is a list of the same fields.
- **Stored settings**: `measure = "prob"`, `R = draws`, `prior`, `draws`,
  `ci`, `mean_threshold`, `bound_threshold`, `adjust`, `groups`, `fits`
  (two-group) or `comparisons` (pairwise).
- **Seed**: when `seed` is given, `.lsa_get_seed()` /
  `.lsa_restore_seed()` scope it to the call.
