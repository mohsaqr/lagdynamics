# Compare Groups' Transition Structures

Permutation test for whether groups produce different LSA transition
structures. For each pair of groups it pools their sequences, repeatedly
reassigns the group label of whole sequences (preserving the original
group sizes), refits each pseudo-group, and builds a permutation
distribution of the per-edge difference in the chosen `measure`. The
two-sided p-value is `(1 + #{ |diff_perm| >= |diff_obs| }) / (1 + R)`
(Phipson & Smyth 2010). A single omnibus test of overall difference is
reported from the same permutations.

## Usage

``` r
compare_lsa(
  x,
  y = NULL,
  R = 1000L,
  measure = c("log_or", "adj_res", "yules_q", "prob", "count", "lift"),
  adjust = "none",
  min_count = 5L,
  parallel = FALSE,
  n_cores = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  Either an `lsa_group` object (from `lsa(..., group = )`) with two or
  more groups, or a single `lsa` fit for the first group.

- y:

  When `x` is a single `lsa` fit, the second group's `lsa` fit. Ignored
  (and must be `NULL`) when `x` is an `lsa_group`.

- R:

  Integer. Number of label permutations per comparison. Default `1000`.

- measure:

  Character. The per-edge quantity compared between groups. Default
  `"log_or"`: the per-cell log odds ratio of the 2x2 transition collapse
  (Haldane-Anscombe corrected on empty cells) – an N-invariant LSA
  effect size, so the comparison reflects behaviour rather than sample
  size. Other options: `"yules_q"` (also N-invariant, but saturates at
  +/-1 on zero cells), `"adj_res"` (adjusted residuals – the LSA *test
  statistic*, which scales with sqrt(N) and is therefore confounded by
  group size; a message warns when groups differ in size), `"prob"`
  (transition probabilities – a raw rate, i.e. the TNA quantity, with no
  independence baseline), `"count"`, or `"lift"` (observed / expected).

- adjust:

  Multiple-comparison correction; any method accepted by
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) (e.g.
  `"holm"`, `"BH"`, `"bonferroni"`). Default `"none"`. For more than two
  groups it is applied across the pooled per-edge p-values of all pairs.

- min_count:

  Integer. Minimum pooled observed count (group a + group b) for a
  transition to be tested. Default `5`. Rarer cells carry an unstable
  odds ratio and a near-degenerate permutation null that produces
  spurious small p-values, so they get `p = NA` and are excluded from
  the multiple-comparison family and the omnibus rather than flagged
  significant. Set `0` to test every cell.

- parallel:

  Logical. Use multi-core resampling. Default `FALSE`.

- n_cores:

  Integer. Worker count when `parallel = TRUE`.

- verbose:

  Logical. Print progress every 100 permutations.

- ...:

  Reserved.

  **NA handling.** Non-estimable cells (structural zeros, zero-margin
  rows in a permuted pseudo-group) carry `NA` in the measure matrix and
  are never coerced to zero. Such cells get `p_perm = NA` rather than a
  spurious significant flag, and the exceedance tally and omnibus
  statistic are computed with `na.rm = TRUE`, matching
  [`permute_lsa()`](https://saqr.me/lagdynamics/reference/permute_lsa.md).

  **Interpretation caveats.** The odds ratio is non-collapsible: the
  per-group log odds ratios are group-specific departure-from-
  independence measures and should not be pooled across groups that have
  different marginal state distributions. As with any LSA, a
  between-group difference can also be driven by subgroup composition
  (Simpson's paradox); confirm with subgroup analysis when a confound is
  plausible.

## Value

For two groups, an object of class `c("lsa_comparison", "list")` with:

- edges:

  Tidy per-edge data frame: `from`, `to`, the measure in each group
  (`<measure>_a`, `<measure>_b`), their difference `diff` (= a - b), the
  permutation p-value `p_perm`, the adjusted p-value `p_adj`, and a
  `significant` flag.

- global:

  Omnibus test list: `statistic` (observed sum of squared edge
  differences), `p_value`, and `R`.

- perm_diff:

  `R x K^2` matrix of permuted edge differences.

- measure, R, adjust, groups:

  Call metadata; `groups` is the length-two character vector of group
  labels (a, b).

- fits:

  The two original fits, named by group.

For more than two groups, an object of class
`c("lsa_comparison_pairwise", "list")` with:

- edges:

  Tidy per-edge data frame across all pairs, prefixed by `group_a`,
  `group_b`; `p_adj` and `significant` reflect the family-wide
  correction.

- global:

  One row per pair: `group_a`, `group_b`, `statistic`, `p_value`, and
  the across-pairs adjusted `p_adj`.

- comparisons:

  Named list of the underlying two-group `lsa_comparison` objects (each
  fit with `adjust = "none"`), for drill-down.

- measure, R, adjust, groups:

  Call metadata; `groups` lists all group labels.

## Details

With exactly two groups a single comparison is returned. With more than
two groups every pairwise comparison is run and the requested `adjust`
correction is applied **once across the whole family** of per-edge
p-values (and separately across the per-pair omnibus tests), giving
family-wise control rather than per-pair control.

## References

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero. *Statistical Applications in Genetics and Molecular Biology*,
9(1), Article 39.

van Borkulo, C. D., et al. (2022). Comparing network structures on three
aspects: A permutation test. *Psychological Methods*.

## See also

[`permute_lsa()`](https://saqr.me/lagdynamics/reference/permute_lsa.md),
[`bootstrap_lsa()`](https://saqr.me/lagdynamics/reference/bootstrap_lsa.md)

## Examples

``` r
# \donttest{
# group_regulation is wide sequences with no grouping column, so
# derive one: sessions whose first regulation act is planning vs not.
grp <- ifelse(group_regulation$T1 == "plan", "starts_plan", "other")
g <- lsa(group_regulation, group = grp)
cmp <- compare_lsa(g, R = 200)
head(cmp$edges)
#>         from    to  log_or_a    log_or_b        diff    p_perm     p_adj
#> 1      adapt adapt -2.934854 -1.22410616 -1.71074827        NA        NA
#> 2   cohesion adapt -1.878259 -2.41906008  0.54080059 0.6019900 0.6019900
#> 3  consensus adapt -1.773051 -1.53109175 -0.24195921 0.6865672 0.6865672
#> 4 coregulate adapt -0.309429 -0.09747994 -0.21194909 0.6318408 0.6318408
#> 5    discuss adapt  1.882091  1.89742114 -0.01533041 0.9502488 0.9502488
#> 6    emotion adapt -2.247226 -2.33556543  0.08833958 1.0000000 1.0000000
#>   significant
#> 1       FALSE
#> 2       FALSE
#> 3       FALSE
#> 4       FALSE
#> 5       FALSE
#> 6       FALSE
cmp$global
#> $statistic
#> [1] 18.33942
#> 
#> $p_value
#> [1] 0.07462687
#> 
#> $R
#> [1] 200
#> 
# }
```
