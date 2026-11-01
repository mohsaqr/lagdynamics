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

Read the per-edge table with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md);
[`print()`](https://rdrr.io/r/base/print.html) shows the group labels,
the number of significant edges and the omnibus test.

## Details

With exactly two groups a single comparison is returned. With more than
two groups every pairwise comparison is run and the requested `adjust`
correction is applied **once across the whole family** of per-edge
p-values (and separately across the per-pair omnibus tests), giving
family-wise control rather than per-pair control.

**NA handling.** Non-estimable cells (structural zeros, zero-margin rows
in a permuted pseudo-group) carry `NA` in the measure matrix and are
never coerced to zero. Such cells get `p_perm = NA` rather than a
spurious significant flag, and the exceedance tally and omnibus
statistic are computed with `na.rm = TRUE`, matching
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md).

**Interpretation caveats.** The odds ratio is non-collapsible: the
per-group log odds ratios are group-specific departure-from-
independence measures and should not be pooled across groups that have
different marginal state distributions. As with any LSA, a between-group
difference can also be driven by subgroup composition (Simpson's
paradox); confirm with subgroup analysis when a confound is plausible.

## References

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology*, 9(1),
Article 39.

van Borkulo, C. D., et al. (2022). Comparing network structures on three
aspects: A permutation test. *Psychological Methods*.

## See also

[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md)

## Examples

``` r
g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
         time = "Time", group = "Achiever")
cmp <- compare_lsa(g, R = 100)
cmp
#> <lsa_comparison>
#>   groups:   High vs Low
#>   measure:  log_or difference (High - Low)
#>   R:        100 label permutations
#>   edges:    43 significant of 78 tested (adjust = none)
#>   omnibus:  statistic = 79.48, p = 0.009901
transitions(cmp, significant = TRUE, sort = "strength")
#>          from         to    log_or_a    log_or_b       diff     p_perm
#> 1    cohesion      adapt -0.79358879 -3.92090141  3.1273126 0.00990099
#> 2     emotion  synthesis -1.61865035 -4.52153508  2.9028847 0.00990099
#> 3        plan  synthesis -2.14522408 -4.95050449  2.8052804 0.00990099
#> 4    cohesion  synthesis -1.35530113 -3.94457485  2.5892737 0.00990099
#> 5   synthesis    monitor -3.22853627 -0.97455001 -2.2539863 0.02970297
#> 6   synthesis       plan -0.58056730 -2.68043829  2.0998710 0.00990099
#> 7       adapt coregulate -3.14559229 -1.17601724 -1.9695750 0.00990099
#> 8        plan      adapt -2.33764937 -4.23038186  1.8927325 0.02970297
#> 9  coregulate      adapt  0.77815349 -1.07661947  1.8547730 0.00990099
#> 10   cohesion   cohesion -0.58689097 -2.33812956  1.7512386 0.00990099
#> 11  synthesis coregulate -1.73380230 -0.31919923 -1.4146031 0.00990099
#> 12    discuss      adapt  1.00995539  2.25071658 -1.2407612 0.00990099
#> 13  synthesis    discuss -1.85303029 -0.62037619 -1.2326541 0.00990099
#> 14       plan coregulate -1.39829837 -2.43430237  1.0360040 0.00990099
#> 15    monitor      adapt -0.01744120 -1.05200843  1.0345672 0.04950495
#> 16    discuss  synthesis  2.76148241  3.77841337 -1.0169310 0.00990099
#> 17   cohesion    monitor -0.98452606 -0.01973429 -0.9647918 0.00990099
#> 18    discuss  consensus  0.84183283 -0.09611980  0.9379526 0.00990099
#> 19  consensus    discuss  0.72233489 -0.14721343  0.8695483 0.00990099
#> 20   cohesion    discuss -1.54566436 -0.70245954 -0.8432048 0.00990099
#> 21    emotion    emotion -0.77996171 -0.05035498 -0.7296067 0.00990099
#> 22   cohesion coregulate  0.09080815  0.77561247 -0.6848043 0.00990099
#> 23 coregulate  consensus -1.17441323 -0.50890395 -0.6655093 0.00990099
#> 24    emotion coregulate -1.32586326 -0.73147590 -0.5943874 0.04950495
#> 25  consensus   cohesion -1.62877101 -2.18856416  0.5597932 0.01980198
#> 26  synthesis  consensus  1.33223346  0.77429226  0.5579412 0.00990099
#> 27 coregulate    monitor  0.98297844  0.44335452  0.5396239 0.01980198
#> 28    discuss   cohesion -0.21681786 -0.73390412  0.5170863 0.00990099
#> 29    emotion    discuss -0.32240670 -0.80671989  0.4843132 0.00990099
#> 30 coregulate    discuss  0.55245810  1.02535891 -0.4729008 0.00990099
#> 31 coregulate       plan  0.21564991 -0.24687495  0.4625249 0.00990099
#> 32    discuss    discuss  0.11514755  0.57379923 -0.4586517 0.00990099
#> 33       plan    emotion  0.68789133  0.25378782  0.4341035 0.00990099
#> 34  consensus    monitor -0.24860635  0.14315978 -0.3917661 0.00990099
#> 35       plan    monitor  0.88978802  0.51554999  0.3742380 0.00990099
#> 36   cohesion       plan -0.54489516 -0.90097330  0.3560781 0.02970297
#> 37       plan       plan  0.66801871  1.02175484 -0.3537361 0.00990099
#> 38       plan   cohesion -1.06300840 -1.40300799  0.3399996 0.03960396
#> 39       plan    discuss -1.26146322 -0.97661999 -0.2848432 0.00990099
#> 40    monitor       plan -0.01892834 -0.29879500  0.2798667 0.03960396
#> 41       plan  consensus  0.15190086  0.42448522 -0.2725844 0.00990099
#> 42  consensus  consensus -1.72572997 -1.45713398 -0.2685960 0.00990099
#> 43  consensus       plan  0.96085946  1.07782451 -0.1169651 0.01980198
#>         p_adj significant
#> 1  0.00990099        TRUE
#> 2  0.00990099        TRUE
#> 3  0.00990099        TRUE
#> 4  0.00990099        TRUE
#> 5  0.02970297        TRUE
#> 6  0.00990099        TRUE
#> 7  0.00990099        TRUE
#> 8  0.02970297        TRUE
#> 9  0.00990099        TRUE
#> 10 0.00990099        TRUE
#> 11 0.00990099        TRUE
#> 12 0.00990099        TRUE
#> 13 0.00990099        TRUE
#> 14 0.00990099        TRUE
#> 15 0.04950495        TRUE
#> 16 0.00990099        TRUE
#> 17 0.00990099        TRUE
#> 18 0.00990099        TRUE
#> 19 0.00990099        TRUE
#> 20 0.00990099        TRUE
#> 21 0.00990099        TRUE
#> 22 0.00990099        TRUE
#> 23 0.00990099        TRUE
#> 24 0.04950495        TRUE
#> 25 0.01980198        TRUE
#> 26 0.00990099        TRUE
#> 27 0.01980198        TRUE
#> 28 0.00990099        TRUE
#> 29 0.00990099        TRUE
#> 30 0.00990099        TRUE
#> 31 0.00990099        TRUE
#> 32 0.00990099        TRUE
#> 33 0.00990099        TRUE
#> 34 0.00990099        TRUE
#> 35 0.00990099        TRUE
#> 36 0.02970297        TRUE
#> 37 0.00990099        TRUE
#> 38 0.03960396        TRUE
#> 39 0.00990099        TRUE
#> 40 0.03960396        TRUE
#> 41 0.00990099        TRUE
#> 42 0.00990099        TRUE
#> 43 0.01980198        TRUE
# \donttest{
# A full analysis uses more label permutations
compare_lsa(g, R = 1000, adjust = "BH")
#> <lsa_comparison>
#>   groups:   High vs Low
#>   measure:  log_or difference (High - Low)
#>   R:        1000 label permutations
#>   edges:    38 significant of 78 tested (adjust = BH)
#>   omnibus:  statistic = 79.48, p = 0.000999
# }
```
