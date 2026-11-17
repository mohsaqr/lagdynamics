# Read the tables of a group comparison

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the tables of a group comparison from
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
or
[`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
as a plain data frame, and
[`print()`](https://rdrr.io/r/base/print.html) shows a short summary of
the comparison. A comparison holds two tables, chosen with `what`: the
comparison of every transition and the omnibus result for the structure
as a whole. The per-transition table is the one
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
reads; use
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
when you want to keep only the significant transitions or order them by
the size of the difference.

## Usage

``` r
# S3 method for class 'lsa_comparison'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("edges", "global"),
  ...
)

# S3 method for class 'lsa_comparison_pairwise'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("edges", "global"),
  ...
)

# S3 method for class 'lsa_comparison'
print(x, ...)

# S3 method for class 'lsa_comparison_pairwise'
print(x, ...)
```

## Arguments

- x:

  A comparison from
  [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
  or
  [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
  (class `lsa_comparison` or `lsa_comparison_pairwise`).

- row.names, optional:

  Not used; present for consistency with the generic.

- what:

  `"edges"` (default): one row per transition. `"global"`: the omnibus
  result, one row per pair of groups.

- ...:

  Not used.

## Value

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
data frame.

- `what = "edges"`: one row per transition (and per pair of groups when
  more than two groups were compared, with leading `group_a`, `group_b`
  columns); the columns are listed under *Value* in
  [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
  and
  [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md).

- `what = "global"`: one row per pair of groups with `group_a`,
  `group_b` and, for
  [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md),
  the omnibus `statistic` (sum of squared differences), its permutation
  `p_value` and either the number of permutations `R` (two groups) or
  the across-pairs adjusted `p_adj` (more than two groups); for
  [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md),
  `n_credible`, the number of credibly different transitions.

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## See also

[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md),
[`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md),
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md).

## Examples

``` r
g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
         time = "Time", group = "Achiever")
cmp <- compare_lsa(g, R = 50, adjust = "BH")
cmp
#> <lsa_comparison>
#>   groups:   High vs Low
#>   measure:  log_or difference (High - Low)
#>   R:        50 label permutations
#>   edges:    37 significant of 78 tested (adjust = BH)
#>   omnibus:  statistic = 79.48, p = 0.01961
as.data.frame(cmp, what = "global")
#>   group_a group_b statistic    p_value  R
#> 1    High     Low  79.47723 0.01960784 50
head(as.data.frame(cmp))
#>         from    to   log_or_a  log_or_b       diff     p_perm      p_adj
#> 1      adapt adapt -1.1833106 -3.167522  1.9842113         NA         NA
#> 2   cohesion adapt -0.7935888 -3.920901  3.1273126 0.01960784 0.04133545
#> 3  consensus adapt -1.2188354 -1.966624  0.7477882 0.09803922 0.15294118
#> 4 coregulate adapt  0.7781535 -1.076619  1.8547730 0.01960784 0.04133545
#> 5    discuss adapt  1.0099554  2.250717 -1.2407612 0.01960784 0.04133545
#> 6    emotion adapt -1.3519489 -3.102996  1.7510475 0.23529412 0.32198142
#>   significant
#> 1       FALSE
#> 2        TRUE
#> 3       FALSE
#> 4        TRUE
#> 5        TRUE
#> 6       FALSE
```
