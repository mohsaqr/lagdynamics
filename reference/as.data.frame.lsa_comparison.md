# Tidy a Group Comparison

Returns the per-edge comparison table, the same table
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
reads, as a plain data frame.
[`print()`](https://rdrr.io/r/base/print.html) shows the header summary
(groups, measure, significant edges, omnibus test).

## Usage

``` r
# S3 method for class 'lsa_comparison'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_comparison_pairwise'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_comparison'
print(x, ...)

# S3 method for class 'lsa_comparison_pairwise'
print(x, ...)
```

## Arguments

- x:

  An `lsa_comparison` or `lsa_comparison_pairwise` object.

- row.names, optional, ...:

  Standard
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  arguments (unused; present for method consistency).

## Value

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): the
per-edge data frame described under *Value* in
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
(one row per transition; for a pairwise object one row per transition
and pair, with leading `group_a`, `group_b`).
[`print()`](https://rdrr.io/r/base/print.html): `x`, invisibly.

## Examples

``` r
g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
         time = "Time", group = "Achiever")
cmp <- compare_lsa(g, R = 50)
cmp
#> <lsa_comparison>
#>   groups:   High vs Low
#>   measure:  log_or difference (High - Low)
#>   R:        50 label permutations
#>   edges:    40 significant of 78 tested (adjust = none)
#>   omnibus:  statistic = 79.48, p = 0.01961
head(as.data.frame(cmp))
#>         from    to   log_or_a  log_or_b       diff     p_perm      p_adj
#> 1      adapt adapt -1.1833106 -3.167522  1.9842113         NA         NA
#> 2   cohesion adapt -0.7935888 -3.920901  3.1273126 0.01960784 0.01960784
#> 3  consensus adapt -1.2188354 -1.966624  0.7477882 0.09803922 0.09803922
#> 4 coregulate adapt  0.7781535 -1.076619  1.8547730 0.01960784 0.01960784
#> 5    discuss adapt  1.0099554  2.250717 -1.2407612 0.01960784 0.01960784
#> 6    emotion adapt -1.3519489 -3.102996  1.7510475 0.23529412 0.23529412
#>   significant
#> 1       FALSE
#> 2        TRUE
#> 3       FALSE
#> 4        TRUE
#> 5        TRUE
#> 6       FALSE
```
