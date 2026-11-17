# Read the tables of a split-half reliability result

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the tables of a
[`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md)
result, which holds the split-half correlation of every random split:
these correlations, one row per split, or a one-row summary of their
distribution, chosen with `what`. For a grouped result the tables of the
groups are stacked.

## Usage

``` r
# S3 method for class 'lsa_reliability'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("replicates", "summary"),
  ...
)

# S3 method for class 'lsa_reliability_group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `lsa_reliability` or `lsa_reliability_group` object from
  [`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md).

- row.names, optional:

  Not used; present for consistency with the generic.

- what:

  `"replicates"` (default): one row per split. `"summary"`: one row
  summarizing all splits.

- ...:

  For a grouped result, passed on to each group (use it for `what`);
  otherwise not used.

## Value

A data frame.

- `what = "replicates"`: one row per split, with `replicate` and
  `correlation`; a failed replicate (a half with no transitions, or no
  variation to correlate) is kept as `NA`.

- `what = "summary"`: one row with `weights`, `method`, `R`, `n_failed`,
  `mean`, `sd`, `ci_low` and `ci_high` (2.5% and 97.5% quantiles of the
  finite correlations).

A grouped result gains a leading `group` column.

## See also

[`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md).

## Examples

``` r
rel <- reliability_lsa(lsa(engagement), R = 20)
as.data.frame(rel) |> head()
#>   replicate correlation
#> 1         1   0.9813133
#> 2         2   0.9683861
#> 3         3   0.9885036
#> 4         4   0.9873595
#> 5         5   0.9729948
#> 6         6   0.9741100
as.data.frame(rel, what = "summary")
#>   weights  method  R n_failed      mean         sd    ci_low   ci_high
#> 1    prob pearson 20        0 0.9651927 0.03308002 0.8891692 0.9939467

g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
         time = "Time", group = "Achiever")
grel <- reliability_lsa(g, R = 10)
as.data.frame(grel, what = "summary")
#>   group weights  method  R n_failed      mean          sd    ci_low   ci_high
#> 1  High    prob pearson 10        0 0.9880429 0.002878287 0.9842305 0.9924178
#> 2   Low    prob pearson 10        0 0.9861624 0.004595129 0.9775524 0.9919682
```
