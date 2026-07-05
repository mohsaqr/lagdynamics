# Tidy the per-replicate split-half correlations

Tidy the per-replicate split-half correlations

## Usage

``` r
# S3 method for class 'lsa_reliability'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_reliability_group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `lsa_reliability` (or `lsa_reliability_group`) object.

- row.names, optional, ...:

  Ignored (method signature compatibility).

## Value

A `data.frame`, one row per replicate, with columns `replicate` and
`correlation` (a grouped object gains a leading `group` column). `NA`
correlations from degenerate splits are kept.

## Examples

``` r
# \donttest{
fit <- lsa(group_regulation)
rel <- reliability_lsa(fit, R = 200)
head(as.data.frame(rel))
#>   replicate correlation
#> 1         1   0.9884478
#> 2         2   0.9931805
#> 3         3   0.9939218
#> 4         4   0.9916258
#> 5         5   0.9936414
#> 6         6   0.9922792
# }
```
