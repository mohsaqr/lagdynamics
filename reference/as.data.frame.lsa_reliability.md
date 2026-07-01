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
