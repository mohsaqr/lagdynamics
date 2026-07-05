# Tidy the Canonical Sequence Object

Returns the canonical `lsa_data` as a tidy data frame: one row per event
(`seq_id`, within-sequence `index`, `state`) for event-level input, or
one row per `from`/`to`/`count` cell for transition-matrix input.

## Usage

``` r
# S3 method for class 'lsa_data'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `lsa_data` object from
  [`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md).

- row.names, optional, ...:

  Standard
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  arguments (unused; present for method consistency).

## Value

A tidy data frame.

## Examples

``` r
d <- lsa_data(c("a", "b", "a", "c", "b"))
as.data.frame(d)
#>   seq_id index state
#> 1      1     1     a
#> 2      1     2     b
#> 3      1     3     a
#> 4      1     4     c
#> 5      1     5     b
```
