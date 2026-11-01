# Tidy the Canonical Sequence Object

Returns the canonical `lsa_data` as a tidy data frame: one row per event
(`seq_id`, within-sequence `index`, `state`) for event-level input, or
one row per `from`/`to`/`count` cell for transition-matrix input.

## Usage

``` r
# S3 method for class 'lsa_data'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_data'
print(x, ...)
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

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): a
`data.frame` with one row per event and columns `seq_id`, `index`,
`state` (event-level input), or one row per cell and columns `from`,
`to`, `count` (transition-matrix input).
[`print()`](https://rdrr.io/r/base/print.html): `x`, invisibly.

## Details

[`print()`](https://rdrr.io/r/base/print.html) shows the input source,
the number of states, sequences, events and lag-1 transitions.

## Examples

``` r
d <- lsa_data(c("a", "b", "a", "c", "b"))
d
#> <lsa_data>
#>   source:        events
#>   states (K):    3
#>   labels:        a, b, c
#>   sequences:     1
#>   events:        5
#>   transitions:   4 (lag 1, summed across sequences)
as.data.frame(d)
#>   seq_id index state
#> 1      1     1     a
#> 2      1     2     b
#> 3      1     3     a
#> 4      1     4     c
#> 5      1     5     b
```
