# Read the events of canonical sequence data as a table

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the canonical sequence data from
[`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md),
every event with its sequence and position, as a tidy table, which shows
how the input was read: which rows became sequences, where missing
values were dropped, and which labels were assigned. For input given as
a transition-count matrix, it returns one row per cell instead.
[`print()`](https://rdrr.io/r/base/print.html) shows the input source,
the number of states, sequences, events and lag-1 transitions.

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

  Ignored; present for consistency with the generic.

## Value

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): a
`data.frame` with one row per event and the columns `seq_id`, `index`
(position within the sequence) and `state`; or, for transition-matrix
input, one row per cell and the columns `from`, `to` and `count`.
[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## See also

[`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md),
[`lsa_transitions()`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md)

## Examples

``` r
d <- lsa_data(list(c("a", "b", NA, "a"), c("b", "c")))
d
#> <lsa_data>
#>   source:        events
#>   states (K):    3
#>   labels:        a, b, c
#>   sequences:     2
#>   events:        5
#>   transitions:   3 (lag 1, summed across sequences)
as.data.frame(d)
#>   seq_id index state
#> 1      1     1     a
#> 2      1     2     b
#> 3      1     3     a
#> 4      2     1     b
#> 5      2     2     c
```
