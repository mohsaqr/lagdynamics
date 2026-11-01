# Tidy Transition Counts at a Given Lag

Computes the `K x K` transition count matrix from canonical lag
sequential data, optionally returning a tidy long-format edge table
alongside the matrix.

## Usage

``` r
lsa_transitions(x, lag = 1)

# S3 method for class 'lsa_transitions'
print(x, ...)

# S3 method for class 'lsa_transitions'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  Either an
  [lsa_data](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md)
  object or any input accepted by
  [`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md)
  (which will be coerced).

- lag:

  Integer. The lag at which to count transitions; default `1`. A
  positive lag counts successors (`from` at `t`, `to` at `t + lag`), a
  negative lag counts predecessors, and `0` pairs each event with itself
  (a degenerate diagonal). Must be a single finite whole number.

- row.names, optional, ...:

  Ignored by the [`print()`](https://rdrr.io/r/base/print.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods; present for consistency with the generics.

## Value

An object of class `c("lsa_transitions", "list")` with elements:

- obs:

  The `K x K` observed transition count matrix with `dimnames` set to
  the labels.

- row_totals:

  Length-`K` vector `rowSums(obs)`.

- col_totals:

  Length-`K` vector `colSums(obs)`.

- n_transitions:

  Scalar `sum(obs)`.

- lag:

  The lag used.

- labels:

  Character vector of state labels.

- edges:

  Tidy long-format data.frame with one row per `(from, to)` cell
  containing columns `from`, `to`, `lag`, `count`, `row_total`,
  `col_total`, `n_transitions`.

[`print()`](https://rdrr.io/r/base/print.html) shows the lag, the number
of states and transitions and the count matrix, and returns `x`
invisibly.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the `edges` table: one row per `from -> to` cell with the columns listed
above.

## Details

Transitions are counted within sequences only; no transition spans a
sequence boundary. For input that was supplied as a pre-computed
transition matrix (`source = "transitions"` on the `lsa_data` object),
the input matrix is returned at lag 1 and an error is raised for any
other lag.

## See also

[`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md),
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)

## Examples

``` r
d <- lsa_data(c("a", "b", "a", "c", "b", "a"))
tx <- lsa_transitions(d, lag = 1)
tx
#> <lsa_transitions>
#>   lag:           1
#>   states (K):    3
#>   transitions:   5
#>   obs:
#>   a b c
#> a 0 1 1
#> b 2 0 0
#> c 0 1 0
as.data.frame(tx)
#>   from to lag count row_total col_total n_transitions
#> 1    a  a   1     0         2         2             5
#> 2    b  a   1     2         2         2             5
#> 3    c  a   1     0         1         2             5
#> 4    a  b   1     1         2         2             5
#> 5    b  b   1     0         2         2             5
#> 6    c  b   1     1         1         2             5
#> 7    a  c   1     1         2         1             5
#> 8    b  c   1     0         2         1             5
#> 9    c  c   1     0         1         1             5
```
