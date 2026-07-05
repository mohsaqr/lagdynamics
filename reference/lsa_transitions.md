# Tidy Transition Counts at a Given Lag

Computes the `K x K` transition count matrix from canonical lag
sequential data, optionally returning a tidy long-format edge table
alongside the matrix.

## Usage

``` r
lsa_transitions(x, lag = 1)
```

## Arguments

- x:

  Either an
  [lsa_data](https://saqr.me/lagdynamics/reference/lsa_data.md) object
  or any input accepted by
  [`lsa_data()`](https://saqr.me/lagdynamics/reference/lsa_data.md)
  (which will be coerced).

- lag:

  Integer. The lag at which to count transitions; default `1`. A
  positive lag counts successors (`from` at `t`, `to` at `t + lag`), a
  negative lag counts predecessors, and `0` pairs each event with itself
  (a degenerate diagonal). Must be a single finite whole number.

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

## Details

Transitions are counted within sequences only; no transition spans a
sequence boundary. For input that was supplied as a pre-computed
transition matrix (`source = "transitions"` on the `lsa_data` object),
the input matrix is returned at lag 1 and an error is raised for any
other lag.

## See also

[`lsa_data()`](https://saqr.me/lagdynamics/reference/lsa_data.md),
[`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md)

## Examples

``` r
d <- lsa_data(c("a", "b", "a", "c", "b", "a"))
tx <- lsa_transitions(d, lag = 1)
tx$obs
#>   a b c
#> a 0 1 1
#> b 2 0 0
#> c 0 1 0
head(tx$edges)
#>   from to lag count row_total col_total n_transitions
#> 1    a  a   1     0         2         2             5
#> 2    b  a   1     2         2         2             5
#> 3    c  a   1     0         1         2             5
#> 4    a  b   1     1         2         2             5
#> 5    b  b   1     0         2         2             5
#> 6    c  b   1     1         1         2             5
```
