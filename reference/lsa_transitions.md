# Count transitions at a given lag

`lsa_transitions()` counts the transitions at a given lag and returns
the table of observed transitions with its row and column totals: for
each pair of states, how often the second occurs `lag` events after the
first. This table is the starting point of every lag sequential
analysis. At lag 1 it records immediate succession; at larger lags it
records what tends to follow a state two or more events later, and at
negative lags what tends to precede it. Transitions are counted within
sequences only; no transition spans the boundary between two sequences.

The count involves no model or test;
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) adds
the expected counts, residuals and the other statistics. Use
`lsa_transitions()` to inspect the raw counts, or to check how many
transitions a lag leaves before fitting.

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

  An
  [`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md)
  object or any input it accepts.

- lag:

  Integer lag (default `1`). Positive lags count successors, negative
  lags predecessors. Input given as a transition-count matrix supports
  `lag = 1` only.

- row.names, optional, ...:

  Ignored; present for consistency with the generics.

## Value

An object of class `c("lsa_transitions", "list")`.
[`print()`](https://rdrr.io/r/base/print.html) shows the lag, the number
of states and transitions and the count matrix, and returns `x`
invisibly.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per `from -> to` pair with the columns `from`, `to`, `lag`,
`count`, `row_total`, `col_total` and `n_transitions`.

## Conditions

Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`:
a lag other than 1 for a transition-count matrix, or input that
[`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md)
rejects.

## See also

[`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md),
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
[`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)

## Examples

``` r
tx <- lsa_transitions(engagement)
tx
#> <lsa_transitions>
#>   lag:           1
#>   states (K):    3
#>   transitions:   1734
#>   obs:
#>            Active Average Disengaged
#> Active        459     176         23
#> Average       153     458        140
#> Disengaged     39     129        157
as.data.frame(tx)
#>         from         to lag count row_total col_total n_transitions
#> 1     Active     Active   1   459       658       651          1734
#> 2    Average     Active   1   153       751       651          1734
#> 3 Disengaged     Active   1    39       325       651          1734
#> 4     Active    Average   1   176       658       763          1734
#> 5    Average    Average   1   458       751       763          1734
#> 6 Disengaged    Average   1   129       325       763          1734
#> 7     Active Disengaged   1    23       658       320          1734
#> 8    Average Disengaged   1   140       751       320          1734
#> 9 Disengaged Disengaged   1   157       325       320          1734

# What precedes each state
lsa_transitions(engagement, lag = -1)
#> <lsa_transitions>
#>   lag:           -1
#>   states (K):    3
#>   transitions:   1734
#>   obs:
#>            Active Average Disengaged
#> Active        459     153         39
#> Average       176     458        129
#> Disengaged     23     140        157
```
