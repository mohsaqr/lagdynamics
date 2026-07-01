# Lag Sequential Analysis Across Several Lags

Fits [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)
at each requested lag and returns the fits together, so you can compare
a transition's strength across lags (a *lag profile*). Each element is
an ordinary `lsa` fit.

## Usage

``` r
lsa_lags(data, lags = 1:3, ...)
```

## Arguments

- data:

  Sequence input (any form accepted by
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)).

- lags:

  Integer vector of lags. May include negative lags (predecessors) and
  `0`. Default `1:3`.

- ...:

  Passed to
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)
  (e.g. `engine`, `alpha`, `structural_zeros`).

## Value

An object of class `c("lsa_lags", "list")`: a named list of `lsa` fits
(names `"lag1"`, `"lag2"`, ...), with a `lags` attribute.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on it
row-binds
[`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md)
of every fit (each already carries its `lag` column) into one tidy long
frame with the same columns as
[`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md).

## See also

[`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)

## Examples

``` r
prof <- lsa_lags(engagement, lags = 1:3)
prof
#> <lsa_lags>
#>   engine: classical
#>   lags:   1, 2, 3
#>     lag 1    7 of 9 transitions significant
#>     lag 2    8 of 9 transitions significant
#>     lag 3    8 of 9 transitions significant
# Track one transition across lags with the dedicated verb:
lag_profile(engagement, from = "Active", to = "Average", lags = 1:3)
#>   lag   from      to count      prob    adj_res            p significant
#> 1   1 Active Average   176 0.2674772 -11.319132 1.055115e-29        TRUE
#> 2   2 Active Average   174 0.2824675  -9.949353 2.538293e-23        TRUE
#> 3   3 Active Average   175 0.3075571  -8.341561 7.331604e-17        TRUE
```
