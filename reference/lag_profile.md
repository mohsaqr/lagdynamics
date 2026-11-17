# Trace one transition across lags

`lag_profile()` traces one `from -> to` transition across lags, one row
per lag with the count, the transition probability, the adjusted
residual and its p-value. This lag profile shows how the association
between two states changes with the distance between them: whether a
transition that is over-represented at lag 1 persists, decays or
reverses at lags 2, 3 and beyond (Sackett, 1979). For example, it shows
whether disengagement is still unlikely several time points after an
active one, or only at the next.

`lag_profile()` fits the lags itself from raw data, or reads an existing
[`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)
result.

## Usage

``` r
lag_profile(x, from, to, lags = 1:3, ...)
```

## Arguments

- x:

  Sequences in any form accepted by
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or
  an
  [`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)
  result.

- from, to:

  State names of the transition.

- lags:

  Integer vector of lags (default `1:3`); ignored when `x` is an
  `lsa_lags` result.

- ...:

  Passed to
  [`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)
  when `x` is raw data.

## Value

A `data.frame` with one row per lag and the columns `lag`, `from`, `to`,
`count`, `prob`, `adj_res`, `p` and `significant`.

## Conditions

Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`
when `from` or `to` is not a state of the data.

## References

Sackett, G. P. (1979). The lag sequential analysis of contingency and
cyclicity in behavioral interaction research. In J. D. Osofsky (Ed.),
*Handbook of infant development* (pp. 623-649). Wiley.

## See also

[`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md),
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)

## Examples

``` r
lag_profile(engagement, from = "Active", to = "Disengaged", lags = 1:4)
#>   lag   from         to count       prob    adj_res            p significant
#> 1   1 Active Disengaged    23 0.03495441 -12.556925 3.642124e-36        TRUE
#> 2   2 Active Disengaged    33 0.05357143 -10.648529 1.771422e-26        TRUE
#> 3   3 Active Disengaged    35 0.06151142  -9.416031 4.684613e-21        TRUE
#> 4   4 Active Disengaged    39 0.07617188  -7.937923 2.055951e-15        TRUE
```
