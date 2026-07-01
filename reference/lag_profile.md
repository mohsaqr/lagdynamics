# Lag Profile of a Single Transition

How one `from -> to` transition behaves across lags, as a tidy
one-row-per-lag data frame. A clean shortcut for "track this transition
over lags 1, 2, 3, ...".

## Usage

``` r
lag_profile(x, from, to, lags = 1:3, ...)
```

## Arguments

- x:

  Sequence input (any form accepted by
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)) or
  an existing
  [`lsa_lags()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_lags.md)
  object.

- from, to:

  State labels of the transition to profile.

- lags:

  Integer vector of lags. Default `1:3`. Ignored when `x` is already an
  `lsa_lags` object.

- ...:

  Passed to
  [`lsa_lags()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_lags.md)
  when `x` is raw data.

## Value

A tidy `data.frame`, one row per lag, with columns `lag`, `from`, `to`,
`count`, `prob`, `adj_res`, `p`, and `significant`.

## See also

[`lsa_lags()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_lags.md)

## Examples

``` r
lag_profile(group_regulation, "plan", "consensus", lags = 1:3)
#>   lag from        to count      prob   adj_res            p significant
#> 1   1 plan consensus  1788 0.2904012  8.526657 1.506377e-17        TRUE
#> 2   2 plan consensus  1316 0.2323036 -3.619001 2.957420e-04        TRUE
#> 3   3 plan consensus  1261 0.2433893 -1.135550 2.561448e-01       FALSE
```
