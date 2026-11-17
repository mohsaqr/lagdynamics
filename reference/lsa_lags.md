# Fit lag sequential analysis models at several lags

`lsa_lags()` fits
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) once
per requested lag, with the same engine and options, and keeps the fits
together. The result is a lag profile in the sense of Sackett (1979),
who proposed examining a transition across a range of lags to separate
short-range from longer-range dependence and to detect cycles: a
sequential contingency need not be immediate, a state may tend to recur
two or three events later, and an association that is strong at lag 1
may fade with distance.

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
on the result stacks the per-lag transition tables into one table with a
`lag` column, with the same filtering and sorting arguments as for a
single fit. To follow one transition across lags, use
[`lag_profile()`](https://pak.dynasite.org/lagdynamics/reference/lag_profile.md).

## Usage

``` r
lsa_lags(data, lags = 1:3, ...)

# S3 method for class 'lsa_lags'
print(x, ...)

# S3 method for class 'lsa_lags'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- data:

  Sequences in any form accepted by
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

- lags:

  Integer vector of lags (default `1:3`). Negative lags look back to
  predecessors.

- ...:

  Passed to
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) (for
  example `engine`, `alpha`, `loops`). Ignored by the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods.

- x:

  An `lsa_lags` object.

- row.names:

  Optional row names for the stacked data frame.

- optional:

  Ignored; present for consistency with the generic.

## Value

An object of class `c("lsa_lags", "list")`, one `lsa` fit per lag.
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
and [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
return a `data.frame` with one row per transition and lag and the
columns of
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
for a single fit. [`print()`](https://rdrr.io/r/base/print.html) lists
the number of significant transitions at each lag and returns `x`
invisibly.

## References

Sackett, G. P. (1979). The lag sequential analysis of contingency and
cyclicity in behavioral interaction research. In J. D. Osofsky (Ed.),
*Handbook of infant development* (pp. 623-649). Wiley.

## See also

[`lag_profile()`](https://pak.dynasite.org/lagdynamics/reference/lag_profile.md)
for one transition across lags;
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md);
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md).

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
transitions(prof, direction = "over", sort = "strength")
#>          from         to lag count  expected      prob  prob_col   adj_res
#> 1      Active     Active   1   459 247.03460 0.6975684 0.7050691 21.662811
#> 2  Disengaged Disengaged   1   157  59.97693 0.4830769 0.4906250 15.390466
#> 3     Average    Average   1   458 330.45732 0.6098535 0.6002621 12.452615
#> 4      Active     Active   2   409 232.48908 0.6639610 0.6760331 18.697346
#> 5  Disengaged Disengaged   2   144  53.36868 0.4965517 0.4881356 15.175127
#> 6     Average    Average   2   420 305.67124 0.6025825 0.5974395 11.608507
#> 7      Active     Active   3   359 213.85753 0.6309315 0.6480144 16.032879
#> 8  Disengaged Disengaged   3    98  48.72456 0.3629630 0.3684211  8.627944
#> 9     Average    Average   3   346 281.74355 0.5448819 0.5290520  6.802981
#> 10    Average Disengaged   3   133 114.59294 0.2094488 0.5000000  2.517610
#>                p   yules_q       kappa    kappa_z      kappa_p     lift sign
#> 1  4.602765e-104 0.8278779  0.44383295 19.1877564 4.684307e-82 1.858039 over
#> 2   1.896612e-53 0.7544102  0.31355560 13.5556144 7.339504e-42 2.617673 over
#> 3   1.353374e-35 0.5530510  0.22748944  9.8348080 7.971952e-23 1.385958 over
#> 4   5.203079e-78 0.7771326  0.32986430 14.2606711 3.847367e-46 1.759222 over
#> 5   5.167777e-52 0.7671697  0.26825723 11.5972787 4.253909e-31 2.698212 over
#> 6   3.730729e-31 0.5389456  0.14479817  6.2599045 3.852132e-10 1.374025 over
#> 7   7.531131e-58 0.7231612  0.21589564  9.3335858 1.023484e-20 1.678688 over
#> 8   6.246447e-18 0.5568914  0.10797072  4.6677830 3.044672e-06 2.011306 over
#> 9   1.024760e-11 0.3472710 -0.02110386 -0.7017496 4.828353e-01 1.228067 over
#> 10  1.181541e-02 0.1688742 -0.13582949 -2.4886197 1.282400e-02 1.160630 over
#>    significant
#> 1         TRUE
#> 2         TRUE
#> 3         TRUE
#> 4         TRUE
#> 5         TRUE
#> 6         TRUE
#> 7         TRUE
#> 8         TRUE
#> 9         TRUE
#> 10        TRUE
```
