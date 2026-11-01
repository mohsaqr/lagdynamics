# Lag Sequential Analysis Across Several Lags

Fits [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) at
each requested lag and returns the fits together, so you can compare a
transition's strength across lags (a *lag profile*). Each element is an
ordinary `lsa` fit.

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

  Sequence input (any form accepted by
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)).

- lags:

  Integer vector of lags. May include negative lags (predecessors) and
  `0`. Default `1:3`.

- ...:

  Passed to
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) (e.g.
  `engine`, `alpha`, `structural_zeros`). Ignored by the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods.

- x:

  An `lsa_lags` object ([`print()`](https://rdrr.io/r/base/print.html)
  and [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods).

- row.names:

  Optional row names for the stacked data frame.

- optional:

  Ignored; present for consistency with the generic.

## Value

An object of class `c("lsa_lags", "list")`: a named list of `lsa` fits
(names `"lag1"`, `"lag2"`, ...), with a `lags` attribute.
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
(with all its filtering arguments) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) row-bind
the per-lag transition tables (each carries its `lag` column) into one
tidy frame, one row per transition and lag, with the same columns as
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md).
[`print()`](https://rdrr.io/r/base/print.html) lists the number of
significant transitions per lag and returns `x` invisibly.

## See also

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)

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
transitions(prof, significant = TRUE, sort = "strength")
#>          from         to lag count  expected       prob   prob_col    adj_res
#> 1      Active     Active   1   459 247.03460 0.69756839 0.70506912  21.662811
#> 2  Disengaged Disengaged   1   157  59.97693 0.48307692 0.49062500  15.390466
#> 3     Average     Active   1   153 281.94983 0.20372836 0.23502304 -12.906038
#> 4      Active Disengaged   1    23 121.43022 0.03495441 0.07187500 -12.556925
#> 5     Average    Average   1   458 330.45732 0.60985353 0.60026212  12.452615
#> 6      Active    Average   1   176 289.53518 0.26747720 0.23066841 -11.319132
#> 7  Disengaged     Active   1    39 122.01557 0.12000000 0.05990783 -10.549486
#> 8      Active     Active   2   409 232.48908 0.66396104 0.67603306  18.697346
#> 9  Disengaged Disengaged   2   144  53.36868 0.49655172 0.48813559  15.175127
#> 10    Average    Average   2   420 305.67124 0.60258250 0.59743954  11.608507
#> 11    Average     Active   2   159 263.05989 0.22812052 0.26280992 -10.815838
#> 12     Active Disengaged   2    33 113.36245 0.05357143 0.11186441 -10.648529
#> 13     Active    Average   2   174 270.14847 0.28246753 0.24751067  -9.949353
#> 14 Disengaged     Active   2    37 109.45103 0.12758621 0.06115702  -9.697740
#> 15 Disengaged    Average   2   109 127.18029 0.37586207 0.15504979  -2.377227
#> 16     Active     Active   3   359 213.85753 0.63093146 0.64801444  16.032879
#> 17     Active Disengaged   3    35 102.68250 0.06151142 0.13157895  -9.416031
#> 18    Average     Active   3   156 238.66350 0.24566929 0.28158845  -8.977251
#> 19 Disengaged     Active   3    39 101.47897 0.14444444 0.07039711  -8.686327
#> 20 Disengaged Disengaged   3    98  48.72456 0.36296296 0.36842105   8.627944
#> 21     Active    Average   3   175 252.45997 0.30755712 0.26758410  -8.341561
#> 22    Average    Average   3   346 281.74355 0.54488189 0.52905199   6.802981
#> 23    Average Disengaged   3   133 114.59294 0.20944882 0.50000000   2.517610
#>                p    yules_q       kappa     kappa_z      kappa_p      lift
#> 1  4.602765e-104  0.8278779  0.44383295  19.1877564 4.684307e-82 1.8580393
#> 2   1.896612e-53  0.7544102  0.31355560  13.5556144 7.339504e-42 2.6176731
#> 3   4.161672e-38 -0.6010580 -0.49940424 -14.6962621 6.811736e-49 0.5426497
#> 4   3.642124e-36 -0.8264798 -0.82717329 -13.4120579 5.139269e-41 0.1894092
#> 5   1.353374e-35  0.5530510  0.22748944   9.8348080 7.971952e-23 1.3859581
#> 6   1.055115e-29 -0.5335259 -0.42415128 -12.4817491 9.390092e-36 0.6078709
#> 7   5.107575e-26 -0.6983917 -0.70694602 -11.4626536 2.031922e-30 0.3196313
#> 8   5.203079e-78  0.7771326  0.32986430  14.2606711 3.847367e-46 1.7592224
#> 9   5.167777e-52  0.7671697  0.26825723  11.5972787 4.253909e-31 2.6982116
#> 10  3.730729e-31  0.5389456  0.14479817   6.2599045 3.852132e-10 1.3740252
#> 11  2.896321e-27 -0.5327831 -0.47977303 -14.1185631 2.918823e-45 0.6044251
#> 12  1.771422e-26 -0.7291576 -0.75203125 -12.1936802 3.359020e-34 0.2911017
#> 13  2.538293e-23 -0.4916143 -0.43069502 -12.6743155 8.207202e-37 0.6440903
#> 14  3.082515e-22 -0.6781087 -0.72197443 -11.7063291 1.182910e-31 0.3380507
#> 15  1.744333e-02 -0.1567797 -0.29177004  -5.3457070 9.006481e-08 0.8570511
#> 16  7.531131e-58  0.7231612  0.21589564   9.3335858 1.023484e-20 1.6786877
#> 17  4.684613e-21 -0.6789253 -0.73700284 -11.9500047 6.492159e-33 0.3408565
#> 18  2.776155e-19 -0.4696536 -0.48958864 -14.4074126 4.647992e-47 0.6536400
#> 19  3.743475e-18 -0.6314902 -0.70694602 -11.4626536 2.031922e-30 0.3843161
#> 20  6.246447e-18  0.5568914  0.10797072   4.6677830 3.044672e-06 2.0113060
#> 21  7.331604e-17 -0.4336742 -0.42742315 -12.5780323 2.788902e-36 0.6931792
#> 22  1.024760e-11  0.3472710 -0.02110386  -0.7017496 4.828353e-01 1.2280671
#> 23  1.181541e-02  0.1688742 -0.13582949  -2.4886197 1.282400e-02 1.1606299
#>     sign significant
#> 1   over        TRUE
#> 2   over        TRUE
#> 3  under        TRUE
#> 4  under        TRUE
#> 5   over        TRUE
#> 6  under        TRUE
#> 7  under        TRUE
#> 8   over        TRUE
#> 9   over        TRUE
#> 10  over        TRUE
#> 11 under        TRUE
#> 12 under        TRUE
#> 13 under        TRUE
#> 14 under        TRUE
#> 15 under        TRUE
#> 16  over        TRUE
#> 17 under        TRUE
#> 18 under        TRUE
#> 19 under        TRUE
#> 20  over        TRUE
#> 21 under        TRUE
#> 22  over        TRUE
#> 23  over        TRUE
# Track one transition across lags with the dedicated verb:
lag_profile(engagement, from = "Active", to = "Average", lags = 1:3)
#>   lag   from      to count      prob    adj_res            p significant
#> 1   1 Active Average   176 0.2674772 -11.319132 1.055115e-29        TRUE
#> 2   2 Active Average   174 0.2824675  -9.949353 2.538293e-23        TRUE
#> 3   3 Active Average   175 0.3075571  -8.341561 7.331604e-17        TRUE
```
