# Read the transitions of a fit or an inference result

`transitions()` returns the transitions of a lag sequential analysis as
a table, one row per ordered pair of states, from a fit and from every
result built on a fit. The transitions are the substance of the
analysis: for each ordered pair of states, how often the second follows
the first at the chosen lag, how often independence would predict, and
whether the difference is larger than chance. A positive adjusted
residual marks an over-represented transition, one the process prefers;
a negative residual marks an avoided one. Yule's Q, kappa and lift
express the same association on bounded or ratio scales.

On a fit from
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), a
grouped fit or a multi-lag fit from
[`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md),
it returns the estimated statistics of each transition. On the result of
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
or
[`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md),
it returns that method's per-transition quantities. The arguments narrow
and order the table: `significant = TRUE` keeps the transitions the
method flags (significant for a fit, a permutation test or a group
comparison; stable under resampling for the bootstrap, certainty and
stability results), `direction` keeps over- or under-represented
transitions, `min_count` drops rare ones, and `sort` puts the strongest
first.

Use
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md) and
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
for the other tables of a fit, and
[`lag_profile()`](https://pak.dynasite.org/lagdynamics/reference/lag_profile.md)
to follow one transition across lags.

## Usage

``` r
transitions(fit, ...)

# S3 method for class 'lsa'
transitions(
  fit,
  significant = FALSE,
  direction = c("any", "over", "under"),
  min_count = NULL,
  alpha = NULL,
  sort = c("none", "strength", "count", "prob"),
  ...
)

# S3 method for class 'lsa_group'
transitions(
  fit,
  significant = FALSE,
  direction = c("any", "over", "under"),
  min_count = NULL,
  alpha = NULL,
  sort = c("none", "strength", "count", "prob"),
  ...
)

# S3 method for class 'lsa_bootstrap'
transitions(fit, significant = FALSE, sort = c("none", "strength"), ...)

# S3 method for class 'lsa_certainty'
transitions(fit, significant = FALSE, sort = c("none", "strength"), ...)

# S3 method for class 'lsa_stability'
transitions(fit, significant = FALSE, sort = c("none", "strength"), ...)

# S3 method for class 'lsa_permutation'
transitions(fit, significant = FALSE, sort = c("none", "strength"), ...)

# S3 method for class 'lsa_comparison'
transitions(fit, significant = FALSE, sort = c("none", "strength"), ...)

# S3 method for class 'lsa_comparison_pairwise'
transitions(fit, significant = FALSE, sort = c("none", "strength"), ...)

# S3 method for class 'lsa_lags'
transitions(
  fit,
  significant = FALSE,
  direction = c("any", "over", "under"),
  min_count = NULL,
  alpha = NULL,
  sort = c("none", "strength", "count", "prob"),
  ...
)
```

## Arguments

- fit:

  An `lsa` fit, a grouped fit, an `lsa_lags` result, or the result of
  [`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
  [`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md),
  [`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md),
  [`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
  [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
  or
  [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md).

- ...:

  Arguments passed to the method.

- significant:

  Keep only the transitions the method flags as significant or stable
  (default `FALSE`).

- direction:

  For fits only: `"any"` (default), `"over"` (significant with a
  positive residual) or `"under"` (significant with a negative
  residual). `"over"` and `"under"` imply `significant = TRUE`.

- min_count:

  For fits only: keep transitions observed at least this many times.
  Default `NULL` (no filter).

- alpha:

  For fits only: significance level. Default `NULL` uses the `alpha` of
  the fit.

- sort:

  Row order: `"none"` (default, table order) or `"strength"` (largest
  absolute effect first: the adjusted residual for a fit, a bootstrap,
  certainty or permutation result; the stability proportion for
  [`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md);
  the group difference for a comparison). Fits also accept `"count"` and
  `"prob"` (largest first).

## Value

A `data.frame` with one row per transition. For a fit the columns are
`from`, `to`, `lag`, `count`, `expected`, `prob` (row conditional,
\\P(\text{to} \mid \text{from})\\), `prob_col` (column conditional,
\\P(\text{from} \mid \text{to})\\), `adj_res`, `p`, `yules_q`, `kappa`,
`kappa_z`, `kappa_p`, `lift` (observed over expected), `sign` (`"over"`,
`"under"` or `"expected"`, comparing the count with the expected count)
and `significant`; engines that compute further statistics append them
(the two-cell engine adds `odds_ratio`, `log_or` and `log_or_se`). For
an inference result the columns are `from`, `to` and that method's
quantities, described on its help page; a pairwise comparison adds
`group_a` and `group_b`. A grouped fit gains a leading `group` column,
an `lsa_lags` result one row per transition and lag.

## See also

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md),
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md),
[`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md),
[`lag_profile()`](https://pak.dynasite.org/lagdynamics/reference/lag_profile.md)

## Examples

``` r
fit <- lsa(engagement)
transitions(fit)
#>         from         to lag count  expected       prob   prob_col     adj_res
#> 1     Active     Active   1   459 247.03460 0.69756839 0.70506912  21.6628107
#> 2    Average     Active   1   153 281.94983 0.20372836 0.23502304 -12.9060384
#> 3 Disengaged     Active   1    39 122.01557 0.12000000 0.05990783 -10.5494860
#> 4     Active    Average   1   176 289.53518 0.26747720 0.23066841 -11.3191319
#> 5    Average    Average   1   458 330.45732 0.60985353 0.60026212  12.4526148
#> 6 Disengaged    Average   1   129 143.00750 0.39692308 0.16906946  -1.7364613
#> 7     Active Disengaged   1    23 121.43022 0.03495441 0.07187500 -12.5569254
#> 8    Average Disengaged   1   140 138.59285 0.18641811 0.43750000   0.1757996
#> 9 Disengaged Disengaged   1   157  59.97693 0.48307692 0.49062500  15.3904661
#>               p     yules_q       kappa    kappa_z      kappa_p      lift  sign
#> 1 4.602765e-104  0.82787794  0.44383295  19.187756 4.684307e-82 1.8580393  over
#> 2  4.161672e-38 -0.60105804 -0.49940424 -14.696262 6.811736e-49 0.5426497 under
#> 3  5.107575e-26 -0.69839167 -0.70694602 -11.462654 2.031922e-30 0.3196313 under
#> 4  1.055115e-29 -0.53352593 -0.42415128 -12.481749 9.390092e-36 0.6078709 under
#> 5  1.353374e-35  0.55305097  0.22748944   9.834808 7.971952e-23 1.3859581  over
#> 6  8.248229e-02 -0.10831746 -0.16181958  -2.964801 3.028788e-03 0.9020506 under
#> 7  3.642124e-36 -0.82647982 -0.82717329 -13.412058 5.139269e-41 0.1894092 under
#> 8  8.604514e-01  0.01097122 -0.09034684  -1.655303 9.786312e-02 1.0101531  over
#> 9  1.896612e-53  0.75441020  0.31355560  13.555614 7.339504e-42 2.6176731  over
#>   significant
#> 1        TRUE
#> 2        TRUE
#> 3        TRUE
#> 4        TRUE
#> 5        TRUE
#> 6       FALSE
#> 7        TRUE
#> 8       FALSE
#> 9        TRUE
transitions(fit, significant = TRUE, sort = "strength")
#>         from         to lag count  expected       prob   prob_col   adj_res
#> 1     Active     Active   1   459 247.03460 0.69756839 0.70506912  21.66281
#> 2 Disengaged Disengaged   1   157  59.97693 0.48307692 0.49062500  15.39047
#> 3    Average     Active   1   153 281.94983 0.20372836 0.23502304 -12.90604
#> 4     Active Disengaged   1    23 121.43022 0.03495441 0.07187500 -12.55693
#> 5    Average    Average   1   458 330.45732 0.60985353 0.60026212  12.45261
#> 6     Active    Average   1   176 289.53518 0.26747720 0.23066841 -11.31913
#> 7 Disengaged     Active   1    39 122.01557 0.12000000 0.05990783 -10.54949
#>               p    yules_q      kappa    kappa_z      kappa_p      lift  sign
#> 1 4.602765e-104  0.8278779  0.4438330  19.187756 4.684307e-82 1.8580393  over
#> 2  1.896612e-53  0.7544102  0.3135556  13.555614 7.339504e-42 2.6176731  over
#> 3  4.161672e-38 -0.6010580 -0.4994042 -14.696262 6.811736e-49 0.5426497 under
#> 4  3.642124e-36 -0.8264798 -0.8271733 -13.412058 5.139269e-41 0.1894092 under
#> 5  1.353374e-35  0.5530510  0.2274894   9.834808 7.971952e-23 1.3859581  over
#> 6  1.055115e-29 -0.5335259 -0.4241513 -12.481749 9.390092e-36 0.6078709 under
#> 7  5.107575e-26 -0.6983917 -0.7069460 -11.462654 2.031922e-30 0.3196313 under
#>   significant
#> 1        TRUE
#> 2        TRUE
#> 3        TRUE
#> 4        TRUE
#> 5        TRUE
#> 6        TRUE
#> 7        TRUE
transitions(fit, direction = "under")
#>         from         to lag count expected       prob   prob_col   adj_res
#> 1    Average     Active   1   153 281.9498 0.20372836 0.23502304 -12.90604
#> 2 Disengaged     Active   1    39 122.0156 0.12000000 0.05990783 -10.54949
#> 3     Active    Average   1   176 289.5352 0.26747720 0.23066841 -11.31913
#> 4     Active Disengaged   1    23 121.4302 0.03495441 0.07187500 -12.55693
#>              p    yules_q      kappa   kappa_z      kappa_p      lift  sign
#> 1 4.161672e-38 -0.6010580 -0.4994042 -14.69626 6.811736e-49 0.5426497 under
#> 2 5.107575e-26 -0.6983917 -0.7069460 -11.46265 2.031922e-30 0.3196313 under
#> 3 1.055115e-29 -0.5335259 -0.4241513 -12.48175 9.390092e-36 0.6078709 under
#> 4 3.642124e-36 -0.8264798 -0.8271733 -13.41206 5.139269e-41 0.1894092 under
#>   significant
#> 1        TRUE
#> 2        TRUE
#> 3        TRUE
#> 4        TRUE
transitions(fit, min_count = 150, sort = "count")
#>         from         to lag count  expected      prob  prob_col   adj_res
#> 1     Active     Active   1   459 247.03460 0.6975684 0.7050691  21.66281
#> 2    Average    Average   1   458 330.45732 0.6098535 0.6002621  12.45261
#> 3     Active    Average   1   176 289.53518 0.2674772 0.2306684 -11.31913
#> 4 Disengaged Disengaged   1   157  59.97693 0.4830769 0.4906250  15.39047
#> 5    Average     Active   1   153 281.94983 0.2037284 0.2350230 -12.90604
#>               p    yules_q      kappa    kappa_z      kappa_p      lift  sign
#> 1 4.602765e-104  0.8278779  0.4438330  19.187756 4.684307e-82 1.8580393  over
#> 2  1.353374e-35  0.5530510  0.2274894   9.834808 7.971952e-23 1.3859581  over
#> 3  1.055115e-29 -0.5335259 -0.4241513 -12.481749 9.390092e-36 0.6078709 under
#> 4  1.896612e-53  0.7544102  0.3135556  13.555614 7.339504e-42 2.6176731  over
#> 5  4.161672e-38 -0.6010580 -0.4994042 -14.696262 6.811736e-49 0.5426497 under
#>   significant
#> 1        TRUE
#> 2        TRUE
#> 3        TRUE
#> 4        TRUE
#> 5        TRUE

# The same verb reads an inference result
transitions(certainty_lsa(fit), significant = TRUE)
#>         from         to observed prob_observed prob_mean    prob_se prob_ci_low
#> 1     Active     Active      459     0.6975684 0.6967400 0.01788572   0.6611341
#> 2    Average     Active      153     0.2037284 0.2039867 0.01467978   0.1759752
#> 3     Active    Average      176     0.2674772 0.2676270 0.01722641   0.2345470
#> 4    Average    Average      458     0.6098535 0.6093023 0.01777441   0.5742005
#> 5 Disengaged    Average      129     0.3969231 0.3966309 0.02703206   0.3442858
#> 6    Average Disengaged      140     0.1864181 0.1867110 0.01419599   0.1596920
#> 7 Disengaged Disengaged      157     0.4830769 0.4823890 0.02761180   0.4284082
#>   prob_ci_high      p_value stable adj_res_observed adj_res_stable
#> 1    0.7312156 5.560385e-20   TRUE       21.6628107           TRUE
#> 2    0.2334883 6.063196e-04   TRUE      -12.9060384           TRUE
#> 3    0.3020418 1.205186e-04   TRUE      -11.3191319           TRUE
#> 4    0.6438538 3.141676e-17   TRUE       12.4526148           TRUE
#> 5    0.4501757 2.232689e-04   TRUE       -1.7364613           TRUE
#> 6    0.2153070 1.168633e-03   TRUE        0.1757996           TRUE
#> 7    0.5365741 9.440665e-06   TRUE       15.3904661           TRUE
```
