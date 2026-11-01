# Bootstrap Confidence Intervals for an LSA Fit

Non-parametric bootstrap for any LSA fit produced by
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).
Resamples the underlying sequence data (whole sequences when more than
one is available; geometric-block stationary bootstrap on events
otherwise), refits the engine on each resample using the immutable
recipe stored in `fit$params`, and aggregates per-edge statistics into a
tidy data frame.

## Usage

``` r
bootstrap_lsa(
  fit,
  R = 1000L,
  level = c("auto", "sequence", "event"),
  block_length = NULL,
  level_alpha = 0.95,
  indices = NULL,
  parallel = FALSE,
  n_cores = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- fit:

  An `lsa` object returned by
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

- R:

  Integer. Number of bootstrap replicates. Default `1000`.

- level:

  Character. Resampling unit: `"auto"` (default) picks `"sequence"` when
  the fit has more than one sequence and `"event"` otherwise;
  `"sequence"` resamples whole sequences with replacement; `"event"`
  runs a stationary block bootstrap on the event stream.

- block_length:

  For event-level bootstrap, mean geometric block length. Default `NULL`
  -\> `ceiling(sqrt(T))`.

- level_alpha:

  Numeric. Confidence level for percentile intervals. Default `0.95`.

- indices:

  Optional integer matrix of replay indices, one row per resample (row
  `b` is used for resample `b`). For sequence-level bootstrap it must be
  `R x S` with each entry a sequence index in `1..S`. For event-level
  bootstrap it is `R x T` with each entry an event position in `1..T`
  (the fully expanded positions, not block starts). When supplied,
  replaces internal RNG and enables bit-identical reproducibility across
  runs. Dimensions and ranges are validated. See Details.

- parallel:

  Logical. Use multi-core resampling. Default `FALSE`. Requires base R
  only (`parallel` package).

- n_cores:

  Integer. Worker count when `parallel = TRUE`. Default `NULL` -\>
  `parallel::detectCores() - 1`.

- verbose:

  Logical. Print progress every 100 replicates. Default `FALSE`.

- ...:

  Reserved for future use.

## Value

An object of class `c("lsa_bootstrap", "list")`. Read the per-edge table
with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md);
see
[lsa-inference-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)
for [`print()`](https://rdrr.io/r/base/print.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).
Elements:

- edges:

  Tidy data frame, one row per `from -> to` cell: `from`, `to`,
  `observed` (observed count); `count_mean`, `count_se`, `count_ci_low`,
  `count_ci_high`; `adj_res_observed`, `adj_res_mean`, `adj_res_se`,
  `adj_res_ci_low`, `adj_res_ci_high`, the two-sided bootstrap p-value
  `adj_res_p_boot`, and `adj_res_stable` (`TRUE` when the residual
  interval excludes zero); `prob_observed`, `prob_mean`, `prob_ci_low`,
  `prob_ci_high`; and `yules_q_observed`, `yules_q_mean`,
  `yules_q_ci_low`, `yules_q_ci_high`.

- boot_obs, boot_adj_res, boot_prob, boot_yulesq:

  `R x K^2` numeric matrices of the per-replicate counts, adjusted
  residuals, transition probabilities and Yule's Q (cells in
  `as.vector(obs)` order).

- R, level, level_alpha, indices_used:

  Recipe metadata; `indices_used` is the replay matrix accepted by
  `indices`.

- fit:

  The original fit.

## Details

**Sequence-level resampling (default for multi-sequence input).** Each
resample draws `S` sequence indices with replacement from `seq_len(S)`
and rebuilds the multi-sequence input as the corresponding list of event
vectors. Preserves within-sequence structure.

**Event-level resampling (single-sequence input).** Implements the
stationary block bootstrap of Politis & Romano (1994). Block length is
geometric with mean `block_length`; resampled blocks wrap around the
event stream and are concatenated until total length equals the original
`T`.

**Reproducibility hook.** Supply `indices` as an `R x S` integer matrix
of sequence indices (sequence-level) or an `R x T` matrix of event
positions (event-level) to deterministically replay the bootstrap across
sessions, processes, or languages. The event-level matrix holds the
fully expanded resampled positions, i.e. the same form produced
internally, so a captured `indices_used` can be fed straight back in.

**NA handling.** Per-cell summary statistics (`mean`, `se`, `ci_low`,
`ci_high`, `p_boot`) are computed with `na.rm = TRUE`, so replicates
that produced `NA` for a given cell (for example structural-zero cells,
or cells whose row marginal collapsed to zero in the resampled data) are
excluded from that cell's summary. The summary therefore reflects only
the finite replicates; cells whose every replicate was `NA` come back as
`NA` themselves.

## References

Efron, B. (1979). Bootstrap methods: another look at the jackknife.
*Annals of Statistics*, 7(1), 1-26.

Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
*Journal of the American Statistical Association*, 89(428), 1303-1313.

## See also

[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md)

## Examples

``` r
fit <- lsa(engagement, engine = "classical")
bs <- bootstrap_lsa(fit, R = 100)
bs
#> <lsa_bootstrap>
#>   engine:        classical
#>   level:         sequence
#>   replicates:    100
#>   CI level:      95%
#>   stable edges:  7 of 9
transitions(bs, significant = TRUE)
#>         from         to observed count_mean count_se count_ci_low count_ci_high
#> 1     Active     Active      459     471.68 56.69106      361.375       577.075
#> 2    Average     Active      153     152.52 16.08279      125.950       187.100
#> 3 Disengaged     Active       39      38.21  6.15670       27.000        49.000
#> 4     Active    Average      176     175.02 17.54013      148.475       212.525
#> 5    Average    Average      458     452.10 42.07605      380.475       543.775
#> 6     Active Disengaged       23      23.23  5.18127       14.000        33.000
#> 7 Disengaged Disengaged      157     159.24 28.15653      106.850       209.675
#>   adj_res_observed adj_res_mean adj_res_se adj_res_ci_low adj_res_ci_high
#> 1         21.66281     21.82386   1.497096       19.22208       24.481665
#> 2        -12.90604    -12.95315   1.418268      -15.55977      -10.181116
#> 3        -10.54949    -10.78628   1.211023      -12.77849       -8.086579
#> 4        -11.31913    -11.46380   1.532969      -14.34037       -8.733901
#> 5         12.45261     12.67136   1.273073       10.40148       14.742772
#> 6        -12.55693    -12.67236   1.033695      -14.50446      -10.666906
#> 7         15.39047     15.69518   2.013464       11.21309       19.008666
#>   adj_res_p_boot adj_res_stable prob_observed prob_mean prob_ci_low
#> 1              0           TRUE    0.69756839 0.7025653  0.64724584
#> 2              0           TRUE    0.20372836 0.2063853  0.16721114
#> 3              0           TRUE    0.12000000 0.1189835  0.08077802
#> 4              0           TRUE    0.26747720 0.2625375  0.21975016
#> 5              0           TRUE    0.60985353 0.6092515  0.56502846
#> 6              0           TRUE    0.03495441 0.0348972  0.02141810
#> 7              0           TRUE    0.48307692 0.4881223  0.38415666
#>   prob_ci_high yules_q_observed yules_q_mean yules_q_ci_low yules_q_ci_high
#> 1    0.7531284        0.8278779    0.8292671      0.7702192       0.8769173
#> 2    0.2556696       -0.6010580   -0.6018424     -0.6912605      -0.4917702
#> 3    0.1658154       -0.6983917   -0.7078359     -0.8071431      -0.5827928
#> 4    0.3200289       -0.5335259   -0.5379198     -0.6443172      -0.4288457
#> 5    0.6501866        0.5530510    0.5611605      0.4687221       0.6368029
#> 6    0.0520404       -0.8264798   -0.8273677     -0.8898637      -0.7596529
#> 7    0.5679351        0.7544102    0.7602738      0.6499925       0.8369990
# \donttest{
# A full analysis uses more replicates
transitions(bootstrap_lsa(fit, R = 1000), sort = "strength")
#>         from         to observed count_mean  count_se count_ci_low
#> 1     Active     Active      459    461.943 50.815332      368.000
#> 2 Disengaged Disengaged      157    156.757 27.544012      103.000
#> 3    Average     Active      153    153.544 14.470503      125.000
#> 4     Active Disengaged       23     22.735  5.211135       13.975
#> 5    Average    Average      458    457.412 39.761807      377.975
#> 6     Active    Average      176    177.250 15.588602      147.975
#> 7 Disengaged     Active       39     39.171  5.988792       28.000
#> 8 Disengaged    Average      129    128.166 13.926823      102.000
#> 9    Average Disengaged      140    139.568 14.569064      111.000
#>   count_ci_high adj_res_observed adj_res_mean adj_res_se adj_res_ci_low
#> 1       561.000       21.6628107   21.5737454   1.504985      18.632735
#> 2       215.000       15.3904661   15.3313800   1.962549      11.375243
#> 3       182.000      -12.9060384  -12.8743075   1.466178     -15.679731
#> 4        33.025      -12.5569254  -12.5716615   1.064613     -14.627083
#> 5       537.050       12.4526148   12.4036141   1.392297       9.759276
#> 6       209.000      -11.3191319  -11.2380812   1.602112     -14.357223
#> 7        51.000      -10.5494860  -10.4982963   1.099916     -12.658793
#> 8       156.000       -1.7364613   -1.7303818   1.449272      -4.492032
#> 9       169.000        0.1757996    0.2468326   1.442574      -2.416885
#>   adj_res_ci_high adj_res_p_boot adj_res_stable prob_observed  prob_mean
#> 1       24.404272          0.000           TRUE    0.69756839 0.69652158
#> 2       18.897963          0.000           TRUE    0.48307692 0.48086637
#> 3       -9.822118          0.000           TRUE    0.20372836 0.20517119
#> 4      -10.469285          0.000           TRUE    0.03495441 0.03451026
#> 5       15.112067          0.000           TRUE    0.60985353 0.60872363
#> 6       -8.017311          0.000           TRUE    0.26747720 0.26896816
#> 7       -8.393433          0.000           TRUE    0.12000000 0.12181236
#> 8        1.067600          0.232          FALSE    0.39692308 0.39732127
#> 9        2.973862          0.858          FALSE    0.18641811 0.18610518
#>   prob_ci_low prob_ci_high yules_q_observed yules_q_mean yules_q_ci_low
#> 1  0.63862833   0.74620366       0.82787794   0.82441161      0.7610255
#> 2  0.38519187   0.56403636       0.75441020   0.75072766      0.6456307
#> 3  0.16451416   0.24819182      -0.60105804  -0.59828201     -0.6993383
#> 4  0.02026974   0.05128738      -0.82647982  -0.82725466     -0.9021004
#> 5  0.56210618   0.65518845       0.55305097   0.55012505      0.4455137
#> 6  0.22143278   0.32460753      -0.53352593  -0.52832161     -0.6452188
#> 7  0.08308331   0.16247847      -0.69839167  -0.69468403     -0.7904546
#> 8  0.32073274   0.48152428      -0.10831746  -0.10670202     -0.2705137
#> 9  0.15425792   0.22256878       0.01097122   0.01668867     -0.1434465
#>   yules_q_ci_high
#> 1      0.87900558
#> 2      0.83256817
#> 3     -0.48163552
#> 4     -0.73651347
#> 5      0.64775636
#> 6     -0.39390622
#> 7     -0.58396174
#> 8      0.07048928
#> 9      0.19375740
# }
```
