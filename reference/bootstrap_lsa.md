# Estimate bootstrap confidence intervals for transitions

`bootstrap_lsa()` estimates bootstrap confidence intervals for the
transitions of a fit from
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md): it
resamples the sequences behind the fit and refits it with the same
engine, lag, alternative, significance level and structural zeros. The
non-parametric bootstrap (Efron, 1979) estimates how much the count,
transition probability, adjusted residual and Yule's Q of each
transition, all estimated from one sample of sequences, would vary over
repeated samples from the same population: it draws many resamples of
the data, refits the model on each, and uses the spread of the refitted
values as the sampling distribution. A transition whose residual
interval excludes zero is over- or under-represented in a way that
survives resampling, not only in the sample at hand.

With two or more sequences, whole sequences are drawn with replacement,
which keeps the dependence between events of one sequence intact. A
single sequence is resampled with the stationary block bootstrap
(Politis & Romano, 1994), which draws blocks of consecutive events of
random length. For every transition it returns the bootstrap mean,
standard error and percentile interval of the count, the transition
probability, the adjusted residual and Yule's Q, a two-sided bootstrap
p-value for the residual, and a flag for transitions whose residual
interval excludes zero.

For an interval on transition probabilities without resampling, use
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md);
it is exact under a Dirichlet-Multinomial model but ignores dependence
between the transitions of one sequence. For a p-value against the null
hypothesis of no sequential structure, use
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md);
for the share of subsamples in which a transition stays significant, use
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md).

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

  An `lsa` fit from
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
  built from sequences (not from a transition matrix).

- R:

  Number of bootstrap resamples. Default `1000`.

- level:

  Resampling unit: `"auto"` (default; `"sequence"` when the fit has two
  or more sequences, `"event"` otherwise), `"sequence"` (whole sequences
  with replacement) or `"event"` (stationary block bootstrap of the
  event stream). A single sequence always uses `"event"`, with a
  warning.

- block_length:

  Mean block length for `level = "event"`. Default `NULL`, the square
  root of the number of events, rounded up.

- level_alpha:

  Confidence level of the percentile intervals. Default `0.95`.

- indices:

  Optional matrix of resampling indices to replay, one row per resample:
  `R` x number of sequences (sequence indices) for `level = "sequence"`,
  or `R` x number of events (event positions) for `level = "event"`.
  Default `NULL` draws them at random.

- parallel:

  Run the resamples on several cores. Default `FALSE`.

- n_cores:

  Number of cores when `parallel = TRUE`. Default `NULL`, all cores but
  one.

- verbose:

  Report progress every 100 resamples. Default `FALSE`.

- ...:

  Not used.

## Value

An object of class `c("lsa_bootstrap", "list")`. Read it with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
which returns a data frame with one row per transition (`from -> to`
pair) and the columns:

- `from`, `to`: source and target state.

- `observed`: the observed count; `count_mean`, `count_se`,
  `count_ci_low`, `count_ci_high`: its bootstrap mean, standard error
  and interval.

- `adj_res_observed`, `adj_res_mean`, `adj_res_se`, `adj_res_ci_low`,
  `adj_res_ci_high`: the same for the adjusted residual;
  `adj_res_p_boot`: two-sided bootstrap p-value of the residual;
  `adj_res_stable`: `TRUE` when the residual interval excludes zero.

- `prob_observed`, `prob_mean`, `prob_ci_low`, `prob_ci_high`: the
  transition probability.

- `yules_q_observed`, `yules_q_mean`, `yules_q_ci_low`,
  `yules_q_ci_high`: Yule's Q.

`transitions(x, significant = TRUE)` keeps the transitions with
`adj_res_stable = TRUE`; `sort = "strength"` orders by the absolute
observed residual. [`print()`](https://rdrr.io/r/base/print.html) shows
the settings and the number of stable transitions (see
[lsa-inference-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md));
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)
draws the intervals. The object also stores the per-resample matrices,
the settings, the resampling indices used (which can be passed back as
`indices`) and the original fit.

## Conditions

Errors carry the class `lagdynamics_error` and one of:

- `lagdynamics_needs_events`: the fit was built from a transition
  matrix, so there are no sequences to resample.

- `lagdynamics_bad_input`: `indices` is not a matrix of whole numbers
  with at least `R` rows, the expected number of columns and values in
  range, or `block_length` is not a single number of at least 1.

## References

Efron, B. (1979). Bootstrap methods: Another look at the jackknife. *The
Annals of Statistics, 7*(1), 1-26.

Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
*Journal of the American Statistical Association, 89*(428), 1303-1313.

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
to read the result;
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)
to draw the intervals;
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
for analytic intervals;
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md)
and
[`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md)
for other checks of a fit;
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
to compare groups.

## Examples

``` r
fit <- lsa(engagement)
bs <- bootstrap_lsa(fit, R = 100)
bs
#> <lsa_bootstrap>
#>   engine:        classical
#>   level:         sequence
#>   replicates:    100
#>   CI level:      95%
#>   stable edges:  7 of 9
transitions(bs, significant = TRUE, sort = "strength")
#>         from         to observed count_mean  count_se count_ci_low
#> 1     Active     Active      459     460.81 44.461695      372.425
#> 2 Disengaged Disengaged      157     156.55 26.090122      111.425
#> 3    Average     Active      153     153.83 14.189949      127.375
#> 4     Active Disengaged       23      22.79  5.238349       14.000
#> 5    Average    Average      458     456.49 41.086874      389.425
#> 6     Active    Average      176     177.76 14.832479      148.950
#> 7 Disengaged     Active       39      39.35  6.733566       28.000
#>   count_ci_high adj_res_observed adj_res_mean adj_res_se adj_res_ci_low
#> 1       528.200         21.66281     21.49281   1.557722      19.181097
#> 2       210.000         15.39047     15.37902   1.727643      12.141444
#> 3       180.775        -12.90604    -12.83189   1.470068     -15.684310
#> 4        33.525        -12.55693    -12.56243   1.051631     -14.619661
#> 5       533.775         12.45261     12.38897   1.226563       9.832774
#> 6       206.000        -11.31913    -11.17702   1.618545     -13.551526
#> 7        53.525        -10.54949    -10.46740   1.195836     -12.468922
#>   adj_res_ci_high adj_res_p_boot adj_res_stable prob_observed  prob_mean
#> 1       24.574148              0           TRUE    0.69756839 0.69574828
#> 2       18.832847              0           TRUE    0.48307692 0.48204320
#> 3      -10.066614              0           TRUE    0.20372836 0.20598847
#> 4      -10.702942              0           TRUE    0.03495441 0.03451675
#> 5       14.521364              0           TRUE    0.60985353 0.60852965
#> 6       -7.502885              0           TRUE    0.26747720 0.26973498
#> 7       -8.105727              0           TRUE    0.12000000 0.12249664
#>   prob_ci_low prob_ci_high yules_q_observed yules_q_mean yules_q_ci_low
#> 1  0.64513014   0.74382975        0.8278779    0.8227804      0.7640340
#> 2  0.40384735   0.55598670        0.7544102    0.7533705      0.6678447
#> 3  0.16563725   0.24642036       -0.6010580   -0.5967991     -0.7043651
#> 4  0.02266560   0.05086872       -0.8264798   -0.8268602     -0.8948834
#> 5  0.56241815   0.64560328        0.5530510    0.5503139      0.4464577
#> 6  0.22476687   0.31584554       -0.5335259   -0.5258960     -0.6219965
#> 7  0.08693525   0.16456610       -0.6983917   -0.6927039     -0.7873195
#>   yules_q_ci_high
#> 1       0.8805156
#> 2       0.8184751
#> 3      -0.4893245
#> 4      -0.7482625
#> 5       0.6293966
#> 6      -0.3709085
#> 7      -0.5712618
# \donttest{
# A full analysis uses more resamples
bootstrap_lsa(fit, R = 1000) |>
  transitions(sort = "strength") |>
  head()
#>         from         to observed count_mean  count_se count_ci_low
#> 1     Active     Active      459    460.627 53.905395      356.975
#> 2 Disengaged Disengaged      157    155.143 28.226216      102.975
#> 3    Average     Active      153    152.771 15.077729      124.000
#> 4     Active Disengaged       23     23.153  5.345106       13.000
#> 5    Average    Average      458    457.798 39.026738      383.975
#> 6     Active    Average      176    175.905 16.191916      146.000
#>   count_ci_high adj_res_observed adj_res_mean adj_res_se adj_res_ci_low
#> 1       562.000         21.66281     21.59385   1.581340      18.390932
#> 2       210.000         15.39047     15.07333   1.958248      11.192064
#> 3       181.000        -12.90604    -12.92917   1.531189     -15.784908
#> 4        34.000        -12.55693    -12.47731   1.095469     -14.591847
#> 5       536.025         12.45261     12.37141   1.374770       9.665542
#> 6       207.000        -11.31913    -11.31600   1.663754     -14.366648
#>   adj_res_ci_high adj_res_p_boot adj_res_stable prob_observed  prob_mean
#> 1       24.510267              0           TRUE    0.69756839 0.69672987
#> 2       18.617808              0           TRUE    0.48307692 0.47631437
#> 3      -10.012508              0           TRUE    0.20372836 0.20395045
#> 4      -10.355861              0           TRUE    0.03495441 0.03527769
#> 5       15.073422              0           TRUE    0.60985353 0.60870161
#> 6       -8.119518              0           TRUE    0.26747720 0.26799245
#>   prob_ci_low prob_ci_high yules_q_observed yules_q_mean yules_q_ci_low
#> 1  0.63251813   0.75000000        0.8278779    0.8250796      0.7587794
#> 2  0.37927881   0.55946189        0.7544102    0.7441168      0.6391659
#> 3  0.16235907   0.24660255       -0.6010580   -0.6006064     -0.7010812
#> 4  0.01959984   0.05280939       -0.8264798   -0.8228775     -0.9036668
#> 5  0.56339644   0.65454602        0.5530510    0.5491329      0.4453994
#> 6  0.21732816   0.32775989       -0.5335259   -0.5318388     -0.6436483
#>   yules_q_ci_high
#> 1       0.8805664
#> 2       0.8284375
#> 3      -0.4888701
#> 4      -0.7277132
#> 5       0.6447136
#> 6      -0.4001190
# }
```
