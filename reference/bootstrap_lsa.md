# Bootstrap Confidence Intervals for an LSA Fit

Non-parametric bootstrap for any LSA fit produced by
[`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md).
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
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md).

- R:

  Integer. Number of bootstrap replicates. Default `1000`.

- level:

  Character. Resampling unit: `"sequence"` (resample whole sequences
  with replacement, default when fit has more than one sequence),
  `"event"` (stationary block bootstrap on the event stream, used
  automatically for single-sequence input), or `"auto"` (default; pick
  based on `fit$data$n_sequences`).

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

An object of class `c("lsa_bootstrap", "list")` with:

- edges:

  Tidy per-edge data frame with observed + bootstrap `mean`, `se`,
  `ci_low`, `ci_high`, `p_boot`, and `stable` for `count`, `adj_res`,
  `prob`, and `yules_q`.

- boot_obs:

  `R x K^2` numeric matrix: cell-wise observed count from each replicate
  (flattened in `as.vector(obs)` order).

- boot_adj_res:

  `R x K^2` matrix of adjusted residuals.

- R, level, level_alpha, indices_used:

  Recipe metadata.

- fit:

  Reference to the original fit (for \$params / labels).

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

[`permute_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/permute_lsa.md),
[`stability_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/stability_lsa.md)

## Examples

``` r
# \donttest{
fit <- lsa(engagement, engine = "classical")
bs <- bootstrap_lsa(fit, R = 200)
head(bs$edges)
#>         from      to observed count_mean  count_se count_ci_low count_ci_high
#> 1     Active  Active      459    469.300 52.308997      376.975       569.150
#> 2    Average  Active      153    152.005 16.286844      122.000       186.025
#> 3 Disengaged  Active       39     38.315  6.084479       27.975        49.000
#> 4     Active Average      176    175.075 17.506621      141.875       213.025
#> 5    Average Average      458    450.865 42.901509      369.875       538.050
#> 6 Disengaged Average      129    127.660 14.242847      101.975       154.025
#>   adj_res_observed adj_res_mean adj_res_se adj_res_ci_low adj_res_ci_high
#> 1        21.662811     21.83015   1.521667      18.806855      24.6240543
#> 2       -12.906038    -12.94027   1.496602     -15.530175     -10.0943842
#> 3       -10.549486    -10.77165   1.218591     -12.757864      -8.3070225
#> 4       -11.319132    -11.37662   1.572194     -14.382033      -8.0782759
#> 5        12.452615     12.52827   1.402022       9.788089      15.0397178
#> 6        -1.736461     -1.66657   1.468575      -4.454881       0.9783829
#>   adj_res_p_boot adj_res_stable prob_observed prob_mean prob_ci_low
#> 1           0.00           TRUE     0.6975684 0.7023435  0.65231852
#> 2           0.00           TRUE     0.2037284 0.2057721  0.16464788
#> 3           0.00           TRUE     0.1200000 0.1188032  0.08244323
#> 4           0.00           TRUE     0.2674772 0.2634867  0.20960424
#> 5           0.00           TRUE     0.6098535 0.6073409  0.54974454
#> 6           0.25          FALSE     0.3969231 0.3943909  0.31606167
#>   prob_ci_high yules_q_observed yules_q_mean yules_q_ci_low yules_q_ci_high
#> 1    0.7555365        0.8278779    0.8295468      0.7604949      0.87912442
#> 2    0.2549351       -0.6010580   -0.6015555     -0.6952226     -0.49358512
#> 3    0.1605894       -0.6983917   -0.7068677     -0.8023661     -0.58578503
#> 4    0.3158819       -0.5335259   -0.5347029     -0.6436904     -0.40695051
#> 5    0.6525239        0.5530510    0.5557749      0.4484348      0.64459557
#> 6    0.4722248       -0.1083175   -0.1028465     -0.2764733      0.06399258
# }
```
