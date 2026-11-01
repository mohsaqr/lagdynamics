# Permutation Test for an LSA Fit

Empirical null-distribution p-values for every cell of an LSA transition
matrix. Repeatedly shuffles the input event vector (within sequence
boundaries) and recomputes the engine's residual matrix, producing a
permutation distribution for each cell. The two-sided p-value is
`(1 + #{ |stat_perm| >= |stat_obs| }) / (1 + R)` (Phipson & Smyth 2010).

## Usage

``` r
permute_lsa(
  fit,
  R = 1000L,
  within_sequence = TRUE,
  shuffles = NULL,
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

  Integer. Number of permutations. Default `1000`.

- within_sequence:

  Logical. When `TRUE` (default for multi-sequence input), each sequence
  is shuffled independently; when `FALSE`, the whole event stream is
  shuffled across sequence boundaries.

- shuffles:

  Optional list of length `R`, each element an integer permutation of
  `seq_len(n_events)`. When supplied, replaces internal RNG.

- parallel:

  Logical. Use multi-core resampling. Default `FALSE`.

- n_cores:

  Integer. Worker count when `parallel = TRUE`.

- verbose:

  Logical. Print progress every 100 replicates.

- ...:

  Reserved.

## Value

An object of class `c("lsa_permutation", "list")`. Read the per-edge
table with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md);
see
[lsa-inference-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)
for [`print()`](https://rdrr.io/r/base/print.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).
Elements:

- edges:

  Tidy data frame, one row per `from -> to` cell: `from`, `to`,
  `observed_count`, `observed_adj_res`, the empirical permutation
  p-value `p_perm`, and `significant` (`p_perm < alpha` of the fit).

- perm_adj_res:

  `R x K^2` numeric matrix of permuted residuals (cells in
  `as.vector(adj_res)` order).

- R, within_sequence:

  Recipe metadata.

- fit:

  Reference to the original fit.

## Details

**NA handling.** The exceedance count that drives `p_perm` is computed
with `na.rm = TRUE`, so replicates that produced `NA` for a cell
(structural-zero cells, zero-margin cells in the permuted table) are
excluded from that cell's tally rather than counted as either an
exceedance or a non-exceedance.

## References

Castellan, N. J. (1992). Shuffling arrays: appearances may be deceiving.
*Behavior Research Methods, Instruments, & Computers*, 24(1), 72-77.

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology*, 9(1),
Article 39.

## See also

[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md)

## Examples

``` r
fit <- lsa(engagement, engine = "classical")
pm <- permute_lsa(fit, R = 100)
pm
#> <lsa_permutation>
#>   engine:        classical
#>   replicates:    100
#>   within seq:    TRUE
#>   significant edges (p_perm < 0.050): 7 of 9
transitions(pm, significant = TRUE)
#>         from         to observed_count observed_adj_res     p_perm significant
#> 1     Active     Active            459        21.662811 0.00990099        TRUE
#> 2    Average     Active            153       -12.906038 0.00990099        TRUE
#> 3     Active    Average            176       -11.319132 0.02970297        TRUE
#> 4    Average    Average            458        12.452615 0.00990099        TRUE
#> 5 Disengaged    Average            129        -1.736461 0.01980198        TRUE
#> 6     Active Disengaged             23       -12.556925 0.00990099        TRUE
#> 7 Disengaged Disengaged            157        15.390466 0.00990099        TRUE
# \donttest{
# A full analysis uses more permutations
transitions(permute_lsa(fit, R = 1000), sort = "strength")
#>         from         to observed_count observed_adj_res      p_perm significant
#> 1     Active     Active            459       21.6628107 0.000999001        TRUE
#> 2 Disengaged Disengaged            157       15.3904661 0.000999001        TRUE
#> 3    Average     Active            153      -12.9060384 0.000999001        TRUE
#> 4     Active Disengaged             23      -12.5569254 0.000999001        TRUE
#> 5    Average    Average            458       12.4526148 0.000999001        TRUE
#> 6     Active    Average            176      -11.3191319 0.013986014        TRUE
#> 7 Disengaged     Active             39      -10.5494860 0.085914086       FALSE
#> 8 Disengaged    Average            129       -1.7364613 0.026973027        TRUE
#> 9    Average Disengaged            140        0.1757996 0.820179820       FALSE
# }
```
