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
  [`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md).

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

  **NA handling.** The exceedance count that drives `p_perm` is computed
  with `na.rm = TRUE`, so replicates that produced `NA` for a cell
  (structural-zero cells, zero-margin cells in the permuted table) are
  excluded from that cell's tally rather than counted as either an
  exceedance or a non-exceedance.

## Value

An object of class `c("lsa_permutation", "list")` with:

- edges:

  Tidy per-edge data frame with observed counts and residuals, the
  empirical permutation p-value `p_perm`, and a `significant` flag at
  the recipe's alpha threshold.

- perm_adj_res:

  `R x K^2` numeric matrix of permuted residuals (cells in
  `as.vector(adj_res)` order).

- R, within_sequence:

  Recipe metadata.

- fit:

  Reference to the original fit.

## References

Castellan, N. J. (1992). Shuffling arrays: appearances may be deceiving.
*Behavior Research Methods, Instruments, & Computers*, 24(1), 72-77.

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology*, 9(1),
Article 39.

## See also

[`bootstrap_lsa()`](https://saqr.me/lagdynamics/reference/bootstrap_lsa.md),
[`stability_lsa()`](https://saqr.me/lagdynamics/reference/stability_lsa.md)

## Examples

``` r
# \donttest{
fit <- lsa(engagement, engine = "classical")
pm <- permute_lsa(fit, R = 200)
head(pm$edges)
#>         from      to observed_count observed_adj_res      p_perm significant
#> 1     Active  Active            459        21.662811 0.004975124        TRUE
#> 2    Average  Active            153       -12.906038 0.004975124        TRUE
#> 3 Disengaged  Active             39       -10.549486 0.094527363       FALSE
#> 4     Active Average            176       -11.319132 0.014925373        TRUE
#> 5    Average Average            458        12.452615 0.004975124        TRUE
#> 6 Disengaged Average            129        -1.736461 0.034825871        TRUE
# }
```
