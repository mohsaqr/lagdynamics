# Test transitions against no sequential structure by permutation

`permute_lsa()` tests every transition of a fit from
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) against
no sequential structure by permutation: whether the transition is over-
or under-represented beyond what the same events in a random order would
give, keeping the frequency of every state fixed. The adjusted residual
of a transition is otherwise tested against a normal reference
distribution, which holds only approximately, for large counts and
independent transitions. A permutation test builds the reference
distribution from the data instead: if the order of events carried no
information, shuffling the events would produce residuals as large as
the observed ones.

The events of the fit are shuffled many times, by default within each
sequence so that every sequence keeps its own length and state
frequencies; the model is refitted with the same engine and settings,
and the adjusted residual of every transition is recorded. The two-sided
permutation p-value of a transition is the share of shuffles whose
residual is at least as large in absolute value as the observed one,
counted so that it is never zero (Phipson & Smyth, 2010). A transition
is significant when its p-value is below the `alpha` of the fit.

The p-values are not corrected for the number of transitions tested. For
intervals on the residuals and other statistics, use
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md);
for the share of subsamples in which a transition stays significant, use
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md).
To test whether groups differ, use
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md).

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

  An `lsa` fit from
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
  built from sequences (not from a transition matrix).

- R:

  Number of permutations. Default `1000`.

- within_sequence:

  Shuffle events within each sequence (`TRUE`, default) or across the
  whole event stream (`FALSE`). A single sequence is always shuffled as
  a whole, with a warning when `within_sequence = TRUE`.

- shuffles:

  Optional list of at least `R` permutations of `seq_len(n_events)` to
  replay instead of random shuffles. Default `NULL`.

- parallel:

  Run the permutations on several cores. Default `FALSE`.

- n_cores:

  Number of cores when `parallel = TRUE`. Default `NULL`, all cores but
  one.

- verbose:

  Report progress every 100 permutations. Default `FALSE`.

- ...:

  Not used.

## Value

An object of class `c("lsa_permutation", "list")`. Read it with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
which returns a data frame with one row per transition (`from -> to`
pair) and the columns `from`, `to`, `observed_count`,
`observed_adj_res`, `p_perm` (permutation p-value; `NA` for transitions
that cannot be estimated, such as structural zeros) and `significant`
(`p_perm` below the fit's `alpha`). `transitions(x, significant = TRUE)`
keeps the significant transitions; `sort = "strength"` orders by the
absolute observed residual.
[`print()`](https://rdrr.io/r/base/print.html) shows the settings and
the number of significant transitions (see
[lsa-inference-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)).
The object also stores the permuted residuals, the settings and the
original fit.

## Conditions

Errors carry the class `lagdynamics_error` and one of:

- `lagdynamics_needs_events`: the fit was built from a transition
  matrix, so there are no events to shuffle.

- `lagdynamics_bad_input`: `shuffles` is not a list of at least `R`
  permutations of `seq_len(n_events)`.

## References

Castellan, N. J. (1992). Shuffling arrays: Appearances may be deceiving.
*Behavior Research Methods, Instruments, & Computers, 24*(1), 72-77.

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: Calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology, 9*(1),
Article 39.

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
to read the result;
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md)
and
[`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md)
for other checks of a fit;
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
to compare groups.

## Examples

``` r
fit <- lsa(engagement)
pm <- permute_lsa(fit, R = 100)
pm
#> <lsa_permutation>
#>   engine:        classical
#>   replicates:    100
#>   within seq:    TRUE
#>   significant edges (p_perm < 0.050): 7 of 9
transitions(pm, significant = TRUE, sort = "strength")
#>         from         to observed_count observed_adj_res     p_perm significant
#> 1     Active     Active            459        21.662811 0.00990099        TRUE
#> 2 Disengaged Disengaged            157        15.390466 0.00990099        TRUE
#> 3    Average     Active            153       -12.906038 0.00990099        TRUE
#> 4     Active Disengaged             23       -12.556925 0.00990099        TRUE
#> 5    Average    Average            458        12.452615 0.00990099        TRUE
#> 6     Active    Average            176       -11.319132 0.02970297        TRUE
#> 7 Disengaged    Average            129        -1.736461 0.00990099        TRUE
# \donttest{
# A full analysis uses more permutations
permute_lsa(fit, R = 1000) |>
  transitions(sort = "strength")
#>         from         to observed_count observed_adj_res      p_perm significant
#> 1     Active     Active            459       21.6628107 0.000999001        TRUE
#> 2 Disengaged Disengaged            157       15.3904661 0.000999001        TRUE
#> 3    Average     Active            153      -12.9060384 0.000999001        TRUE
#> 4     Active Disengaged             23      -12.5569254 0.000999001        TRUE
#> 5    Average    Average            458       12.4526148 0.000999001        TRUE
#> 6     Active    Average            176      -11.3191319 0.006993007        TRUE
#> 7 Disengaged     Active             39      -10.5494860 0.068931069       FALSE
#> 8 Disengaged    Average            129       -1.7364613 0.022977023        TRUE
#> 9    Average Disengaged            140        0.1757996 0.828171828       FALSE
# }
```
