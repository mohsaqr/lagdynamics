# Assess the case-drop stability of significant transitions

`stability_lsa()` assesses the case-drop stability of the significant
transitions of a fit from
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md):
whether a transition's significance survives when part of the data is
left out, since a transition found significant in the full data may owe
its significance to a few sequences. It follows the case-dropping
approach that Epskamp, Borsboom and Fried (2018) describe for the
stability of psychological networks: subsets of the cases are drawn
without replacement, the model is refitted on each, and the share of
subsets in which a transition is significant measures how much its
significance depends on the particular cases sampled.

`R` subsets are drawn, each keeping a share `proportion` of the
sequences (or, for a single sequence, a run of consecutive events of
that share of its length); the model is refitted with the same engine
and settings, and whether each transition is significant at the fit's
`alpha` is recorded. A transition's stability is the share of subsets in
which it was significant; it is flagged stable when that share reaches
`min_stable`. Replicates in which a subset contained no transitions are
reported as failed and count as not significant.

Stability concerns the significance decision, in either direction, not
the size of the effect. For intervals on the statistics use
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md);
for a p-value against no sequential structure use
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md);
for the agreement of the whole network between halves of the data use
[`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md).

## Usage

``` r
stability_lsa(
  fit,
  R = 500L,
  proportion = 0.8,
  min_stable = 0.95,
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

  Number of subsets. Default `500`.

- proportion:

  Share of the sequences kept in each subset, between 0 and 1. Default
  `0.8`.

- min_stable:

  Share of subsets in which a transition must be significant to be
  flagged stable, greater than 0 and at most 1. Default `0.95`.

- parallel:

  Run the subsets on several cores. Default `FALSE`.

- n_cores:

  Number of cores when `parallel = TRUE`. Default `NULL`, all cores but
  one.

- verbose:

  Report progress every 100 subsets. Default `FALSE`.

- ...:

  Not used.

## Value

An object of class `c("lsa_stability", "list")`. Read it with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
which returns a data frame with one row per transition (`from -> to`
pair) and the columns `from`, `to`, `observed_sig` (significant in the
full fit), `stability` (share of subsets in which it was significant)
and `stable` (`stability >= min_stable`).
`transitions(x, significant = TRUE)` keeps the stable transitions;
`sort = "strength"` orders by stability.
[`print()`](https://rdrr.io/r/base/print.html) shows the settings, the
number of failed replicates and the number of stable transitions (see
[lsa-inference-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)).
The object also stores the per-subset significance decisions, the
settings and the original fit.

## Conditions

Errors carry the class `lagdynamics_error` and one of:

- `lagdynamics_needs_events`: the fit was built from a transition
  matrix, so there are no sequences to drop.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods, 50*(1), 195-212.

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
to read the result;
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)
and
[`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md)
for other checks of a fit.

## Examples

``` r
fit <- lsa(engagement)
st <- stability_lsa(fit, R = 100)
st
#> <lsa_stability>
#>   engine:        classical
#>   replicates:    100
#>   proportion:    80%
#>   min stable:    95%
#>   failed:        0 of 100 replicates
#>   stable edges:  7 of 9 (>= 95% across replicates)
transitions(st, significant = TRUE)
#>         from         to observed_sig stability stable
#> 1     Active     Active         TRUE         1   TRUE
#> 2    Average     Active         TRUE         1   TRUE
#> 3 Disengaged     Active         TRUE         1   TRUE
#> 4     Active    Average         TRUE         1   TRUE
#> 5    Average    Average         TRUE         1   TRUE
#> 6     Active Disengaged         TRUE         1   TRUE
#> 7 Disengaged Disengaged         TRUE         1   TRUE
# \donttest{
# A full analysis uses more subsets
stability_lsa(fit, R = 500, proportion = 0.5) |>
  transitions(sort = "strength")
#>         from         to observed_sig stability stable
#> 1     Active     Active         TRUE     1.000   TRUE
#> 2    Average     Active         TRUE     1.000   TRUE
#> 3 Disengaged     Active         TRUE     1.000   TRUE
#> 4     Active    Average         TRUE     1.000   TRUE
#> 5    Average    Average         TRUE     1.000   TRUE
#> 6     Active Disengaged         TRUE     1.000   TRUE
#> 7 Disengaged Disengaged         TRUE     1.000   TRUE
#> 8 Disengaged    Average        FALSE     0.248  FALSE
#> 9    Average Disengaged        FALSE     0.058  FALSE
# }
```
