# Bayesian Comparison of Group Transition Structures (Dirichlet-Multinomial)

Closed-form Bayesian alternative to
[`compare_lsa()`](https://saqr.me/lagdynamics/reference/compare_lsa.md)
for comparing the transition structures of two (or, pairwise, more)
groups. Each state's outgoing transitions are modelled as
Dirichlet-Multinomial with a Jeffreys prior, so each transition
probability is marginally Beta. The per-edge posterior mean difference
`prob_a - prob_b` is exact; a credible interval, the probability of
direction `pd`, and a two-sided Bayesian p-equivalent `2 * (1 - pd)`
come from a Monte Carlo draw on the Beta marginals.

## Usage

``` r
bayes_compare_lsa(
  x,
  y = NULL,
  prior = 0.5,
  draws = 10000L,
  ci = 0.95,
  mean_threshold = 0.01,
  bound_threshold = 0.001,
  adjust = "none",
  seed = NULL
)
```

## Arguments

- x:

  An `lsa_group` (two or more groups), or a single `lsa` fit for the
  first group.

- y:

  The second group's `lsa` fit when `x` is a single fit; otherwise
  `NULL`.

- prior:

  Numeric \> 0. Dirichlet prior concentration added to every cell.
  Default `0.5` (Jeffreys). Use `1` for a uniform (Laplace) prior.

- draws:

  Integer. Monte Carlo posterior draws for the credible intervals.
  Default `10000`.

- ci:

  Numeric in (0, 1). Credible-interval mass. Default `0.95`.

- mean_threshold, bound_threshold:

  An edge is flagged credibly different only if its credible interval
  excludes zero, `|posterior mean diff|` exceeds `mean_threshold`
  (default `0.01`), and the credible bound nearest zero exceeds
  `bound_threshold` (default `0.001`). The thresholds guard against
  differences that are detectable but negligibly small.

- adjust:

  Multiple-comparison correction applied to the two-sided Bayesian p
  across edges (and family-wide across pairs); any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"none"`.

- seed:

  Optional integer for reproducible credible intervals.

## Value

For two groups, class `c("lsa_bayes", "lsa_comparison", "list")` with an
`edges` data frame (`from`, `to`, `prob_a`, `prob_b`, `diff`, `ci_low`,
`ci_high`, `pd`, `effect_size`, `p_value`, `p_adj`, `significant`), the
two `fits`, and the Bayesian settings. For more than two groups, an
all-pairwise
`c("lsa_bayes_pairwise", "lsa_comparison_pairwise", "list")`.

## Details

This complements
[`compare_lsa()`](https://saqr.me/lagdynamics/reference/compare_lsa.md):
the permutation test asks whether a difference is more extreme than
chance; the Bayesian comparison asks for the plausible range of the true
difference and how precisely it is estimated. An edge whose source state
is rarely visited gets a wide credible interval even when its
row-normalised probability looks decisive.

The result carries class `c("lsa_bayes", "lsa_comparison")` (and the
pairwise object `c("lsa_bayes_pairwise", "lsa_comparison_pairwise")`),
so
[plot()](https://saqr.me/lagdynamics/reference/plot.lsa_comparison.md)
(barrel / heatmap) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) work as
for a permutation comparison.

## References

Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.

## See also

[`compare_lsa()`](https://saqr.me/lagdynamics/reference/compare_lsa.md),
[`certainty_lsa()`](https://saqr.me/lagdynamics/reference/certainty_lsa.md)

## Examples

``` r
# \donttest{
g <- lsa(group_regulation,
         group = ifelse(group_regulation$T1 == "plan", "p", "o"))
bc <- bayes_compare_lsa(g, seed = 1)
head(as.data.frame(bc))
#>         from    to      prob_a      prob_b          diff       ci_low
#> 1      adapt adapt 0.001166861 0.005586592 -0.0044197310 -0.027133070
#> 2   cohesion adapt 0.003850193 0.001814882  0.0020353105 -0.005554945
#> 3  consensus adapt 0.004742088 0.005548442 -0.0008063537 -0.006101002
#> 4 coregulate adapt 0.016252683 0.018651363 -0.0023986798 -0.020029418
#> 5    discuss adapt 0.071332032 0.072164948 -0.0008329166 -0.023384136
#> 6    emotion adapt 0.002758328 0.003064351 -0.0003060234 -0.007038335
#>       ci_high     pd effect_size p_value  p_adj significant
#> 1 0.003861497 0.7288 -0.56057932  0.5424 0.5424       FALSE
#> 2 0.007019141 0.8324  0.67017951  0.3352 0.3352       FALSE
#> 3 0.003256711 0.6010 -0.34274817  0.7980 0.7980       FALSE
#> 4 0.010899959 0.5969 -0.30189745  0.8062 0.8062       FALSE
#> 5 0.020201560 0.5158 -0.07450965  0.9684 0.9684       FALSE
#> 6 0.003725192 0.5305 -0.11245333  0.9390 0.9390       FALSE
# }
```
