# Analytic Certainty of Transition Edges (Dirichlet-Multinomial)

Closed-form Bayesian alternative to
[`bootstrap_lsa()`](https://saqr.me/lagdynamics/reference/bootstrap_lsa.md)
for the transition-probability edges of an `lsa` fit. Each state's
outgoing transitions are modelled as Dirichlet-Multinomial: with a
Jeffreys prior the posterior for a row is `Dirichlet(count + prior)`, so
each edge probability is marginally `Beta(a, b)` and its posterior mean,
standard deviation, credible interval and stability decision are
available analytically. No resampling, so it runs in microseconds.

## Usage

``` r
certainty_lsa(
  fit,
  prior = 0.5,
  level_alpha = 0.95,
  inference = c("stability", "threshold"),
  consistency_range = c(0.75, 1.25),
  edge_threshold = NULL
)
```

## Arguments

- fit:

  An `lsa` fit from
  [`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md), or an
  `lsa_group`.

- prior:

  Numeric \> 0. Dirichlet prior concentration added to every cell.
  Default `0.5` (the Jeffreys prior).

- level_alpha:

  Numeric in (0, 1). Credible-interval level. Default `0.95` (a 95%
  interval), matching
  [`bootstrap_lsa()`](https://saqr.me/lagdynamics/reference/bootstrap_lsa.md).

- inference:

  `"stability"` (default) flags an edge whose posterior keeps it within
  a multiplicative `consistency_range` of its observed probability;
  `"threshold"` flags an edge whose posterior mass lies above
  `edge_threshold`.

- consistency_range:

  Length-2 multiplicative bounds for stability inference. Default
  `c(0.75, 1.25)`.

- edge_threshold:

  Numeric or `NULL`. Fixed threshold for `inference = "threshold"`;
  `NULL` uses the 0.10 quantile of non-zero edge probabilities.

## Value

An object of class `c("lsa_certainty", "lsa_bootstrap", "list")` with an
`edges` data frame (`from`, `to`, `prob_observed`, `prob_mean`,
`prob_se`, `prob_ci_low`, `prob_ci_high`, `p_value`, `stable`, plus
`adj_res_observed`/`adj_res_stable` for plotting), the posterior
matrices (`mean`, `sd`, `ci_lower`, `ci_upper`), and call metadata
(`prior`, `level_alpha`, `inference`, ...). For an `lsa_group`, a named
list of these (class `lsa_certainty_group`).

## Details

The result carries class `c("lsa_certainty", "lsa_bootstrap")` and an
`edges` table with the columns
[`plot_forest()`](https://saqr.me/lagdynamics/reference/plot_forest.md)
and [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
expect, so it is a drop-in for a bootstrap result (use
`metric = "prob"`).

**Certainty vs bootstrap.** Both answer "how precisely is this edge
pinned down?". They agree on homogeneous data. The Dirichlet posterior
treats transitions as independent, so on strongly heterogeneous data (a
mixture of latent classes with long sequences) it reports *more*
certainty than the sequence bootstrap – prefer
[`bootstrap_lsa()`](https://saqr.me/lagdynamics/reference/bootstrap_lsa.md)
then.

## References

Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.

## See also

[`bootstrap_lsa()`](https://saqr.me/lagdynamics/reference/bootstrap_lsa.md),
[`stability_lsa()`](https://saqr.me/lagdynamics/reference/stability_lsa.md),
[`plot_forest()`](https://saqr.me/lagdynamics/reference/plot_forest.md)

## Examples

``` r
# \donttest{
fit  <- lsa(engagement)
cert <- certainty_lsa(fit)
cert
#> <lsa_certainty>  (analytic Dirichlet-Multinomial)
#>   engine:        classical
#>   prior:         Dirichlet(0.50)
#>   CI level:      95%  |  inference: stability
#>   certain edges: 7 of 9
head(as.data.frame(cert))
#>         from      to observed prob_observed prob_mean    prob_se prob_ci_low
#> 1     Active  Active      459     0.6975684 0.6967400 0.01788572  0.66113410
#> 2    Average  Active      153     0.2037284 0.2039867 0.01467978  0.17597523
#> 3 Disengaged  Active       39     0.1200000 0.1209801 0.01801983  0.08793233
#> 4     Active Average      176     0.2674772 0.2676270 0.01722641  0.23454702
#> 5    Average Average      458     0.6098535 0.6093023 0.01777441  0.57420054
#> 6 Disengaged Average      129     0.3969231 0.3966309 0.02703206  0.34428578
#>   prob_ci_high      p_value stable adj_res_observed adj_res_stable
#> 1    0.7312156 5.560385e-20   TRUE        21.662811           TRUE
#> 2    0.2334883 6.063196e-04   TRUE       -12.906038           TRUE
#> 3    0.1584213 9.448149e-02  FALSE       -10.549486          FALSE
#> 4    0.3020418 1.205186e-04   TRUE       -11.319132           TRUE
#> 5    0.6438538 3.141020e-17   TRUE        12.452615           TRUE
#> 6    0.4501757 2.232689e-04   TRUE        -1.736461           TRUE
# }
```
