# Analytic Certainty of Transition Edges (Dirichlet-Multinomial)

Closed-form Bayesian alternative to
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md)
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

# S3 method for class 'lsa_certainty'
print(x, ...)

# S3 method for class 'lsa_certainty_group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_certainty_group'
print(x, ...)
```

## Arguments

- fit:

  An `lsa` fit from
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or
  an `lsa_group`.

- prior:

  Numeric \> 0. Dirichlet prior concentration added to every cell.
  Default `0.5` (the Jeffreys prior).

- level_alpha:

  Numeric in (0, 1). Credible-interval level. Default `0.95` (a 95%
  interval), matching
  [`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md).

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

- x:

  An `lsa_certainty` or `lsa_certainty_group` object
  ([`print()`](https://rdrr.io/r/base/print.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods).

- row.names, optional, ...:

  Ignored by the methods; present for consistency with the generics.

## Value

An object of class `c("lsa_certainty", "lsa_bootstrap", "list")` with an
`edges` data frame (`from`, `to`, `prob_observed`, `prob_mean`,
`prob_se`, `prob_ci_low`, `prob_ci_high`, `p_value`, `stable`, plus
`adj_res_observed`/`adj_res_stable` for plotting), the posterior
matrices (`mean`, `sd`, `ci_lower`, `ci_upper`), and call metadata
(`prior`, `level_alpha`, `inference`, ...). Read the per-edge table with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md).
For an `lsa_group`, a named list of these, one per group (class
`lsa_certainty_group`), whose
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) stacks
the group tables with a leading `group` column.
[`print()`](https://rdrr.io/r/base/print.html) shows a header summary
(the number of certain edges) and returns `x` invisibly.

## Details

The result carries class `c("lsa_certainty", "lsa_bootstrap")` and an
`edges` table with the columns
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)
and
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
expect, so it is a drop-in for a bootstrap result (use
`metric = "prob"`).

**Certainty vs bootstrap.** Both answer "how precisely is this edge
pinned down?". They agree on homogeneous data. The Dirichlet posterior
treats transitions as independent, so on strongly heterogeneous data (a
mixture of latent classes with long sequences) it reports *more*
certainty than the sequence bootstrap – prefer
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md)
then.

## References

Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.

## See also

[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md),
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)

## Examples

``` r
fit  <- lsa(engagement)
cert <- certainty_lsa(fit)
cert
#> <lsa_certainty>  (analytic Dirichlet-Multinomial)
#>   engine:        classical
#>   prior:         Dirichlet(0.50)
#>   CI level:      95%  |  inference: stability
#>   certain edges: 7 of 9
transitions(cert, significant = TRUE)
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
#> 4    0.6438538 3.141020e-17   TRUE       12.4526148           TRUE
#> 5    0.4501757 2.232689e-04   TRUE       -1.7364613           TRUE
#> 6    0.2153070 1.168633e-03   TRUE        0.1757996           TRUE
#> 7    0.5365741 9.440665e-06   TRUE       15.3904661           TRUE

# One certainty table per group
g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
         time = "Time", group = "Achiever")
gcert <- certainty_lsa(g)
gcert
#> <lsa_certainty_group>
#>   groups: 2 (High, Low)
#>     - High:        40 certain edges
#>     - Low:         43 certain edges
head(as.data.frame(gcert))
#>   group       from    to observed prob_observed   prob_mean     prob_se
#> 1  High      adapt adapt        0   0.000000000 0.003436426 0.004834897
#> 2  High   cohesion adapt        5   0.005330490 0.005835544 0.002479698
#> 3  High  consensus adapt       14   0.004127358 0.004269101 0.001118560
#> 4  High coregulate adapt       20   0.022371365 0.022815804 0.004978580
#> 5  High    discuss adapt       48   0.023964054 0.024159402 0.003426073
#> 6  High    emotion adapt        5   0.003225806 0.003538115 0.001505502
#>    prob_ci_low prob_ci_high    p_value stable adj_res_observed adj_res_stable
#> 1 3.392283e-06  0.017203890 1.00000000  FALSE        -1.273187          FALSE
#> 2 2.029216e-03  0.011601136 0.56804423  FALSE        -1.784125          FALSE
#> 3 2.364900e-03  0.006723496 0.34397924  FALSE        -4.589083          FALSE
#> 4 1.410149e-02  0.033535629 0.25336217  FALSE         3.273701          FALSE
#> 5 1.790475e-02  0.031310858 0.07902483  FALSE         5.884023          FALSE
#> 6 1.229139e-03  0.007040395 0.56880929  FALSE        -3.194153          FALSE
```
