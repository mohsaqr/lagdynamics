# Compute the analytic certainty of transition probabilities

`certainty_lsa()` computes the analytic certainty of every transition
probability of a fit from
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md): how
precisely the data pin down each transition probability, and which
transitions are estimated precisely enough to be called certain. Under a
Dirichlet-Multinomial model the outgoing transitions of each state are a
multinomial draw whose probabilities have a Dirichlet posterior, so each
transition probability has a Beta posterior with a closed-form mean,
standard deviation and credible interval. A transition probability
estimated from counts is uncertain in proportion to how often its source
state was visited: a probability of 0.5 from 4 departures is far less
certain than the same value from 400.

The posterior is computed with a Jeffreys prior by default and needs no
resampling. Each transition gets its posterior mean, standard deviation
and credible interval, and a posterior tail probability that decides
whether it is certain. With `inference = "stability"` a transition is
certain when its posterior stays within a band around the observed
probability (by default 75% to 125% of it); with
`inference = "threshold"` it is certain when its posterior lies above a
minimum probability. A grouped fit gives one result per group.

The model treats the transitions as independent draws. When sequences
differ strongly from one another (for example, a mixture of latent
groups with long sequences) it reports more certainty than the data
support; use
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
which resamples whole sequences, in that case. To compare two groups'
transition probabilities with the same model, use
[`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md).

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
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or a
  grouped fit from `lsa(..., group = )`.

- prior:

  Dirichlet concentration added to every cell. Default `0.5` (the
  Jeffreys prior).

- level_alpha:

  Credible-interval level. Default `0.95`. The certainty decision uses
  `1 - level_alpha` as its threshold.

- inference:

  `"stability"` (default): a transition is certain when little posterior
  mass lies outside `consistency_range` times its observed probability.
  `"threshold"`: certain when little posterior mass lies below
  `edge_threshold`.

- consistency_range:

  Two multipliers of the observed probability that bound the band for
  `inference = "stability"`. Default `c(0.75, 1.25)`.

- edge_threshold:

  Minimum probability for `inference = "threshold"`. Default `NULL`, the
  10% quantile of the non-zero transition probabilities.

- x:

  An `lsa_certainty` or `lsa_certainty_group` object, for the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods.

- ...:

  Not used.

- row.names, optional:

  Not used; present for consistency with the generic.

## Value

An object of class `c("lsa_certainty", "lsa_bootstrap", "list")`. Read
it with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
which returns a data frame with one row per transition (`from -> to`
pair) and the columns `from`, `to`, `observed` (count), `prob_observed`,
`prob_mean` (posterior mean), `prob_se` (posterior standard deviation),
`prob_ci_low`, `prob_ci_high` (credible interval), `p_value` (posterior
mass outside the band, or below the threshold; the upper tail is
computed directly with `pbeta(lower.tail = FALSE)`, so tiny masses keep
their relative accuracy), `stable` (certain), and `adj_res_observed`,
`adj_res_stable` (the fit's adjusted residual and a copy of `stable`,
used by
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)).
Transitions that cannot be estimated (a source state never left, or a
structural zero) are `NA` and never certain.
`transitions(x, significant = TRUE)` keeps the certain transitions.

For a grouped fit, an object of class `lsa_certainty_group`, a list with
one `lsa_certainty` per group;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) stacks
their tables with a leading `group` column.
[`print()`](https://rdrr.io/r/base/print.html) shows the settings and
the number of certain transitions and returns `x` invisibly.

## References

DiDonato, A. R., & Morris, A. H. (1992). Algorithm 708: Significant
digit computation of the incomplete beta function ratios. *ACM
Transactions on Mathematical Software, 18*(3), 360-373.
[doi:10.1145/131766.131776](https://doi.org/10.1145/131766.131776)

Johnston, L., & Jendoubi, T. (2026). How delivery mode reshapes resource
engagement: A Bayesian differential network analysis. TNA Workshop 2026.

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
to read the result;
[`plot.lsa_certainty()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_certainty.md)
and
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)
to draw it;
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md)
for resampling intervals;
[`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
to compare groups under the same model.

## Examples

``` r
fit <- lsa(engagement)
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
#> 4    0.6438538 3.141676e-17   TRUE       12.4526148           TRUE
#> 5    0.4501757 2.232689e-04   TRUE       -1.7364613           TRUE
#> 6    0.2153070 1.168633e-03   TRUE        0.1757996           TRUE
#> 7    0.5365741 9.440665e-06   TRUE       15.3904661           TRUE

# One table per group
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
