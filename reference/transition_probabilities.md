# Transition-Probability Matrix of an LSA Fit

The row-stochastic transition-probability matrix \\P(\text{to} \mid
\text{from})\\: each row holds the distribution over next states out of
a state, so the entries in a row sum to 1 (up to structural zeros). This
is the transition-probability matrix a Transition Network Analysis
reads; pair it with
[`initial()`](https://saqr.me/lagdynamics/reference/initial.md) for the
initial-state probabilities.

## Usage

``` r
transition_probabilities(fit)

# S3 method for class 'lsa'
transition_probabilities(fit)

# S3 method for class 'lsa_group'
transition_probabilities(fit)
```

## Arguments

- fit:

  An `lsa` fit from
  [`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md).

## Value

A square numeric matrix with `from` states in rows and `to` states in
columns. For an `lsa_group`, a named list of such matrices, one per
group.

## See also

[`initial()`](https://saqr.me/lagdynamics/reference/initial.md)
(initial-state probabilities),
[`transitions()`](https://saqr.me/lagdynamics/reference/transitions.md)
(the tidy per-edge view),
[`nodes()`](https://saqr.me/lagdynamics/reference/nodes.md)

## Examples

``` r
transition_probabilities(lsa(group_regulation))
#>                   adapt   cohesion  consensus coregulate    discuss    emotion
#> adapt      0.0000000000 0.27308448 0.47740668 0.02161100 0.05893910 0.11984283
#> cohesion   0.0029498525 0.02713864 0.49793510 0.11917404 0.05958702 0.11563422
#> consensus  0.0047400853 0.01485227 0.08200348 0.18770738 0.18802338 0.07268131
#> coregulate 0.0162436548 0.03604061 0.13451777 0.02335025 0.27360406 0.17208122
#> discuss    0.0713743356 0.04758289 0.32118451 0.08428246 0.19488737 0.10579600
#> emotion    0.0024673951 0.32534367 0.32040888 0.03419105 0.10186817 0.07684173
#> monitor    0.0111653873 0.05582694 0.15910677 0.05792045 0.37543615 0.09071877
#> plan       0.0009745006 0.02517460 0.29040117 0.01721618 0.06789021 0.14682475
#> synthesis  0.2346625767 0.03374233 0.46625767 0.04447853 0.06288344 0.07055215
#>               monitor       plan   synthesis
#> adapt      0.03339882 0.01571709 0.000000000
#> cohesion   0.03303835 0.14100295 0.003539823
#> consensus  0.04661084 0.39579712 0.007584137
#> coregulate 0.08629442 0.23908629 0.018781726
#> discuss    0.02227284 0.01164262 0.140976968
#> emotion    0.03630596 0.09975326 0.002819880
#> monitor    0.01814375 0.21563154 0.016050244
#> plan       0.07552379 0.37420822 0.001786584
#> synthesis  0.01226994 0.07515337 0.000000000
```
