# Initial-State Distribution of an LSA Fit (Tidy)

The proportion of sequences starting in each state, as a tidy
`data.frame`. These are the initial-state probabilities (init P) that
complement the transition-probability matrix from
[`transition_probabilities()`](https://mohsaqr.github.io/lagdynamics/reference/transition_probabilities.md).

## Usage

``` r
initial(fit)

# S3 method for class 'lsa'
initial(fit)

# S3 method for class 'lsa_group'
initial(fit)
```

## Arguments

- fit:

  An `lsa` fit from
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md).

## Value

A `data.frame` with columns `state`, `init_prob`; zero rows when the fit
came from a transition matrix (no initial states).

## See also

[`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md),
[`nodes()`](https://mohsaqr.github.io/lagdynamics/reference/nodes.md)

## Examples

``` r
initial(lsa(group_regulation))
#>        state init_prob
#> 1      adapt    0.0115
#> 2   cohesion    0.0605
#> 3  consensus    0.2140
#> 4 coregulate    0.0190
#> 5    discuss    0.1755
#> 6    emotion    0.1515
#> 7    monitor    0.1440
#> 8       plan    0.2045
#> 9  synthesis    0.0195
```
