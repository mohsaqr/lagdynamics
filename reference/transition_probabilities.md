# Read the transition probabilities of a fit

`transition_probabilities()` returns the transition probabilities of a
fit, one row per ordered pair of states. The transition probability
\\P(\text{to} \mid \text{from})\\ is the share of the transitions
leaving a state that go to each other state, so the probabilities out of
one state sum to 1. They describe where the process goes next,
regardless of chance expectation; the adjusted residuals of
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
say which of these moves are more or less frequent than independence
predicts. Together with the initial-state probabilities of
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md),
they are the parameters of the first-order Markov model read by
Transition Network Analysis.

A cell declared a structural zero has probability 0, as does every cell
of a state that is never left.

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
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or a
  grouped fit from `lsa(..., group = )`.

## Value

A `data.frame` with one row per `from -> to` pair (`K^2` rows for `K`
states, ordered by `from` and then `to`) and the columns `from`, `to`
and `prob`. A grouped fit stacks the groups under a leading `group`
column.

## See also

[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
for the initial-state probabilities;
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
which also carries `prob`;
[`lsa_to_tna()`](https://pak.dynasite.org/lagdynamics/reference/lsa_to_tna.md).

## Examples

``` r
transition_probabilities(lsa(engagement))
#>         from         to       prob
#> 1     Active     Active 0.69756839
#> 2     Active    Average 0.26747720
#> 3     Active Disengaged 0.03495441
#> 4    Average     Active 0.20372836
#> 5    Average    Average 0.60985353
#> 6    Average Disengaged 0.18641811
#> 7 Disengaged     Active 0.12000000
#> 8 Disengaged    Average 0.39692308
#> 9 Disengaged Disengaged 0.48307692
```
