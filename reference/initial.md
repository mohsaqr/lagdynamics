# Read the initial-state distribution of a fit

`initial()` returns the initial-state distribution of a fit, one row per
state: the proportion of sequences that begin in each state. It
describes how episodes open, which the transition table cannot show, and
it completes the description of the process as a first-order Markov
chain together with
[`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md).

A fit built from a transition-count matrix has no sequences and so no
initial states; the table then has no rows.

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
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or a
  grouped fit from `lsa(..., group = )`.

## Value

A `data.frame` with one row per state and the columns `state` and
`init_prob` (summing to 1). A grouped fit gains a leading `group`
column.

## See also

[`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md),
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)

## Examples

``` r
initial(lsa(engagement))
#>        state init_prob
#> 1     Active 0.3823529
#> 2    Average 0.3676471
#> 3 Disengaged 0.2500000
```
