# Read the states of a lag sequential analysis fit

`nodes()` returns the states of a lag sequential analysis fit, one row
per state with its total outgoing and incoming transition counts. In a
lag transition network the states are the nodes, and how often a state
is left and how often it is entered describe its place in the process: a
state with many outgoing but few incoming transitions tends to open
episodes, one with the reverse tends to close them. The marginal totals
also set the expected counts against which every transition is tested.

The shares of all transitions and the initial-state probability are in
`as.data.frame(summary(fit), what = "nodes")`; the per-transition table
is
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md).

## Usage

``` r
nodes(fit)

# S3 method for class 'lsa'
nodes(fit)

# S3 method for class 'lsa_group'
nodes(fit)
```

## Arguments

- fit:

  An `lsa` fit from
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or a
  grouped fit from `lsa(..., group = )`.

## Value

A `data.frame` with one row per state and the columns `state` (matching
the `from` and `to` names of
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)),
`outgoing` and `incoming`. A grouped fit gives one row per group and
state, with a leading `group` column.

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md),
[lsa-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
for [`summary()`](https://rdrr.io/r/base/summary.html).

## Examples

``` r
nodes(lsa(engagement))
#>        state outgoing incoming
#> 1     Active      658      651
#> 2    Average      751      763
#> 3 Disengaged      325      320
```
