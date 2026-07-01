# Nodes of an LSA Fit (Tidy)

The states and their incoming / outgoing transition totals, as a tidy
`data.frame` (one row per state).

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
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md).

## Value

A `data.frame` with columns `state` (the state name, matching the
`from`/`to` endpoints of
[`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md)),
`outgoing`, and `incoming` (its total out- and in-transition counts).

## See also

[`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md),
[`tests()`](https://mohsaqr.github.io/lagdynamics/reference/tests.md)

## Examples

``` r
nodes(lsa(group_regulation))
#>        state outgoing incoming
#> 1      adapt      509      531
#> 2   cohesion     1695     1718
#> 3  consensus     6329     6369
#> 4 coregulate     1970     2095
#> 5    discuss     3951     3916
#> 6    emotion     2837     2772
#> 7    monitor     1433     1228
#> 8       plan     6157     6214
#> 9  synthesis      652      690
```
