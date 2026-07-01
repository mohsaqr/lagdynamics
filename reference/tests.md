# Tablewise Independence Tests of an LSA Fit (Tidy)

The global tests of independence (likelihood-ratio G^2 and Pearson
chi-square) as a one-row-per-test `data.frame`.

## Usage

``` r
tests(fit)

# S3 method for class 'lsa'
tests(fit)

# S3 method for class 'lsa_group'
tests(fit)
```

## Arguments

- fit:

  An `lsa` fit from
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md).

## Value

A `data.frame` with columns `test` (`"lrx2"` / `"x2"`), `statistic`,
`df`, `p`. Tests the engine did not compute are omitted.

## See also

[`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md),
[`nodes()`](https://mohsaqr.github.io/lagdynamics/reference/nodes.md)

## Examples

``` r
tests(lsa(group_regulation))
#>   test statistic df p
#> 1 lrx2  13203.77 64 0
#> 2   x2  15435.26 64 0
```
