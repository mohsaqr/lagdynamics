# Read the tablewise tests of independence of a fit

`tests()` returns the tablewise tests of independence of a fit, one row
per test. They test whether the table as a whole departs from
independence: whether the state at one time carries any information
about the state `lag` events later (Bakeman & Quera, 1995). The
likelihood-ratio statistic \\G^2\\ and the Pearson \\\chi^2\\ statistic
test this omnibus hypothesis against a chi-square distribution, and are
read before individual transitions are interpreted.

The tests are computed by the classical engine; the other built-in
engines test cells only, and `tests()` then returns a table with no
rows.

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
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or a
  grouped fit from `lsa(..., group = )`.

## Value

A `data.frame` with one row per test and the columns `test` (`"lrx2"`
for the likelihood-ratio test, `"x2"` for Pearson's test), `statistic`,
`df` and `p`. A grouped fit gains a leading `group` column.

## References

Bakeman, R., & Quera, V. (1995). *Analyzing interaction: Sequential
analysis with SDIS and GSEQ*. Cambridge University Press.

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
for the cell-wise tests;
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md);
[lsa-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
for [`summary()`](https://rdrr.io/r/base/summary.html).

## Examples

``` r
tests(lsa(engagement))
#>   test statistic df             p
#> 1 lrx2  618.3053  4 1.691520e-132
#> 2   x2  629.2022  4 7.406567e-135
```
