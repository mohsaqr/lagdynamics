# Fit expected counts with structural zeros by iterative proportional fitting

`lsa_ipf()` fits the expected counts of a count table with structural
zeros by iterative proportional fitting. When some transitions cannot
occur by construction (structural zeros, such as self-transitions in a
coding scheme that merges repeated codes), the quasi-independence model
fits expected counts that are zero in the forbidden cells and reproduce
the observed row and column totals; it has no closed form and is found
by iterative proportional fitting (Wickens, 1989), which rescales the
rows and columns of the table in turn until both sets of totals match.
Under ordinary independence, the expected count of a transition is the
product of its row and column totals divided by the grand total, a
formula that assigns expected counts to impossible cells.

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) calls
`lsa_ipf()` when `loops = FALSE` or `structural_zeros` is given, so it
is needed directly only to inspect the expected table or to use it
outside
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

## Usage

``` r
lsa_ipf(obs, structure = NULL, tol = 1e-08, max_iter = 200L)
```

## Arguments

- obs:

  Square numeric matrix of observed transition counts.

- structure:

  Square 0/1 matrix of the same size: `0` marks a structural zero, `1`
  an estimable cell. Default `NULL`: every cell is estimable, which
  gives the ordinary independence table.

- tol:

  Convergence tolerance on the largest difference between fitted and
  observed totals (default `1e-8`).

- max_iter:

  Maximum number of row-and-column passes (default `200`).

## Value

A list with the fitted `K x K` expected-count matrix (`fit`), the number
of passes used (`iterations`), whether the totals converged within
`max_iter` (`converged`), and the largest remaining difference between
fitted and observed totals (`max_margin_diff`).

## Conditions

Errors carry the class `lagdynamics_error` and one of:

- `lagdynamics_bad_input`: `obs` contains missing, infinite or negative
  values.

- `lagdynamics_ipf_infeasible`: a row or column has observed transitions
  but every one of its cells is a structural zero.

## References

Wickens, T. D. (1989). *Multiway contingency tables analysis for the
social sciences* (pp. 107-112). Lawrence Erlbaum.

## See also

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) and its
`loops` and `structural_zeros` arguments.

## Examples

``` r
obs <- matrix(c(0, 4, 6,
                3, 0, 5,
                7, 2, 0), nrow = 3, byrow = TRUE)
# Every cell estimable: the ordinary independence table
lsa_ipf(obs)
#> $fit
#>          [,1]     [,2]     [,3]
#> [1,] 3.703704 2.222222 4.074074
#> [2,] 2.962963 1.777778 3.259259
#> [3,] 3.333333 2.000000 3.666667
#> 
#> $iterations
#> [1] 1
#> 
#> $converged
#> [1] TRUE
#> 
#> $max_margin_diff
#> [1] 0
#> 

# Self-transitions forbidden: zero expected counts on the diagonal,
# observed row and column totals kept
lsa_ipf(obs, structure = 1 - diag(3))
#> $fit
#>          [,1]     [,2]     [,3]
#> [1,] 0.000000 3.104423 6.895577
#> [2,] 3.895577 0.000000 4.104423
#> [3,] 6.104423 2.895577 0.000000
#> 
#> $iterations
#> [1] 20
#> 
#> $converged
#> [1] TRUE
#> 
#> $max_margin_diff
#> [1] 4.908683e-09
#> 
```
