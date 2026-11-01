# Iterative Proportional Fitting for Two-Way Tables with Structural Zeros

Fits expected cell frequencies under the row + column independence model
when some cells are constrained to be exactly zero (structural zeros).
The fitted table has the same row and column marginals as `obs`,
satisfies `E[i, j] = 0` wherever `structure[i, j] = 0`, and is the
maximum likelihood estimate under the independence model restricted to
the non-zero pattern.

## Usage

``` r
lsa_ipf(obs, structure = NULL, tol = 1e-08, max_iter = 200L)
```

## Arguments

- obs:

  Numeric `K x K` matrix of observed counts.

- structure:

  Numeric `K x K` 0/1 matrix. **Default `matrix(1, K, K)`: every cell,
  including the diagonal, is estimable – self-transitions and every
  observed cell are kept.** A `0` marks a cell as a structural zero
  (forbidden), a `1` marks it as estimable. Pass an explicit pattern
  (e.g. `1 - diag(K)`) only when you want to *opt out* of specific cells
  because the coding scheme makes them impossible by construction.

- tol:

  Numeric. Convergence tolerance on marginal differences. Default
  `1e-8`.

- max_iter:

  Integer. Maximum number of row+column scaling passes. Default `200L`.

## Value

A plain list (no class) with elements:

- fit:

  The fitted `K x K` expected-frequency matrix.

- iterations:

  Number of row+column passes used.

- converged:

  Logical scalar: whether convergence was reached within `max_iter`.

- max_margin_diff:

  Maximum absolute difference between observed and fitted marginals at
  termination.

## Details

Implementation follows Wickens (1989), pp. 107-112: alternately scale
rows then columns of an initialized expected table until row and column
marginals converge to those of `obs` within `tol`.

## References

Wickens, T. D. (1989). *Multiway contingency tables analysis for the
social sciences*, pp. 107-112. Lawrence Erlbaum.

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

# Diagonal forbidden: the empty self-transition cells become structural
# zeros and the fit keeps the observed row and column totals
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
