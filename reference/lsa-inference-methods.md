# Print and Tidy Resampling Results

Methods for the results of
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)
and
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md).
[`print()`](https://rdrr.io/r/base/print.html) shows a header summary:
the engine, the number of replicates, the resampling settings, and how
many edges are stable or significant.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the per-edge table described under *Value* of the constructor;
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
returns the same table and adds `significant =` and `sort =` arguments.

## Usage

``` r
# S3 method for class 'lsa_bootstrap'
print(x, ...)

# S3 method for class 'lsa_permutation'
print(x, ...)

# S3 method for class 'lsa_stability'
print(x, ...)

# S3 method for class 'lsa_bootstrap'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_permutation'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_stability'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `lsa_bootstrap`, `lsa_permutation` or `lsa_stability` object.

- row.names, optional, ...:

  Ignored; present for consistency with the generics.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
`data.frame` with one row per `from -> to` cell: the bootstrap,
permutation or stability columns listed in
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)
and
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md).

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)

## Examples

``` r
fit <- lsa(engagement)
st <- stability_lsa(fit, R = 20)
st
#> <lsa_stability>
#>   engine:        classical
#>   replicates:    20
#>   proportion:    80%
#>   min stable:    95%
#>   stable edges:  7 of 9 (>= 95% across replicates)
head(as.data.frame(st))
#>         from      to observed_sig stability stable
#> 1     Active  Active         TRUE      1.00   TRUE
#> 2    Average  Active         TRUE      1.00   TRUE
#> 3 Disengaged  Active         TRUE      1.00   TRUE
#> 4     Active Average         TRUE      1.00   TRUE
#> 5    Average Average         TRUE      1.00   TRUE
#> 6 Disengaged Average        FALSE      0.15  FALSE
```
