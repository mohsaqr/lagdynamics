# Print and read resampling results

[`print()`](https://rdrr.io/r/base/print.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) read the
results of
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)
and
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md).
[`print()`](https://rdrr.io/r/base/print.html) shows a short summary:
the engine, the number of replicates, the resampling settings and how
many transitions are stable or significant (with the number of failed
replicates for a stability result).
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the per-transition table of the result.

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
returns the same table and adds the arguments `significant` (keep the
transitions the method flags) and `sort` (order by the size of the
effect); it is the usual way to read these results. Results of
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
use the bootstrap methods.

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

- ...:

  Not used.

- row.names, optional:

  Not used; present for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
data frame with one row per transition (`from -> to` pair); its columns
are listed under *Value* in
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)
and
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md).

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md).

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
#>   failed:        0 of 20 replicates
#>   stable edges:  7 of 9 (>= 95% across replicates)
head(as.data.frame(st))
#>         from      to observed_sig stability stable
#> 1     Active  Active         TRUE       1.0   TRUE
#> 2    Average  Active         TRUE       1.0   TRUE
#> 3 Disengaged  Active         TRUE       1.0   TRUE
#> 4     Active Average         TRUE       1.0   TRUE
#> 5    Average Average         TRUE       1.0   TRUE
#> 6 Disengaged Average        FALSE       0.1  FALSE
```
