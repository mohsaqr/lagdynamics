# List All Registered LSA Engines

List All Registered LSA Engines

## Usage

``` r
list_lsa_engines()
```

## Value

A data.frame with columns `name`, `description`, `requires`.

## See also

[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md),
[`get_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/get_lsa_engine.md)

## Examples

``` r
list_lsa_engines()
#>                    name
#> 1         bidirectional
#> 2             classical
#> 3 nonparallel_dominance
#> 4    parallel_dominance
#> 5              two_cell
#>                                                            description requires
#> 1 Sackett's bidirectional / matched-pair test on the symmetrized table         
#> 2                    Bakeman & Quera classical lag sequential analysis         
#> 3       Sackett's non-parallel-dominance (observed-SE + binomial) test         
#> 4                      Sackett's parallel-dominance (expected-SE) test         
#> 5                  2x2 cell test (odds ratio, log-OR Wald z, Yule's Q)         
```
