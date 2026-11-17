# List the registered engines

`list_lsa_engines()` lists the engines available in the session, the
five built-in ones and any added with
[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md),
with a one-line description of each. Each engine of
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
implements one test of the transition table. Use the list to choose the
`engine` argument of
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

## Usage

``` r
list_lsa_engines()
```

## Value

A `data.frame` with one row per registered engine, in alphabetical
order, and the columns `name`, `description` and `requires`
(comma-separated package names, `""` when none).

## See also

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
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
