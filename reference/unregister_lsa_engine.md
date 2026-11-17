# Remove a registered engine

`unregister_lsa_engine()` removes an engine from the registry for the
rest of the session, so that `lsa(engine = name)` no longer finds it. It
is mainly used to clean up an engine added with
[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md).
A removed built-in engine returns when the package is loaded again.

## Usage

``` r
unregister_lsa_engine(name)
```

## Arguments

- name:

  Name of the engine.

## Value

`NULL`, invisibly.

## Conditions

Errors carry the class `lagdynamics_error` and
`lagdynamics_unknown_engine` when no engine is registered under `name`.

## See also

[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md),
[`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md)

## Examples

``` r
placeholder <- function(transitions, ...) {
  obs <- transitions$obs
  expected <- outer(rowSums(obs), colSums(obs)) / sum(obs)
  list(obs = obs, exp = expected, prob = obs / rowSums(obs),
       adj_res = obs * 0, p = obs * 0 + 1)
}
register_lsa_engine("placeholder", placeholder,
                    description = "Placeholder engine")
unregister_lsa_engine("placeholder")
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
