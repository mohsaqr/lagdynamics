# Retrieve a registered engine

`get_lsa_engine()` returns the registry entry of one engine: its name,
its function, its description and the packages it requires. It is the
lookup [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
performs before fitting, and it is useful for inspecting a built-in
engine's function or for checking that a registered engine can run in
the current session. To see which engines exist, use
[`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md).

## Usage

``` r
get_lsa_engine(name)
```

## Arguments

- name:

  Name of the engine.

## Value

A list with the elements `name`, `fn`, `description` and `requires`, as
given to
[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md).

## Conditions

Errors carry the class `lagdynamics_error` and one of:

- `lagdynamics_unknown_engine`: no engine is registered under `name`.

- `lagdynamics_missing_package`: a package the engine requires is not
  installed.

## See also

[`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md),
[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md)

## Examples

``` r
str(get_lsa_engine("classical"), give.attr = FALSE)
#> List of 4
#>  $ name       : chr "classical"
#>  $ fn         :function (transitions, structural_zeros = NULL, alternative = c("two.sided", 
#>     "greater", "less"), n_events = NULL, ...)  
#>  $ description: chr "Bakeman & Quera classical lag sequential analysis"
#>  $ requires   : chr(0) 
```
