# Retrieve a Registered LSA Engine

Retrieve a Registered LSA Engine

## Usage

``` r
get_lsa_engine(name)
```

## Arguments

- name:

  Character scalar. The engine's identifier.

## Value

The registry entry: a list with elements `name`, `fn`, `description`,
`requires`.

## See also

[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md),
[`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md)

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
