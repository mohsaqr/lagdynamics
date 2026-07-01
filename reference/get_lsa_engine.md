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

[`register_lsa_engine()`](https://mohsaqr.github.io/lagdynamics/reference/register_lsa_engine.md),
[`list_lsa_engines()`](https://mohsaqr.github.io/lagdynamics/reference/list_lsa_engines.md)
