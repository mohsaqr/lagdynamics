# Remove a Registered LSA Engine

Remove a Registered LSA Engine

## Usage

``` r
unregister_lsa_engine(name)
```

## Arguments

- name:

  Character scalar. The engine's identifier.

## Value

Invisibly `NULL`.

## See also

[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md)

## Examples

``` r
my_engine <- function(transitions, ...) {
  get_lsa_engine("classical")$fn(transitions, ...)
}
register_lsa_engine("temporary_engine", my_engine, "Temporary alias")
unregister_lsa_engine("temporary_engine")
```
