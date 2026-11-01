# Register a Lag Sequential Analysis Engine

Adds a new engine to the lagdynamics registry so it can be referenced by
name via `lsa(..., engine = "<name>")`. Built-in engines (`"classical"`,
`"two_cell"`, `"bidirectional"`, `"parallel_dominance"`,
`"nonparallel_dominance"`) are registered automatically when the package
loads.

## Usage

``` r
register_lsa_engine(name, fn, description, requires = character())
```

## Arguments

- name:

  Character scalar. The engine's identifier as used in
  `lsa(engine = name)`.

- fn:

  A function.
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) calls
  it with `transitions` (the `lsa_transitions` object produced by
  [`lsa_transitions()`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md)),
  `structural_zeros` (`NULL` or a `K x K` 0/1 matrix), `alternative`,
  `n_events` (`NULL` for transition-matrix input), and any engine
  parameters supplied through `params` or `...` of
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), so
  it should accept `transitions` and `...`. Must return a named list
  with at least the matrix elements `obs`, `exp`, `prob`, `adj_res`, and
  `p` (each `K x K`). Additional `K x K` matrices are preserved as
  engine-specific edge statistics.

- description:

  Character scalar. One-line human-readable description shown by
  [`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md).

- requires:

  Character vector. Names of packages the engine depends on. Empty by
  default.

## Value

Invisibly returns `name` (a character scalar).

## See also

[`get_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/get_lsa_engine.md),
[`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md),
[`unregister_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/unregister_lsa_engine.md),
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)

## Examples

``` r
my_engine <- function(transitions, ...) {
  get_lsa_engine("classical")$fn(transitions, ...)
}
register_lsa_engine("my_classical", my_engine, "Classical test alias")
fit <- lsa(engagement, engine = "my_classical")
unregister_lsa_engine("my_classical")
fit
#> Lag Sequential Analysis  -  my_classical  (lag 1, directed)
#>   3 states | 1734 transitions | 1870 events | 136 sequences
#>   states: Active, Average, Disengaged
#>   independence: G² = 618.3, df = 4, p <2e-16
#> 
#>   Significant transitions (p < 0.05): 7 of 9
#>   strongest over-represented (of 3):
#>     Active -> Active          z =  +21.7  ***
#>     Disengaged -> Disengaged  z =  +15.4  ***
#>     Average -> Average        z =  +12.5  ***
#> 
#>   Initial states:
#>     Active     0.382  ████████████████████████
#>     Average    0.368  ███████████████████████
#>     Disengaged 0.250  ████████████████
```
