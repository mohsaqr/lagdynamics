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

  A function. Must accept a `transitions` argument (a tidy transition
  table produced by
  [`lsa_transitions()`](https://saqr.me/lagdynamics/reference/lsa_transitions.md))
  and arbitrary named `...` arguments forwarded from
  `lsa(params = list(...))`. Must return a named list with at least the
  matrix elements `obs`, `exp`, `prob`, `adj_res`, and `p` (each
  `K x K`). See the built-in `.engine_classical` for the full contract.

- description:

  Character scalar. One-line human-readable description shown by
  [`list_lsa_engines()`](https://saqr.me/lagdynamics/reference/list_lsa_engines.md).

- requires:

  Character vector. Names of packages the engine depends on. Empty by
  default.

## Value

Invisibly returns `name`.

## See also

[`get_lsa_engine()`](https://saqr.me/lagdynamics/reference/get_lsa_engine.md),
[`list_lsa_engines()`](https://saqr.me/lagdynamics/reference/list_lsa_engines.md),
[`unregister_lsa_engine()`](https://saqr.me/lagdynamics/reference/unregister_lsa_engine.md),
[`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md)

## Examples

``` r
if (FALSE) { # \dontrun{
my_engine <- function(transitions, ...) {
  # ... compute and return a list with obs, exp, prob, adj_res, p
}
register_lsa_engine("my_engine", my_engine, "Custom LSA variant")
} # }
```
