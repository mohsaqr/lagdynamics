# Implementation notes removed from `?transfer_entropy` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/transfer_entropy.R); it is kept here for maintainers.

- **State-flow mode** (`.te_state_network()`): binary occupancy
  decomposition. For each ordered pair of distinct states (s, r), the
  source series is the indicator `x == s` and the target `x == r`; `NA`
  stays `NA`. Self pairs are excluded. With K = 3 the non-target
  indicators are redundant, which is why pairs share values.
- **History**: the target's last `history` values are combined into a
  composite symbol; triples (source past, target history, target future)
  are pooled across sequences.
- **Surrogate test** (`.te_surrogate()`): shuffles the source column
  within the pooled triples; `te_effective = te - mean(surrogates)`,
  `p = (1 + #{surrogate >= observed}) / (R + 1)`.
- **Bivariate labels** come from `substitute(x)` / `substitute(y)`.
- **seed**: validated with a plain `stopifnot()`; RNG restored on exit.
