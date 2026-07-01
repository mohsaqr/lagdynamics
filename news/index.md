# Changelog

## lagdynamics 0.1.0 (in development)

### First public release

`lagdynamics` is a from-scratch, clean-room implementation of lag
sequential analysis for categorical event sequences. The package is
independent of any prior R or non-R LSA implementation; every numerical
method is implemented from primary literature (Bakeman & Quera 1995;
Sackett 1979; Wickens 1989; Christensen 1997; Haberman 1979). See
`inst/REFERENCES.md` for the formula-by-formula citation map.

#### Design

- Unified [`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md)
  constructor with a pluggable engine registry.
- Five built-in engines: `classical`, `two_cell`, `bidirectional`,
  `parallel_dominance`, `nonparallel_dominance`.
- Convenience wrappers:
  [`lsa_classical()`](https://saqr.me/lagdynamics/reference/lsa.md),
  [`lsa_two_cell()`](https://saqr.me/lagdynamics/reference/lsa.md),
  [`lsa_bidirectional()`](https://saqr.me/lagdynamics/reference/lsa.md),
  [`lsa_parallel_dominance()`](https://saqr.me/lagdynamics/reference/lsa.md),
  [`lsa_nonparallel_dominance()`](https://saqr.me/lagdynamics/reference/lsa.md).
- Sequence-level bootstrap
  ([`bootstrap_lsa()`](https://saqr.me/lagdynamics/reference/bootstrap_lsa.md)),
  permutation
  ([`permute_lsa()`](https://saqr.me/lagdynamics/reference/permute_lsa.md)),
  case-drop stability
  ([`stability_lsa()`](https://saqr.me/lagdynamics/reference/stability_lsa.md)),
  and split-half reliability
  ([`reliability_lsa()`](https://saqr.me/lagdynamics/reference/reliability_lsa.md)).
- Multi-group fits via `lsa(data, group = ...)`, returning an
  `lsa_group` object with grouped
  [`transitions()`](https://saqr.me/lagdynamics/reference/transitions.md),
  [`nodes()`](https://saqr.me/lagdynamics/reference/nodes.md),
  [`tests()`](https://saqr.me/lagdynamics/reference/tests.md), and
  [`initial()`](https://saqr.me/lagdynamics/reference/initial.md)
  methods.
- Tidy edge tables and S3 objects of class `c("lsa", "cograph_network")`
  for seamless integration with the `cograph` plotting layer.
- Native transition and initial probabilities:
  [`transition_probabilities()`](https://saqr.me/lagdynamics/reference/transition_probabilities.md)
  returns the row-stochastic P(to \| from) matrix and
  [`initial()`](https://saqr.me/lagdynamics/reference/initial.md) the
  initial-state probabilities.
- Between-group comparison via
  [`compare_lsa()`](https://saqr.me/lagdynamics/reference/compare_lsa.md)
  and
  [`bayes_compare_lsa()`](https://saqr.me/lagdynamics/reference/bayes_compare_lsa.md).
- Reproducibility hooks (`indices=`, `shuffles=`) for bit-identical
  cross-language verification.

#### Experimental

- [`transfer_entropy()`](https://saqr.me/lagdynamics/reference/transfer_entropy.md):
  directed Schreiber transfer entropy for categorical sequences. Two
  modes (a directed state-flow network, and bivariate full-alphabet
  between two series), with effective (surrogate-debiased) and 0-1
  normalised variants, and boundary-safe pooling that never lags across
  sequences. Sign-blind by design; pair with
  [`transitions()`](https://saqr.me/lagdynamics/reference/transitions.md)
  Yule’s Q for direction of effect. Validated in `equivonly/` against
  `infotheo::condinformation`, an independent direct double-sum oracle,
  and exact analytic cases.

#### Bundled data

`ai_long`, `engagement`, `group_regulation`, and `group_regulation_long`
(a long-format event log with a recorded achievement grouping) ship with
the package, so every example, test, and vignette runs without any
external data source.

#### Dependencies

Lean runtime: base R only (`Imports: grDevices, grid, stats, utils`).
The plotting packages (`ggplot2`, `cograph`) are the only soft
`Suggests`, used when present.
