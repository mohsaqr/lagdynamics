# lagseq

> Modern, tidy lag sequential analysis for categorical event sequences.

`lagseq` provides a unified, pipe-friendly interface for lag sequential
analysis (LSA), in the style of [Nestimate](https://github.com/mohsaqr/Nestimate)
and [cograph](https://github.com/sonsoleslp/cograph). A single `lsa()`
constructor with a pluggable engine registry exposes the classical and
extended LSA family — classical, two-cell, bidirectional,
parallel-dominance, and non-parallel-dominance — and returns a tidy edge
table ready for transition-network visualization.

## Design principles

- **Clean room.** Every numerical method is implemented from primary
  literature (Bakeman & Quera 1995; Sackett 1979; Wickens 1989;
  Christensen 1997; Haberman 1979). No code is derived from any prior R
  implementation. See `inst/REFERENCES.md` for the formula-by-formula
  citation map.
- **Best-in-class equivalence.** Engines are cross-validated against
  multiple independent oracles: published worked examples from Bakeman &
  Quera (1995), Wampold (1982, 1984), and Sackett (1979); plus base-R
  primitives `stats::loglin()` (for iterative proportional fitting of
  expected frequencies with structural zeros), `chisq.test()$stdres`
  (for Haberman adjusted residuals when no structural zeros are
  present), and `pchisq()` (for likelihood-ratio tests). See
  `inst/EQUIVALENCE.md` for the oracle table and pass rates.
- **Stable S3 class.** `lsa()` returns an object of class
  `c("lsa", "cograph_network")` with named slots
  (`$obs`, `$exp`, `$prob`, `$adj_res`, `$p`, `$yules_q`, `$kappa`,
  `$edges`, `$nodes`, `$weights`, `$params`, `$meta`). Inheriting the
  `cograph_network` class means `cograph::splot()` works on a fit
  unchanged.
- **Recipe pattern.** Configuration is snapshotted on the fit in
  `$params`. Bootstrap, permutation, and stability inference all read
  from that single snapshot to prevent config drift.
- **Minimal-dependency policy.** Runtime imports are base R only
  (`stats`, `utils`). `cograph` is a soft, `Suggests`-level dependency
  used only when present for visualization.

## Quick start

```r
library(lagseq)

seq <- c("Question", "Explain", "Agree",
         "Question", "Explain", "Elaborate",
         "Agree", "Question", "Explain")

fit <- lsa(seq, lag = 1, engine = "classical")
fit

# Tidy outputs
fit$edges                   # one row per transition, with residuals + p
fit$nodes                   # one row per state, with in/out totals

# Filter helpers
significant_transitions(fit, alpha = 0.05)
overrepresented_transitions(fit)
underrepresented_transitions(fit)
common_transitions(fit, min_count = 2)

# Inference
boot <- bootstrap_lsa(fit, R = 1000)
perm <- permute_lsa(fit, R = 1000)
stab <- stability_lsa(fit, R = 500)
rel  <- reliability_lsa(fit, R = 100)

# TNA / igraph interop (require `tna` / `igraph` to be installed).
# lagseq converts; the downstream package does the analysis.
net   <- lsa_to_tna(fit, weights = "prob")
cents <- tna::centralities(net)
g     <- igraph::as.igraph(fit)

# Visualization
# lsa fits inherit class "cograph_network", so the cograph package
# can render them once installed. Native plot methods are on the
# roadmap (see Status).
```

## Engines

| Engine | Method | Reference |
|---|---|---|
| `classical` | Lag-1 sequential with adjusted residuals, Yule's Q, and unidirectional kappa | Bakeman & Quera (1995) |
| `two_cell` | 2×2 single-transition test (Yule's Q + odds ratio) | Bakeman & Gottman (1997) ch. 7 |
| `bidirectional` | Matched-pair / mutual-influence test | Sackett (1979) |
| `parallel_dominance` | Parallel-dominance analysis | Sackett (1979) |
| `nonparallel_dominance` | Non-parallel-dominance analysis | Sackett (1979) |

Users can register custom engines via `register_lsa_engine()`.

## Status

v0.1.0 — under active development. See `lsa_plan.md` for the build
roadmap and `HANDOFF.md` for the current session state.

**Implemented:** classical / two-cell / bidirectional / parallel-
and non-parallel-dominance engines; structural-zero handling via IPF
with rank-based quasi-independence df (classical engine); bootstrap
(sequence-level + stationary block), permutation, case-drop stability,
and split-half reliability inference; significance / over- / under- /
common-transition filter helpers; TNA / igraph bridge
(`lsa_to_tna.lsa()`, `as.igraph.lsa()`) that converts a fit into the native
object of the broader `tna` / `igraph` toolkits (centralities, pruning,
communities, etc.) without making either a hard dependency. Multi-group
fits (`lsa(group = )` → `lsa_group`) are supported across the filter,
reliability, and bridge layers.

**Roadmap:** native plot methods (`plot.lsa()`), between-group
comparison (`compare_lsa()`, `group_lsa()`), stationarity tests
(`stationarity_lsa()`).

## License

MIT © 2026 Mohammed Saqr
