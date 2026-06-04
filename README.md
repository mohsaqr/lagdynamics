# lagseq

> Modern, tidy lag sequential analysis for categorical event sequences.

`lagseq` provides a unified, pipe-friendly interface for lag sequential
analysis (LSA). A single `lsa()`
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
- **Minimal-dependency policy.** Runtime imports are base packages only
  (`grid`, `stats`, `utils`). The plotting and interop packages
  (`ggplot2`, `cograph`, `tna`, `Nestimate`, `igraph`, `TraMineR`) are
  soft `Suggests`, used only when present.

## Quick start

```r
library(lagseq)

seq <- c("Question", "Explain", "Agree",
         "Question", "Explain", "Elaborate",
         "Agree", "Question", "Explain")

fit <- lsa(seq, lag = 1, engine = "classical")
fit

# Reading the fit — one verb per result, each returns a tidy data.frame
transitions(fit)                          # one row per transition
transitions(fit, significant = TRUE)      # significant only
transitions(fit, direction = "over")      # over-represented
transitions(fit, min_count = 2)           # frequently observed
nodes(fit)                                # one row per state
tests(fit)                                # tablewise independence tests
initial(fit)                              # initial-state distribution

# Inference
boot <- bootstrap_lsa(fit, R = 1000)
perm <- permute_lsa(fit, R = 1000)
stab <- stability_lsa(fit, R = 500)
rel  <- reliability_lsa(fit, R = 100)

# Plotting — one verb, pick the view with `type`
plot(fit)                                 # residual heatmap (default)
plot(fit, type = "network")               # transition network
plot(fit, type = "chord")                 # chord diagram
plot(fit, type = "sunburst")              # polar sunburst
plot(bootstrap_lsa(fit))                  # circular bootstrap CI forest

# Interop — convert a fit to another toolkit's native object
# (lagseq converts; the downstream package does the analysis).
net   <- lsa_to_tna(fit, weights = "prob")  # -> tna object
cents <- tna::centralities(net)
g     <- igraph::as.igraph(fit)             # -> igraph object
```

## Compatibility with `tna` and `Nestimate`

`lagseq` is a converter, not a competing analyser: it interoperates with
the `tna` and `Nestimate` ecosystems at both ends, and keeps **zero
exported-name overlap** with either.

**Input** — `lsa()` accepts the sequence objects those packages already
use, so a fit drops straight into an existing pipeline:

```r
lsa(tna_object)                  # a `tna` model / sequence object
lsa(nestimate_data)              # a `Nestimate` `nestimate_data` object
lsa(stslist)                     # a `TraMineR` state-sequence object
lsa(long_log, actor = , action = , time = )   # or a raw long event log
```

**Output** — a fit is a directed weighted network, so it converts in one
call to the native object those toolkits analyse (centralities, pruning,
communities, …):

```r
lsa_to_tna(fit, weights = "prob")    # -> `tna` object (also works on a grouped fit)
igraph::as.igraph(fit)               # -> `igraph` object
```

`tna`, `Nestimate`, `igraph`, and `TraMineR` stay optional (`Suggests`);
the converters and ingestors error informatively only if the relevant
package is missing.

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

v0.1.0.

**Implemented:** classical / two-cell / bidirectional / parallel-
and non-parallel-dominance engines; structural-zero handling via IPF
with rank-based quasi-independence df (classical engine); bootstrap
(sequence-level + stationary block), permutation, case-drop stability,
and split-half reliability inference; a single tidy `transitions()` verb
(with `significant` / `direction` / `min_count` selectors) plus
`nodes()`, `tests()`, and `initial()`; plotting via one `plot(fit, type = )`
verb — heatmap, network, chord, sunburst, and a circular bootstrap
forest, with grouped fits drawn one panel per group; and the `tna` /
`igraph` bridge (`lsa_to_tna()`, `as.igraph()`). Multi-group fits
(`lsa(group = )` → `lsa_group`) are supported across every reading,
plotting, inference, and bridge layer.

**Roadmap:** between-group comparison (`compare_lsa()`) and
stationarity tests (`stationarity_lsa()`).

## License

MIT © 2026 Mohammed Saqr
