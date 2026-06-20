# lagseq

> Modern, tidy lag sequential analysis for categorical event sequences.

<!-- badges: start -->
[![r-universe status](https://mohsaqr.r-universe.dev/badges/lagseq)](https://mohsaqr.r-universe.dev/lagseq)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`lagseq` provides a modern, tidy, pipe-friendly interface for lag
sequential analysis (LSA). A single `lsa()` constructor with a pluggable
engine registry exposes the classical and extended LSA family — classical,
two-cell, bidirectional, parallel-dominance, and non-parallel-dominance —
and every result is read through a verb that returns a tidy
one-row-per-observation `data.frame`.

It is the lag-sequential member of the **Dynalytics** framework: every
edge is a *tested departure from independence*, backed by a confirmatory
testing battery — bootstrap and analytic certainty for edges, split-half
reliability for the whole network, case-drop stability, permutation tests,
and permutation- or Bayesian-based group comparison. Fits visualize
through a single `plot()` verb and interoperate with the `tna` and
`Nestimate` ecosystems at both ends.

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
  (`ggplot2`, `cograph`, `tna`, `Nestimate`, `TraMineR`) are
  soft `Suggests`, used only when present.

## Installation

Install the latest build from [r-universe](https://mohsaqr.r-universe.dev/lagseq)
(recommended — pre-built, no compiler required):

```r
install.packages("lagseq", repos = "https://mohsaqr.r-universe.dev")
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("mohsaqr/lagseq")
```

`lagseq` needs only base R at runtime (`grid`, `stats`, `utils`). The
plotting and interop packages (`ggplot2`, `cograph`, `tna`, `Nestimate`,
`TraMineR`) are optional `Suggests` — install only the ones you need.

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

# Confirmatory testing — quantify the evidence behind each claim
boot <- bootstrap_lsa(fit, R = 1000)      # edge uncertainty (resampling)
cert <- certainty_lsa(fit)                # edge uncertainty (analytic, Bayesian)
perm <- permute_lsa(fit, R = 1000)        # more than chance
stab <- stability_lsa(fit, R = 500)       # significant edges under case-dropping
rel  <- reliability_lsa(fit, R = 100)     # whole-network split-half reliability

# Groups and comparison
g   <- lsa(sequences, group = group_labels)  # one fit per group -> lsa_group
cmp <- compare_lsa(g)                         # permutation test of the difference
bc  <- bayes_compare_lsa(g)                   # Bayesian group comparison

# Every result is tidy: as.data.frame() returns a one-row-per-edge frame
as.data.frame(cmp)

# Plotting — one verb, pick the view with `type`
plot(fit)                                 # residual heatmap (default)
plot(fit, type = "network")               # residual network (blue = more than chance)
plot(fit, type = "network", weights = "prob")  # transition network (a TNA model)
plot(fit, type = "chord")                 # chord diagram
plot(fit, type = "sunburst")              # polar sunburst
plot(bootstrap_lsa(fit))                  # circular bootstrap CI forest
plot(cmp)                                 # back-to-back group-comparison barrel

# Interop — convert a fit to the tna ecosystem's native object
# (lagseq converts; tna does the analysis).
net   <- lsa_to_tna(fit, weights = "prob")  # -> tna object
cents <- tna::centralities(net)
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
call to the `tna` object those toolkits analyse (centralities, pruning,
communities, …):

```r
lsa_to_tna(fit, weights = "prob")    # -> `tna` object (also works on a grouped fit)
```

`tna`, `Nestimate`, and `TraMineR` stay optional (`Suggests`);
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

## Vignettes

| Vignette | Topic |
|---|---|
| `vignette("lagseq")` | Get started: the method, why lagseq, and a hands-on tour |
| `vignette("workflow")` | A complete analysis from sequences to a group comparison |
| `vignette("confirmatory")` | The confirmatory testing battery: matching claims to evidence |
| `vignette("interop")` | Interoperability with `tna` and `Nestimate` objects |
| `vignette("plotting")` | The full plotting gallery |

## Status

v0.1.0.

**Implemented:** classical / two-cell / bidirectional / parallel-
and non-parallel-dominance engines; multi-lag analysis (`lsa_lags()`,
`lag_profile()`) and structural-zero handling (`loops = FALSE` or an
explicit 0/1 matrix, IPF with rank-based quasi-independence df).
Confirmatory battery: bootstrap (sequence-level + stationary block) and
analytic Dirichlet-Multinomial certainty (`certainty_lsa()`) for edges,
permutation (`permute_lsa()`), case-drop stability, split-half
reliability, and group comparison by permutation (`compare_lsa()`, two-
group and all-pairwise) or Bayesian Dirichlet-Multinomial
(`bayes_compare_lsa()`). Tidy reading: `transitions()` (with
`significant` / `direction` / `min_count` selectors), `nodes()`,
`tests()`, `initial()`, and `as.data.frame()` on every result object.
Plotting via one `plot(fit, type = )` verb — heatmap, residual network,
TNA transition network, chord, sunburst, bootstrap/certainty forests, and
group-comparison barrels — with grouped fits drawn one panel per group.
Interoperability: ingests `tna` / `Nestimate` / `TraMineR` objects and
converts out via `lsa_to_tna()`. Multi-group fits
(`lsa(group = )` → `lsa_group`) are supported across every reading,
plotting, inference, and bridge layer.

**Roadmap:** stationarity tests (`stationarity_lsa()`).

## License

MIT © 2026 Mohammed Saqr
