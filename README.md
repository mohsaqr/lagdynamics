# lagdynamics

> Modern, tidy lag sequential analysis for categorical event sequences.

<!-- badges: start -->
[![R-CMD-check](https://github.com/mohsaqr/lagdynamics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mohsaqr/lagdynamics/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![r-universe](https://mohsaqr.r-universe.dev/badges/lagdynamics)](https://mohsaqr.r-universe.dev/lagdynamics)
<!-- badges: end -->

Lag sequential analysis (LSA) is a statistical method for detecting temporal
contingencies in categorical, time-ordered data. For every ordered pair of
states, it tests whether one follows another at a given lag more — or less —
often than expected under sequential independence, quantifying each contingency
with the adjusted (standardized) residual of observed against expected
transition frequencies. It is widely used to study how categorical processes
unfold over time — in behavioral and dyadic interaction, conversation and group
discourse, psychotherapy process research, human–computer interaction, and
learning analytics — wherever the order of events, not merely their frequency,
carries the meaning.

`lagdynamics` implements LSA as a modern, tidy workflow. A single `lsa()`
constructor fits the classical and extended engine family — classical,
two-cell, bidirectional, and parallel / non-parallel dominance — and every
result is read through a verb that returns a one-row-per-observation
`data.frame`. Its distinguishing commitment is evidence: each transition is
treated as a *tested departure from independence*, backed by a confirmatory
battery of bootstrap intervals, analytic Bayesian certainty, split-half
reliability, case-drop stability, and permutation tests.

Beyond the classical test, the package adds multi-lag analysis, structural-zero
constraints under quasi-independence, and grouped fits that estimate one model
per group in a shared state space, with formal group comparison by permutation
or Bayesian Dirichlet–Multinomial contrast. It is designed to interoperate
across the sequence-analysis ecosystem: `lsa()` ingests long event logs, wide
matrices, lists of sequences, and objects from `tna`, `Nestimate`, and
`TraMineR`, while fitted objects carry a `cograph` network class and expose
transition and initial probabilities for downstream tooling. A single `plot()`
verb renders any fit as a residual heatmap, transition network, chord diagram,
polar sunburst, or uncertainty forest.

## Installation

```r
# CRAN (after acceptance)
install.packages("lagdynamics")

# r-universe (pre-built binaries, no compiler needed)
install.packages("lagdynamics",
                 repos = c("https://mohsaqr.r-universe.dev",
                           "https://cloud.r-project.org"))

# development version
# install.packages("remotes")
remotes::install_github("mohsaqr/lagdynamics")
```

The analytical core depends only on base R; plotting adds `ggplot2` and
`cograph`, which install automatically. Interoperability packages (`tna`,
`Nestimate`) stay optional `Suggests`.

## At a glance

| Task | Functions |
|---|---|
| Fit a model | `lsa()`, `lsa_lags()`, `lag_profile()` |
| Read results (tidy) | `transitions()`, `nodes()`, `tests()`, `initial()`, `as.data.frame()` |
| Weigh the evidence | `bootstrap_lsa()`, `certainty_lsa()`, `permute_lsa()`, `stability_lsa()`, `reliability_lsa()` |
| Compare groups | `compare_lsa()`, `bayes_compare_lsa()` |
| Plot | `plot(fit, type = )`, `plot_transitions()`, `plot_chords()`, `plot_polar()`, `plot_forest()` |
| Probabilities | `transition_probabilities()`, `initial()` |
| Engines | `lsa(engine = )`, `register_lsa_engine()`, `list_lsa_engines()` |

## Example

```r
library(lagdynamics)

seq <- c("Question", "Explain", "Agree",
         "Question", "Explain", "Elaborate",
         "Agree", "Question", "Explain")

fit <- lsa(seq, engine = "classical")
fit

transitions(fit, significant = TRUE)   # tidy table of tested transitions
plot(fit, type = "network")            # residual transition network

# quantify the evidence behind each edge
boot <- bootstrap_lsa(fit, R = 1000)
plot(boot)                             # circular bootstrap CI forest
```

## Vignettes

| Vignette | Topic |
|---|---|
| `vignette("lagdynamics")` | Get started: the method, why lagdynamics, and a hands-on tour |
| `vignette("interop")` | Interoperability with wide data, long logs, `tna`, `Nestimate`, and `cograph` |
| `vignette("workflow")` | A complete analysis from sequences to a group comparison |
| `vignette("confirmatory")` | The confirmatory testing battery: matching claims to evidence |
| `vignette("lag-transition-networks")` | Lag transition networks |
| `vignette("plotting")` | The full plotting gallery |

## License

MIT © 2026 Mohammed Saqr
