# lagdynamics

> Modern, tidy lag sequential analysis for categorical event sequences.

Lag sequential analysis (LSA) is a statistical method for detecting
temporal contingencies in categorical, time-ordered data. For every
ordered pair of states, it tests whether one follows another at a given
lag more — or less — often than expected under sequential independence,
quantifying each contingency with the adjusted (standardized) residual
of observed against expected transition frequencies. It is widely used
to study how categorical processes unfold over time — in behavioral and
dyadic interaction, conversation and group discourse, psychotherapy
process research, human–computer interaction, and learning analytics —
wherever the order of events, not merely their frequency, carries the
meaning.

`lagdynamics` implements LSA as a modern, tidy workflow. A single
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
constructor fits the classical and extended engine family — classical,
two-cell, bidirectional, and parallel / non-parallel dominance — and
every result is read through a verb that returns a
one-row-per-observation `data.frame`. Its distinguishing commitment is
evidence: each transition is treated as a *tested departure from
independence*, backed by a confirmatory battery of bootstrap intervals,
analytic Bayesian certainty, split-half reliability, case-drop
stability, and permutation tests.

Beyond the classical test, the package adds multi-lag analysis,
structural-zero constraints under quasi-independence, and grouped fits
that estimate one model per group in a shared state space, with formal
group comparison by permutation or Bayesian Dirichlet–Multinomial
contrast. It is designed to interoperate across the sequence-analysis
ecosystem:
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) ingests
long event logs, wide matrices, lists of sequences, and objects from
`tna`, `Nestimate`, and `TraMineR`, while fitted objects carry a
`cograph` network class and expose transition and initial probabilities
for downstream tooling. A single
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) verb renders
any fit as a residual heatmap, transition network, chord diagram, polar
sunburst, or uncertainty forest.

## Installation

``` r

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
`cograph`, which install automatically. Interoperability packages
(`tna`, `Nestimate`) stay optional `Suggests`.

## At a glance

| Task | Functions |
|----|----|
| Fit a model | [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), [`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md), [`lag_profile()`](https://pak.dynasite.org/lagdynamics/reference/lag_profile.md) |
| Read results (tidy) | [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md), [`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md), [`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md), [`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md), [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) |
| Weigh the evidence | [`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md), [`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md), [`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md), [`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md), [`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md) |
| Compare groups | [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md), [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md) |
| Plot | `plot(fit, type = )`, [`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md), [`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md), [`plot_polar()`](https://pak.dynasite.org/lagdynamics/reference/plot_polar.md), [`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md) |
| Probabilities | [`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md), [`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md) |
| Engines | `lsa(engine = )`, [`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md), [`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md) |

## Example

``` r

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
|----|----|
| [`vignette("lagdynamics")`](https://pak.dynasite.org/lagdynamics/articles/lagdynamics.md) | Get started: the method, why lagdynamics, and a hands-on tour |
| [`vignette("interop")`](https://pak.dynasite.org/lagdynamics/articles/interop.md) | Interoperability with wide data, long logs, `tna`, `Nestimate`, and `cograph` |
| [`vignette("workflow")`](https://pak.dynasite.org/lagdynamics/articles/workflow.md) | A complete analysis from sequences to a group comparison |
| [`vignette("confirmatory")`](https://pak.dynasite.org/lagdynamics/articles/confirmatory.md) | The confirmatory testing battery: matching claims to evidence |
| [`vignette("lag-transition-networks")`](https://pak.dynasite.org/lagdynamics/articles/lag-transition-networks.md) | Lag transition networks |
| [`vignette("plotting")`](https://pak.dynasite.org/lagdynamics/articles/plotting.md) | The full plotting gallery |

## License

MIT © 2026 Mohammed Saqr
