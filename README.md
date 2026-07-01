# lagdynamics

> Modern, tidy lag sequential analysis for categorical event sequences.

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![r-universe](https://mohsaqr.r-universe.dev/badges/lagdynamics)](https://mohsaqr.r-universe.dev/lagdynamics)
[![Docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mohsaqr.github.io/lagdynamics/)
<!-- badges: end -->

`lagdynamics` provides a modern, tidy, pipe-friendly interface for lag
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
through a single `plot()` verb and expose their transition and
initial probabilities as tidy objects for downstream network tooling.

## Design principles

- **Clean room.** Every numerical method is implemented from primary
  literature (Bakeman & Quera 1995; Sackett 1979; Wickens 1989;
  Christensen 1997; Haberman 1979). No code is derived from any prior R
  implementation. See `inst/REFERENCES.md` for the formula-by-formula
  citation map.
- **Best-in-class equivalence.** Engines are cross-validated against
  independent hand-formula identities and base-R primitives:
  `stats::loglin()` (for iterative proportional fitting of expected
  frequencies with structural zeros), `chisq.test()$stdres` (for
  Haberman adjusted residuals when no structural zeros are present),
  `pchisq()`, `pnorm()`, and `binom.test()`.
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
  (`grDevices`, `grid`, `stats`, `utils`). The plotting packages
  (`ggplot2`, `cograph`) are soft `Suggests`, used only when present.

## Installation

From [r-universe](https://mohsaqr.r-universe.dev/lagdynamics) (pre-built
binaries, no compiler needed):

```r
install.packages("lagdynamics",
                 repos = c("https://mohsaqr.r-universe.dev",
                           "https://cloud.r-project.org"))
```

Or the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("mohsaqr/lagdynamics")
```

`lagdynamics` needs only base R at runtime (`grDevices`, `grid`, `stats`,
`utils`). The plotting packages (`ggplot2`, `cograph`) are optional
`Suggests` — install only the ones you need.

## Quick start

```r
library(lagdynamics)

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
plot(fit, type = "network", weights = "prob")  # transition (TNA-style) network
plot(fit, type = "chord")                 # chord diagram
plot(fit, type = "sunburst")              # polar sunburst
plot(bootstrap_lsa(fit))                  # circular bootstrap CI forest
plot(cmp)                                 # back-to-back group-comparison barrel

# Transition and initial probabilities as native objects
transition_probabilities(fit)             # row-stochastic P(to | from) matrix
initial(fit)                              # initial-state probabilities (tidy)
```

## Transition and initial probabilities

`lsa()` computes the quantities a Transition Network Analysis reads and
exposes them natively — no other package required:

```r
transition_probabilities(fit)    # row-stochastic P(to | from) matrix
initial(fit)                     # initial-state probabilities (tidy data.frame)
```

`lsa()` also reads sequences straight out of common inputs, so a fit
drops into an existing pipeline without reshaping:

```r
lsa(long_log, actor = , action = , time = )   # a raw long event log
lsa(wide_matrix)                               # rows = sequences, cols = time
lsa(list_of_sequences)                         # one vector per sequence
```

It additionally recognises sequence-bearing objects from sibling
state-sequence and transition-model packages when they happen to be
installed; none of them is a dependency of `lagdynamics`.

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
| `vignette("lagdynamics")` | Get started: the method, why lagdynamics, and a hands-on tour |
| `vignette("workflow")` | A complete analysis from sequences to a group comparison |
| `vignette("confirmatory")` | The confirmatory testing battery: matching claims to evidence |
| `vignette("lag-transition-networks")` | Lag transition networks |
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
TNA-style transition network, chord, sunburst, bootstrap/certainty
forests, and group-comparison barrels — with grouped fits drawn one panel
per group. Transition and initial probabilities are exposed natively via
`transition_probabilities()` and `initial()`. Multi-group fits
(`lsa(group = )` → `lsa_group`) are supported across every reading,
plotting, and inference layer.

**Roadmap:** stationarity tests (`stationarity_lsa()`).

## License

MIT © 2026 Mohammed Saqr
