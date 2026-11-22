# Changelog

## lagdynamics 0.35.3

- With structural zeros (`loops = FALSE` or `structural_zeros`) on data
  that hold transitions in the forbidden cells, the quasi-independence
  model is now fitted to the possible cells, with the margins of those
  cells (Bakeman & Quera, 1995, *Psychological Bulletin*; Bishop,
  Fienberg & Holland, 1975): the forbidden counts are left out of the
  expected counts, the adjusted residuals and the tablewise tests. The
  margins used to keep them, so the expected counts did not sum to the
  counts being fitted and the tablewise G^2 could be negative:
  `lsa(engagement, loops = FALSE)` reported G^2 = -1184.633 with p = 1
  and now reports 5.994 on 1 df (p = 0.014), the deviance of the Poisson
  quasi-independence model `count ~ from + to` fitted by
  [`glm()`](https://rdrr.io/r/stats/glm.html) to the six possible cells.
  [`lsa_ipf()`](https://pak.dynasite.org/lagdynamics/reference/lsa_ipf.md)
  likewise reproduces the totals of the possible cells. Tables whose
  forbidden cells are empty are unchanged. The warning
  `lagdynamics_structural_zero_counts` is still raised and now says that
  the counts are left out.
- [`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
  computes the upper tail of the Beta posterior with
  `pbeta(lower.tail = FALSE)` instead of `1 - pbeta()`, which lost
  relative accuracy for masses below about 1e-16 (Average -\> Average in
  `engagement`: 3.14102e-17, now 3.14168e-17). No certainty decision
  changes.

## lagdynamics 0.35.2

- [`print()`](https://rdrr.io/r/base/print.html) of a fit counts
  significant transitions out of the cells that were tested. With
  `loops = FALSE` the self-transitions are structural zeros without a
  p-value; the header said “of 81” on a nine-state table where 72 cells
  are tested, and now says “of 72”.

## lagdynamics 0.35.0

CRAN release: 2026-11-18

This release collects all changes since the CRAN release 0.32.

### Changes to results

- [`summary()`](https://rdrr.io/r/base/summary.html) on a fit now
  returns a `summary.lsa` object instead of printing the fit and its
  matrices. The object holds two tidy tables, a one-row overview and a
  per-state activity table; it prints when shown, and
  `as.data.frame(summary(fit), what = "nodes")` returns either table.
  Grouped fits are unchanged (one row per group).
- [`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md)
  now returns a tidy data frame with columns `from`, `to` and `prob`
  (one row per pair of states; grouped fits gain a leading `group`
  column) instead of a matrix or a list of matrices.
- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on
  [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
  and
  [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
  results gains `what = "global"` for the omnibus test, and on
  [`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md)
  results `what = "summary"` for a one-row summary of the replicates.

### Reliability and stability

- [`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md)
  and
  [`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md)
  now treat only the documented degenerate case (a resample with no
  transitions) as a failed replicate; any other error is reported to the
  caller. The number of failed replicates is stored in the result
  (`n_failed`) and shown by
  [`print()`](https://rdrr.io/r/base/print.html).

### Errors and output

- User-facing errors are classed conditions: `lagdynamics_bad_input`,
  `lagdynamics_needs_events`, `lagdynamics_no_transitions`,
  `lagdynamics_unknown_engine`, `lagdynamics_unsupported_option`,
  `lagdynamics_incompatible_fits`, `lagdynamics_missing_package` and
  `lagdynamics_ipf_infeasible`, each also of class `lagdynamics_error`.
- Console output and plot subtitles use ASCII replacements for the block
  glyph, the superscript in G^2 and the middle dot when the session is
  not UTF-8.
- `transfer_entropy(seed = )` and `bayes_compare_lsa(seed = )` restore
  the caller’s random number stream on exit.

### Documentation

- New documentation structure. Vignettes: lag sequential analysis, lag
  transition networks, evaluating transitions and comparing groups.
  Guides (on the package website): a complete workflow, event logs and
  interoperability, engines, lags and structural zeros, plotting, and
  transfer entropy. Case studies: course engagement, regulation in high
  and low achieving groups, and AI behaviour in programming sessions.

## lagdynamics 0.32

CRAN release: 2026-07-21

### New features

- [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
  now reads inference results, not only fitted models. It accepts the
  output of
  [`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
  [`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md),
  [`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md),
  [`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
  [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md),
  [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
  and
  [`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md),
  and returns the same tidy one-row-per-transition data frame.
  `significant = TRUE` keeps the transitions the method itself flags,
  and `sort = "strength"` orders them by that method’s signed effect.
  Previously the per-transition table of an inference result was
  reachable only through
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  because the print methods show a header summary alone.
- The
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
  generic is now `transitions(fit, ...)`, so each class declares only
  the arguments backed by columns it actually has. Existing calls on a
  fit are unaffected.

### Documentation

- Vignettes restore the user’s
  [`options()`](https://rdrr.io/r/base/options.html) after changing
  `digits`, as CRAN requires.
- Vignettes read every result with
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
  rather than coercing it with
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and
  select rows with the `significant` argument rather than subsetting.
- Removed the citations, reference lists, DOIs, documentation URLs and
  linked author metadata that 0.31 added to the shipped vignettes.

## lagdynamics 0.31

### CRAN candidate

- Prepared the package for first CRAN submission.
- Moved `ggplot2` and `cograph` from Suggests to Imports so the plotting
  surface works out of the box; plotting examples now run during checks.
  The analytical core still depends only on base R.
- Added a dedicated interoperability vignette covering wide data, long
  event logs, `tna`, `Nestimate`, `cograph`, and
  [`lsa_to_tna()`](https://pak.dynasite.org/lagdynamics/reference/lsa_to_tna.md).
- Added linked author metadata and Dynalytics framework links to all
  shipped vignettes.
- Restored
  [`lsa_to_tna()`](https://pak.dynasite.org/lagdynamics/reference/lsa_to_tna.md)
  for handing an `lsa` fit to `tna` tooling.
- Added ingestion of
  [`Nestimate::build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.html)
  netobjects through their prepared sequence data.
- Stored the bundled `engagement` data as a data frame so
  `tna::tna(engagement)` works directly.
- Removed dead package-site URLs from CRAN-visible metadata.
- Added future CRAN installation instructions to the README.
- Included `NEWS.md` in the source package.

## lagdynamics 0.3.0

### Interoperability and documentation

- Added integration tests for `cograph`, `Nestimate`, and the Dynalytics
  evidence surface.
- Added native TNA-style aliases: `weights = "tna"` and
  `weights = "relative"` now map to transition probabilities.
- Updated plotting documentation and vignettes to use `weights = "tna"`
  for probability-weighted transition networks.
- Added and reorganised vignettes:
  - `intro`: conceptual overview and package map.
  - `lagdynamics`: concise quick start.
  - `workflow`: complete applied workflow.
  - `interop`: interoperability with sibling packages.
  - `lag-transition-networks`: transition-network interpretation.
  - `confirmatory`: evidence and uncertainty workflow.
  - `plotting`: plot gallery.
- Removed public documentation references to unexported internals.
- Made long-format input more flexible: `action` is the only mandatory
  long-format column, with optional `actor`, `session`, `time`, and
  `order`.
- Added warnings for single-sequence bootstrap and permutation cases
  where the requested procedure has limited inferential meaning.

## lagdynamics 0.2.0

### Confirmatory workflow and group comparison

- Added the Dynalytics-style confirmatory evidence battery:
  [`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md),
  [`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
  [`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md),
  [`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md),
  and
  [`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md).
- Added group comparison with
  [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
  and Bayesian group comparison with
  [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md).
- Added grouped
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) fits
  through `group = ...`, with grouped methods for
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
  [`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
  [`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md),
  [`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md),
  plotting, reliability, and comparison workflows.
- Added tidy
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  for inference and comparison result objects.
- Added the unified plotting surface: residual heatmaps, residual
  networks, TNA probability networks, chord diagrams, sunbursts,
  uncertainty forests, and group-comparison plots.
- Added native transition and initial probabilities:
  [`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md)
  and
  [`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md).
- Added bundled long-format data for examples and tests.

## lagdynamics 0.1.0

### Initial implementation

- Created a from-scratch, clean-room implementation of lag sequential
  analysis for categorical event sequences.
- Added the unified
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
  constructor and canonical sequence handling through
  [`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md)
  and
  [`lsa_transitions()`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md).
- Added five built-in engines: `classical`, `two_cell`, `bidirectional`,
  `parallel_dominance`, and `nonparallel_dominance`.
- Added convenience wrappers:
  [`lsa_classical()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
  [`lsa_two_cell()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
  [`lsa_bidirectional()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
  [`lsa_parallel_dominance()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
  and
  [`lsa_nonparallel_dominance()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).
- Added the pluggable engine registry:
  [`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md),
  [`get_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/get_lsa_engine.md),
  [`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md),
  and
  [`unregister_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/unregister_lsa_engine.md).
- Added tidy reading verbs:
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
  [`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
  [`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md),
  [`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md),
  and [`summary()`](https://rdrr.io/r/base/summary.html).
- Added multi-lag helpers with
  [`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)
  and
  [`lag_profile()`](https://pak.dynasite.org/lagdynamics/reference/lag_profile.md).
- Added structural-zero handling through `loops = FALSE` and arbitrary
  structural-zero matrices.
- Added experimental
  [`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md)
  for directed categorical information-flow analysis.
- Kept runtime dependencies minimal: only base R packages are imported
  (`grDevices`, `grid`, `stats`, and `utils`).
