# Changelog

## lagdynamics 0.31

### CRAN candidate

- Prepared the package for first CRAN submission.
- Moved `ggplot2` and `cograph` from Suggests to Imports so the plotting
  surface works out of the box; plotting examples now run during checks.
  The analytical core still depends only on base R.
- Added a dedicated interoperability vignette covering wide data, long
  event logs, `tna`, `Nestimate`, `cograph`, and
  [`lsa_to_tna()`](https://saqr.me/lagdynamics/reference/lsa_to_tna.md).
- Added linked author metadata and Dynalytics framework links to all
  shipped vignettes.
- Restored
  [`lsa_to_tna()`](https://saqr.me/lagdynamics/reference/lsa_to_tna.md)
  for handing an `lsa` fit to `tna` tooling.
- Added ingestion of
  [`Nestimate::build_network()`](https://saqr.me/Nestimate/reference/build_network.html)
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
  [`certainty_lsa()`](https://saqr.me/lagdynamics/reference/certainty_lsa.md),
  [`bootstrap_lsa()`](https://saqr.me/lagdynamics/reference/bootstrap_lsa.md),
  [`reliability_lsa()`](https://saqr.me/lagdynamics/reference/reliability_lsa.md),
  [`stability_lsa()`](https://saqr.me/lagdynamics/reference/stability_lsa.md),
  and
  [`permute_lsa()`](https://saqr.me/lagdynamics/reference/permute_lsa.md).
- Added group comparison with
  [`compare_lsa()`](https://saqr.me/lagdynamics/reference/compare_lsa.md)
  and Bayesian group comparison with
  [`bayes_compare_lsa()`](https://saqr.me/lagdynamics/reference/bayes_compare_lsa.md).
- Added grouped [`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md)
  fits through `group = ...`, with grouped methods for
  [`transitions()`](https://saqr.me/lagdynamics/reference/transitions.md),
  [`nodes()`](https://saqr.me/lagdynamics/reference/nodes.md),
  [`tests()`](https://saqr.me/lagdynamics/reference/tests.md),
  [`initial()`](https://saqr.me/lagdynamics/reference/initial.md),
  plotting, reliability, and comparison workflows.
- Added tidy
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  for inference and comparison result objects.
- Added the unified plotting surface: residual heatmaps, residual
  networks, TNA probability networks, chord diagrams, sunbursts,
  uncertainty forests, and group-comparison plots.
- Added native transition and initial probabilities:
  [`transition_probabilities()`](https://saqr.me/lagdynamics/reference/transition_probabilities.md)
  and [`initial()`](https://saqr.me/lagdynamics/reference/initial.md).
- Added bundled long-format data for examples and tests.

## lagdynamics 0.1.0

### Initial implementation

- Created a from-scratch, clean-room implementation of lag sequential
  analysis for categorical event sequences.
- Added the unified
  [`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md) constructor
  and canonical sequence handling through
  [`lsa_data()`](https://saqr.me/lagdynamics/reference/lsa_data.md) and
  [`lsa_transitions()`](https://saqr.me/lagdynamics/reference/lsa_transitions.md).
- Added five built-in engines: `classical`, `two_cell`, `bidirectional`,
  `parallel_dominance`, and `nonparallel_dominance`.
- Added convenience wrappers:
  [`lsa_classical()`](https://saqr.me/lagdynamics/reference/lsa.md),
  [`lsa_two_cell()`](https://saqr.me/lagdynamics/reference/lsa.md),
  [`lsa_bidirectional()`](https://saqr.me/lagdynamics/reference/lsa.md),
  [`lsa_parallel_dominance()`](https://saqr.me/lagdynamics/reference/lsa.md),
  and
  [`lsa_nonparallel_dominance()`](https://saqr.me/lagdynamics/reference/lsa.md).
- Added the pluggable engine registry:
  [`register_lsa_engine()`](https://saqr.me/lagdynamics/reference/register_lsa_engine.md),
  [`get_lsa_engine()`](https://saqr.me/lagdynamics/reference/get_lsa_engine.md),
  [`list_lsa_engines()`](https://saqr.me/lagdynamics/reference/list_lsa_engines.md),
  and
  [`unregister_lsa_engine()`](https://saqr.me/lagdynamics/reference/unregister_lsa_engine.md).
- Added tidy reading verbs:
  [`transitions()`](https://saqr.me/lagdynamics/reference/transitions.md),
  [`nodes()`](https://saqr.me/lagdynamics/reference/nodes.md),
  [`tests()`](https://saqr.me/lagdynamics/reference/tests.md),
  [`initial()`](https://saqr.me/lagdynamics/reference/initial.md), and
  [`summary()`](https://rdrr.io/r/base/summary.html).
- Added multi-lag helpers with
  [`lsa_lags()`](https://saqr.me/lagdynamics/reference/lsa_lags.md) and
  [`lag_profile()`](https://saqr.me/lagdynamics/reference/lag_profile.md).
- Added structural-zero handling through `loops = FALSE` and arbitrary
  structural-zero matrices.
- Added experimental
  [`transfer_entropy()`](https://saqr.me/lagdynamics/reference/transfer_entropy.md)
  for directed categorical information-flow analysis.
- Kept runtime dependencies minimal: only base R packages are imported
  (`grDevices`, `grid`, `stats`, and `utils`).
