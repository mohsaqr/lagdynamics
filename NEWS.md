# lagdynamics 0.31

## CRAN candidate

- Prepared the package for first CRAN submission.
- Moved `ggplot2` and `cograph` from Suggests to Imports so the plotting
  surface works out of the box; plotting examples now run during checks.
  The analytical core still depends only on base R.
- Added a dedicated interoperability vignette covering wide data, long
  event logs, `tna`, `Nestimate`, `cograph`, and `lsa_to_tna()`.
- Added linked author metadata and Dynalytics framework links to all
  shipped vignettes.
- Restored `lsa_to_tna()` for handing an `lsa` fit to `tna` tooling.
- Added ingestion of `Nestimate::build_network()` netobjects through
  their prepared sequence data.
- Stored the bundled `engagement` data as a data frame so
  `tna::tna(engagement)` works directly.
- Removed dead package-site URLs from CRAN-visible metadata.
- Added future CRAN installation instructions to the README.
- Included `NEWS.md` in the source package.

# lagdynamics 0.3.0

## Interoperability and documentation

- Added integration tests for `cograph`, `Nestimate`, and the Dynalytics
  evidence surface.
- Added native TNA-style aliases: `weights = "tna"` and
  `weights = "relative"` now map to transition probabilities.
- Updated plotting documentation and vignettes to use
  `weights = "tna"` for probability-weighted transition networks.
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

# lagdynamics 0.2.0

## Confirmatory workflow and group comparison

- Added the Dynalytics-style confirmatory evidence battery:
  `certainty_lsa()`, `bootstrap_lsa()`, `reliability_lsa()`,
  `stability_lsa()`, and `permute_lsa()`.
- Added group comparison with `compare_lsa()` and Bayesian group
  comparison with `bayes_compare_lsa()`.
- Added grouped `lsa()` fits through `group = ...`, with grouped
  methods for `transitions()`, `nodes()`, `tests()`, `initial()`,
  plotting, reliability, and comparison workflows.
- Added tidy `as.data.frame()` methods for inference and comparison
  result objects.
- Added the unified plotting surface:
  residual heatmaps, residual networks, TNA probability networks, chord
  diagrams, sunbursts, uncertainty forests, and group-comparison plots.
- Added native transition and initial probabilities:
  `transition_probabilities()` and `initial()`.
- Added bundled long-format data for examples and tests.

# lagdynamics 0.1.0

## Initial implementation

- Created a from-scratch, clean-room implementation of lag sequential
  analysis for categorical event sequences.
- Added the unified `lsa()` constructor and canonical sequence handling
  through `lsa_data()` and `lsa_transitions()`.
- Added five built-in engines:
  `classical`, `two_cell`, `bidirectional`, `parallel_dominance`, and
  `nonparallel_dominance`.
- Added convenience wrappers:
  `lsa_classical()`, `lsa_two_cell()`, `lsa_bidirectional()`,
  `lsa_parallel_dominance()`, and `lsa_nonparallel_dominance()`.
- Added the pluggable engine registry:
  `register_lsa_engine()`, `get_lsa_engine()`, `list_lsa_engines()`,
  and `unregister_lsa_engine()`.
- Added tidy reading verbs:
  `transitions()`, `nodes()`, `tests()`, `initial()`, and `summary()`.
- Added multi-lag helpers with `lsa_lags()` and `lag_profile()`.
- Added structural-zero handling through `loops = FALSE` and arbitrary
  structural-zero matrices.
- Added experimental `transfer_entropy()` for directed categorical
  information-flow analysis.
- Kept runtime dependencies minimal: only base R packages are imported
  (`grDevices`, `grid`, `stats`, and `utils`).
