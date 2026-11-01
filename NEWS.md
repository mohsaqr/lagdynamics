# lagdynamics 0.33

## Documentation

- Every exported S3 method now has a help page with a `Value` section:
  `print()`, `summary()`, `as.data.frame()` and `nobs()` for fits (new
  topic `?"lsa-methods"`), and the `print()` / `as.data.frame()` methods
  of every inference, comparison, lag and data object, each documented
  on its constructor's page.
- `?lsa` documents the engine wrappers (`lsa_classical()` and siblings),
  adds method references, and warns that combining `actor` and `session`
  overflows R's integer range once distinct actors x distinct sessions
  exceeds about 2.1e9, silently dropping sequences; the fix (pass a
  globally unique session id as `actor`) is given there and in the
  interoperability vignette.
- Examples read results with `transitions()`, `summary()` and the print
  methods instead of reaching into objects with `$`, and group examples
  use the recorded `Achiever` column of `group_regulation_long` instead
  of a derived grouping vector. Resampling examples now run a small
  unwrapped version before the `\donttest{}` full-size call.
- Corrected return-value descriptions: the column names of
  `bootstrap_lsa()` edges, the conditional columns of
  `transfer_entropy()`, the engine-function contract of
  `register_lsa_engine()`, and the size of the `engagement` data (136
  non-empty sequences). Details that were attached to the wrong roxygen
  tag (`permute_lsa()`, `compare_lsa()`, `reliability_lsa()`) now sit in
  Details.
- README: CRAN install line, a bundled-data example, and the missing
  `intro` vignette. Vignettes gain interpretation of the printed results
  in the quick start, the permutation-null section of the confirmatory
  vignette and the comparison step of the workflow, and correct several
  stale statements (engine count, `transition_probabilities()` returns a
  matrix, state capitalization, time-gap splitting of long logs).

# lagdynamics 0.32

## New features

- `transitions()` now reads inference results, not only fitted models. It
  accepts the output of `bootstrap_lsa()`, `certainty_lsa()`,
  `stability_lsa()`, `permute_lsa()`, `compare_lsa()`,
  `bayes_compare_lsa()` and `lsa_lags()`, and returns the same tidy
  one-row-per-transition data frame. `significant = TRUE` keeps the
  transitions the method itself flags, and `sort = "strength"` orders them
  by that method's signed effect. Previously the per-transition table of an
  inference result was reachable only through `as.data.frame()`, because the
  print methods show a header summary alone.
- The `transitions()` generic is now `transitions(fit, ...)`, so each class
  declares only the arguments backed by columns it actually has. Existing
  calls on a fit are unaffected.

## Documentation

- Vignettes restore the user's `options()` after changing `digits`, as CRAN
  requires.
- Vignettes read every result with `transitions()` rather than coercing it
  with `as.data.frame()`, and select rows with the `significant` argument
  rather than subsetting.
- Removed the citations, reference lists, DOIs, documentation URLs and
  linked author metadata that 0.31 added to the shipped vignettes.

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
