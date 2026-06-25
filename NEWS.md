# lagdynamics 0.1.0 (in development)

## First public release

`lagdynamics` is a from-scratch, clean-room implementation of lag sequential
analysis for categorical event sequences. The package is independent
of any prior R or non-R LSA implementation; every numerical method is
implemented from primary literature (Bakeman & Quera 1995; Sackett
1979; Wickens 1989; Christensen 1997; Haberman 1979). See
`inst/REFERENCES.md` for the formula-by-formula citation map and
`inst/EQUIVALENCE.md` for the oracle table and pass rates.

### Design

- Unified `lsa()` constructor with a pluggable engine registry, in the
  style of `Nestimate::build_network()`.
- Five built-in engines: `classical`, `two_cell`, `bidirectional`,
  `parallel_dominance`, `nonparallel_dominance`.
- Convenience wrappers: `lsa_classical()`, `lsa_two_cell()`,
  `lsa_bidirectional()`, `lsa_parallel_dominance()`,
  `lsa_nonparallel_dominance()`.
- Sequence-level bootstrap (`bootstrap_lsa()`), permutation
  (`permute_lsa()`), case-drop stability (`stability_lsa()`), and
  split-half reliability (`reliability_lsa()`).
- Multi-group fits via `lsa(data, group = ...)`, returning an
  `lsa_group` object with grouped `transitions()`, `nodes()`,
  `tests()`, and `initial()` methods.
- Tidy edge tables and S3 objects of class
  `c("lsa", "cograph_network")` for seamless integration with the
  `cograph` plotting layer.
- Network interop: `lsa_to_tna()` converts a fit to a `tna` object and
  `as.igraph()` to an `igraph` graph (lagdynamics converts; the downstream
  package does the analysis).
- Reproducibility hooks (`indices=`, `shuffles=`) for bit-identical
  cross-language verification.

Planned for a future release: between-group comparison
(`compare_lsa()`) and stationarity testing (`stationarity_lsa()`).

### Dependencies

Lean runtime: base R plus `grid` (`Imports: grid, stats, utils`). The
plotting and interop packages (`ggplot2`, `cograph`, `tna`, `Nestimate`,
`igraph`, `TraMineR`) are soft `Suggests`, used only when present.
