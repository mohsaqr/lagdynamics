# lagseq 0.1.0 (in development)

## First public release

`lagseq` is a from-scratch, clean-room implementation of lag sequential
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
  (`permute_lsa()`), and case-drop stability (`stability_lsa()`).
- Group, stationarity, and comparison verbs: `group_lsa()`,
  `stationarity_lsa()`, `compare_lsa()`.
- Tidy edge tables and S3 objects of class
  `c("lsa", "cograph_network")` for seamless integration with the
  `cograph` plotting layer.
- Pipe-friendly aesthetic setters in the style of cograph's `sn_*`
  chain: `ln_layout()`, `ln_theme()`, `ln_palette()`, `ln_edges()`,
  `ln_nodes()`.
- Reproducibility hooks (`indices=`, `shuffles=`) for bit-identical
  cross-language verification.

### Dependencies

Single-runtime-dep policy: base R plus `ggplot2` and `grid`. All other
packages are in `Suggests:` and used only when present.
