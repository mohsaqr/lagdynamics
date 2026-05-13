.onLoad <- function(libname, pkgname) {
  .register_builtin_engines()
  # .register_builtin_layouts()    # Step 7
  # .register_builtin_themes()     # Step 7
  # .register_builtin_palettes()   # Step 7
  invisible(NULL)
}

.register_builtin_engines <- function() {
  register_lsa_engine(
    "classical", .engine_classical,
    "Bakeman & Quera classical lag sequential analysis"
  )
  register_lsa_engine(
    "two_cell", .engine_two_cell,
    "2x2 cell test (odds ratio, log-OR Wald z, Yule's Q)"
  )
  register_lsa_engine(
    "bidirectional", .engine_bidirectional,
    "Sackett's bidirectional / matched-pair test on the symmetrized table"
  )
  register_lsa_engine(
    "parallel_dominance", .engine_parallel_dominance,
    "Sackett's parallel-dominance (expected-SE) test"
  )
  register_lsa_engine(
    "nonparallel_dominance", .engine_nonparallel_dominance,
    "Sackett's non-parallel-dominance (observed-SE + binomial) test"
  )
  invisible(NULL)
}
