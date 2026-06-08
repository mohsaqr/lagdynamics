#' lagseq: Modern Lag Sequential Analysis with Tidy Transition Networks
#'
#' A unified, pipe-friendly interface for lag sequential analysis (LSA) of
#' categorical event sequences. Provides classical, two-cell,
#' bidirectional, parallel-dominance, and non-parallel-dominance engines
#' via a single [lsa()] constructor with a pluggable engine registry,
#' sequence-level bootstrap and permutation inference, case-drop
#' stability, and tidy edge tables ready for transition-network
#' visualization.
#'
#' @section Design principles:
#'
#' All numerical methods are implemented from primary literature and
#' cross-validated against published worked examples and base-R primitives
#' ([stats::loglin()], [stats::chisq.test()], [stats::pchisq()]). See the
#' `inst/REFERENCES.md` file shipped with the package for the
#' formula-by-formula citation map and `inst/EQUIVALENCE.md` for the
#' oracle table and pass rates.
#'
#' The package returns S3 objects of class `c("lsa", "cograph_network")`,
#' which means fits inherit compatibility with the `cograph` package's
#' visualization layer when it is installed.
#'
#' @section Main functions:
#'
#' Construction: [lsa()], [lsa_data()], [lsa_transitions()]
#'
#' Convenience wrappers: `lsa_classical()`, `lsa_two_cell()`,
#' `lsa_bidirectional()`, `lsa_parallel_dominance()`,
#' `lsa_nonparallel_dominance()`
#'
#' Engine registry: `register_lsa_engine()`, `get_lsa_engine()`,
#' `list_lsa_engines()`
#'
#' Inference: [bootstrap_lsa()], [permute_lsa()], [stability_lsa()],
#' [reliability_lsa()]
#'
#' Reading a fit (tidy data frames): [transitions()], [nodes()],
#' [tests()], [initial()], and [summary()].
#'
#' TNA / igraph bridge: [lsa_to_tna()], [as.igraph.lsa()]. These
#' convert an `lsa` fit into the native object of a downstream network
#' package; the analysis itself (centralities via `tna::centralities()`,
#' centrality stability via `tna::estimate_cs()`, pruning via
#' `tna::prune()`, communities via `tna::communities()`) is performed by
#' that package. lagseq stays the *converter*, not the analyser, and is
#' not coupled to `tna` at install time.
#'
#' @section Roadmap (not yet implemented):
#'
#' A dedicated between-group network comparison verb (`compare_lsa()`)
#' and stationarity tests (`stationarity_lsa()`) are planned but not
#' yet exported. Until then, multi-group fits are built directly with
#' `lsa(data, group = ...)`, and between-group comparison is possible
#' by routing each group's fit through [lsa_to_tna()] and using
#' `tna::compare()`.
#'
#' @references
#'
#' Bakeman, R., & Quera, V. (1995). \emph{Analyzing interaction:
#' Sequential analysis}. Cambridge University Press.
#'
#' Bakeman, R., & Gottman, J. M. (1997). \emph{Observing interaction: An
#' introduction to sequential analysis} (2nd ed.). Cambridge University
#' Press.
#'
#' Christensen, R. (1997). \emph{Log-linear models and logistic
#' regression} (2nd ed.). Springer.
#'
#' Haberman, S. J. (1979). \emph{Analysis of qualitative data: Volume 2,
#' New developments}. Academic Press.
#'
#' Sackett, G. P. (1979). The lag sequential analysis of contingency and
#' cyclicity in behavioral interaction research. In J. D. Osofsky (Ed.),
#' \emph{Handbook of infant development} (pp. 623-649). Wiley.
#'
#' Wampold, B. E. (1982). Sequential analysis: Stationarity, statistical
#' procedures, and applications. \emph{Psychological Bulletin, 92}(2),
#' 380-387.
#'
#' Wampold, B. E. (1984). Tests of dominance in sequential categorical
#' data. \emph{Psychological Bulletin, 96}(2), 424-429.
#'
#' Wickens, T. D. (1989). \emph{Multiway contingency tables analysis for
#' the social sciences}. Lawrence Erlbaum.
#'
#' @keywords internal
"_PACKAGE"
