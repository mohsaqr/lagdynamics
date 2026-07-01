#' lagdynamics: Lag Sequential Analysis, Dynamics, and Lag Transition Networks
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
#' Numerical methods are implemented from primary literature and
#' cross-validated against base-R primitives and hand-formula identities
#' ([stats::loglin()], [stats::chisq.test()], [stats::pchisq()]). See the
#' `inst/REFERENCES.md` file shipped with the package for the
#' formula-by-formula citation map.
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
#' [reliability_lsa()], and analytic certainty via [certainty_lsa()].
#'
#' Group comparison: `compare_lsa()` and `bayes_compare_lsa()` for
#' between-group differences in transition structure.
#'
#' Reading a fit (tidy data frames): [transitions()], [nodes()],
#' [tests()], [initial()], and [summary()].
#'
#' Transition and initial probabilities: [transition_probabilities()]
#' returns the row-stochastic transition-probability matrix and
#' [initial()] the initial-state probabilities.
#'
#' @section Roadmap (not yet implemented):
#'
#' Stationarity tests (`stationarity_lsa()`) are planned but not yet
#' exported. Multi-group fits are built with `lsa(data, group = ...)`;
#' between-group differences are then tested with `compare_lsa()` or
#' `bayes_compare_lsa()`.
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
