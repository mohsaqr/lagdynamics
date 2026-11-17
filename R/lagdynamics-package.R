#' lagdynamics: Lag Sequential Analysis, Dynamics, and Lag Transition Networks
#'
#' @description
#' `lagdynamics` fits lag sequential analysis models, reads them as lag
#' transition networks, and assesses, compares and plots them. Lag sequential
#' analysis (LSA; Sackett, 1979; Bakeman & Gottman, 1997) tests whether, in a
#' stream of coded events, one behaviour follows another more or less often than
#' chance would predict.
#'
#' The observed transitions at a given lag are counted into a contingency table
#' of antecedent by consequent states, an expected table is derived under the
#' hypothesis that the two are independent, and each cell is tested with an
#' adjusted residual (Haberman, 1979). The cells that depart from independence
#' form a lag transition network: a directed graph whose edges are the over- and
#' under-represented sequential contingencies.
#'
#' The model is fitted with one constructor, [lsa()], which accepts event
#' sequences, wide sequence tables, transition-count matrices, long event logs
#' and sequence objects from \pkg{tna}, \pkg{Nestimate} and \pkg{TraMineR}.
#' Results are read with verbs that return tidy data frames, assessed with
#' resampling and analytic inference, compared across groups, and drawn as
#' heatmaps, networks, chord diagrams and polar displays.
#'
#' The documentation site, with tutorials and case studies, is at
#' <https://pak.dynasite.org/lagdynamics/>.
#'
#' @section Workflow:
#' 1. **Fit** a model with [lsa()]: choose the lag, the engine (the
#'    classical test of Bakeman & Quera, 1995, or the two-cell,
#'    bidirectional and dominance tests of Sackett, 1979, and Wampold,
#'    1984), and any structurally forbidden transitions. A `group`
#'    argument fits one model per group on a shared state set.
#' 2. **Read** the fit with [transitions()], [nodes()], [tests()],
#'    [initial()], [transition_probabilities()] and `summary()`.
#' 3. **Evaluate** how far the transitions can be trusted:
#'    [bootstrap_lsa()], [certainty_lsa()], [permute_lsa()],
#'    [stability_lsa()] and [reliability_lsa()].
#' 4. **Compare** groups with [compare_lsa()] (permutation test) or
#'    [bayes_compare_lsa()] (Dirichlet-multinomial posterior).
#' 5. **Plot** the fit or a comparison: [plot.lsa()], [plot_transitions()],
#'    [plot_chords()], [plot_polar()], [plot_forest()] and
#'    [plot.lsa_comparison()].
#' 6. **Go beyond LSA** with multi-lag profiles ([lsa_lags()],
#'    [lag_profile()]), directed transfer entropy ([transfer_entropy()]) and
#'    Transition Network Analysis through [lsa_to_tna()].
#'
#' @section Functions by task:
#' * Fitting: [lsa()] and its engine wrappers `lsa_classical()`,
#'   `lsa_two_cell()`, `lsa_bidirectional()`, `lsa_parallel_dominance()`,
#'   `lsa_nonparallel_dominance()`; [lsa_lags()].
#' * Input: [lsa_data()] (canonical sequences), [lsa_transitions()]
#'   (transition counts at a lag), [lsa_ipf()] (expected counts with
#'   structural zeros).
#' * Reading: [transitions()], [nodes()], [tests()], [initial()],
#'   [transition_probabilities()], [lag_profile()], and the `print()`,
#'   `summary()` and `as.data.frame()` methods in [lsa-methods].
#' * Inference: [bootstrap_lsa()], [certainty_lsa()], [permute_lsa()],
#'   [stability_lsa()], [reliability_lsa()]; every result is read with
#'   [transitions()] or `as.data.frame()`.
#' * Group comparison: [compare_lsa()], [bayes_compare_lsa()].
#' * Plotting: [plot.lsa()], [plot_transitions()], [plot_chords()],
#'   [plot_polar()], [plot_forest()], [plot.lsa_comparison()].
#' * Engines: [list_lsa_engines()], [get_lsa_engine()],
#'   [register_lsa_engine()], [unregister_lsa_engine()].
#' * Other: [transfer_entropy()], [lsa_to_tna()].
#'
#' @section Data:
#' [engagement] (engagement states of 138 students over 15 time points),
#' [group_regulation] and [group_regulation_long] (collaborative regulation
#' actions, as a wide sequence table and as a long event log), and
#' [ai_long] (coded AI behaviours in human-AI programming sessions).
#'
#' @section Formulas and sources:
#' Each formula the package implements is listed with its primary source,
#' page and equation number in the file `REFERENCES.md` installed with the
#' package; its location is given by
#' `system.file("REFERENCES.md", package = "lagdynamics")`.
#'
#' Fits are objects of class `c("lsa", "cograph_network")`, so the network
#' functions of \pkg{cograph} accept them directly.
#'
#' @references
#' Bakeman, R., & Gottman, J. M. (1997). \emph{Observing interaction: An
#' introduction to sequential analysis} (2nd ed.). Cambridge University
#' Press.
#'
#' Bakeman, R., & Quera, V. (1995). \emph{Analyzing interaction:
#' Sequential analysis with SDIS and GSEQ}. Cambridge University Press.
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
#' Wampold, B. E. (1984). Tests of dominance in sequential categorical
#' data. \emph{Psychological Bulletin, 96}(2), 424-429.
#'
#' Wickens, T. D. (1989). \emph{Multiway contingency tables analysis for
#' the social sciences}. Lawrence Erlbaum.
#'
#' @seealso [lsa()] to start; <https://pak.dynasite.org/lagdynamics/> for
#'   the tutorials.
#'
#' @keywords internal
"_PACKAGE"
