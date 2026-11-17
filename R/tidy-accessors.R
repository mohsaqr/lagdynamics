# Reading verbs: nodes(fit), tests(fit), initial(fit). Each takes the
# fit first and returns a tidy data.frame, matching transitions(fit).
# Pure base R; collision-free names (no overlap with tna / Nestimate).

#' Read the states of a lag sequential analysis fit
#'
#' @description
#' `nodes()` returns the states of a lag sequential analysis fit, one row per
#' state with its total outgoing and incoming transition counts. In a lag
#' transition network the states are the nodes, and how often a state is left
#' and how often it is entered describe its place in the process: a state with
#' many outgoing but few incoming transitions tends to open episodes, one with
#' the reverse tends to close them. The marginal totals also set the expected
#' counts against which every transition is tested.
#'
#' The shares of all transitions and the initial-state probability are in
#' `as.data.frame(summary(fit), what = "nodes")`; the per-transition table is
#' [transitions()].
#'
#' @param fit An `lsa` fit from [lsa()], or a grouped fit from
#'   `lsa(..., group = )`.
#'
#' @return A `data.frame` with one row per state and the columns `state`
#'   (matching the `from` and `to` names of [transitions()]), `outgoing`
#'   and `incoming`. A grouped fit gives one row per group and state, with
#'   a leading `group` column.
#'
#' @examples
#' nodes(lsa(engagement))
#'
#' @seealso [transitions()], [initial()], [tests()], [lsa-methods] for
#'   `summary()`.
#'
#' @export
nodes <- function(fit) UseMethod("nodes")

#' @rdname nodes
#' @export
nodes.lsa <- function(fit) {
  n <- fit$nodes
  data.frame(state = n$label,
             outgoing = n$outgoing,
             incoming = n$incoming,
             stringsAsFactors = FALSE, row.names = NULL)
}

#' @rdname nodes
#' @export
nodes.lsa_group <- function(fit) .bind_group_edges(fit, nodes)

#' Read the tablewise tests of independence of a fit
#'
#' @description
#' `tests()` returns the tablewise tests of independence of a fit, one row per
#' test. They test whether the table as a whole departs from independence:
#' whether the state at one time carries any information about the state `lag`
#' events later (Bakeman & Quera, 1995). The likelihood-ratio statistic
#' \eqn{G^2} and the Pearson \eqn{\chi^2} statistic test this omnibus hypothesis
#' against a chi-square distribution, and are read before individual transitions
#' are interpreted.
#'
#' The tests are computed by the classical engine; the other built-in engines
#' test cells only, and `tests()` then returns a table with no rows.
#'
#' @param fit An `lsa` fit from [lsa()], or a grouped fit from
#'   `lsa(..., group = )`.
#'
#' @return A `data.frame` with one row per test and the columns `test`
#'   (`"lrx2"` for the likelihood-ratio test, `"x2"` for Pearson's test),
#'   `statistic`, `df` and `p`. A grouped fit gains a leading `group`
#'   column.
#'
#' @examples
#' tests(lsa(engagement))
#'
#' @references
#' Bakeman, R., & Quera, V. (1995). \emph{Analyzing interaction:
#' Sequential analysis with SDIS and GSEQ}. Cambridge University Press.
#'
#' @seealso [transitions()] for the cell-wise tests; [nodes()];
#'   [lsa-methods] for `summary()`.
#'
#' @export
tests <- function(fit) UseMethod("tests")

#' @rdname tests
#' @export
tests.lsa <- function(fit) {
  rows <- list(lrx2 = fit$lrx2, x2 = fit$x2)
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) {
    return(data.frame(test = character(0L), statistic = numeric(0L),
                      df = integer(0L), p = numeric(0L)))
  }
  data.frame(
    test = names(rows),
    statistic = vapply(rows, function(r) as.numeric(r$statistic), numeric(1L)),
    df = vapply(rows, function(r) as.integer(r$df), integer(1L)),
    p = vapply(rows, function(r) as.numeric(r$p), numeric(1L)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @rdname tests
#' @export
tests.lsa_group <- function(fit) .bind_group_edges(fit, tests)

#' Read the transition probabilities of a fit
#'
#' @description
#' `transition_probabilities()` returns the transition probabilities of a fit,
#' one row per ordered pair of states. The transition probability
#' \eqn{P(\text{to} \mid \text{from})} is the share of the transitions leaving a
#' state that go to each other state, so the probabilities out of one state sum
#' to 1. They describe where the process goes next, regardless of chance
#' expectation; the adjusted residuals of [transitions()] say which of these
#' moves are more or less frequent than independence predicts. Together with the
#' initial-state probabilities of [initial()], they are the parameters of the
#' first-order Markov model read by Transition Network Analysis.
#'
#' A cell declared a structural zero has probability 0, as does every cell of a
#' state that is never left.
#'
#' @param fit An `lsa` fit from [lsa()], or a grouped fit from
#'   `lsa(..., group = )`.
#'
#' @return A `data.frame` with one row per `from -> to` pair (`K^2` rows
#'   for `K` states, ordered by `from` and then `to`) and the columns
#'   `from`, `to` and `prob`. A grouped fit stacks the groups under a
#'   leading `group` column.
#'
#' @examples
#' transition_probabilities(lsa(engagement))
#'
#' @seealso [initial()] for the initial-state probabilities;
#'   [transitions()], which also carries `prob`; [lsa_to_tna()].
#'
#' @export
transition_probabilities <- function(fit) UseMethod("transition_probabilities")

#' @rdname transition_probabilities
#' @export
transition_probabilities.lsa <- function(fit) {
  P <- .lsa_weight_matrix(fit, "prob")
  labels <- rownames(P)
  K <- length(labels)
  # Row-major (from, then to): P is stored column-major, so read t(P).
  data.frame(
    from = rep(labels, each = K),
    to   = rep(labels, times = K),
    prob = as.vector(t(P)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @rdname transition_probabilities
#' @export
transition_probabilities.lsa_group <- function(fit) {
  .bind_group_edges(fit, transition_probabilities)
}

#' Read the initial-state distribution of a fit
#'
#' @description
#' `initial()` returns the initial-state distribution of a fit, one row per
#' state: the proportion of sequences that begin in each state. It describes how
#' episodes open, which the transition table cannot show, and it completes the
#' description of the process as a first-order Markov chain together with
#' [transition_probabilities()].
#'
#' A fit built from a transition-count matrix has no sequences and so no initial
#' states; the table then has no rows.
#'
#' @param fit An `lsa` fit from [lsa()], or a grouped fit from
#'   `lsa(..., group = )`.
#'
#' @return A `data.frame` with one row per state and the columns `state`
#'   and `init_prob` (summing to 1). A grouped fit gains a leading `group`
#'   column.
#'
#' @examples
#' initial(lsa(engagement))
#'
#' @seealso [transition_probabilities()], [nodes()], [transitions()]
#'
#' @export
initial <- function(fit) UseMethod("initial")

#' @rdname initial
#' @export
initial.lsa <- function(fit) {
  if (is.null(fit$inits)) {
    return(data.frame(state = character(0L), init_prob = numeric(0L)))
  }
  data.frame(state = names(fit$inits),
             init_prob = as.numeric(fit$inits),
             stringsAsFactors = FALSE, row.names = NULL)
}

#' @rdname initial
#' @export
initial.lsa_group <- function(fit) .bind_group_edges(fit, initial)
