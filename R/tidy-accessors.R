# Reading verbs: nodes(fit), tests(fit), initial(fit). Each takes the
# fit first and returns a tidy data.frame, matching transitions(fit).
# Pure base R; collision-free names (no overlap with tna / Nestimate).

#' Nodes of an LSA Fit (Tidy)
#'
#' The states and their incoming / outgoing transition totals, as a
#' tidy `data.frame` (one row per state).
#'
#' @param fit An `lsa` fit from [lsa()].
#'
#' @return A `data.frame` with columns `state` (the state name, matching
#'   the `from`/`to` endpoints of [transitions()]), `outgoing`, and
#'   `incoming` (its total out- and in-transition counts).
#'
#' @examples
#' nodes(lsa(group_regulation))
#'
#' @seealso [transitions()], [tests()]
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

#' Tablewise Independence Tests of an LSA Fit (Tidy)
#'
#' The global tests of independence (likelihood-ratio G^2 and Pearson
#' chi-square) as a one-row-per-test `data.frame`.
#'
#' @param fit An `lsa` fit from [lsa()].
#'
#' @return A `data.frame` with columns `test` (`"lrx2"` / `"x2"`),
#'   `statistic`, `df`, `p`. Tests the engine did not compute are
#'   omitted.
#'
#' @examples
#' tests(lsa(group_regulation))
#'
#' @seealso [transitions()], [nodes()]
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

#' Transition-Probability Matrix of an LSA Fit
#'
#' The row-stochastic transition-probability matrix \eqn{P(\text{to} \mid
#' \text{from})}: each row holds the distribution over next states out of
#' a state, so the entries in a row sum to 1 (up to structural zeros).
#' This is the transition-probability matrix a Transition Network Analysis
#' reads; pair it with [initial()] for the initial-state probabilities.
#'
#' @param fit An `lsa` fit from [lsa()].
#'
#' @return A square numeric matrix with `from` states in rows and `to`
#'   states in columns. For an `lsa_group`, a named list of such matrices,
#'   one per group.
#'
#' @examples
#' transition_probabilities(lsa(group_regulation))
#'
#' @seealso [initial()] (initial-state probabilities), [transitions()]
#'   (the tidy per-edge view), [nodes()]
#'
#' @export
transition_probabilities <- function(fit) UseMethod("transition_probabilities")

#' @rdname transition_probabilities
#' @export
transition_probabilities.lsa <- function(fit) {
  .lsa_weight_matrix(fit, "prob")
}

#' @rdname transition_probabilities
#' @export
transition_probabilities.lsa_group <- function(fit) {
  stats::setNames(lapply(fit, transition_probabilities), names(fit))
}

#' Initial-State Distribution of an LSA Fit (Tidy)
#'
#' The proportion of sequences starting in each state, as a tidy
#' `data.frame`. These are the initial-state probabilities (init P) that
#' complement the transition-probability matrix from
#' [transition_probabilities()].
#'
#' @param fit An `lsa` fit from [lsa()].
#'
#' @return A `data.frame` with columns `state`, `init_prob`; zero rows
#'   when the fit came from a transition matrix (no initial states).
#'
#' @examples
#' initial(lsa(group_regulation))
#'
#' @seealso [transitions()], [nodes()]
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
