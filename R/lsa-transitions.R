# Vectorized transition counting. The mathematical contract is in
# `inst/REFERENCES.md` §1. The implementation uses tabulate() on a
# row-major linearized pair index — no for-loops.

#' Count transitions at a given lag
#'
#' @description
#' `lsa_transitions()` counts the transitions at a given lag and returns the
#' table of observed transitions with its row and column totals: for each pair
#' of states, how often the second occurs `lag` events after the first. This
#' table is the starting point of every lag sequential analysis. At lag 1 it
#' records immediate succession; at larger lags it records what tends to follow
#' a state two or more events later, and at negative lags what tends to precede
#' it. Transitions are counted within sequences only; no transition spans the
#' boundary between two sequences.
#'
#' The count involves no model or test; [lsa()] adds the expected counts,
#' residuals and the other statistics. Use `lsa_transitions()` to inspect the
#' raw counts, or to check how many transitions a lag leaves before fitting.
#'
#' @param x An [lsa_data()] object or any input it accepts.
#' @param lag Integer lag (default `1`). Positive lags count successors,
#'   negative lags predecessors. Input given as a transition-count matrix
#'   supports `lag = 1` only.
#' @param row.names,optional,... Ignored; present for consistency with the
#'   generics.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`:
#' a lag other than 1 for a transition-count matrix, or input that
#' [lsa_data()] rejects.
#'
#' @return An object of class `c("lsa_transitions", "list")`. `print()`
#'   shows the lag, the number of states and transitions and the count
#'   matrix, and returns `x` invisibly. `as.data.frame()` returns one row
#'   per `from -> to` pair with the columns `from`, `to`, `lag`, `count`,
#'   `row_total`, `col_total` and `n_transitions`.
#'
#' @examples
#' tx <- lsa_transitions(engagement)
#' tx
#' as.data.frame(tx)
#'
#' # What precedes each state
#' lsa_transitions(engagement, lag = -1)
#'
#' @seealso [lsa_data()], [lsa()], [lsa_lags()]
#'
#' @export
lsa_transitions <- function(x, lag = 1) {
  if (!inherits(x, "lsa_data")) x <- lsa_data(x)
  stopifnot(
    "`lag` must be a single whole number" =
      is.numeric(lag) && length(lag) == 1L && is.finite(lag) &&
      lag == round(lag)
  )
  lag <- as.integer(lag)

  if (identical(x$source, "transitions")) {
    if (lag != 1L) {
      .lsa_stop("Pre-computed transition matrix input only supports lag = 1.",
           call. = FALSE)
    }
    obs <- x$obs_input
    return(.assemble_lsa_transitions(obs, lag = lag, labels = x$labels))
  }

  obs <- .count_transitions(events = x$events, seq_id = x$seq_id,
                            K = x$n_states, lag = lag)
  dimnames(obs) <- list(x$labels, x$labels)
  .assemble_lsa_transitions(obs, lag = lag, labels = x$labels)
}

# Vectorized count of transitions at a given lag, summed over
# sequences. Implements §1 of inst/REFERENCES.md. `lag` may be any
# integer: positive counts successors (from at t, to at t + lag),
# negative counts predecessors (to at t + lag, i.e. |lag| steps before
# from), and 0 pairs each event with itself (the degenerate diagonal --
# meaningful only for concurrent-code data, which lagdynamics does not model).
.count_transitions <- function(events, seq_id, K, lag) {
  n <- length(events)
  if (n < abs(lag) + 1L) {
    return(matrix(0L, K, K))
  }
  if (lag >= 0L) {
    from_idx <- seq_len(n - lag)
  } else {
    from_idx <- seq.int(1L - lag, n)
  }
  to_idx   <- from_idx + lag
  # Only retain transitions where from and to belong to the same seq.
  keep <- seq_id[from_idx] == seq_id[to_idx]
  if (!any(keep)) return(matrix(0L, K, K))
  from_ev <- events[from_idx][keep]
  to_ev   <- events[to_idx][keep]
  pair    <- (from_ev - 1L) * K + to_ev
  counts  <- tabulate(pair, nbins = K * K)
  matrix(counts, nrow = K, ncol = K, byrow = TRUE)
}

# Build the tidy edge frame and wrap into an S3 object.
.assemble_lsa_transitions <- function(obs, lag, labels) {
  K <- nrow(obs)
  rt <- rowSums(obs)
  ct <- colSums(obs)
  n  <- sum(obs)
  grid <- expand.grid(from = labels, to = labels,
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  # expand.grid varies `from` fastest; vec(obs) = c(col1, col2, ...) is
  # column-major. We want the order that matches expand.grid, which is
  # exactly column-major when `from` is the row index. So as.vector(obs)
  # in default column-major order matches grid$from cycling fastest.
  edges <- data.frame(
    from = grid$from,
    to = grid$to,
    lag = lag,
    count = as.vector(obs),
    row_total = rt[grid$from],
    col_total = ct[grid$to],
    n_transitions = n,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  structure(
    list(
      obs = obs,
      row_totals = rt,
      col_totals = ct,
      n_transitions = n,
      lag = lag,
      labels = labels,
      edges = edges
    ),
    class = c("lsa_transitions", "list")
  )
}

#' @rdname lsa_transitions
#' @export
print.lsa_transitions <- function(x, ...) {
  cat("<lsa_transitions>\n")
  cat(sprintf("  lag:           %d\n", x$lag))
  cat(sprintf("  states (K):    %d\n", length(x$labels)))
  cat(sprintf("  transitions:   %d\n", x$n_transitions))
  cat("  obs:\n")
  print(x$obs)
  invisible(x)
}
