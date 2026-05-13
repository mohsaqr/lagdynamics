# Vectorized transition counting. The mathematical contract is in
# `inst/REFERENCES.md` §1. The implementation uses tabulate() on a
# row-major linearized pair index — no for-loops.

#' Tidy Transition Counts at a Given Lag
#'
#' Computes the `K x K` transition count matrix from canonical lag
#' sequential data, optionally returning a tidy long-format edge table
#' alongside the matrix.
#'
#' Transitions are counted within sequences only; no transition spans a
#' sequence boundary. For input that was supplied as a pre-computed
#' transition matrix (`source = "transitions"` on the `lsa_data`
#' object), the input matrix is returned at lag 1 and an error is
#' raised for any other lag.
#'
#' @param x Either an [lsa_data] object or any input accepted by
#'   [lsa_data()] (which will be coerced).
#' @param lag Positive integer. The lag at which to count transitions.
#'   Default `1`.
#'
#' @return An object of class `c("lsa_transitions", "list")` with
#'   elements:
#' \describe{
#'   \item{obs}{The `K x K` observed transition count matrix with
#'     `dimnames` set to the labels.}
#'   \item{row_totals}{Length-`K` vector `rowSums(obs)`.}
#'   \item{col_totals}{Length-`K` vector `colSums(obs)`.}
#'   \item{n_transitions}{Scalar `sum(obs)`.}
#'   \item{lag}{The lag used.}
#'   \item{labels}{Character vector of state labels.}
#'   \item{edges}{Tidy long-format data.frame with one row per
#'     `(from, to)` cell containing columns `from`, `to`, `lag`,
#'     `count`, `row_total`, `col_total`, `n_transitions`.}
#' }
#'
#' @examples
#' d <- lsa_data(c("a", "b", "a", "c", "b", "a"))
#' tx <- lsa_transitions(d, lag = 1)
#' tx$obs
#' head(tx$edges)
#'
#' @seealso [lsa_data()], [lsa()]
#'
#' @export
lsa_transitions <- function(x, lag = 1) {
  if (!inherits(x, "lsa_data")) x <- lsa_data(x)
  stopifnot(
    is.numeric(lag), length(lag) == 1L, lag >= 1, lag == round(lag)
  )
  lag <- as.integer(lag)

  if (identical(x$source, "transitions")) {
    if (lag != 1L) {
      stop("Pre-computed transition matrix input only supports lag = 1.",
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
# sequences. Implements §1 of inst/REFERENCES.md.
.count_transitions <- function(events, seq_id, K, lag) {
  n <- length(events)
  if (n < lag + 1L) {
    return(matrix(0L, K, K))
  }
  from_idx <- seq_len(n - lag)
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
