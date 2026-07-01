# Internal: pull a named weight matrix out of an lsa fit and coerce
# non-finite cells to 0 so downstream consumers (plotting, transition
# probabilities) see a well-defined adjacency. Structural-zero cells are
# forced to 0 even when the raw weight type would otherwise carry a value
# there. Shared by plot_transitions(), plot_chords(), plot_polar(), and
# transition_probabilities().

.lsa_weight_matrix <- function(x, weights, positive_residuals_only = FALSE) {
  W <- switch(weights,
    prob    = x$prob,
    count   = x$obs,
    adj_res = x$adj_res,
    yules_q = x$yules_q,
    lift    = matrix(x$edges$lift,
                     nrow = nrow(x$obs),
                     ncol = ncol(x$obs),
                     dimnames = dimnames(x$obs))
  )
  if (weights == "adj_res" && isTRUE(positive_residuals_only)) {
    W[!is.finite(W) | W < 0] <- 0
  } else {
    W[!is.finite(W)] <- 0
  }
  # Honour structural-zero declarations even when the raw weight type
  # would otherwise expose a nonzero value at a forbidden cell (e.g.
  # prob and obs still carry the row-normalised counts on the diagonal
  # even when the diagonal is structurally forbidden).
  sz <- x$params$structural_zeros
  if (!is.null(sz)) W[sz == 0] <- 0
  W
}
