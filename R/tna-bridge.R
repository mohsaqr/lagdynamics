# Adapter layer that converts an lsa fit into objects consumed by
# downstream network packages. Soft dependencies: `tna` and `igraph`
# are declared in Suggests so the adapter functions error informatively
# when the target package is not installed but lagseq itself loads
# without them.
#
# Weight semantics. TNA and FTNA answer different questions: a
# probability-weighted network models *transition tendency* (where
# does the chain go next?), while a count-weighted network models
# *volume* (how many transitions of this type occurred?). We expose
# four weight choices instead of collapsing them:
#
#   prob    : row-normalised transition probabilities (Standard TNA)
#   count   : raw observed transition counts (Frequency TNA / FTNA)
#   adj_res : Haberman / Christensen adjusted residuals (residual TNA)
#   lift    : observed / expected, association strength
#
# Non-finite cells (NA from structural zeros or zero-margin rows) are
# coerced to 0 in the output weight matrix so downstream consumers see
# a well-defined adjacency.

#' Convert an lsa Fit to a tna Network
#'
#' Convert an `lsa` fit to a `tna`-class network object usable by the
#' `tna` package's centrality, pruning, community, and bootstrap
#' routines. Requires the `tna` package (declared in Suggests).
#'
#' The function is deliberately **not** named `as_tna`: the `Nestimate`
#' package exports an `as_tna` generic, and lagseq avoids overlapping
#' export names with sibling packages so the two can be loaded together
#' without masking each other. `lsa_to_tna()` is lagseq's own,
#' collision-free converter.
#'
#' @param x An `lsa` fit from [lsa()], or an `lsa_group` from
#'   `lsa(..., group = )`.
#' @param weights Character. Which matrix to expose as the tna edge
#'   weights. One of `"prob"` (row-normalised probabilities, default,
#'   matches Standard TNA), `"count"` (raw observed counts, matches
#'   Frequency TNA / FTNA), `"adj_res"` (adjusted residuals, residual
#'   network), or `"lift"` (observed / expected, association strength).
#'   `tna` requires non-negative weights, so the `"adj_res"` choice
#'   always clips negative residuals to `0` (returning the
#'   over-representation network). To work with signed residuals, use
#'   [as.igraph.lsa()] or read `fit$adj_res` directly.
#' @param ... Method-specific arguments.
#'
#' @return For an `lsa` fit, a `tna` object with `weights`, `inits`,
#'   `labels`, `data` slots and `type`/`scaling`/`class` attributes
#'   (`type` is `"frequency"` for counts, `"relative"` otherwise). For
#'   an `lsa_group`, a `group_tna` (named list of `tna` objects).
#'
#' @examples
#' \dontrun{
#' fit <- lsa(engagement, engine = "classical")
#' net <- lsa_to_tna(fit, weights = "prob")
#' tna::centralities(net)
#' tna::prune(net, method = "threshold", threshold = 0.05)
#' }
#'
#' @export
lsa_to_tna <- function(x, ...) UseMethod("lsa_to_tna")

#' @rdname lsa_to_tna
#' @export
lsa_to_tna.lsa <- function(x,
                           weights = c("prob", "count", "adj_res", "lift"),
                           ...) {
  weights <- match.arg(weights)
  if (!requireNamespace("tna", quietly = TRUE)) {
    stop("Package 'tna' is required for lsa_to_tna(). ",
         "Install with install.packages('tna').", call. = FALSE)
  }
  W <- .lsa_weight_matrix(x, weights,
                          positive_residuals_only = TRUE)
  # Construct the tna object directly rather than via build_model().
  # build_model(type = "relative") row-normalises the input, which is
  # idempotent for prob weights but destroys the scale of adj_res and
  # lift. Manual construction matches the published tna contract
  # (weights/inits/labels/data slots, type/scaling/class attributes)
  # and is consumed by tna::centralities(), tna::prune() etc.
  #
  # When the fit was built from event sequences we rebuild tna's
  # `tna_seq_data` matrix and attach it as `$data`, so tna's
  # sequence-based verbs (bootstrap, permutation_test, estimate_cs) and
  # the `inits` it derives all work. NOTE: those verbs re-estimate from
  # the sequences using the object's `type`, so they reflect a Standard
  # (relative) or Frequency network -- they match weights = "prob" /
  # "count" but NOT the residual / lift networks, whose scale lives only
  # in `$weights` and is consumed by the non-resampling verbs.
  # Transition-matrix fits have no sequences, so `$data` stays NULL.
  type <- if (weights == "count") "frequency" else "relative"
  seqdata <- .lsa_seqdata_matrix(x)
  inits <- x$inits
  structure(
    list(weights = W, inits = inits,
         labels = rownames(W), data = seqdata),
    type = type,
    scaling = character(0L),
    class = "tna"
  )
}

#' @rdname lsa_to_tna
#' @export
lsa_to_tna.lsa_group <- function(x,
                                 weights = c("prob", "count", "adj_res",
                                              "lift"),
                                 ...) {
  weights <- match.arg(weights)
  if (!requireNamespace("tna", quietly = TRUE)) {
    stop("Package 'tna' is required for lsa_to_tna(). ",
         "Install with install.packages('tna').", call. = FALSE)
  }
  # Each per-group fit becomes a tna; wrapped as a group_tna, the
  # multi-group container consumed by tna's grouped verbs
  # (tna::centralities(), tna::cliques(), tna::communities(), ...).
  nets <- lapply(x, function(f) lsa_to_tna(f, weights = weights, ...))
  structure(nets, levels = names(x), class = "group_tna")
}

#' Convert an lsa Fit to an igraph Graph
#'
#' Wraps the fit's chosen weight matrix in an `igraph` graph. Requires
#' the `igraph` package (declared in Suggests).
#'
#' @param x An `lsa` object returned by [lsa()].
#' @param weights Character. As in [lsa_to_tna()].
#' @param positive_residuals_only Logical. As in [lsa_to_tna()].
#' @param mode Character. Passed to [igraph::graph_from_adjacency_matrix()].
#'   Default `"directed"`; use `"undirected"` for bidirectional fits.
#' @param ... Unused.
#'
#' @return An `igraph` graph with `weight` edge attribute.
#'
#' @examples
#' \dontrun{
#' fit <- lsa(engagement, engine = "classical")
#' g <- as.igraph(fit, weights = "prob")
#' igraph::betweenness(g)
#' }
#'
#' @exportS3Method igraph::as.igraph
as.igraph.lsa <- function(x,
                          weights = c("prob", "count", "adj_res", "lift"),
                          positive_residuals_only = TRUE,
                          mode = "directed",
                          ...) {
  stopifnot(inherits(x, "lsa"))
  weights <- match.arg(weights)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for as.igraph.lsa(). ",
         "Install with install.packages('igraph').", call. = FALSE)
  }
  W <- .lsa_weight_matrix(x, weights, positive_residuals_only)
  igraph::graph_from_adjacency_matrix(W, mode = mode, weighted = TRUE,
                                       diag = TRUE)
}

#' Convert a Grouped lsa Fit to a List of igraph Graphs
#'
#' Applies [as.igraph.lsa()] to each per-group fit. Unlike the `tna`
#' container, `igraph` has no native multi-graph object, so the result
#' is a plain named list of graphs (one per group), which composes with
#' `lapply()` for batched igraph analysis.
#'
#' @param x An `lsa_group` object returned by `lsa(..., group = )`.
#' @param weights Character. As in [as.igraph.lsa()].
#' @param positive_residuals_only Logical. As in [as.igraph.lsa()].
#' @param mode Character. As in [as.igraph.lsa()].
#' @param ... Passed to [as.igraph.lsa()].
#'
#' @return A named list of `igraph` graphs, one per group.
#'
#' @seealso [as.igraph.lsa()]
#' @exportS3Method igraph::as.igraph
as.igraph.lsa_group <- function(x,
                                weights = c("prob", "count", "adj_res",
                                             "lift"),
                                positive_residuals_only = TRUE,
                                mode = "directed",
                                ...) {
  weights <- match.arg(weights)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for as.igraph.lsa_group(). ",
         "Install with install.packages('igraph').", call. = FALSE)
  }
  graphs <- lapply(x, function(f) {
    igraph::as.igraph(f, weights = weights,
                      positive_residuals_only = positive_residuals_only,
                      mode = mode, ...)
  })
  names(graphs) <- names(x)
  graphs
}

# --- internal -----------------------------------------------------------

# Rebuild tna's `tna_seq_data` matrix (rows = sequences, cols = time,
# integer codes 1..K, NA-padded) from an lsa fit's event-level data,
# with the `alphabet`/`labels`/`colors` attributes tna expects. Returns
# NULL when the fit was built from a pre-computed transition matrix
# (no sequences to recover).
.lsa_seqdata_matrix <- function(x) {
  d <- x$data
  if (!identical(d$source, "events") || is.null(d$events)) return(NULL)
  per <- split(d$events,
               factor(d$seq_id, levels = seq_len(d$n_sequences)))
  maxlen <- max(lengths(per))
  m <- t(vapply(per, function(s) {
    c(s, rep(NA_integer_, maxlen - length(s)))
  }, integer(maxlen)))
  dimnames(m) <- NULL
  structure(m,
            class = c("tna_seq_data", "matrix", "array"),
            alphabet = d$labels,
            labels = d$labels,
            colors = grDevices::rainbow(length(d$labels)))
}

# Pull the chosen weight matrix from an lsa fit and coerce non-finite
# cells to 0 so downstream packages see a well-defined adjacency.
.lsa_weight_matrix <- function(x, weights, positive_residuals_only) {
  W <- switch(weights,
    prob    = x$prob,
    count   = x$obs,
    adj_res = x$adj_res,
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
  # prob and obs still carry the row-normalised counts on the
  # diagonal even when the diagonal is structurally forbidden).
  sz <- x$params$structural_zeros
  if (!is.null(sz)) W[sz == 0] <- 0
  W
}
