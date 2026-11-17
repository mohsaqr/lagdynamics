# Adapter layer that converts an lsa fit into objects consumed by tna.
# tna remains optional: lagdynamics loads without it, and this bridge
# errors informatively if the user calls it without tna installed.

#' Convert a lag sequential analysis fit to a transition network for tna
#'
#' @description
#' `lsa_to_tna()` converts a lag sequential analysis fit to a transition network
#' for the \pkg{tna} package: a `tna` network from an `lsa` fit, or a
#' `group_tna` from a grouped fit, so that a network found by LSA can be
#' analysed with the network tools of \pkg{tna}. Lag sequential analysis and
#' Transition Network Analysis (TNA) read the same transitions from two angles.
#' LSA tests which transitions are more or less frequent than independence
#' predicts; TNA treats the transition probabilities as a weighted network and
#' studies its structure with centralities, communities, cliques and
#' resampling-based pruning.
#'
#' The edge weights are chosen by `weights`: the transition probabilities (the
#' TNA model), the counts, the positive adjusted residuals (the
#' over-representation network) or the lift. With probabilities or counts the
#' network also carries the initial-state probabilities and the sequences, so
#' the resampling functions of \pkg{tna} can be used on it. With residuals or
#' lift it carries neither, because \pkg{tna} would re-estimate probabilities
#' from the sequences; resample with [bootstrap_lsa()] instead.
#'
#' @param x An `lsa` fit from [lsa()], or a grouped fit from
#'   `lsa(..., group = )`.
#' @param weights Edge weights: `"prob"` (default, transition
#'   probabilities), `"count"` (observed counts), `"adj_res"` (adjusted
#'   residuals, negative values set to 0) or `"lift"` (observed over
#'   expected counts).
#' @param ... Passed between methods.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and
#' `lagdynamics_missing_package` when \pkg{tna} is not installed. With
#' `weights = "adj_res"` or `"lift"` and a fit that has sequences, a
#' warning of class `lagdynamics_tna_no_data` notes that the network
#' carries no sequence data.
#'
#' @return A `tna` object (for an `lsa` fit) or a `group_tna` object (one
#'   `tna` per group, for a grouped fit), read with the functions of the
#'   \pkg{tna} package.
#'
#' @examplesIf requireNamespace("tna", quietly = TRUE)
#' fit <- lsa(engagement)
#' net <- lsa_to_tna(fit, weights = "prob")
#' tna::centralities(net)
#'
#' @seealso [transition_probabilities()] and [initial()] for the same
#'   parameters as tidy tables; [lsa()].
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
    .lsa_stop(class = "missing_package", "Package 'tna' is required for lsa_to_tna(). ",
         "Install with install.packages('tna').", call. = FALSE)
  }

  W <- .lsa_weight_matrix(x, weights, positive_residuals_only = TRUE)
  type <- if (weights == "count") "frequency" else "relative"
  resampleable <- weights %in% c("prob", "count")
  seqdata <- if (resampleable) .lsa_seqdata_matrix(x) else NULL

  if (!resampleable && !is.null(.lsa_seqdata_matrix(x))) {
    warning(warningCondition(paste0(
      "lsa_to_tna(weights = '", weights, "') returns no sequence data: ",
      "tna's sequence-resampling verbs would re-estimate a probability ",
      "network, not the '", weights, "' scale of the weights. Use ",
      "weights = 'prob' or 'count' for resampling, or resample with ",
      "lagdynamics instead."),
      class = "lagdynamics_tna_no_data", call = NULL))
  }

  structure(
    list(weights = W, inits = if (resampleable) x$inits else NULL,
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
    .lsa_stop(class = "missing_package", "Package 'tna' is required for lsa_to_tna(). ",
         "Install with install.packages('tna').", call. = FALSE)
  }
  nets <- lapply(x, function(f) lsa_to_tna(f, weights = weights, ...))
  structure(nets, levels = names(x), class = "group_tna")
}

# Okabe-Ito palette (colour-blind safe), recycled for the state colours
# that tna expects on its sequence data.
.lsa_okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                    "#D55E00", "#CC79A7", "#999999", "#000000")

# Rebuild tna's sequence-data matrix from an lsa fit. Returns NULL for
# transition-matrix fits because no event-level sequences are available.
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
            colors = rep_len(.lsa_okabe_ito, length(d$labels)))
}
