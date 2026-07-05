# Adapter layer that converts an lsa fit into objects consumed by tna.
# tna remains optional: lagdynamics loads without it, and this bridge
# errors informatively if the user calls it without tna installed.

#' Convert an lsa Fit to a tna Network
#'
#' Convert an `lsa` fit to a `tna`-class network object usable by the
#' `tna` package's centrality, pruning, community, and bootstrap routines.
#'
#' This function is deliberately named `lsa_to_tna()`, not `as_tna()`,
#' to avoid overlapping export names with sibling packages.
#'
#' @param x An `lsa` fit from [lsa()], or an `lsa_group` from
#'   `lsa(..., group = )`.
#' @param weights Character. Which matrix to expose as the tna edge
#'   weights. One of `"prob"` (row-normalised probabilities, default),
#'   `"count"` (raw observed counts), `"adj_res"` (adjusted residuals,
#'   over-representation network), or `"lift"` (observed / expected
#'   association strength).
#' @param ... Method-specific arguments.
#'
#' @return For an `lsa` fit, a `tna` object with `weights`, `inits`,
#'   `labels`, and, when available, `data` slots. For an `lsa_group`, a
#'   `group_tna` object.
#'
#' @examplesIf requireNamespace("tna", quietly = TRUE)
#' fit <- lsa(engagement)
#' net <- lsa_to_tna(fit, weights = "prob")
#' tna::centralities(net)
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

  W <- .lsa_weight_matrix(x, weights, positive_residuals_only = TRUE)
  type <- if (weights == "count") "frequency" else "relative"
  resampleable <- weights %in% c("prob", "count")
  seqdata <- if (resampleable) .lsa_seqdata_matrix(x) else NULL

  if (!resampleable && !is.null(.lsa_seqdata_matrix(x))) {
    warning("lsa_to_tna(weights = '", weights, "') omits the $data slot: ",
            "tna's sequence-resampling verbs would re-estimate a ",
            "probability network, not the '", weights, "' scale in ",
            "$weights. Use weights = 'prob' or 'count' for resampling, ",
            "or resample with lagdynamics instead.",
            call. = FALSE)
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
    stop("Package 'tna' is required for lsa_to_tna(). ",
         "Install with install.packages('tna').", call. = FALSE)
  }
  nets <- lapply(x, function(f) lsa_to_tna(f, weights = weights, ...))
  structure(nets, levels = names(x), class = "group_tna")
}

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
            colors = grDevices::rainbow(length(d$labels)))
}
