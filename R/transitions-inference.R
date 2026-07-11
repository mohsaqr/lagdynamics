# transitions() methods for the inference result classes. Every inference
# verb produces a per-transition table; transitions() is the one way to read
# it, exactly as transitions(fit) reads a fitted model. Without these methods
# the only route to the table is as.data.frame(), which forces the caller to
# assemble the view by hand.
#
# Each class names two columns: the boolean the method uses to decide an
# edge (`significant = TRUE` keeps those rows) and the signed effect the
# method estimates (`sort = "strength"` orders by its magnitude).

# Read an inference result's per-edge table, then narrow and order it. The
# per-class column names arrive as `flag` and `strength`; everything else is
# shared. as.data.frame() is the class's own tidy view, so this stays a thin
# filter over it.
.inference_transitions <- function(x, significant, sort, flag, strength) {
  stopifnot(is.logical(significant), length(significant) == 1L,
            !is.na(significant))
  sort <- match.arg(sort, c("none", "strength"))
  e <- as.data.frame(x)
  if (isTRUE(significant)) {
    keep <- e[[flag]]
    e <- e[!is.na(keep) & keep, , drop = FALSE]
  }
  if (identical(sort, "strength") && nrow(e) > 0L) {
    key <- abs(e[[strength]])
    key[!is.finite(key)] <- -Inf
    e <- e[order(key, decreasing = TRUE), , drop = FALSE]
  }
  rownames(e) <- NULL
  e
}

#' @rdname transitions
#' @export
transitions.lsa_bootstrap <- function(fit, significant = FALSE,
                                      sort = c("none", "strength"), ...) {
  .inference_transitions(fit, significant, sort,
                         flag = "adj_res_stable",
                         strength = "adj_res_observed")
}

#' @rdname transitions
#' @export
transitions.lsa_certainty <- function(fit, significant = FALSE,
                                      sort = c("none", "strength"), ...) {
  .inference_transitions(fit, significant, sort,
                         flag = "stable",
                         strength = "adj_res_observed")
}

#' @rdname transitions
#' @export
transitions.lsa_stability <- function(fit, significant = FALSE,
                                      sort = c("none", "strength"), ...) {
  .inference_transitions(fit, significant, sort,
                         flag = "stable",
                         strength = "stability")
}

#' @rdname transitions
#' @export
transitions.lsa_permutation <- function(fit, significant = FALSE,
                                        sort = c("none", "strength"), ...) {
  .inference_transitions(fit, significant, sort,
                         flag = "significant",
                         strength = "observed_adj_res")
}

#' @rdname transitions
#' @export
transitions.lsa_comparison <- function(fit, significant = FALSE,
                                       sort = c("none", "strength"), ...) {
  .inference_transitions(fit, significant, sort,
                         flag = "significant",
                         strength = "diff")
}

#' @rdname transitions
#' @export
transitions.lsa_comparison_pairwise <- function(fit, significant = FALSE,
                                                sort = c("none", "strength"),
                                                ...) {
  .inference_transitions(fit, significant, sort,
                         flag = "significant",
                         strength = "diff")
}

#' @rdname transitions
#' @export
transitions.lsa_lags <- function(fit, significant = FALSE,
                                 direction = c("any", "over", "under"),
                                 min_count = NULL, alpha = NULL,
                                 sort = c("none", "strength", "count", "prob"),
                                 ...) {
  direction <- match.arg(direction)
  sort <- match.arg(sort)
  # An lsa_lags is a named list of lsa fits, one per lag. Each fit carries
  # its own alpha, so the filters apply per lag and the results stack.
  out <- do.call(rbind, lapply(unclass(fit), transitions,
                               significant = significant,
                               direction = direction,
                               min_count = min_count,
                               alpha = alpha,
                               sort = sort))
  rownames(out) <- NULL
  out
}
