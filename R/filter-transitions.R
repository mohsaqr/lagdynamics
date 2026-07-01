# The transitions() verb: the one predictable way to read an lsa fit's
# edges as a tidy data.frame, with optional filters. Replaces the old
# significant_/overrepresented_/underrepresented_/common_transitions
# family with a single verb whose arguments select the subset.
#
# transitions(fit) returns every transition; the arguments narrow it.
# On a grouped fit the per-group results are row-bound with a leading
# `group` column. Base data.frame throughout (no tibble dependency).

# Row-bind a per-group result into one data.frame with a leading
# `group` column. Empty groups are dropped; an all-empty result keeps
# the correct columns so the shape is stable.
.bind_group_edges <- function(x, fun, ...) {
  parts <- lapply(x, fun, ...)
  pieces <- Map(function(df, g) {
    if (nrow(df) == 0L) return(NULL)
    cbind(group = rep(g, nrow(df)), df, stringsAsFactors = FALSE)
  }, parts, names(x))
  pieces <- pieces[!vapply(pieces, is.null, logical(1L))]
  if (length(pieces) == 0L) {
    out <- cbind(group = character(0L),
                 parts[[1L]][0, , drop = FALSE],
                 stringsAsFactors = FALSE)
  } else {
    out <- do.call(rbind, pieces)
  }
  rownames(out) <- NULL
  out
}

#' Transitions of an LSA Fit (Tidy)
#'
#' The canonical way to read a fit's transitions as a tidy
#' one-row-per-transition `data.frame`. `transitions(fit)` returns every
#' transition; the arguments narrow it.
#'
#' @param fit An `lsa` fit from [lsa()], or a grouped `lsa_group`.
#' @param significant Logical. Keep only transitions whose adjusted-
#'   residual p-value is below `alpha`. Default `FALSE` (keep all).
#' @param direction One of `"any"` (default), `"over"`
#'   (over-represented: significant with a positive residual), or
#'   `"under"` (under-represented: significant with a negative
#'   residual). Selecting a direction implies `significant = TRUE`.
#' @param min_count Optional integer. Keep only transitions observed at
#'   least this many times. Default `NULL` (no count filter).
#' @param alpha Significance threshold. Default `NULL`, which uses the
#'   alpha recorded on the fit (`fit$params$alpha`).
#' @param sort Row ordering. `"none"` (default) keeps the matrix
#'   (column-major) order; `"strength"` orders by `|adj_res|`, `"count"` by
#'   observed count, `"prob"` by transition probability -- each descending,
#'   so the table reads strongest-first.
#'
#' @return A `data.frame`, one row per transition, with columns `from`,
#'   `to` (the source and target **state names**), `lag`, `count`,
#'   `expected`, `prob` (row-conditional), `prob_col` (column-
#'   conditional), `adj_res`, `p`, `yules_q`, `kappa`, `kappa_z`,
#'   `kappa_p`, `lift`, `sign`, `significant`. Engines that compute extra
#'   per-cell statistics append them as further columns (e.g. the two-cell
#'   engine adds `odds_ratio`, `log_or`, `log_or_se`). A grouped fit gains a
#'   leading `group` column. Row names are reset.
#'
#' @examples
#' fit <- lsa(group_regulation)
#' transitions(fit)                       # all transitions
#' transitions(fit, significant = TRUE)   # significant ones
#' transitions(fit, direction = "over")   # over-represented
#' transitions(fit, min_count = 500)      # frequently observed
#'
#' @seealso [lsa()], [nodes()], [tests()]
#'
#' @export
transitions <- function(fit, significant = FALSE,
                        direction = c("any", "over", "under"),
                        min_count = NULL, alpha = NULL,
                        sort = c("none", "strength", "count", "prob")) {
  UseMethod("transitions")
}

#' @rdname transitions
#' @export
transitions.lsa <- function(fit, significant = FALSE,
                            direction = c("any", "over", "under"),
                            min_count = NULL, alpha = NULL,
                            sort = c("none", "strength", "count", "prob")) {
  direction <- match.arg(direction)
  sort <- match.arg(sort)
  if (is.null(alpha)) alpha <- fit$params$alpha
  stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1)
  e <- fit$edges
  if (isTRUE(significant) || direction != "any") {
    e <- e[is.finite(e$p) & e$p < alpha, , drop = FALSE]
  }
  if (direction == "over") {
    e <- e[is.finite(e$adj_res) & e$adj_res > 0, , drop = FALSE]
  } else if (direction == "under") {
    e <- e[is.finite(e$adj_res) & e$adj_res < 0, , drop = FALSE]
  }
  if (!is.null(min_count)) {
    stopifnot(is.numeric(min_count), length(min_count) == 1L,
              is.finite(min_count), min_count >= 1)
    e <- e[is.finite(e$count) & e$count >= min_count, , drop = FALSE]
  }
  # Project a tidy, name-keyed view. The endpoints are the state names
  # (`from`/`to`); the integer ids, duplicate `weight`/`edge`, and other
  # cograph-protocol fields stay on `fit$edges`, not in the accessor.
  keep <- c("from_label", "to_label", "lag", "count", "expected", "prob",
            "prob_col", "adj_res", "p", "yules_q", "kappa", "kappa_z",
            "kappa_p", "lift", "sign", "significant")
  # Engine-specific statistic columns (e.g. odds_ratio / log_or / log_or_se
  # from the two-cell engine) are surfaced too: anything in the edge table
  # beyond the standard schema and the cograph-protocol id/duplicate fields.
  extra_cols <- setdiff(names(e), c(keep, "from", "to", "edge", "weight"))
  keep <- c(keep, extra_cols)
  e <- e[, keep, drop = FALSE]
  names(e)[1:2] <- c("from", "to")
  # Optional ordering, strongest first, so the table reads top-down. Default
  # "none" preserves the matrix (column-major) order.
  if (sort != "none" && nrow(e) > 0L) {
    key <- switch(sort,
                  strength = abs(e$adj_res),
                  count    = e$count,
                  prob     = e$prob)
    key[!is.finite(key)] <- -Inf
    e <- e[order(key, decreasing = TRUE), , drop = FALSE]
  }
  rownames(e) <- NULL
  e
}

#' @rdname transitions
#' @export
transitions.lsa_group <- function(fit, significant = FALSE,
                                  direction = c("any", "over", "under"),
                                  min_count = NULL, alpha = NULL,
                                  sort = c("none", "strength", "count",
                                           "prob")) {
  direction <- match.arg(direction)
  sort <- match.arg(sort)
  .bind_group_edges(fit, transitions, significant = significant,
                    direction = direction, min_count = min_count,
                    alpha = alpha, sort = sort)
}
