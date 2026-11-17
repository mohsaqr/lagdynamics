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

#' Read the transitions of a fit or an inference result
#'
#' @description
#' `transitions()` returns the transitions of a lag sequential analysis as a
#' table, one row per ordered pair of states, from a fit and from every result
#' built on a fit. The transitions are the substance of the analysis: for each
#' ordered pair of states, how often the second follows the first at the chosen
#' lag, how often independence would predict, and whether the difference is
#' larger than chance. A positive adjusted residual marks an over-represented
#' transition, one the process prefers; a negative residual marks an avoided
#' one. Yule's Q, kappa and lift express the same association on bounded or
#' ratio scales.
#'
#' On a fit from [lsa()], a grouped fit or a multi-lag fit from [lsa_lags()], it
#' returns the estimated statistics of each transition. On the result of
#' [bootstrap_lsa()], [certainty_lsa()], [stability_lsa()], [permute_lsa()],
#' [compare_lsa()] or [bayes_compare_lsa()], it returns that method's
#' per-transition quantities. The arguments narrow and order the table:
#' `significant = TRUE` keeps the transitions the method flags (significant for
#' a fit, a permutation test or a group comparison; stable under resampling for
#' the bootstrap, certainty and stability results), `direction` keeps over- or
#' under-represented transitions, `min_count` drops rare ones, and `sort` puts
#' the strongest first.
#'
#' Use [nodes()], [tests()] and [initial()] for the other tables of a fit,
#' and [lag_profile()] to follow one transition across lags.
#'
#' @param fit An `lsa` fit, a grouped fit, an `lsa_lags` result, or the
#'   result of [bootstrap_lsa()], [certainty_lsa()], [stability_lsa()],
#'   [permute_lsa()], [compare_lsa()] or [bayes_compare_lsa()].
#' @param significant Keep only the transitions the method flags as
#'   significant or stable (default `FALSE`).
#' @param direction For fits only: `"any"` (default), `"over"` (significant
#'   with a positive residual) or `"under"` (significant with a negative
#'   residual). `"over"` and `"under"` imply `significant = TRUE`.
#' @param min_count For fits only: keep transitions observed at least this
#'   many times. Default `NULL` (no filter).
#' @param alpha For fits only: significance level. Default `NULL` uses the
#'   `alpha` of the fit.
#' @param sort Row order: `"none"` (default, table order) or `"strength"`
#'   (largest absolute effect first: the adjusted residual for a fit, a
#'   bootstrap, certainty or permutation result; the stability proportion
#'   for [stability_lsa()]; the group difference for a comparison). Fits
#'   also accept `"count"` and `"prob"` (largest first).
#' @param ... Arguments passed to the method.
#'
#' @return A `data.frame` with one row per transition. For a fit the
#'   columns are `from`, `to`, `lag`, `count`, `expected`, `prob` (row
#'   conditional, \eqn{P(\text{to} \mid \text{from})}), `prob_col` (column
#'   conditional, \eqn{P(\text{from} \mid \text{to})}), `adj_res`, `p`,
#'   `yules_q`, `kappa`, `kappa_z`, `kappa_p`, `lift` (observed over
#'   expected), `sign` (`"over"`, `"under"` or `"expected"`, comparing
#'   the count with the expected count) and `significant`; engines
#'   that compute further statistics append them (the two-cell engine adds
#'   `odds_ratio`, `log_or` and `log_or_se`). For an inference result the
#'   columns are `from`, `to` and that method's quantities, described on
#'   its help page; a pairwise comparison adds `group_a` and `group_b`. A
#'   grouped fit gains a leading `group` column, an `lsa_lags` result one
#'   row per transition and lag.
#'
#' @examples
#' fit <- lsa(engagement)
#' transitions(fit)
#' transitions(fit, significant = TRUE, sort = "strength")
#' transitions(fit, direction = "under")
#' transitions(fit, min_count = 150, sort = "count")
#'
#' # The same verb reads an inference result
#' transitions(certainty_lsa(fit), significant = TRUE)
#'
#' @seealso [lsa()], [nodes()], [tests()], [initial()],
#'   [transition_probabilities()], [lag_profile()]
#'
#' @export
transitions <- function(fit, ...) {
  UseMethod("transitions")
}

#' @rdname transitions
#' @export
transitions.lsa <- function(fit, significant = FALSE,
                            direction = c("any", "over", "under"),
                            min_count = NULL, alpha = NULL,
                            sort = c("none", "strength", "count", "prob"),
                            ...) {
  direction <- match.arg(direction)
  sort <- match.arg(sort)
  if (is.null(alpha)) alpha <- fit$params$alpha
  stopifnot(
    "`alpha` must be a single number in (0, 1)" =
      is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1
  )
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
    stopifnot(
      "`min_count` must be a single number >= 1" =
        is.numeric(min_count) && length(min_count) == 1L &&
        is.finite(min_count) && min_count >= 1
    )
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
                                           "prob"),
                                  ...) {
  direction <- match.arg(direction)
  sort <- match.arg(sort)
  .bind_group_edges(fit, transitions, significant = significant,
                    direction = direction, min_count = min_count,
                    alpha = alpha, sort = sort)
}
