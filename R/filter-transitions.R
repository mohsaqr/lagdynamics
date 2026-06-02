# Tidy filter helpers over an lsa fit's edge table. Each helper returns
# a data.frame with the same columns as fit$edges, so they compose with
# base R subsetting and other tidy operations.
#
# Each is an S3 generic: the `.lsa` method filters a single fit; the
# `.lsa_group` method maps that filter over the per-group fits and
# row-binds the results with a leading `group` column. The grouped
# output stays a plain base data.frame (not a tibble) to match the
# package's base-R house style.

# Row-bind a per-group filter into one long data.frame with a leading
# `group` column. Groups that contribute no rows are dropped; if every
# group is empty, a zero-row frame with the correct columns is returned
# so callers can rely on a stable shape.
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

#' Filter Transitions by Significance
#'
#' Returns the subset of `fit$edges` whose adjusted-residual p-value
#' falls below `alpha`. Cells with non-finite p-values (structural
#' zeros, zero-margin rows, or non-estimable cells) are excluded.
#'
#' @param fit An `lsa` object returned by [lsa()].
#' @param alpha Numeric scalar in `(0, 1)`. Significance threshold.
#'   Default `NULL`, which uses the alpha snapshotted on the fit
#'   (`fit$params$alpha`).
#'
#' @return A data.frame with the same columns as `fit$edges`,
#'   containing only the significant rows. Row order is preserved.
#'
#' @examples
#' fit <- lsa(engagement, engine = "classical")
#' significant_transitions(fit)
#' significant_transitions(fit, alpha = 0.01)
#'
#' For a grouped fit (`lsa_group`), the filter is applied within each
#' group and the results are row-bound into one data.frame with a
#' leading `group` column.
#'
#' @seealso [overrepresented_transitions()],
#'   [underrepresented_transitions()], [common_transitions()]
#'
#' @export
significant_transitions <- function(fit, ...) {
  UseMethod("significant_transitions")
}

#' @rdname significant_transitions
#' @export
significant_transitions.lsa <- function(fit, alpha = NULL, ...) {
  if (is.null(alpha)) alpha <- fit$params$alpha
  stopifnot(is.numeric(alpha), length(alpha) == 1L,
            alpha > 0, alpha < 1)
  e <- fit$edges
  e[is.finite(e$p) & e$p < alpha, , drop = FALSE]
}

#' @rdname significant_transitions
#' @export
significant_transitions.lsa_group <- function(fit, alpha = NULL, ...) {
  .bind_group_edges(fit, significant_transitions, alpha = alpha, ...)
}

#' Filter Overrepresented (Positive-Residual) Transitions
#'
#' Returns significant transitions where the observed count exceeds
#' the expected count under independence (i.e. positive adjusted
#' residual). These are the cells where the focal transition occurs
#' more often than chance.
#'
#' @inheritParams significant_transitions
#'
#' @return A data.frame with the same columns as `fit$edges`.
#'
#' @examples
#' fit <- lsa(engagement, engine = "classical")
#' overrepresented_transitions(fit)
#'
#' @seealso [significant_transitions()],
#'   [underrepresented_transitions()]
#'
#' @export
overrepresented_transitions <- function(fit, ...) {
  UseMethod("overrepresented_transitions")
}

#' @rdname overrepresented_transitions
#' @export
overrepresented_transitions.lsa <- function(fit, alpha = NULL, ...) {
  e <- significant_transitions(fit, alpha = alpha)
  e[is.finite(e$adj_res) & e$adj_res > 0, , drop = FALSE]
}

#' @rdname overrepresented_transitions
#' @export
overrepresented_transitions.lsa_group <- function(fit, alpha = NULL,
                                                   ...) {
  .bind_group_edges(fit, overrepresented_transitions, alpha = alpha, ...)
}

#' Filter Underrepresented (Negative-Residual) Transitions
#'
#' Returns significant transitions where the observed count is below
#' the expected count under independence (negative adjusted residual).
#' These are the cells where the focal transition is actively avoided
#' relative to chance.
#'
#' @inheritParams significant_transitions
#'
#' @return A data.frame with the same columns as `fit$edges`.
#'
#' @examples
#' fit <- lsa(engagement, engine = "classical")
#' underrepresented_transitions(fit)
#'
#' @seealso [significant_transitions()],
#'   [overrepresented_transitions()]
#'
#' @export
underrepresented_transitions <- function(fit, ...) {
  UseMethod("underrepresented_transitions")
}

#' @rdname underrepresented_transitions
#' @export
underrepresented_transitions.lsa <- function(fit, alpha = NULL, ...) {
  e <- significant_transitions(fit, alpha = alpha)
  e[is.finite(e$adj_res) & e$adj_res < 0, , drop = FALSE]
}

#' @rdname underrepresented_transitions
#' @export
underrepresented_transitions.lsa_group <- function(fit, alpha = NULL,
                                                    ...) {
  .bind_group_edges(fit, underrepresented_transitions, alpha = alpha,
                    ...)
}

#' Filter Frequently-Observed Transitions
#'
#' Returns transitions whose observed count is at least `min_count`.
#' This is a descriptive volume filter independent of statistical
#' significance, useful for trimming long-tail noise from a transition
#' network before visualisation.
#'
#' @param fit An `lsa` object returned by [lsa()].
#' @param min_count Integer scalar `>= 1`. Minimum observed count.
#'   Default `1L` (keep every transition that was observed at least
#'   once).
#'
#' @return A data.frame with the same columns as `fit$edges`.
#'
#' @examples
#' fit <- lsa(engagement, engine = "classical")
#' common_transitions(fit, min_count = 3)
#'
#' @seealso [significant_transitions()]
#'
#' @export
common_transitions <- function(fit, ...) {
  UseMethod("common_transitions")
}

#' @rdname common_transitions
#' @export
common_transitions.lsa <- function(fit, min_count = 1L, ...) {
  stopifnot(is.numeric(min_count), length(min_count) == 1L,
            min_count >= 1L)
  e <- fit$edges
  e[is.finite(e$count) & e$count >= min_count, , drop = FALSE]
}

#' @rdname common_transitions
#' @export
common_transitions.lsa_group <- function(fit, min_count = 1L, ...) {
  .bind_group_edges(fit, common_transitions, min_count = min_count, ...)
}
