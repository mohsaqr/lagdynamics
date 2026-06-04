# Multi-lag analysis: fit the same data at several lags in one call,
# producing a lag profile. Returns a named list of `lsa`
# fits (one per lag); as.data.frame() stacks their edge tables with the
# lag column so a transition's residual can be tracked across lags.

#' Lag Sequential Analysis Across Several Lags
#'
#' Fits [lsa()] at each requested lag and returns the fits together, so
#' you can compare a transition's strength across lags (a *lag
#' profile*). Each element is an ordinary `lsa` fit.
#'
#' @param data Sequence input (any form accepted by [lsa()]).
#' @param lags Integer vector of lags. May include negative lags
#'   (predecessors) and `0`. Default `1:3`.
#' @param ... Passed to [lsa()] (e.g. `engine`, `alpha`,
#'   `structural_zeros`).
#'
#' @return An object of class `c("lsa_lags", "list")`: a named list of
#'   `lsa` fits (names `"lag1"`, `"lag2"`, ...), with a `lags` attribute.
#'   [as.data.frame()] on it row-binds [transitions()] of every fit (each
#'   already carries its `lag` column) into one tidy long frame with the
#'   same columns as `transitions()`.
#'
#' @examples
#' prof <- lsa_lags(engagement, lags = 1:3)
#' prof
#' # Track one transition across lags:
#' d <- as.data.frame(prof)
#' d[d$from == "Active" & d$to == "Average",
#'   c("lag", "count", "adj_res", "p")]
#'
#' @seealso [lsa()]
#'
#' @export
lsa_lags <- function(data, lags = 1:3, ...) {
  stopifnot(is.numeric(lags), length(lags) >= 1L,
            all(lags == round(lags)))
  lags <- as.integer(lags)
  fits <- lapply(lags, function(L) lsa(data, lag = L, ...))
  names(fits) <- paste0("lag", lags)
  structure(fits, lags = lags, class = c("lsa_lags", "list"))
}

#' @export
print.lsa_lags <- function(x, ...) {
  lags <- attr(x, "lags")
  cat("<lsa_lags>\n")
  cat(sprintf("  engine: %s\n", x[[1L]]$method))
  cat(sprintf("  lags:   %s\n", paste(lags, collapse = ", ")))
  sig <- vapply(x, function(f) sum(f$edges$significant, na.rm = TRUE),
                integer(1L))
  ne <- vapply(x, function(f) nrow(f$edges), integer(1L))
  cat(sprintf("    lag %-3d  %d of %d transitions significant\n",
              lags, sig, ne), sep = "")
  invisible(x)
}

#' @export
as.data.frame.lsa_lags <- function(x, row.names = NULL, optional = FALSE,
                                   ...) {
  out <- do.call(rbind, lapply(unclass(x), transitions))
  rownames(out) <- row.names
  out
}

#' Lag Profile of a Single Transition
#'
#' How one `from -> to` transition behaves across lags, as a tidy
#' one-row-per-lag data frame. A clean shortcut for "track this
#' transition over lags 1, 2, 3, ...".
#'
#' @param x Sequence input (any form accepted by [lsa()]) or an existing
#'   [lsa_lags()] object.
#' @param from,to State labels of the transition to profile.
#' @param lags Integer vector of lags. Default `1:3`. Ignored when `x`
#'   is already an `lsa_lags` object.
#' @param ... Passed to [lsa_lags()] when `x` is raw data.
#'
#' @return A tidy `data.frame`, one row per lag, with columns `lag`,
#'   `from`, `to`, `count`, `prob`, `adj_res`, `p`, and `significant`.
#'
#' @examples
#' lag_profile(group_regulation, "plan", "consensus", lags = 1:3)
#'
#' @seealso [lsa_lags()]
#'
#' @export
lag_profile <- function(x, from, to, lags = 1:3, ...) {
  fits <- if (inherits(x, "lsa_lags")) x else lsa_lags(x, lags = lags, ...)
  d <- as.data.frame(fits)
  if (!from %in% d$from || !to %in% d$to) {
    stop(sprintf("Transition %s -> %s not found; states are: %s",
                 from, to,
                 paste(sort(unique(d$from)), collapse = ", ")),
         call. = FALSE)
  }
  out <- d[d$from == from & d$to == to,
           c("lag", "from", "to", "count", "prob",
             "adj_res", "p", "significant")]
  rownames(out) <- NULL
  out
}
