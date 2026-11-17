# Multi-lag analysis: fit the same data at several lags in one call,
# producing a lag profile. Returns a named list of `lsa`
# fits (one per lag); as.data.frame() stacks their edge tables with the
# lag column so a transition's residual can be tracked across lags.

#' Fit lag sequential analysis models at several lags
#'
#' @description
#' `lsa_lags()` fits [lsa()] once per requested lag, with the same engine and
#' options, and keeps the fits together. The result is a lag profile in the
#' sense of Sackett (1979), who proposed examining a transition across a range
#' of lags to separate short-range from longer-range dependence and to detect
#' cycles: a sequential contingency need not be immediate, a state may tend to
#' recur two or three events later, and an association that is strong at lag 1
#' may fade with distance.
#'
#' [transitions()] on the result stacks the per-lag transition tables into one
#' table with a `lag` column, with the same filtering and sorting arguments as
#' for a single fit. To follow one transition across lags, use [lag_profile()].
#'
#' @param data Sequences in any form accepted by [lsa()].
#' @param lags Integer vector of lags (default `1:3`). Negative lags look
#'   back to predecessors.
#' @param ... Passed to [lsa()] (for example `engine`, `alpha`, `loops`).
#'   Ignored by the `print()` and `as.data.frame()` methods.
#' @param x An `lsa_lags` object.
#' @param row.names Optional row names for the stacked data frame.
#' @param optional Ignored; present for consistency with the generic.
#'
#' @return An object of class `c("lsa_lags", "list")`, one `lsa` fit per
#'   lag. [transitions()] and `as.data.frame()` return a `data.frame` with
#'   one row per transition and lag and the columns of [transitions()] for
#'   a single fit. `print()` lists the number of significant transitions at
#'   each lag and returns `x` invisibly.
#'
#' @examples
#' prof <- lsa_lags(engagement, lags = 1:3)
#' prof
#' transitions(prof, direction = "over", sort = "strength")
#'
#' @references
#' Sackett, G. P. (1979). The lag sequential analysis of contingency and
#' cyclicity in behavioral interaction research. In J. D. Osofsky (Ed.),
#' \emph{Handbook of infant development} (pp. 623-649). Wiley.
#'
#' @seealso [lag_profile()] for one transition across lags; [lsa()];
#'   [transitions()].
#'
#' @export
lsa_lags <- function(data, lags = 1:3, ...) {
  stopifnot(
    "`lags` must be a vector of whole numbers" =
      is.numeric(lags) && length(lags) >= 1L && all(is.finite(lags)) &&
      all(lags == round(lags))
  )
  lags <- as.integer(lags)
  fits <- lapply(lags, function(L) lsa(data, lag = L, ...))
  names(fits) <- paste0("lag", lags)
  structure(fits, lags = lags, class = c("lsa_lags", "list"))
}

#' @rdname lsa_lags
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

#' @rdname lsa_lags
#' @export
as.data.frame.lsa_lags <- function(x, row.names = NULL, optional = FALSE,
                                   ...) {
  out <- do.call(rbind, lapply(unclass(x), transitions))
  rownames(out) <- row.names
  out
}

#' Trace one transition across lags
#'
#' @description
#' `lag_profile()` traces one `from -> to` transition across lags, one row per
#' lag with the count, the transition probability, the adjusted residual and its
#' p-value. This lag profile shows how the association between two states
#' changes with the distance between them: whether a transition that is
#' over-represented at lag 1 persists, decays or reverses at lags 2, 3 and
#' beyond (Sackett, 1979). For example, it shows whether disengagement is still
#' unlikely several time points after an active one, or only at the next.
#'
#' `lag_profile()` fits the lags itself from raw data, or reads an existing
#' [lsa_lags()] result.
#'
#' @param x Sequences in any form accepted by [lsa()], or an [lsa_lags()]
#'   result.
#' @param from,to State names of the transition.
#' @param lags Integer vector of lags (default `1:3`); ignored when `x` is
#'   an `lsa_lags` result.
#' @param ... Passed to [lsa_lags()] when `x` is raw data.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`
#' when `from` or `to` is not a state of the data.
#'
#' @return A `data.frame` with one row per lag and the columns `lag`,
#'   `from`, `to`, `count`, `prob`, `adj_res`, `p` and `significant`.
#'
#' @examples
#' lag_profile(engagement, from = "Active", to = "Disengaged", lags = 1:4)
#'
#' @references
#' Sackett, G. P. (1979). The lag sequential analysis of contingency and
#' cyclicity in behavioral interaction research. In J. D. Osofsky (Ed.),
#' \emph{Handbook of infant development} (pp. 623-649). Wiley.
#'
#' @seealso [lsa_lags()], [transitions()]
#'
#' @export
lag_profile <- function(x, from, to, lags = 1:3, ...) {
  fits <- if (inherits(x, "lsa_lags")) x else lsa_lags(x, lags = lags, ...)
  d <- as.data.frame(fits)
  if (!from %in% d$from || !to %in% d$to) {
    .lsa_stop(sprintf("Transition %s -> %s not found; states are: %s",
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
