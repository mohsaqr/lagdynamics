# S3 methods for the lsa class. Print is a compact, visual model summary
# (typed header, key stats, an initial-state bar chart, and the
# strongest transitions).
# summary() returns a "summary.lsa" object (a one-row overview and a
# per-state table) with its own print method; as.data.frame() returns the
# tidy edge table for downstream tooling.

# Horizontal bar, scaled so the largest value fills `width`. Vectorised
# over `frac` (already in [0, 1] relative to the max). Uses the full-block
# glyph in a UTF-8 session and "#" otherwise.
.lsa_bar <- function(frac, width = 24L) {
  n <- as.integer(round(frac * width))
  glyph <- .lsa_glyph("\u2588", "#")
  vapply(n, function(k) strrep(glyph, max(0L, k)), character(1L))
}

# Significance stars from a p-value (vectorised).
.lsa_stars <- function(p) {
  ifelse(!is.finite(p), "   ",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "** ",
  ifelse(p < 0.05,  "*  ", "   "))))
}

#' @rdname lsa-methods
#' @export
print.lsa <- function(x, ...) {
  labels <- rownames(x$obs)
  K <- length(labels)
  alpha <- x$params$alpha
  e <- x$edges
  nsig <- sum(e$significant, na.rm = TRUE)
  dir <- if (isTRUE(x$directed)) "directed" else "undirected"

  # --- header ---
  cat(sprintf("Lag Sequential Analysis  -  %s  (lag %d, %s)\n",
              x$method, x$params$lag, dir))
  if (identical(x$data$source, "events")) {
    cat(sprintf("  %d states | %d transitions | %d events | %d sequences\n",
                K, sum(x$obs), x$data$n_events, x$data$n_sequences))
  } else {
    cat(sprintf("  %d states | %d transitions (from transition matrix)\n",
                K, sum(x$obs)))
  }
  cat(sprintf("  states: %s\n",
              paste(utils::head(labels, 12), collapse = ", ")))
  lr <- x$lrx2
  if (!is.null(lr) && is.finite(lr$statistic)) {
    cat(sprintf("  independence: %s = %.1f, df = %d, p %s\n",
                .lsa_glyph("G\u00b2", "G^2"), lr$statistic, lr$df,
                format.pval(lr$p, digits = 3, eps = 2.2e-16)))
  }

  # --- strongest significant over-represented transitions ---
  cat(sprintf("\n  Significant transitions (p < %g): %d of %d\n",
              alpha, nsig, nrow(e)))
  sig <- e[is.finite(e$p) & e$significant & is.finite(e$adj_res) &
             e$adj_res > 0, , drop = FALSE]
  if (nrow(sig)) {
    sig <- sig[order(-sig$adj_res), , drop = FALSE]
    top <- utils::head(sig, 5L)
    w <- max(nchar(top$edge))
    cat(sprintf("  strongest over-represented (of %d):\n", nrow(sig)))
    cat(sprintf("    %-*s  z = %+6.1f  %s\n",
                w, top$edge, top$adj_res, .lsa_stars(top$p)),
        sep = "")
    if (nrow(sig) > nrow(top)) {
      cat(sprintf("    ... and %d more\n", nrow(sig) - nrow(top)))
    }
  }

  # --- initial-state distribution bar chart ---
  if (!is.null(x$inits) && any(x$inits > 0)) {
    cat("\n  Initial states:\n")
    ord <- order(-x$inits)
    nm <- names(x$inits)[ord]
    pv <- as.numeric(x$inits)[ord]
    wlab <- max(nchar(nm))
    cat(sprintf("    %-*s %5.3f  %s\n",
                wlab, nm, pv, .lsa_bar(pv / max(pv))),
        sep = "")
  }
  invisible(x)
}

#' @rdname lsa-methods
#' @export
summary.lsa <- function(object, ...) {
  structure(
    list(overview = .lsa_glance(object),
         nodes = .lsa_node_summary(object)),
    class = "summary.lsa"
  )
}

# Per-state profile of a fit: transition totals, their shares, and the
# initial-state probability. One row per state, in label order.
.lsa_node_summary <- function(fit) {
  obs <- fit$obs
  total <- sum(obs)
  labels <- rownames(obs)
  init <- if (is.null(fit$inits)) rep(NA_real_, length(labels))
          else as.numeric(fit$inits[labels])
  data.frame(
    state     = labels,
    outgoing  = as.numeric(rowSums(obs)),
    incoming  = as.numeric(colSums(obs)),
    out_share = if (total > 0) as.numeric(rowSums(obs)) / total
                else rep(NA_real_, length(labels)),
    in_share  = if (total > 0) as.numeric(colSums(obs)) / total
                else rep(NA_real_, length(labels)),
    init_prob = init,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @rdname lsa-methods
#' @export
print.summary.lsa <- function(x, ...) {
  o <- x$overview
  cat(sprintf("Lag sequential analysis: %s engine, lag %d\n",
              o$engine, o$lag))
  sizes <- if (is.na(o$n_sequences)) {
    sprintf("  %d states, %d transitions (from a transition matrix)\n",
            o$n_states, o$n_transitions)
  } else {
    sprintf("  %d states, %d transitions, %d events, %d sequences\n",
            o$n_states, o$n_transitions, o$n_events, o$n_sequences)
  }
  cat(sizes)
  cat(sprintf("  significant transitions (p < %g): %d of %d\n",
              o$alpha, o$n_significant, o$n_states^2))
  if (is.finite(o$lrx2)) {
    cat(sprintf("  likelihood-ratio test: %s = %.2f, df = %d, p %s\n",
                .lsa_glyph("G\u00b2", "G^2"), o$lrx2, o$lrx2_df,
                format.pval(o$lrx2_p, digits = 3, eps = 2.2e-16)))
  }
  if (is.finite(o$x2)) {
    cat(sprintf("  Pearson test: %s = %.2f, df = %d, p %s\n",
                .lsa_glyph("\u03c7\u00b2", "X^2"), o$x2, o$x2_df,
                format.pval(o$x2_p, digits = 3, eps = 2.2e-16)))
  }

  n <- x$nodes
  scale <- max(c(n$out_share, n$in_share), na.rm = TRUE)
  wlab <- max(nchar(n$state), 5L)
  cat("\n  State activity (share of transitions out and in):\n")
  # Pad the bars by character count: sprintf("%-12s") counts bytes, and
  # the UTF-8 block glyph is three bytes wide.
  out_bar <- .lsa_bar(n$out_share / scale, 12L)
  out_bar <- paste0(out_bar, strrep(" ", 12L - nchar(out_bar)))
  cat(sprintf("    %-*s  out %5.3f %s  in %5.3f %s\n",
              wlab, n$state,
              n$out_share, out_bar,
              n$in_share, .lsa_bar(n$in_share / scale, 12L)),
      sep = "")
  invisible(x)
}

#' @rdname lsa-methods
#' @param what For `as.data.frame()` on a summary: `"overview"` (default,
#'   the one-row overview of the fit) or `"nodes"` (one row per state).
#' @export
as.data.frame.summary.lsa <- function(x, row.names = NULL,
                                      optional = FALSE,
                                      what = c("overview", "nodes"), ...) {
  what <- match.arg(what)
  out <- x[[what]]
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

#' @rdname lsa-methods
#' @export
as.data.frame.lsa <- function(x, row.names = NULL, optional = FALSE,
                              ...) {
  out <- x$edges
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

#' @rdname lsa_transitions
#' @export
as.data.frame.lsa_transitions <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  x$edges
}

#' @rdname lsa-methods
#' @export
as.data.frame.lsa_group <- function(x, row.names = NULL, optional = FALSE,
                                    ...) {
  .grouped_df(x)
}

# Stack a named list of result objects into one tidy data frame, prefixing
# a `group` column carrying each element's name. Each element is coerced
# with its own as.data.frame() method, so the same helper tidies a list of
# fits, certainty results, reliability results, etc.
.grouped_df <- function(x, ...) {
  pieces <- Map(function(df, g) {
    cbind(group = rep(g, nrow(df)), df, stringsAsFactors = FALSE)
  }, lapply(x, as.data.frame, ...), names(x))
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

# One-row tidy overview of a fit (engine, sizes, significance counts,
# tablewise tests). The `overview` table of summary.lsa(), and one row
# per group of summary.lsa_group().
.lsa_glance <- function(fit) {
  e <- fit$edges
  d <- fit$data
  scal <- function(v) if (is.null(v)) NA_real_ else as.numeric(v)
  data.frame(
    engine        = fit$method,
    lag           = fit$params$lag,
    n_states      = nrow(fit$obs),
    n_sequences   = if (is.null(d$n_sequences)) NA_integer_
                    else as.integer(d$n_sequences),
    n_events      = if (is.null(d$n_events)) NA_integer_
                    else as.integer(d$n_events),
    n_transitions = sum(fit$obs),
    n_significant = sum(e$significant, na.rm = TRUE),
    alpha         = fit$params$alpha,
    lrx2          = scal(fit$lrx2$statistic),
    lrx2_df       = if (is.null(fit$lrx2)) NA_integer_
                    else as.integer(fit$lrx2$df),
    lrx2_p        = scal(fit$lrx2$p),
    x2            = scal(fit$x2$statistic),
    x2_df         = if (is.null(fit$x2)) NA_integer_
                    else as.integer(fit$x2$df),
    x2_p          = scal(fit$x2$p),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @rdname lsa-methods
#' @export
summary.lsa_group <- function(object, ...) {
  # Grouped summary: one tidy row per group (no rich per-group dump).
  parts <- Map(function(f, g) cbind(group = g, .lsa_glance(f),
                                    stringsAsFactors = FALSE),
               unclass(object), names(object))
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' @rdname lsa-methods
#' @importFrom stats nobs
#' @export
nobs.lsa <- function(object, ...) sum(object$obs)
