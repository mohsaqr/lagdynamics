# S3 methods for the lsa class. Print is a compact, visual model summary
# (typed header, key stats, an initial-state bar chart, and the
# strongest transitions).
# summary() adds the full matrices; as.data.frame() returns the tidy
# edge table for downstream tooling.

# Horizontal unicode bar, scaled so the largest value fills `width`.
# Vectorised over `frac` (already in [0, 1] relative to the max).
.lsa_bar <- function(frac, width = 24L) {
  n <- as.integer(round(frac * width))
  vapply(n, function(k) {
    strrep("\u2588", max(0L, k))
  }, character(1L))
}

# Significance stars from a p-value (vectorised).
.lsa_stars <- function(p) {
  ifelse(!is.finite(p), "   ",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "** ",
  ifelse(p < 0.05,  "*  ", "   "))))
}

#' @export
print.lsa <- function(x, ...) {
  labels <- rownames(x$obs)
  K <- length(labels)
  alpha <- x$params$alpha
  e <- x$edges
  nsig <- sum(e$significant, na.rm = TRUE)
  dir <- if (isTRUE(x$directed)) "directed" else "undirected"

  # --- header ---
  cat(sprintf("Lag Sequential Analysis  \u2014  %s  (lag %d, %s)\n",
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
    cat(sprintf("  independence: G\u00b2 = %.1f, df = %d, p %s\n",
                lr$statistic, lr$df,
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

#' @export
summary.lsa <- function(object, ...) {
  cat("Lag Sequential Analysis\n")
  cat("=======================\n")
  print(object)

  # Node-activity bar chart: each state's share of outgoing and incoming
  # transitions. Complements the initial-state bars shown by print() and
  # gives a quick source/target profile per state.
  labels <- rownames(object$obs)
  out_share <- rowSums(object$obs) / sum(object$obs)
  in_share  <- colSums(object$obs) / sum(object$obs)
  scale <- max(out_share, in_share)
  wlab <- max(nchar(labels))
  cat("\nNode activity (share of transitions):\n")
  cat(sprintf("    %-*s  out %5.3f %-12s  in %5.3f %s\n",
              wlab, labels,
              out_share, .lsa_bar(out_share / scale, 12L),
              in_share,  .lsa_bar(in_share / scale, 12L)),
      sep = "")

  cat("\nObserved counts (obs):\n")
  print(object$obs)
  cat("\nExpected counts (exp):\n")
  print(round(object$exp, 3))
  cat("\nTransitional probabilities (prob):\n")
  print(round(object$prob, 3))
  cat("\nAdjusted residuals (adj_res):\n")
  print(round(object$adj_res, 3))
  # Return the tidy one-row summary (invisibly): summary() prints the
  # human view above, but `s <- summary(fit)` captures tidy data.
  invisible(.lsa_glance(object))
}

#' @export
as.data.frame.lsa <- function(x, row.names = NULL, optional = FALSE,
                              ...) {
  out <- x$edges
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

#' @export
as.data.frame.lsa_transitions <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  x$edges
}

#' @export
as.data.frame.lsa_group <- function(x, row.names = NULL, optional = FALSE,
                                    ...) {
  .grouped_df(x)
}

# Stack a named list of result objects into one tidy data frame, prefixing
# a `group` column carrying each element's name. Each element is coerced
# with its own as.data.frame() method, so the same helper tidies a list of
# fits, certainty results, reliability results, etc.
.grouped_df <- function(x) {
  pieces <- Map(function(df, g) {
    cbind(group = rep(g, nrow(df)), df, stringsAsFactors = FALSE)
  }, lapply(x, as.data.frame), names(x))
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

# One-row tidy overview of a fit (engine, sizes, significance counts,
# tablewise tests). Returned invisibly by summary.lsa(); also the body
# of the public fit_summary() accessor.
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

#' @importFrom stats nobs
#' @export
nobs.lsa <- function(object, ...) sum(object$obs)
