# S3 methods for the lsa class. Print is a compact, visual summary in
# the style of the tna / Nestimate model prints (typed header, key
# stats, an initial-state bar chart, and the strongest transitions).
# summary() adds the full matrices; as.data.frame() returns the tidy
# edge table for downstream tooling.

# Horizontal unicode bar, scaled so the largest value fills `width`.
# Vectorised over `frac` (already in [0, 1] relative to the max).
.lsa_bar <- function(frac, width = 24L) {
  n <- as.integer(round(frac * width))
  vapply(n, function(k) {
    strrep("█", max(0L, k))
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
  cat(sprintf("Lag Sequential Analysis  —  %s  (lag %d, %s)\n",
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
    cat(sprintf("  independence: G² = %.1f, df = %d, p %s\n",
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

  # --- initial-state distribution bar chart (tna / Nestimate style) ---
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
  invisible(object)
}

#' @export
as.data.frame.lsa <- function(x, row.names = NULL, optional = FALSE,
                              ...) {
  out <- x$edges
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

#' @importFrom stats nobs
#' @export
nobs.lsa <- function(object, ...) sum(object$obs)
