# S3 methods for the lsa class. Print is concise; summary shows the
# full picture. as.data.frame returns the tidy edge table for use with
# downstream tooling.

#' @export
print.lsa <- function(x, ...) {
  cat("<lsa>\n")
  cat(sprintf("  engine:        %s\n", x$engine))
  cat(sprintf("  lag:           %d\n", x$params$lag))
  cat(sprintf("  states (K):    %d\n", nrow(x$obs)))
  cat(sprintf("  transitions:   %d\n", sum(x$obs)))
  cat(sprintf("  labels:        %s\n",
              paste(utils::head(rownames(x$obs), 8), collapse = ", ")))
  lr <- x$lrx2
  if (!is.null(lr)) {
    cat(sprintf("  G^2 = %.3f  df = %d  p = %s\n",
                lr$statistic, lr$df, format.pval(lr$p, digits = 4)))
  }
  alpha <- x$params$alpha
  sig <- sum(x$edges$significant, na.rm = TRUE)
  cat(sprintf("  significant edges (p < %.3f): %d of %d\n",
              alpha, sig, nrow(x$edges)))
  invisible(x)
}

#' @export
summary.lsa <- function(object, ...) {
  cat("Lag Sequential Analysis\n")
  cat("=======================\n")
  print(object)
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
