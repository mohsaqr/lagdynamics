# Print methods for the three inference S3 classes.

#' @export
print.lsa_bootstrap <- function(x, ...) {
  cat("<lsa_bootstrap>\n")
  cat(sprintf("  engine:        %s\n", x$fit$engine))
  cat(sprintf("  level:         %s\n", x$level))
  cat(sprintf("  replicates:    %d\n", x$R))
  cat(sprintf("  CI level:      %.0f%%\n", 100 * x$level_alpha))
  n_stable <- sum(x$edges$adj_res_stable, na.rm = TRUE)
  cat(sprintf("  stable edges:  %d of %d\n",
              n_stable, nrow(x$edges)))
  invisible(x)
}

#' @export
print.lsa_permutation <- function(x, ...) {
  cat("<lsa_permutation>\n")
  cat(sprintf("  engine:        %s\n", x$fit$engine))
  cat(sprintf("  replicates:    %d\n", x$R))
  cat(sprintf("  within seq:    %s\n", x$within_sequence))
  alpha <- x$fit$params$alpha
  n_sig <- sum(x$edges$significant, na.rm = TRUE)
  cat(sprintf("  significant edges (p_perm < %.3f): %d of %d\n",
              alpha, n_sig, nrow(x$edges)))
  invisible(x)
}

#' @export
print.lsa_stability <- function(x, ...) {
  cat("<lsa_stability>\n")
  cat(sprintf("  engine:        %s\n", x$fit$engine))
  cat(sprintf("  replicates:    %d\n", x$R))
  cat(sprintf("  proportion:    %.0f%%\n", 100 * x$proportion))
  cat(sprintf("  min stable:    %.0f%%\n", 100 * x$min_stable))
  n_stable <- sum(x$edges$stable, na.rm = TRUE)
  cat(sprintf("  stable edges:  %d of %d (>= %.0f%% across replicates)\n",
              n_stable, nrow(x$edges), 100 * x$min_stable))
  invisible(x)
}

#' @export
as.data.frame.lsa_bootstrap <- function(x, row.names = NULL,
                                        optional = FALSE, ...) {
  x$edges
}

#' @export
as.data.frame.lsa_permutation <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  x$edges
}

#' @export
as.data.frame.lsa_stability <- function(x, row.names = NULL,
                                        optional = FALSE, ...) {
  x$edges
}
