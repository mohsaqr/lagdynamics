# Print methods for the three inference S3 classes.

#' Print and read resampling results
#'
#' @description
#' `print()` and `as.data.frame()` read the results of [bootstrap_lsa()],
#' [permute_lsa()] and [stability_lsa()]. `print()` shows a short summary: the
#' engine, the number of replicates, the resampling settings and how many
#' transitions are stable or significant (with the number of failed replicates
#' for a stability result). `as.data.frame()` returns the per-transition table
#' of the result.
#'
#' [transitions()] returns the same table and adds the arguments
#' `significant` (keep the transitions the method flags) and `sort`
#' (order by the size of the effect); it is the usual way to read these
#' results. Results of [certainty_lsa()] use the bootstrap methods.
#'
#' @param x An `lsa_bootstrap`, `lsa_permutation` or `lsa_stability`
#'   object.
#' @param row.names,optional Not used; present for consistency with the
#'   generic.
#' @param ... Not used.
#'
#' @return `print()` returns `x` invisibly. `as.data.frame()` returns a
#'   data frame with one row per transition (`from -> to` pair); its
#'   columns are listed under *Value* in [bootstrap_lsa()],
#'   [permute_lsa()] and [stability_lsa()].
#'
#' @examples
#' fit <- lsa(engagement)
#' st <- stability_lsa(fit, R = 20)
#' st
#' head(as.data.frame(st))
#'
#' @seealso [transitions()], [bootstrap_lsa()], [permute_lsa()],
#'   [stability_lsa()].
#'
#' @name lsa-inference-methods
NULL

#' @rdname lsa-inference-methods
#' @export
print.lsa_bootstrap <- function(x, ...) {
  cat("<lsa_bootstrap>\n")
  cat(sprintf("  engine:        %s\n", x$fit$method))
  cat(sprintf("  level:         %s\n", x$level))
  cat(sprintf("  replicates:    %d\n", x$R))
  cat(sprintf("  CI level:      %.0f%%\n", 100 * x$level_alpha))
  n_stable <- sum(x$edges$adj_res_stable, na.rm = TRUE)
  cat(sprintf("  stable edges:  %d of %d\n",
              n_stable, nrow(x$edges)))
  invisible(x)
}

#' @rdname lsa-inference-methods
#' @export
print.lsa_permutation <- function(x, ...) {
  cat("<lsa_permutation>\n")
  cat(sprintf("  engine:        %s\n", x$fit$method))
  cat(sprintf("  replicates:    %d\n", x$R))
  cat(sprintf("  within seq:    %s\n", x$within_sequence))
  alpha <- x$fit$params$alpha
  n_sig <- sum(x$edges$significant, na.rm = TRUE)
  cat(sprintf("  significant edges (p_perm < %.3f): %d of %d\n",
              alpha, n_sig, nrow(x$edges)))
  invisible(x)
}

#' @rdname lsa-inference-methods
#' @export
print.lsa_stability <- function(x, ...) {
  cat("<lsa_stability>\n")
  cat(sprintf("  engine:        %s\n", x$fit$method))
  cat(sprintf("  replicates:    %d\n", x$R))
  cat(sprintf("  proportion:    %.0f%%\n", 100 * x$proportion))
  cat(sprintf("  min stable:    %.0f%%\n", 100 * x$min_stable))
  cat(sprintf("  failed:        %d of %d replicates\n",
              x$n_failed %||% 0L, x$R))
  n_stable <- sum(x$edges$stable, na.rm = TRUE)
  cat(sprintf("  stable edges:  %d of %d (>= %.0f%% across replicates)\n",
              n_stable, nrow(x$edges), 100 * x$min_stable))
  invisible(x)
}

#' @rdname lsa-inference-methods
#' @export
as.data.frame.lsa_bootstrap <- function(x, row.names = NULL,
                                        optional = FALSE, ...) {
  x$edges
}

#' @rdname lsa-inference-methods
#' @export
as.data.frame.lsa_permutation <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  x$edges
}

#' @rdname lsa-inference-methods
#' @export
as.data.frame.lsa_stability <- function(x, row.names = NULL,
                                        optional = FALSE, ...) {
  x$edges
}
