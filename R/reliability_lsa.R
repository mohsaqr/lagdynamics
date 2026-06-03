# Split-half network reliability for an lsa fit. Implements the
# random-half resampling described in Epskamp et al. (2018) for
# psychometric networks, adapted to LSA edge matrices.
#
# Each replicate: partition sequences uniformly at random into halves
# H1 and H2, refit lsa on each half, correlate the two edge-weight
# vectors. The distribution of these correlations across R replicates
# is the reliability summary.

#' Split-Half Reliability for an LSA Fit
#'
#' Estimates the reliability of an LSA network by repeated random
#' split-half resampling of sequences. Each replicate draws two
#' disjoint halves of the sequences without replacement, refits the
#' engine on each half, and computes the correlation between the two
#' half-network edge-weight vectors. Returns the distribution of
#' replicate correlations plus a point summary.
#'
#' @param fit An `lsa` object returned by [lsa()]. Must be built from
#'   event-level data (sequences), not from a pre-computed transition
#'   matrix.
#' @param R Integer. Number of split-half replicates. Default `100`.
#' @param weights Character. Which edge matrix to correlate across
#'   halves: `"prob"` (default), `"count"`, or `"adj_res"`.
#' @param method Character. Correlation method: `"pearson"` (default)
#'   or `"spearman"`.
#' @param parallel Logical. Use multi-core resampling. Default `FALSE`.
#' @param n_cores Integer. Worker count when `parallel = TRUE`.
#' @param verbose Logical. Print progress every 100 replicates.
#' @param ... Reserved.
#'
#' @return An object of class `c("lsa_reliability", "list")` with:
#' \describe{
#'   \item{correlations}{Numeric vector of length `R`: the split-half
#'     correlation of each replicate.}
#'   \item{mean, sd}{Mean and standard deviation of the finite
#'     replicate correlations.}
#'   \item{ci_low, ci_high}{Empirical 2.5% and 97.5% quantiles.}
#'   \item{R, weights, method, n_sequences}{Recipe metadata.}
#'   \item{fit}{Reference to the original fit.}
#' }
#'
#' @examples
#' \donttest{
#' fit <- lsa(engagement, engine = "classical")
#' rel <- reliability_lsa(fit, R = 50)
#' rel
#' }
#'
#' @references
#' Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#' psychological networks and their accuracy: A tutorial paper.
#' \emph{Behavior Research Methods, 50}(1), 195-212.
#'
#' For a grouped fit (`lsa_group`), reliability is estimated separately
#' within each group and the per-group `lsa_reliability` objects are
#' returned in an `lsa_reliability_group` container with its own print
#' method.
#'
#' @seealso [bootstrap_lsa()], [stability_lsa()], [permute_lsa()]
#'
#' @export
reliability_lsa <- function(fit, ...) UseMethod("reliability_lsa")

#' @rdname reliability_lsa
#' @export
reliability_lsa.lsa <- function(fit,
                                R = 100L,
                                weights = c("prob", "count", "adj_res"),
                                method = c("pearson", "spearman"),
                                parallel = FALSE,
                                n_cores = NULL,
                                verbose = FALSE,
                                ...) {
  weights <- match.arg(weights)
  method <- match.arg(method)
  stopifnot(is.numeric(R), length(R) == 1L, R >= 1L)
  R <- as.integer(R)

  recipe <- fit$params
  d <- fit$data
  if (identical(d$source, "transitions")) {
    stop("reliability_lsa() requires event-level input. The fit was ",
         "built from a pre-computed transition matrix.", call. = FALSE)
  }
  per_seq <- split(d$events, d$seq_id)
  S <- length(per_seq)
  if (S < 2L) {
    stop("reliability_lsa() needs at least 2 sequences; got ", S, ".",
         call. = FALSE)
  }

  pull_weights <- function(refitted) {
    switch(weights,
      prob    = refitted$prob,
      count   = refitted$obs,
      adj_res = refitted$adj_res
    )
  }

  worker <- function(b) {
    half <- max(1L, S %/% 2L)
    idx_1 <- sample.int(S, half, replace = FALSE)
    idx_2 <- setdiff(seq_len(S), idx_1)
    if (length(idx_2) == 0L) return(NA_real_)

    refit_half <- function(idx) {
      hs <- per_seq[idx]
      ev <- unlist(hs, use.names = FALSE)
      sid <- rep.int(seq_along(idx),
                     times = vapply(hs, length, integer(1)))
      # A half made up only of singleton sequences contributes zero
      # transitions and the engine errors. That is a property of the
      # random split, not of the fit, so the replicate returns NA
      # rather than crashing the whole run.
      tryCatch(
        .refit_from_events(events = ev, seq_id = sid,
                           labels = d$labels, recipe = recipe),
        error = function(e) NULL
      )
    }

    f1 <- refit_half(idx_1)
    f2 <- refit_half(idx_2)
    if (is.null(f1) || is.null(f2)) return(NA_real_)
    w1 <- as.vector(pull_weights(f1))
    w2 <- as.vector(pull_weights(f2))
    ok <- is.finite(w1) & is.finite(w2)
    if (sum(ok) < 2L || stats::sd(w1[ok]) == 0 ||
        stats::sd(w2[ok]) == 0) {
      return(NA_real_)
    }
    suppressWarnings(stats::cor(w1[ok], w2[ok], method = method))
  }

  if (verbose) message("Running ", R, " split-half replicates ...")
  results <- .run_parallel(worker, R = R, parallel = parallel,
                           n_cores = n_cores, verbose = verbose)
  cors <- vapply(results, function(v) {
    if (is.numeric(v) && length(v) == 1L) v else NA_real_
  }, numeric(1))

  finite <- is.finite(cors)
  ci <- if (sum(finite) >= 2L) {
    stats::quantile(cors[finite], c(0.025, 0.975), names = FALSE)
  } else {
    c(NA_real_, NA_real_)
  }
  structure(
    list(
      correlations = cors,
      mean         = if (any(finite)) mean(cors[finite]) else NA_real_,
      sd           = if (sum(finite) >= 2L) stats::sd(cors[finite])
                     else NA_real_,
      ci_low       = ci[1L],
      ci_high      = ci[2L],
      R            = R,
      weights      = weights,
      method       = method,
      n_sequences  = S,
      fit          = fit
    ),
    class = c("lsa_reliability", "list")
  )
}

#' @rdname reliability_lsa
#' @export
reliability_lsa.lsa_group <- function(fit, ...) {
  rels <- lapply(fit, function(f) reliability_lsa(f, ...))
  names(rels) <- names(fit)
  structure(rels, levels = names(fit),
            class = c("lsa_reliability_group", "list"))
}

#' @export
print.lsa_reliability_group <- function(x, ...) {
  cat("<lsa_reliability_group>\n")
  cat(sprintf("  groups: %d\n\n", length(x)))
  for (i in seq_along(x)) {
    cat(sprintf("[%s]\n", names(x)[i]))
    print(x[[i]])
    if (i < length(x)) cat("\n")
  }
  invisible(x)
}

#' @export
print.lsa_reliability <- function(x, ...) {
  cat("<lsa_reliability>\n")
  cat(sprintf("  engine:        %s\n", x$fit$method))
  cat(sprintf("  replicates:    %d\n", x$R))
  cat(sprintf("  weights:       %s\n", x$weights))
  cat(sprintf("  method:        %s\n", x$method))
  cat(sprintf("  n sequences:   %d\n", x$n_sequences))
  cat(sprintf("  split-half r:  %.3f  (sd = %.3f)\n", x$mean, x$sd))
  cat(sprintf("  95%% CI:        [%.3f, %.3f]\n", x$ci_low, x$ci_high))
  n_finite <- sum(is.finite(x$correlations))
  if (n_finite < x$R) {
    cat(sprintf("  note:          %d / %d replicates returned NA\n",
                x$R - n_finite, x$R))
  }
  invisible(x)
}
