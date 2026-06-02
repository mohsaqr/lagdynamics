# Case-drop stability for an LSA fit. Adapted from the Nestimate /
# bootnet "centrality stability" pattern: subsample sequences (or
# events for single-sequence input), refit, and track which edges
# remain significant. An edge is "stable" if it stays significant
# under most subsamples.

#' Case-Drop Stability for an LSA Fit
#'
#' Resamples the underlying sequence data **without** replacement at
#' a specified retention proportion, refits the engine on each
#' subsample, and records which edges remain significant at the
#' recipe's alpha threshold. Returns a per-edge "stability"
#' proportion: the fraction of subsamples in which the edge was
#' significant. Edges with stability >= `min_stable` (default
#' `0.95`) are flagged as robust.
#'
#' @param fit An `lsa` object returned by [lsa()].
#' @param R Integer. Number of subsamples. Default `500`.
#' @param proportion Numeric in (0, 1). Fraction of cases retained
#'   per subsample. Default `0.8`.
#' @param min_stable Numeric in (0, 1). Stability threshold for
#'   the `stable` flag in the output edge frame. Default `0.95`.
#' @param parallel Logical. Use multi-core resampling. Default
#'   `FALSE`.
#' @param n_cores Integer. Worker count when `parallel = TRUE`.
#' @param verbose Logical. Print progress every 100 replicates.
#' @param ... Reserved.
#'
#' @return An object of class `c("lsa_stability", "list")` with:
#' \describe{
#'   \item{edges}{Tidy per-edge data frame with `observed_sig`
#'     (whether the cell was significant in the original fit),
#'     `stability` (fraction of subsamples in which the cell was
#'     significant), and `stable` (`stability >= min_stable`).}
#'   \item{stability_matrix}{`R x K^2` 0/1 matrix recording per-cell
#'     significance across replicates.}
#'   \item{R, proportion, min_stable}{Recipe metadata.}
#'   \item{fit}{Reference to the original fit.}
#' }
#'
#' @examples
#' \donttest{
#' fit <- lsa(engagement, engine = "classical")
#' st <- stability_lsa(fit, R = 100)
#' head(st$edges[order(-st$edges$stability), ])
#' }
#'
#' @seealso [bootstrap_lsa()], [permute_lsa()]
#'
#' @export
stability_lsa <- function(fit,
                          R = 500L,
                          proportion = 0.8,
                          min_stable = 0.95,
                          parallel = FALSE,
                          n_cores = NULL,
                          verbose = FALSE,
                          ...) {
  stopifnot(inherits(fit, "lsa"))
  stopifnot(is.numeric(R), R >= 1L)
  stopifnot(is.numeric(proportion), proportion > 0, proportion < 1)
  stopifnot(is.numeric(min_stable), min_stable > 0, min_stable <= 1)
  R <- as.integer(R)

  recipe <- fit$params
  d <- fit$data
  if (identical(d$source, "transitions")) {
    stop("stability_lsa() requires event-level input.",
         call. = FALSE)
  }

  per_seq <- split(d$events, d$seq_id)
  S <- length(per_seq)
  n_keep <- max(1L, round(S * proportion))

  worker <- function(b) {
    if (S >= 2L) {
      idx <- sample.int(S, n_keep, replace = FALSE)
      new_seqs <- per_seq[idx]
      new_events <- unlist(new_seqs, use.names = FALSE)
      new_seq_id <- rep.int(seq_along(idx),
                            times = vapply(new_seqs, length, integer(1)))
    } else {
      T <- d$n_events
      m <- max(2L, round(T * proportion))
      start <- sample.int(T - m + 1L, 1L)
      new_events <- d$events[start:(start + m - 1L)]
      new_seq_id <- rep.int(1L, length(new_events))
    }
    fb <- .refit_from_events(events = new_events,
                              seq_id = new_seq_id,
                              labels = d$labels, recipe = recipe)
    as.vector(fb$p < recipe$alpha)
  }

  if (verbose) message("Running ", R, " stability replicates ...")
  results <- .run_parallel(worker, R = R, parallel = parallel,
                           n_cores = n_cores, verbose = verbose)

  K <- d$n_states
  stab_mat <- do.call(rbind, results)        # R x K^2 logical
  stab_mat[is.na(stab_mat)] <- FALSE
  storage.mode(stab_mat) <- "logical"
  stability <- colMeans(stab_mat)

  labels <- d$labels
  grid <- expand.grid(from = labels, to = labels,
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  edges <- data.frame(
    from         = grid$from,
    to           = grid$to,
    observed_sig = {
      p_vec <- as.vector(fit$p)
      is.finite(p_vec) & p_vec < recipe$alpha
    },
    stability    = stability,
    stable       = stability >= min_stable,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  structure(
    list(
      edges            = edges,
      stability_matrix = stab_mat,
      R                = R,
      proportion       = proportion,
      min_stable       = min_stable,
      fit              = fit
    ),
    class = c("lsa_stability", "list")
  )
}
