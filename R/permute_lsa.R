# Permutation test for an LSA fit. Implements §14 of
# inst/REFERENCES.md following Castellan (1992) and Phipson & Smyth
# (2010).
#
# Each replicate shuffles the event vector within sequence boundaries
# (or globally if the input was a single sequence), recomputes the
# transition matrix, and tallies the proportion of permuted cells
# that exceed the observed cell's magnitude. The +1 correction
# (Phipson-Smyth 2010) prevents p = 0 artifacts.

#' Permutation Test for an LSA Fit
#'
#' Empirical null-distribution p-values for every cell of an LSA
#' transition matrix. Repeatedly shuffles the input event vector
#' (within sequence boundaries) and recomputes the engine's residual
#' matrix, producing a permutation distribution for each cell. The
#' two-sided p-value is `(1 + #{ |stat_perm| >= |stat_obs| }) / (1 + R)`
#' (Phipson & Smyth 2010).
#'
#' @param fit An `lsa` object returned by [lsa()].
#' @param R Integer. Number of permutations. Default `1000`.
#' @param within_sequence Logical. When `TRUE` (default for
#'   multi-sequence input), each sequence is shuffled independently;
#'   when `FALSE`, the whole event stream is shuffled across
#'   sequence boundaries.
#' @param shuffles Optional list of length `R`, each element an
#'   integer permutation of `seq_len(n_events)`. When supplied,
#'   replaces internal RNG.
#' @param parallel Logical. Use multi-core resampling. Default
#'   `FALSE`.
#' @param n_cores Integer. Worker count when `parallel = TRUE`.
#' @param verbose Logical. Print progress every 100 replicates.
#' @param ... Reserved.
#'
#' **NA handling.** The exceedance count that drives `p_perm` is
#' computed with `na.rm = TRUE`, so replicates that produced `NA`
#' for a cell (structural-zero cells, zero-margin cells in the
#' permuted table) are excluded from that cell's tally rather than
#' counted as either an exceedance or a non-exceedance.
#'
#' @return An object of class `c("lsa_permutation", "list")` with:
#' \describe{
#'   \item{edges}{Tidy per-edge data frame with observed counts and
#'     residuals, the empirical permutation p-value `p_perm`, and a
#'     `significant` flag at the recipe's alpha threshold.}
#'   \item{perm_adj_res}{`R x K^2` numeric matrix of permuted
#'     residuals (cells in `as.vector(adj_res)` order).}
#'   \item{R, within_sequence}{Recipe metadata.}
#'   \item{fit}{Reference to the original fit.}
#' }
#'
#' @examples
#' \donttest{
#' fit <- lsa(engagement, engine = "classical")
#' pm <- permute_lsa(fit, R = 200)
#' head(pm$edges)
#' }
#'
#' @references
#' Castellan, N. J. (1992). Shuffling arrays: appearances may be
#' deceiving. \emph{Behavior Research Methods, Instruments, &
#' Computers}, 24(1), 72-77.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should
#' never be zero: calculating exact p-values when permutations are
#' randomly drawn. \emph{Statistical Applications in Genetics and
#' Molecular Biology}, 9(1), Article 39.
#'
#' @seealso [bootstrap_lsa()], [stability_lsa()]
#'
#' @export
permute_lsa <- function(fit,
                        R = 1000L,
                        within_sequence = TRUE,
                        shuffles = NULL,
                        parallel = FALSE,
                        n_cores = NULL,
                        verbose = FALSE,
                        ...) {
  stopifnot(inherits(fit, "lsa"))
  stopifnot(is.numeric(R), length(R) == 1L, R >= 1L)
  R <- as.integer(R)

  recipe <- fit$params
  d <- fit$data
  if (identical(d$source, "transitions")) {
    stop("permute_lsa() requires event-level input. The fit was ",
         "built from a pre-computed transition matrix.",
         call. = FALSE)
  }

  events <- d$events
  seq_id <- d$seq_id
  labels <- d$labels
  K <- d$n_states
  T <- d$n_events

  # Per-sequence index lists for within-sequence shuffling.
  seq_positions <- split(seq_len(T), seq_id)

  # Pre-validate user-supplied shuffles.
  if (!is.null(shuffles)) {
    stopifnot(is.list(shuffles), length(shuffles) >= R)
    shuffles <- shuffles[seq_len(R)]
  }

  worker <- function(b) {
    if (!is.null(shuffles)) {
      perm <- as.integer(shuffles[[b]])
      new_events <- events[perm]
    } else if (within_sequence && d$n_sequences > 1L) {
      new_events <- .shuffle_within_sequences(events, seq_positions)
    } else {
      new_events <- events[sample.int(T)]
    }
    # Recompute the transition count matrix and adjusted residuals
    # via the registered engine.
    .refit_from_events(events = new_events, seq_id = seq_id,
                       labels = labels, recipe = recipe)
  }

  if (verbose) message("Running ", R, " permutations ...")
  results <- .run_parallel(worker, R = R, parallel = parallel,
                           n_cores = n_cores, verbose = verbose)

  # Stack residuals across replicates.
  perm_adj_res <- matrix(NA_real_, R, K * K)
  for (b in seq_len(R)) {
    perm_adj_res[b, ] <- as.vector(results[[b]]$adj_res)
  }

  obs_adj <- as.vector(fit$adj_res)
  # Phipson-Smyth two-sided p-value with +1 correction.
  abs_obs <- abs(obs_adj)
  exceed <- colSums(abs(perm_adj_res) >= matrix(abs_obs, R, K * K,
                                                 byrow = TRUE),
                     na.rm = TRUE)
  p_perm <- (1 + exceed) / (1 + R)

  alpha <- recipe$alpha
  grid <- expand.grid(from = labels, to = labels,
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  edges <- data.frame(
    from           = grid$from,
    to             = grid$to,
    observed_count = as.vector(fit$obs),
    observed_adj_res = obs_adj,
    p_perm         = p_perm,
    significant    = is.finite(p_perm) & p_perm < alpha,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  structure(
    list(
      edges          = edges,
      perm_adj_res   = perm_adj_res,
      R              = R,
      within_sequence = within_sequence,
      fit            = fit
    ),
    class = c("lsa_permutation", "list")
  )
}

# Independently permute the events inside each sequence's positions.
# seq_positions is a list of integer index vectors into `events`, one
# entry per sequence. Reads and writes stay inside each sp, so
# disjoint-alphabet sequences cannot bleed into one another.
#
# Indexing note: the source positions are sp[sample.int(length(sp))],
# NOT sample.int(length(sp)). The latter (pre-fix) form yields indices
# 1..length(sp), which for sequences after the first reads from the
# wrong region of events and copies one sequence's content into
# another. sample.int() is used (not sample(sp)) to avoid R's
# sample(x) quirk where a length-1 integer x is treated as 1:x.
.shuffle_within_sequences <- function(events, seq_positions) {
  for (sp in seq_positions) {
    if (length(sp) >= 2L) {
      events[sp] <- events[sp[sample.int(length(sp))]]
    }
  }
  events
}
