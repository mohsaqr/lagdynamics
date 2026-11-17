# Case-drop stability for an LSA fit. Follows the "centrality
# stability" pattern: subsample sequences (or events for
# single-sequence input), refit, and track which edges remain
# significant. An edge is "stable" if it stays significant under most
# subsamples.

#' Assess the case-drop stability of significant transitions
#'
#' @description
#' `stability_lsa()` assesses the case-drop stability of the significant
#' transitions of a fit from [lsa()]: whether a transition's significance
#' survives when part of the data is left out, since a transition found
#' significant in the full data may owe its significance to a few sequences. It
#' follows the case-dropping approach that Epskamp, Borsboom and Fried (2018)
#' describe for the stability of psychological networks: subsets of the cases
#' are drawn without replacement, the model is refitted on each, and the share
#' of subsets in which a transition is significant measures how much its
#' significance depends on the particular cases sampled.
#'
#' `R` subsets are drawn, each keeping a share `proportion` of the sequences
#' (or, for a single sequence, a run of consecutive events of that share of its
#' length); the model is refitted with the same engine and settings, and whether
#' each transition is significant at the fit's `alpha` is recorded. A
#' transition's stability is the share of subsets in which it was significant;
#' it is flagged stable when that share reaches `min_stable`. Replicates in
#' which a subset contained no transitions are reported as failed and count as
#' not significant.
#'
#' Stability concerns the significance decision, in either direction,
#' not the size of the effect. For intervals on the statistics use
#' [bootstrap_lsa()]; for a p-value against no sequential structure use
#' [permute_lsa()]; for the agreement of the whole network between halves
#' of the data use [reliability_lsa()].
#'
#' @param fit An `lsa` fit from [lsa()], built from sequences (not from a
#'   transition matrix).
#' @param R Number of subsets. Default `500`.
#' @param proportion Share of the sequences kept in each subset, between
#'   0 and 1. Default `0.8`.
#' @param min_stable Share of subsets in which a transition must be
#'   significant to be flagged stable, greater than 0 and at most 1.
#'   Default `0.95`.
#' @param parallel Run the subsets on several cores. Default `FALSE`.
#' @param n_cores Number of cores when `parallel = TRUE`. Default `NULL`,
#'   all cores but one.
#' @param verbose Report progress every 100 subsets. Default `FALSE`.
#' @param ... Not used.
#'
#' @return An object of class `c("lsa_stability", "list")`. Read it with
#'   [transitions()], which returns a data frame with one row per
#'   transition (`from -> to` pair) and the columns `from`, `to`,
#'   `observed_sig` (significant in the full fit), `stability` (share of
#'   subsets in which it was significant) and `stable`
#'   (`stability >= min_stable`). `transitions(x, significant = TRUE)`
#'   keeps the stable transitions; `sort = "strength"` orders by
#'   stability. `print()` shows the settings, the number of failed
#'   replicates and the number of stable transitions (see
#'   [lsa-inference-methods]). The object also stores the per-subset
#'   significance decisions, the settings and the original fit.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and one of:
#' * `lagdynamics_needs_events`: the fit was built from a transition
#'   matrix, so there are no sequences to drop.
#'
#' @examples
#' fit <- lsa(engagement)
#' st <- stability_lsa(fit, R = 100)
#' st
#' transitions(st, significant = TRUE)
#' \donttest{
#' # A full analysis uses more subsets
#' stability_lsa(fit, R = 500, proportion = 0.5) |>
#'   transitions(sort = "strength")
#' }
#'
#' @references
#' Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#' psychological networks and their accuracy: A tutorial paper.
#' \emph{Behavior Research Methods, 50}(1), 195-212.
#'
#' @seealso [transitions()] to read the result; [bootstrap_lsa()],
#'   [permute_lsa()] and [reliability_lsa()] for other checks of a fit.
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
  stopifnot("`fit` must be an `lsa` fit from lsa()" = inherits(fit, "lsa"))
  stopifnot(
    "`R` must be a single whole number >= 1" =
      is.numeric(R) && length(R) == 1L && is.finite(R) && R >= 1L &&
      R == floor(R)
  )
  stopifnot(
    "`proportion` must be a single number in (0, 1)" =
      is.numeric(proportion) && length(proportion) == 1L &&
      is.finite(proportion) && proportion > 0 && proportion < 1,
    "`min_stable` must be a single number in (0, 1]" =
      is.numeric(min_stable) && length(min_stable) == 1L &&
      is.finite(min_stable) && min_stable > 0 && min_stable <= 1
  )
  R <- as.integer(R)

  recipe <- fit$params
  d <- fit$data
  if (identical(d$source, "transitions")) {
    .lsa_stop(class = "needs_events", "stability_lsa() requires event-level input.",
         call. = FALSE)
  }

  per_seq <- split(d$events, d$seq_id)
  S <- length(per_seq)
  n_keep <- max(1L, round(S * proportion))
  K <- d$n_states

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
    # A subsample of only singleton sequences contributes zero
    # transitions, and the engine signals `lagdynamics_no_transitions`.
    # That is a property of the random case-drop, not of the fit, so the
    # replicate yields an all-NA significance row (counted in `n_failed`)
    # rather than aborting the run. Any other error propagates.
    fb <- tryCatch(
      .refit_from_events(events = new_events,
                         seq_id = new_seq_id,
                         labels = d$labels, recipe = recipe),
      lagdynamics_no_transitions = function(e) NULL
    )
    if (is.null(fb)) return(rep(NA, K * K))
    as.vector(fb$p < recipe$alpha)
  }

  if (verbose) message("Running ", R, " stability replicates ...")
  results <- .run_parallel(worker, R = R, parallel = parallel,
                           n_cores = n_cores, verbose = verbose)

  stab_mat <- do.call(rbind, results)        # R x K^2 logical
  # A failed replicate is an all-NA row. It is counted, and (as before)
  # contributes "not significant" to every cell's stability.
  n_failed <- sum(rowSums(!is.na(stab_mat)) == 0L)
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
      n_failed         = n_failed,
      fit              = fit
    ),
    class = c("lsa_stability", "list")
  )
}
