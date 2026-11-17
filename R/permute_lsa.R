# Permutation test for an LSA fit. Implements §14 of
# inst/REFERENCES.md following Castellan (1992) and Phipson & Smyth
# (2010).
#
# Each replicate shuffles the event vector within sequence boundaries
# (or globally if the input was a single sequence), recomputes the
# transition matrix, and tallies the proportion of permuted cells
# that exceed the observed cell's magnitude. The +1 correction
# (Phipson-Smyth 2010) prevents p = 0 artifacts.

#' Test transitions against no sequential structure by permutation
#'
#' @description
#' `permute_lsa()` tests every transition of a fit from [lsa()] against no
#' sequential structure by permutation: whether the transition is over- or
#' under-represented beyond what the same events in a random order would give,
#' keeping the frequency of every state fixed. The adjusted residual of a
#' transition is otherwise tested against a normal reference distribution, which
#' holds only approximately, for large counts and independent transitions. A
#' permutation test builds the reference distribution from the data instead: if
#' the order of events carried no information, shuffling the events would
#' produce residuals as large as the observed ones.
#'
#' The events of the fit are shuffled many times, by default within each
#' sequence so that every sequence keeps its own length and state frequencies;
#' the model is refitted with the same engine and settings, and the adjusted
#' residual of every transition is recorded. The two-sided permutation p-value
#' of a transition is the share of shuffles whose residual is at least as large
#' in absolute value as the observed one, counted so that it is never zero
#' (Phipson & Smyth, 2010). A transition is significant when its p-value is
#' below the `alpha` of the fit.
#'
#' The p-values are not corrected for the number of transitions tested.
#' For intervals on the residuals and other statistics, use
#' [bootstrap_lsa()]; for the share of subsamples in which a transition
#' stays significant, use [stability_lsa()]. To test whether groups
#' differ, use [compare_lsa()].
#'
#' @param fit An `lsa` fit from [lsa()], built from sequences (not from a
#'   transition matrix).
#' @param R Number of permutations. Default `1000`.
#' @param within_sequence Shuffle events within each sequence (`TRUE`,
#'   default) or across the whole event stream (`FALSE`). A single
#'   sequence is always shuffled as a whole, with a warning when
#'   `within_sequence = TRUE`.
#' @param shuffles Optional list of at least `R` permutations of
#'   `seq_len(n_events)` to replay instead of random shuffles. Default
#'   `NULL`.
#' @param parallel Run the permutations on several cores. Default
#'   `FALSE`.
#' @param n_cores Number of cores when `parallel = TRUE`. Default `NULL`,
#'   all cores but one.
#' @param verbose Report progress every 100 permutations. Default `FALSE`.
#' @param ... Not used.
#'
#' @return An object of class `c("lsa_permutation", "list")`. Read it with
#'   [transitions()], which returns a data frame with one row per
#'   transition (`from -> to` pair) and the columns `from`, `to`,
#'   `observed_count`, `observed_adj_res`, `p_perm` (permutation p-value;
#'   `NA` for transitions that cannot be estimated, such as structural
#'   zeros) and `significant` (`p_perm` below the fit's `alpha`).
#'   `transitions(x, significant = TRUE)` keeps the significant
#'   transitions; `sort = "strength"` orders by the absolute observed
#'   residual. `print()` shows the settings and the number of significant
#'   transitions (see [lsa-inference-methods]). The object also stores the
#'   permuted residuals, the settings and the original fit.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and one of:
#' * `lagdynamics_needs_events`: the fit was built from a transition
#'   matrix, so there are no events to shuffle.
#' * `lagdynamics_bad_input`: `shuffles` is not a list of at least `R`
#'   permutations of `seq_len(n_events)`.
#'
#' @examples
#' fit <- lsa(engagement)
#' pm <- permute_lsa(fit, R = 100)
#' pm
#' transitions(pm, significant = TRUE, sort = "strength")
#' \donttest{
#' # A full analysis uses more permutations
#' permute_lsa(fit, R = 1000) |>
#'   transitions(sort = "strength")
#' }
#'
#' @references
#' Castellan, N. J. (1992). Shuffling arrays: Appearances may be
#' deceiving. \emph{Behavior Research Methods, Instruments, & Computers,
#' 24}(1), 72-77.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never
#' be zero: Calculating exact p-values when permutations are randomly
#' drawn. \emph{Statistical Applications in Genetics and Molecular
#' Biology, 9}(1), Article 39.
#'
#' @seealso [transitions()] to read the result; [bootstrap_lsa()],
#'   [stability_lsa()] and [reliability_lsa()] for other checks of a fit;
#'   [compare_lsa()] to compare groups.
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
  stopifnot("`fit` must be an `lsa` fit from lsa()" = inherits(fit, "lsa"))
  stopifnot(
    "`R` must be a single whole number >= 1" =
      is.numeric(R) && length(R) == 1L && is.finite(R) && R >= 1L &&
      R == floor(R)
  )
  R <- as.integer(R)

  recipe <- fit$params
  d <- fit$data
  if (identical(d$source, "transitions")) {
    .lsa_stop(class = "needs_events", "permute_lsa() requires event-level input. The fit was ",
         "built from a pre-computed transition matrix.",
         call. = FALSE)
  }

  events <- d$events
  seq_id <- d$seq_id
  labels <- d$labels
  K <- d$n_states
  T <- d$n_events

  if (isTRUE(within_sequence) && d$n_sequences < 2L) {
    warning("Only one sequence available; permute_lsa() will shuffle the ",
            "whole event stream. Within-sequence permutation is identical ",
            "to global permutation for a single sequence.", call. = FALSE)
  }

  # Per-sequence index lists for within-sequence shuffling.
  seq_positions <- split(seq_len(T), seq_id)

  # Pre-validate user-supplied shuffles. Each element must be a genuine
  # permutation of seq_len(T): used directly as events[perm], a vector
  # with duplicates, gaps, out-of-range indices, or NA would silently
  # drop/duplicate events or inject NA and still return an object,
  # breaking the documented reproducible-replay contract.
  if (!is.null(shuffles)) {
    if (!is.list(shuffles) || length(shuffles) < R) {
      .lsa_stop(sprintf("`shuffles` must be a list of length >= R (%d).", R),
           call. = FALSE)
    }
    shuffles <- shuffles[seq_len(R)]
    target <- seq_len(T)
    bad <- which(!vapply(shuffles, .is_permutation, logical(1), n = T))
    if (length(bad) > 0L) {
      .lsa_stop(sprintf(
        "`shuffles[[%d]]` is not a permutation of seq_len(%d). Each ",
        bad[1L], T),
        "shuffle must contain every integer 1..n_events exactly once, ",
        "with no NA, duplicates, or out-of-range values.",
        call. = FALSE)
    }
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
  # Phipson-Smyth two-sided p-value with +1 correction. Non-estimable
  # cells (NA observed residual, e.g. structural zeros or zero-margin
  # rows) must NOT be tested: with an NA observed value every
  # comparison is NA, colSums(na.rm = TRUE) collapses to 0 exceedances,
  # and p would spuriously become 1/(R + 1) -- flagging a forbidden
  # transition as significant. Such cells get p = NA. The denominator
  # is the count of FINITE null replicates per cell (not a flat R), so
  # cells whose null is partly non-estimable are still scaled correctly.
  abs_obs <- abs(obs_adj)
  exceed <- colSums(abs(perm_adj_res) >= matrix(abs_obs, R, K * K,
                                                 byrow = TRUE),
                     na.rm = TRUE)
  n_finite <- colSums(is.finite(perm_adj_res))
  p_perm <- (1 + exceed) / (1 + n_finite)
  p_perm[!is.finite(obs_adj) | n_finite == 0L] <- NA_real_

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

# TRUE iff `p` is a permutation of seq_len(n): integer-valued, length n,
# finite, and a bijection onto 1..n. Used to reject malformed replay
# hooks before they reach events[perm].
.is_permutation <- function(p, n) {
  if (length(p) != n) return(FALSE)
  p <- suppressWarnings(as.integer(p))
  if (anyNA(p)) return(FALSE)
  all(sort(p) == seq_len(n))
}
