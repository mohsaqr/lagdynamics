# Sequence-level (default) and event-level bootstrap for LSA fits.
# Implements §13 of inst/REFERENCES.md following Efron (1979) and
# Politis & Romano (1994).
#
# The recipe pattern: every bootstrap resample is refit using the
# config snapshot in fit$params (engine, lag, alternative, alpha,
# structural_zeros). No config drift possible between the original
# fit and its bootstrap.

#' Bootstrap Confidence Intervals for an LSA Fit
#'
#' Non-parametric bootstrap for any LSA fit produced by [lsa()].
#' Resamples the underlying sequence data (whole sequences when more
#' than one is available; geometric-block stationary bootstrap on
#' events otherwise), refits the engine on each resample using the
#' immutable recipe stored in `fit$params`, and aggregates per-edge
#' statistics into a tidy data frame.
#'
#' @param fit An `lsa` object returned by [lsa()].
#' @param R Integer. Number of bootstrap replicates. Default `1000`.
#' @param level Character. Resampling unit: `"sequence"` (resample
#'   whole sequences with replacement, default when fit has more than
#'   one sequence), `"event"` (stationary block bootstrap on the
#'   event stream, used automatically for single-sequence input), or
#'   `"auto"` (default; pick based on `fit$data$n_sequences`).
#' @param block_length For event-level bootstrap, mean geometric
#'   block length. Default `NULL` -> `ceiling(sqrt(T))`.
#' @param level_alpha Numeric. Confidence level for percentile
#'   intervals. Default `0.95`.
#' @param indices Optional `R x S` integer matrix where row `b` lists
#'   the sequence (or event-block start) indices to use on resample
#'   `b`. When supplied, replaces internal RNG and enables
#'   bit-identical reproducibility across runs. See Details.
#' @param parallel Logical. Use multi-core resampling. Default
#'   `FALSE`. Requires base R only (`parallel` package).
#' @param n_cores Integer. Worker count when `parallel = TRUE`.
#'   Default `NULL` -> `parallel::detectCores() - 1`.
#' @param verbose Logical. Print progress every 100 replicates.
#'   Default `FALSE`.
#' @param ... Reserved for future use.
#'
#' @details
#' **Sequence-level resampling (default for multi-sequence input).**
#' Each resample draws `S` sequence indices with replacement from
#' `seq_len(S)` and rebuilds the multi-sequence input as the
#' corresponding list of event vectors. Preserves within-sequence
#' structure.
#'
#' **Event-level resampling (single-sequence input).** Implements
#' the stationary block bootstrap of Politis & Romano (1994). Block
#' length is geometric with mean `block_length`; resampled blocks
#' wrap around the event stream and are concatenated until total
#' length equals the original `T`.
#'
#' **Reproducibility hook.** Supply `indices` as an `R x S` integer
#' matrix (sequence-level) or `R x n_blocks` matrix of start
#' positions (event-level) to deterministically replay the bootstrap
#' across sessions, processes, or languages.
#'
#' @return An object of class `c("lsa_bootstrap", "list")` with:
#' \describe{
#'   \item{edges}{Tidy per-edge data frame with observed + bootstrap
#'     `mean`, `se`, `ci_low`, `ci_high`, `p_boot`, and `stable` for
#'     `count`, `adj_res`, `prob`, and `yules_q`.}
#'   \item{boot_obs}{`R x K^2` numeric matrix: cell-wise observed
#'     count from each replicate (flattened in `as.vector(obs)`
#'     order).}
#'   \item{boot_adj_res}{`R x K^2` matrix of adjusted residuals.}
#'   \item{R, level, level_alpha, indices_used}{Recipe metadata.}
#'   \item{fit}{Reference to the original fit (for $params / labels).}
#' }
#'
#' @examples
#' \donttest{
#' fit <- lsa(engagement, engine = "classical")
#' bs <- bootstrap_lsa(fit, R = 200)
#' head(bs$edges)
#' }
#'
#' @references
#' Efron, B. (1979). Bootstrap methods: another look at the
#' jackknife. \emph{Annals of Statistics}, 7(1), 1-26.
#'
#' Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association}, 89(428),
#' 1303-1313.
#'
#' @seealso [permute_lsa()], [stability_lsa()]
#'
#' @export
bootstrap_lsa <- function(fit,
                          R = 1000L,
                          level = c("auto", "sequence", "event"),
                          block_length = NULL,
                          level_alpha = 0.95,
                          indices = NULL,
                          parallel = FALSE,
                          n_cores = NULL,
                          verbose = FALSE,
                          ...) {
  stopifnot(inherits(fit, "lsa"))
  stopifnot(is.numeric(R), length(R) == 1L, R >= 1L)
  stopifnot(is.numeric(level_alpha), level_alpha > 0, level_alpha < 1)
  level <- match.arg(level)
  R <- as.integer(R)

  recipe <- fit$params
  d <- fit$data
  if (identical(d$source, "transitions")) {
    stop("bootstrap_lsa() requires event-level input. The fit was ",
         "built from a pre-computed transition matrix and has no ",
         "underlying sequences to resample.", call. = FALSE)
  }

  # Decide resampling level.
  level_used <- if (level == "auto") {
    if (d$n_sequences >= 2L) "sequence" else "event"
  } else level
  if (level_used == "sequence" && d$n_sequences < 2L) {
    warning("Only one sequence available; falling back to event-level ",
            "stationary block bootstrap.", call. = FALSE)
    level_used <- "event"
  }

  # Pre-extract per-sequence event vectors as integer codes (faster
  # downstream than rebuilding lsa_data per resample).
  per_seq <- split(d$events, d$seq_id)

  # Pre-compute resample indices.
  if (!is.null(indices)) {
    stopifnot(is.matrix(indices), nrow(indices) >= R)
    indices <- indices[seq_len(R), , drop = FALSE]
    storage.mode(indices) <- "integer"
  } else if (level_used == "sequence") {
    S <- d$n_sequences
    indices <- matrix(sample.int(S, R * S, replace = TRUE), R, S)
  } else {
    bl <- if (is.null(block_length)) ceiling(sqrt(d$n_events)) else
            as.integer(block_length)
    indices <- .stationary_indices(R = R, T = d$n_events,
                                   mean_block = bl)
  }

  # Worker closure: produce a fitted lsa for one replicate.
  worker <- function(b) {
    if (level_used == "sequence") {
      idx <- indices[b, ]
      new_seqs <- per_seq[idx]
      new_events <- unlist(new_seqs, use.names = FALSE)
      new_seq_id <- rep.int(seq_along(idx),
                            times = vapply(new_seqs, length, integer(1)))
    } else {
      # Event-level stationary block bootstrap
      idx <- indices[b, ]
      new_events <- d$events[idx]
      new_seq_id <- rep.int(1L, length(new_events))
    }
    .refit_from_events(events = new_events, seq_id = new_seq_id,
                       labels = d$labels, recipe = recipe)
  }

  if (verbose) message("Running ", R, " bootstrap replicates ...")
  results <- .run_parallel(worker, R = R, parallel = parallel,
                           n_cores = n_cores, verbose = verbose)

  # Aggregate per-cell.
  K <- d$n_states
  cell_n <- K * K
  boot_obs     <- matrix(NA_real_, R, cell_n)
  boot_adj_res <- matrix(NA_real_, R, cell_n)
  boot_prob    <- matrix(NA_real_, R, cell_n)
  boot_yulesq  <- matrix(NA_real_, R, cell_n)
  for (b in seq_len(R)) {
    fb <- results[[b]]
    boot_obs[b, ]     <- as.vector(fb$obs)
    boot_adj_res[b, ] <- as.vector(fb$adj_res)
    boot_prob[b, ]    <- as.vector(fb$prob)
    boot_yulesq[b, ]  <- as.vector(fb$yules_q)
  }

  edges <- .summarize_bootstrap(
    fit = fit,
    boot_obs = boot_obs, boot_adj_res = boot_adj_res,
    boot_prob = boot_prob, boot_yulesq = boot_yulesq,
    level_alpha = level_alpha
  )

  structure(
    list(
      edges         = edges,
      boot_obs      = boot_obs,
      boot_adj_res  = boot_adj_res,
      boot_prob     = boot_prob,
      boot_yulesq   = boot_yulesq,
      R             = R,
      level         = level_used,
      level_alpha   = level_alpha,
      indices_used  = indices,
      fit           = fit
    ),
    class = c("lsa_bootstrap", "list")
  )
}

# --- helpers ----------------------------------------------------------

# Pre-compute indices for stationary block bootstrap.
# Each row is a length-T integer vector of event positions to sample.
# Block lengths are geometric with mean `mean_block`; blocks wrap.
.stationary_indices <- function(R, T, mean_block) {
  out <- matrix(NA_integer_, R, T)
  p <- 1 / mean_block
  for (b in seq_len(R)) {
    pos <- integer(T)
    i <- 1L
    while (i <= T) {
      start <- sample.int(T, 1L)
      blen  <- 1L + stats::rgeom(1L, p)
      end <- min(i + blen - 1L, T)
      n <- end - i + 1L
      # Wrap blocks past T
      idx <- ((start + seq_len(n) - 2L) %% T) + 1L
      pos[i:end] <- idx
      i <- end + 1L
    }
    out[b, ] <- pos
  }
  out
}

# Run a worker function R times, optionally in parallel. Falls back to
# sequential on Windows or when parallel = FALSE.
.run_parallel <- function(worker, R, parallel = FALSE, n_cores = NULL,
                          verbose = FALSE) {
  if (!parallel) {
    out <- vector("list", R)
    for (b in seq_len(R)) {
      out[[b]] <- worker(b)
      if (verbose && b %% 100L == 0L) message("  replicate ", b, "/", R)
    }
    return(out)
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("parallel = TRUE requires the base 'parallel' package.",
         call. = FALSE)
  }
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
  if (.Platform$OS.type != "windows") {
    return(parallel::mclapply(seq_len(R), worker, mc.cores = n_cores))
  }
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::parLapply(cl, seq_len(R), worker)
}

# Refit the engine given new events/seq_id, preserving the recipe.
# Avoids rebuilding lsa_data via the public coercer (much faster).
.refit_from_events <- function(events, seq_id, labels, recipe) {
  K <- length(labels)
  d <- structure(
    list(
      events = events, seq_id = seq_id, labels = labels,
      n_states = K, n_sequences = max(seq_id),
      n_events = length(events),
      transitions_per_seq = NA_integer_,  # not needed downstream
      source = "events", obs_input = NULL
    ),
    class = c("lsa_data", "list")
  )
  tx <- lsa_transitions(d, lag = recipe$lag)
  event_totals_col <- tabulate(events, nbins = K)
  attr(tx, "event_totals_col") <- event_totals_col
  entry <- get_lsa_engine(recipe$engine)
  do.call(entry$fn, c(
    list(transitions = tx,
         structural_zeros = recipe$structural_zeros,
         alternative = recipe$alternative,
         n_events = length(events)),
    recipe$params %||% list()
  ))
}

# Null-coalesce helper.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Summarize per-cell bootstrap distributions into a tidy edges frame.
.summarize_bootstrap <- function(fit, boot_obs, boot_adj_res,
                                 boot_prob, boot_yulesq,
                                 level_alpha) {
  labels <- rownames(fit$obs)
  K <- length(labels)
  grid <- expand.grid(from = labels, to = labels,
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  ci_lo <- (1 - level_alpha) / 2
  ci_hi <- 1 - ci_lo

  summarize_one <- function(boot_mat, obs_vec, label) {
    means <- colMeans(boot_mat, na.rm = TRUE)
    ses   <- apply(boot_mat, 2L, function(x) stats::sd(x, na.rm = TRUE))
    lows  <- apply(boot_mat, 2L,
                   function(x) stats::quantile(x, ci_lo, na.rm = TRUE,
                                                names = FALSE))
    highs <- apply(boot_mat, 2L,
                   function(x) stats::quantile(x, ci_hi, na.rm = TRUE,
                                                names = FALSE))
    list(observed = obs_vec, mean = means, se = ses,
         ci_low = lows, ci_high = highs)
  }

  obs_count   <- as.vector(fit$obs)
  obs_adj     <- as.vector(fit$adj_res)
  obs_prob    <- as.vector(fit$prob)
  obs_yulesq  <- as.vector(fit$yules_q)

  s_count  <- summarize_one(boot_obs,     obs_count,  "count")
  s_adj    <- summarize_one(boot_adj_res, obs_adj,    "adj_res")
  s_prob   <- summarize_one(boot_prob,    obs_prob,   "prob")
  s_yulesq <- summarize_one(boot_yulesq,  obs_yulesq, "yules_q")

  # Two-sided bootstrap p-value: 2 * min(P(stat <= 0), P(stat >= 0))
  p_boot_count   <- 2 * pmin(colMeans(boot_obs     <= 0, na.rm = TRUE),
                              colMeans(boot_obs     >= 0, na.rm = TRUE))
  p_boot_adj     <- 2 * pmin(colMeans(boot_adj_res <= 0, na.rm = TRUE),
                              colMeans(boot_adj_res >= 0, na.rm = TRUE))
  stable_adj <- sign(s_adj$ci_low) == sign(s_adj$ci_high) &
                is.finite(s_adj$ci_low) & is.finite(s_adj$ci_high) &
                s_adj$ci_low != 0 & s_adj$ci_high != 0

  data.frame(
    from         = grid$from,
    to           = grid$to,
    observed     = obs_count,
    count_mean   = s_count$mean,
    count_se     = s_count$se,
    count_ci_low = s_count$ci_low,
    count_ci_high = s_count$ci_high,
    adj_res_observed = obs_adj,
    adj_res_mean = s_adj$mean,
    adj_res_se   = s_adj$se,
    adj_res_ci_low  = s_adj$ci_low,
    adj_res_ci_high = s_adj$ci_high,
    adj_res_p_boot  = p_boot_adj,
    adj_res_stable  = stable_adj,
    prob_observed   = obs_prob,
    prob_mean       = s_prob$mean,
    prob_ci_low     = s_prob$ci_low,
    prob_ci_high    = s_prob$ci_high,
    yules_q_observed = obs_yulesq,
    yules_q_mean    = s_yulesq$mean,
    yules_q_ci_low  = s_yulesq$ci_low,
    yules_q_ci_high = s_yulesq$ci_high,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
