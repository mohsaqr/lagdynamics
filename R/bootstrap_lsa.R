# Sequence-level (default) and event-level bootstrap for LSA fits.
# Implements §13 of inst/REFERENCES.md following Efron (1979) and
# Politis & Romano (1994).
#
# The recipe pattern: every bootstrap resample is refit using the
# config snapshot in fit$params (engine, lag, alternative, alpha,
# structural_zeros). No config drift possible between the original
# fit and its bootstrap.

#' Estimate bootstrap confidence intervals for transitions
#'
#' @description
#' `bootstrap_lsa()` estimates bootstrap confidence intervals for the
#' transitions of a fit from [lsa()]: it resamples the sequences behind the fit
#' and refits it with the same engine, lag, alternative, significance level and
#' structural zeros. The non-parametric bootstrap (Efron, 1979) estimates how
#' much the count, transition probability, adjusted residual and Yule's Q of
#' each transition, all estimated from one sample of sequences, would vary over
#' repeated samples from the same population: it draws many resamples of the
#' data, refits the model on each, and uses the spread of the refitted values as
#' the sampling distribution. A transition whose residual interval excludes zero
#' is over- or under-represented in a way that survives resampling, not only in
#' the sample at hand.
#'
#' With two or more sequences, whole sequences are drawn with replacement, which
#' keeps the dependence between events of one sequence intact. A single sequence
#' is resampled with the stationary block bootstrap (Politis & Romano, 1994),
#' which draws blocks of consecutive events of random length. For every
#' transition it returns the bootstrap mean, standard error and percentile
#' interval of the count, the transition probability, the adjusted residual and
#' Yule's Q, a two-sided bootstrap p-value for the residual, and a flag for
#' transitions whose residual interval excludes zero.
#'
#' For an interval on transition probabilities without resampling, use
#' [certainty_lsa()]; it is exact under a Dirichlet-Multinomial model but
#' ignores dependence between the transitions of one sequence. For a
#' p-value against the null hypothesis of no sequential structure, use
#' [permute_lsa()]; for the share of subsamples in which a transition
#' stays significant, use [stability_lsa()].
#'
#' @param fit An `lsa` fit from [lsa()], built from sequences (not from a
#'   transition matrix).
#' @param R Number of bootstrap resamples. Default `1000`.
#' @param level Resampling unit: `"auto"` (default; `"sequence"` when the
#'   fit has two or more sequences, `"event"` otherwise), `"sequence"`
#'   (whole sequences with replacement) or `"event"` (stationary block
#'   bootstrap of the event stream). A single sequence always uses
#'   `"event"`, with a warning.
#' @param block_length Mean block length for `level = "event"`. Default
#'   `NULL`, the square root of the number of events, rounded up.
#' @param level_alpha Confidence level of the percentile intervals.
#'   Default `0.95`.
#' @param indices Optional matrix of resampling indices to replay, one row
#'   per resample: `R` x number of sequences (sequence indices) for
#'   `level = "sequence"`, or `R` x number of events (event positions) for
#'   `level = "event"`. Default `NULL` draws them at random.
#' @param parallel Run the resamples on several cores. Default `FALSE`.
#' @param n_cores Number of cores when `parallel = TRUE`. Default `NULL`,
#'   all cores but one.
#' @param verbose Report progress every 100 resamples. Default `FALSE`.
#' @param ... Not used.
#'
#' @return An object of class `c("lsa_bootstrap", "list")`. Read it with
#'   [transitions()], which returns a data frame with one row per
#'   transition (`from -> to` pair) and the columns:
#'   * `from`, `to`: source and target state.
#'   * `observed`: the observed count; `count_mean`, `count_se`,
#'     `count_ci_low`, `count_ci_high`: its bootstrap mean, standard
#'     error and interval.
#'   * `adj_res_observed`, `adj_res_mean`, `adj_res_se`, `adj_res_ci_low`,
#'     `adj_res_ci_high`: the same for the adjusted residual;
#'     `adj_res_p_boot`: two-sided bootstrap p-value of the residual;
#'     `adj_res_stable`: `TRUE` when the residual interval excludes zero.
#'   * `prob_observed`, `prob_mean`, `prob_ci_low`, `prob_ci_high`: the
#'     transition probability.
#'   * `yules_q_observed`, `yules_q_mean`, `yules_q_ci_low`,
#'     `yules_q_ci_high`: Yule's Q.
#'
#'   `transitions(x, significant = TRUE)` keeps the transitions with
#'   `adj_res_stable = TRUE`; `sort = "strength"` orders by the absolute
#'   observed residual. `print()` shows the settings and the number of
#'   stable transitions (see [lsa-inference-methods]); [plot_forest()]
#'   draws the intervals. The object also stores the per-resample
#'   matrices, the settings, the resampling indices used (which can be
#'   passed back as `indices`) and the original fit.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and one of:
#' * `lagdynamics_needs_events`: the fit was built from a transition
#'   matrix, so there are no sequences to resample.
#' * `lagdynamics_bad_input`: `indices` is not a matrix of whole numbers
#'   with at least `R` rows, the expected number of columns and values in
#'   range, or `block_length` is not a single number of at least 1.
#'
#' @examples
#' fit <- lsa(engagement)
#' bs <- bootstrap_lsa(fit, R = 100)
#' bs
#' transitions(bs, significant = TRUE, sort = "strength")
#' \donttest{
#' # A full analysis uses more resamples
#' bootstrap_lsa(fit, R = 1000) |>
#'   transitions(sort = "strength") |>
#'   head()
#' }
#'
#' @references
#' Efron, B. (1979). Bootstrap methods: Another look at the jackknife.
#' \emph{The Annals of Statistics, 7}(1), 1-26.
#'
#' Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association, 89}(428),
#' 1303-1313.
#'
#' @seealso [transitions()] to read the result; [plot_forest()] to draw
#'   the intervals; [certainty_lsa()] for analytic intervals;
#'   [permute_lsa()], [stability_lsa()] and [reliability_lsa()] for other
#'   checks of a fit; [compare_lsa()] to compare groups.
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
  stopifnot("`fit` must be an `lsa` fit from lsa()" = inherits(fit, "lsa"))
  stopifnot(
    "`R` must be a single whole number >= 1" =
      is.numeric(R) && length(R) == 1L && is.finite(R) && R >= 1L &&
      R == floor(R)
  )
  stopifnot(
    "`level_alpha` must be a single number in (0, 1)" =
      is.numeric(level_alpha) && length(level_alpha) == 1L &&
      is.finite(level_alpha) && level_alpha > 0 && level_alpha < 1
  )
  level <- match.arg(level)
  R <- as.integer(R)

  recipe <- fit$params
  d <- fit$data
  if (identical(d$source, "transitions")) {
    .lsa_stop(class = "needs_events", "bootstrap_lsa() requires event-level input. The fit was ",
         "built from a pre-computed transition matrix and has no ",
         "underlying sequences to resample.", call. = FALSE)
  }

  # Decide resampling level.
  level_used <- if (level == "auto") {
    if (d$n_sequences >= 2L) "sequence" else "event"
  } else level
  single_sequence_warning <- NULL
  if (d$n_sequences < 2L && level == "auto") {
    single_sequence_warning <- paste0(
      "Only one sequence available; using event-level stationary ",
      "block bootstrap. Sequence-level bootstrap needs at least ",
      "two independent sequences."
    )
  } else if (level_used == "sequence" && d$n_sequences < 2L) {
    single_sequence_warning <- paste0(
      "Only one sequence available; falling back to event-level ",
      "stationary block bootstrap."
    )
    level_used <- "event"
  }

  # Pre-extract per-sequence event vectors as integer codes (faster
  # downstream than rebuilding lsa_data per resample).
  per_seq <- split(d$events, d$seq_id)

  # Pre-compute resample indices.
  if (!is.null(indices)) {
    indices <- .validate_boot_indices(indices, R = R, level = level_used,
                                      S = d$n_sequences, T = d$n_events)
  } else if (level_used == "sequence") {
    S <- d$n_sequences
    indices <- matrix(sample.int(S, R * S, replace = TRUE), R, S)
  } else {
    bl <- .validate_block_length(block_length, T = d$n_events)
    indices <- .stationary_indices(R = R, T = d$n_events,
                                   mean_block = bl)
  }

  if (!is.null(single_sequence_warning)) {
    warning(single_sequence_warning, call. = FALSE)
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

# Validate a user-supplied replay-index matrix. For sequence-level
# resampling each row must list exactly S sequence indices in 1..S; for
# event-level resampling each row is a vector of event positions in
# 1..T (the same expanded form .stationary_indices() produces, NOT
# block-start positions). Without these checks a wrong-width matrix
# silently resamples a different number of units and an out-of-range or
# NA index reaches per_seq[idx] / d$events[idx], breaking the documented
# bit-identical replay contract.
.validate_boot_indices <- function(indices, R, level, S, T) {
  if (!is.matrix(indices) || !is.numeric(indices)) {
    .lsa_stop("`indices` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(indices) < R) {
    .lsa_stop(sprintf("`indices` has %d rows but R = %d are required.",
                 nrow(indices), R), call. = FALSE)
  }
  indices <- indices[seq_len(R), , drop = FALSE]
  if (anyNA(indices) || !all(is.finite(indices)) ||
      !isTRUE(all(indices == floor(indices)))) {
    .lsa_stop("`indices` must contain only finite whole numbers ",
         "(no NA, NaN, Inf, or fractional values).", call. = FALSE)
  }
  hi <- if (level == "sequence") S else T
  unit <- if (level == "sequence") "sequence" else "event"
  if (level == "sequence" && ncol(indices) != S) {
    .lsa_stop(sprintf(
      "`indices` has %d columns but sequence-level replay needs exactly ",
      ncol(indices)),
      sprintf("S = %d (one column per resampled sequence).", S),
      call. = FALSE)
  }
  if (level == "event" && ncol(indices) != T) {
    .lsa_stop(sprintf(
      "`indices` has %d columns but event-level replay needs exactly ",
      ncol(indices)),
      sprintf("T = %d (one column per resampled event).", T),
      call. = FALSE)
  }
  if (min(indices) < 1L || max(indices) > hi) {
    .lsa_stop(sprintf("`indices` values must be in 1..%d (%s positions).",
                 hi, unit), call. = FALSE)
  }
  storage.mode(indices) <- "integer"
  indices
}

# Validate the event-level mean block length. NULL -> ceiling(sqrt(T)).
.validate_block_length <- function(block_length, T) {
  if (is.null(block_length)) return(ceiling(sqrt(T)))
  if (!is.numeric(block_length) || length(block_length) != 1L ||
      !is.finite(block_length) || block_length < 1) {
    .lsa_stop("`block_length` must be a single finite number >= 1.",
         call. = FALSE)
  }
  as.integer(block_length)
}

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
    return(lapply(seq_len(R), function(b) {
      out <- worker(b)
      if (verbose && b %% 100L == 0L) message("  replicate ", b, "/", R)
      out
    }))
  }
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
  if (.Platform$OS.type != "windows") {
    return(parallel::mclapply(seq_len(R), worker, mc.cores = n_cores))
  }
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # Seed the PSOCK workers from the caller's RNG state so that
  # set.seed(); <resample>(parallel = TRUE) is reproducible on Windows
  # (mclapply already inherits the seed on Unix).
  parallel::clusterSetRNGStream(cl)
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

# Null-coalesce helper (base R has `%||%` only from 4.4.0).
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

  # Two-sided bootstrap p-value on adjusted residuals:
  # 2 * min(P(stat <= 0), P(stat >= 0)). Only residuals get a bootstrap
  # p-value; a p-value "around zero" is not meaningful for non-negative
  # counts or probabilities, so those columns report CIs only.
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
