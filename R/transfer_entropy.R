#' Estimate directed transfer entropy between categorical sequences
#'
#' @description
#' `transfer_entropy()` (**experimental**) estimates directed transfer entropy
#' between categorical sequences. Transfer entropy (Schreiber, 2000) measures
#' directed predictive information: how much knowing the past of a source series
#' reduces uncertainty about the future of a target series, beyond what the
#' target's own past already tells. Because it conditions on the target's
#' history, it is not inflated when a state merely tends to persist, which
#' inflates lagged association measures for processes with strong momentum.
#' Transfer entropy is unsigned: a large value means strong directed dependence,
#' which may make the target more or less likely. Read it alongside the signed
#' residuals and Yule's Q of [transitions()].
#'
#' It works in two modes. With one set of sequences (`y = NULL`), it treats each
#' state as an indicator series (in the state or not) and returns the transfer
#' entropy from every state to every other state, a directed information-flow
#' network among the states. With very few states the indicator series of the
#' other states carry nearly the same information, so several pairs can share a
#' value; the bivariate mode avoids this. With two aligned sets of sequences
#' (`x` and `y`), it returns the transfer entropy between the two series in both
#' directions, using their full alphabets.
#'
#' Transfer entropy is in bits. The surrogate test (default) shuffles the
#' source to build a null distribution, which gives a p-value and an
#' effective transfer entropy with the average surrogate value subtracted.
#' For the classical test of sequential dependence use [lsa()].
#'
#' @param x Categorical sequences: a vector (one sequence) or a matrix or
#'   data frame with one sequence per row and one time point per column,
#'   padded with `NA`, as accepted by [lsa()].
#' @param y Optional second set of sequences of the same shape as `x`, for
#'   the bivariate mode. Default `NULL` (state-flow network of `x`).
#' @param lag Prediction horizon, an integer of at least 1 (default `1`).
#' @param history Number of past steps of the target conditioned on, an
#'   integer of at least 1 (default `1`).
#' @param test `"surrogate"` (default) for the surrogate test, or
#'   `"none"`.
#' @param R Number of surrogates (default `199`).
#' @param normalize Add `te_normalised`, transfer entropy as a share of the
#'   target's remaining uncertainty given its own history, in `[0, 1]`
#'   (default `TRUE`).
#' @param seed Optional seed for the surrogate test (default `NULL`). The
#'   session's random number state is restored on exit.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`
#' when `lag`, `history` or `R` is not a positive whole number, or when `x`
#' and `y` are not aligned (different numbers of sequences or different
#' lengths within a sequence).
#'
#' @return A `data.frame` with one row per ordered pair (`from`, `to`) and
#'   the columns `te` (bits), `te_effective` (with `test = "surrogate"`),
#'   `te_normalised` (with `normalize = TRUE`), `p` (with `test =
#'   "surrogate"`) and `n` (number of transitions used), sorted by
#'   decreasing `te`. In the state-flow mode there is one row per ordered
#'   pair of distinct states; in the bivariate mode two rows, named after
#'   the arguments.
#'
#' @references
#' Schreiber, T. (2000). Measuring information transfer. \emph{Physical
#' Review Letters, 85}(2), 461-464.
#'
#' @examples
#' # Information flow among the engagement states
#' transfer_entropy(engagement, test = "none")
#'
#' # With the surrogate test (an analysis would use more surrogates)
#' transfer_entropy(engagement, R = 19, seed = 1)
#'
#' # Between two aligned series
#' mood <- c("calm", "calm", "tense", "tense", "calm", "tense", "tense", "calm")
#' effort <- c("low", "low", "low", "high", "high", "low", "high", "high")
#' transfer_entropy(mood, effort, test = "none")
#'
#' @seealso [lsa()], [transitions()], [lsa_lags()]
#'
#' @export
transfer_entropy <- function(x, y = NULL, lag = 1L, history = 1L,
                             test = c("surrogate", "none"), R = 199L,
                             normalize = TRUE, seed = NULL) {
  test <- match.arg(test)
  .te_check_integer(lag, "lag")
  .te_check_integer(history, "history")
  .te_check_integer(R, "R")
  if (!is.null(seed)) {
    stopifnot("`seed` must be a single number" =
                is.numeric(seed) && length(seed) == 1L)
    old_seed <- .lsa_get_seed()
    on.exit(.lsa_restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }
  lag <- as.integer(lag); history <- as.integer(history); R <- as.integer(R)

  x_list <- .te_as_seq_list(x)
  if (is.null(y)) {
    states <- .te_states(x_list)
    out <- .te_state_network(x_list, states, lag, history, test, R)
  } else {
    y_list <- .te_as_seq_list(y)
    labels <- c(.te_label(substitute(x), "x"), .te_label(substitute(y), "y"))
    out <- .te_bivariate(x_list, y_list, labels, lag, history, test, R)
  }
  if (!normalize) out$te_normalised <- NULL
  if (test == "none") { out$te_effective <- NULL; out$p <- NULL }
  out <- out[order(-out$te), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# --- internals --------------------------------------------------------------

# Split input into a list of per-unit character sequences.
.te_as_seq_list <- function(z) {
  if (is.null(dim(z))) {
    list(as.character(z))
  } else {
    m <- as.matrix(z)
    lapply(seq_len(nrow(m)), function(i) as.character(m[i, ]))
  }
}

.te_states <- function(x_list) {
  s <- sort(unique(stats::na.omit(unlist(x_list))))
  s[s != ""]
}

.te_label <- function(sym, fallback) {
  lab <- tryCatch(deparse(sym), error = function(e) fallback)
  if (length(lab) != 1L || nchar(lab) == 0L || nchar(lab) > 20L) fallback else lab
}

.te_check_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != floor(x)) {
    .lsa_stop(sprintf("`%s` must be a single finite integer >= 1.", name),
         call. = FALSE)
  }
  invisible(TRUE)
}

.te_validate_aligned <- function(x_list, y_list) {
  if (length(x_list) != length(y_list)) {
    .lsa_stop("`x` and `y` must contain the same number of aligned sequences.",
         call. = FALSE)
  }
  lx <- lengths(x_list)
  ly <- lengths(y_list)
  if (!identical(lx, ly)) {
    .lsa_stop("`x` and `y` must have the same length within every aligned sequence.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Aligned (target future, target history, source past) triples for ONE
# sequence. The lag never crosses the sequence boundary; NA / "" are dropped.
.te_triples <- function(src, tgt, lag, history) {
  src[src == ""] <- NA; tgt[tgt == ""] <- NA
  n <- length(tgt)
  t <- seq.int(history, n - lag)
  if (length(t) < 1L) return(NULL)
  yf <- tgt[t + lag]
  xp <- src[t]
  hist_mat <- vapply(seq_len(history) - 1L,
                     function(d) tgt[t - d], character(length(t)))
  dim(hist_mat) <- c(length(t), history)
  ok <- !is.na(yf) & !is.na(xp) & rowSums(is.na(hist_mat)) == 0L
  if (!any(ok)) return(NULL)
  yh <- apply(hist_mat[ok, , drop = FALSE], 1L, paste, collapse = "|")
  data.frame(yf = yf[ok], yh = yh, xp = xp[ok], stringsAsFactors = FALSE)
}

# Plug-in joint Shannon entropy (bits).
.te_entropy <- function(...) {
  p <- prop.table(table(...))
  p <- p[p > 0]
  -sum(p * log2(p))
}

# TE and the target's leftover uncertainty from pooled triples.
.te_value <- function(d) {
  leftover <- .te_entropy(d$yf, d$yh) - .te_entropy(d$yh)            # H(yf | yh)
  cond     <- .te_entropy(d$yf, d$xp, d$yh) - .te_entropy(d$xp, d$yh) # H(yf | xp, yh)
  c(te = leftover - cond, leftover = leftover)
}

# Source-permutation surrogate: breaks the source->target coupling while
# preserving marginals, giving a p-value and bias-corrected effective TE.
.te_surrogate <- function(d, R) {
  observed   <- unname(.te_value(d)["te"])
  surrogates <- replicate(R, unname(.te_value(within(d, xp <- sample(xp)))["te"]))
  list(effective = observed - mean(surrogates),
       p = (1 + sum(surrogates >= observed)) / (R + 1))
}

# One ordered pair: pool triples across units, compute TE (+ surrogate).
.te_one <- function(src_list, tgt_list, lag, history, test, R) {
  trip <- do.call(rbind, Map(.te_triples, src_list, tgt_list,
                             MoreArgs = list(lag = lag, history = history)))
  if (is.null(trip) || nrow(trip) == 0L) {
    return(c(te = NA_real_, te_effective = NA_real_, te_normalised = NA_real_,
             p = NA_real_, n = 0))
  }
  val <- .te_value(trip)
  sur <- if (test == "surrogate") .te_surrogate(trip, R)
         else list(effective = NA_real_, p = NA_real_)
  c(te = unname(val["te"]),
    te_effective = sur$effective,
    te_normalised = if (val["leftover"] > 0) unname(val["te"] / val["leftover"]) else NA_real_,
    p = sur$p,
    n = nrow(trip))
}

# Directed state-flow network via binary occupancy channels.
.te_state_network <- function(x_list, states, lag, history, test, R) {
  grid <- expand.grid(from = states, to = states, stringsAsFactors = FALSE)
  grid <- grid[grid$from != grid$to, , drop = FALSE]
  ind <- function(v, s) ifelse(is.na(v), NA_character_,
                               as.character(as.integer(v == s)))
  one <- function(s, r) {
    .te_one(lapply(x_list, ind, s), lapply(x_list, ind, r),
            lag, history, test, R)
  }
  res <- t(mapply(one, grid$from, grid$to))
  data.frame(grid, res, row.names = NULL)
}

# Bivariate TE in both directions between two aligned series.
.te_bivariate <- function(x_list, y_list, labels, lag, history, test, R) {
  .te_validate_aligned(x_list, y_list)
  rbind(
    data.frame(from = labels[1], to = labels[2],
               t(.te_one(x_list, y_list, lag, history, test, R)),
               row.names = NULL),
    data.frame(from = labels[2], to = labels[1],
               t(.te_one(y_list, x_list, lag, history, test, R)),
               row.names = NULL)
  )
}
