#' Directed transfer entropy for categorical sequences (experimental)
#'
#' @description
#' **Experimental.** Transfer entropy (Schreiber, 2000) measures *directed*
#' predictive coupling: how much knowing the source's past reduces uncertainty
#' about the target's future, **beyond what the target's own past already
#' explains**. Because it conditions on the target's own history, transfer
#' entropy is immune to the autocorrelation confound that inflates plain lagged
#' association when a process has strong momentum.
#'
#' Unlike [lsa()]'s Yule's Q / adjusted residuals, transfer entropy is
#' **sign-blind**: a large value means "strong directed predictive structure",
#' which may be *facilitating* or *suppressing*. Read it alongside the signed
#' measures from [transitions()] to interpret direction of effect.
#'
#' Two modes:
#' * **State-flow network** (default, `y = NULL`): given categorical
#'   sequences, returns directed transfer entropy between every ordered pair of
#'   states, via a binary occupancy decomposition. Note: with very few states
#'   the non-target source channels become redundant; prefer the bivariate mode
#'   or a larger alphabet when that matters.
#' * **Bivariate** (`y` supplied): full-alphabet transfer entropy between two
#'   aligned categorical series, in both directions.
#'
#' @param x Categorical sequence data: a vector (one sequence), or a matrix /
#'   data.frame with one sequence per row and one time-step per column
#'   (`NA`-padded), exactly the shape [lsa()] consumes.
#' @param y Optional second series, same shape as `x`, for bivariate transfer
#'   entropy. When `NULL`, the directed state-flow network of `x` is returned.
#' @param lag Integer >= 1. Prediction horizon: the target's future is taken
#'   `lag` steps ahead. Default `1`.
#' @param history Integer >= 1. Order of the target's own history conditioned
#'   on (and combined into a composite symbol). Default `1`.
#' @param test `"surrogate"` (default) runs a source-permutation null to give a
#'   p-value and the bias-corrected *effective* transfer entropy; `"none"`
#'   skips it.
#' @param R Integer. Number of surrogate permutations. Default `199`.
#' @param normalize Logical. Add `te_normalised`, transfer entropy as a share
#'   of the target's leftover uncertainty `H(future | history)`, in `[0, 1]`.
#'   Default `TRUE`.
#' @param seed Optional integer seed for the surrogate test.
#'
#' @return A tidy `data.frame`, one row per ordered pair, with columns
#'   `from`, `to`, `te` (bits), `te_effective` (surrogate-debiased),
#'   `te_normalised` (0-1, if `normalize`), `p` (surrogate p-value), and
#'   `n` (pooled transitions used). Rows are ordered by descending `te`.
#'
#' @references Schreiber, T. (2000). Measuring information transfer.
#'   *Physical Review Letters*, 85(2), 461-464.
#'
#' @examples
#' # Directed information-flow network over engagement states
#' transfer_entropy(engagement, test = "none")
#'
#' # Bivariate transfer entropy between two aligned series
#' a <- c("calm", "calm", "tense", "tense", "calm", "tense", "tense", "calm")
#' b <- c("low", "low", "low", "high", "high", "low", "high", "high")
#' transfer_entropy(a, b, test = "none")
#' @export
transfer_entropy <- function(x, y = NULL, lag = 1L, history = 1L,
                             test = c("surrogate", "none"), R = 199L,
                             normalize = TRUE, seed = NULL) {
  test <- match.arg(test)
  .te_check_integer(lag, "lag")
  .te_check_integer(history, "history")
  .te_check_integer(R, "R")
  if (!is.null(seed)) set.seed(seed)
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
    stop(sprintf("`%s` must be a single finite integer >= 1.", name),
         call. = FALSE)
  }
  invisible(TRUE)
}

.te_validate_aligned <- function(x_list, y_list) {
  if (length(x_list) != length(y_list)) {
    stop("`x` and `y` must contain the same number of aligned sequences.",
         call. = FALSE)
  }
  lx <- lengths(x_list)
  ly <- lengths(y_list)
  if (!identical(lx, ly)) {
    stop("`x` and `y` must have the same length within every aligned sequence.",
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
