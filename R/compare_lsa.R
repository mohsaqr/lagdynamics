# Group comparison for LSA fits. Confirmatory counterpart to
# permute_lsa(): instead of testing each edge against a no-sequential-
# structure null, it tests whether groups' transition structures differ,
# by permuting the group label of whole sequences.
#
# The unit of permutation is the sequence (not the event), so every
# within-sequence dependency is preserved and only the group<->structure
# association is tested. The per-edge statistic is the difference in the
# chosen measure (adjusted residuals by default); the omnibus statistic
# is the Frobenius norm of the difference matrix (sum of squared edge
# differences), an NCT-style test of overall network invariance.
#
# Two groups -> one lsa_comparison. More than two -> all C(k, 2) pairwise
# comparisons (lsa_comparison_pairwise), with the multiple-comparison
# correction applied once across the whole family of per-edge p-values.

#' Compare transition structures between groups with a permutation test
#'
#' @description
#' `compare_lsa()` tests whether groups of sequences (high and low achievers,
#' courses, conditions) differ in their transition structure, transition by
#' transition and as a whole, with a permutation test. Groups can differ in
#' which transitions they over- or under-produce. The permutation test checks
#' whether an observed difference between groups is larger than would arise if
#' group membership were unrelated to the sequences: the group labels of whole
#' sequences are shuffled many times, the models are refitted, and the observed
#' difference is set against the differences produced by the shuffles. Because
#' whole sequences are permuted, the dependence between events of one sequence
#' is kept and only the link between group and transition structure is tested.
#'
#' For each transition, `compare_lsa()` compares a per-group measure, by default
#' the log odds ratio of the transition against all other transitions, which
#' does not grow with the number of transitions observed, so groups of different
#' sizes can be compared. Each transition gets the group difference, a
#' permutation p-value (Phipson & Smyth, 2010) and an adjusted p-value;
#' transitions too rare to be compared are left untested. An omnibus test of
#' overall difference, the sum of squared differences over the tested
#' transitions, is computed from the same permutations, in the spirit of the
#' network comparison test of van Borkulo et al. (2023). With more than two
#' groups every pair is compared and the `adjust` correction is applied once
#' across the p-values of all pairs.
#'
#' The groups must share their states, engine, lag and structural zeros;
#' fitting them together with `lsa(data, group = )` ensures this. For a
#' Bayesian comparison of transition probabilities that reports the
#' plausible size of each difference without resampling, use
#' [bayes_compare_lsa()]. A difference between groups can also arise from
#' a difference in their composition; check subgroups when such a
#' confound is plausible.
#'
#' @param x A grouped fit from `lsa(data, group = )` with two or more
#'   groups, or the `lsa` fit of the first group.
#' @param y The `lsa` fit of the second group when `x` is a single fit
#'   (the groups are then labelled `"a"` and `"b"`). Default `NULL`.
#' @param R Number of label permutations per comparison. Default `1000`.
#' @param measure Per-transition measure compared between groups:
#'   `"log_or"` (default; log odds ratio of the transition, with 0.5 added
#'   to every cell of its 2 x 2 table when one cell is empty), `"yules_q"`
#'   (also independent of sample size, but bounded at -1 and 1),
#'   `"adj_res"` (adjusted residual; grows with the number of
#'   transitions, so a message is shown when one group has more than 1.5
#'   times the transitions of another), `"prob"` (transition probability),
#'   `"count"` or `"lift"` (observed / expected count).
#' @param adjust Multiple-comparison correction, any method of
#'   [stats::p.adjust()] such as `"BH"`, `"holm"` or `"bonferroni"`.
#'   Default `"none"`.
#' @param min_count Minimum combined count of a transition in the two
#'   groups for it to be tested. Default `5`; `0` tests every transition.
#' @param parallel Run the permutations on several cores. Default `FALSE`.
#' @param n_cores Number of cores when `parallel = TRUE`. Default `NULL`,
#'   all cores but one.
#' @param verbose Report progress. Default `FALSE`.
#' @param ... Not used.
#'
#' @return With two groups, an object of class
#'   `c("lsa_comparison", "list")`; with more than two, an object of class
#'   `c("lsa_comparison_pairwise", "list")`. Read the per-transition table
#'   with [transitions()]: one row per transition (`from -> to` pair), and
#'   for more than two groups one row per transition and pair of groups,
#'   with the columns
#'   * `group_a`, `group_b` (more than two groups only): the pair compared;
#'   * `from`, `to`: source and target state;
#'   * `<measure>_a`, `<measure>_b` (for example `log_or_a`): the measure
#'     in each group;
#'   * `diff`: the difference, group a minus group b;
#'   * `p_perm`: permutation p-value (`NA` when not tested);
#'   * `p_adj`: p-value after `adjust`;
#'   * `significant`: `p_adj` below the `alpha` of the fits.
#'
#'   `transitions(x, significant = TRUE)` keeps the significant
#'   transitions; `sort = "strength"` orders by the absolute difference.
#'   The omnibus test is read with `as.data.frame(x, what = "global")`
#'   (one row per pair of groups). `print()` shows the groups, the number
#'   of significant transitions and the omnibus test (see
#'   [as.data.frame.lsa_comparison()]); [plot()][plot.lsa_comparison]
#'   draws the differences. The object also stores the permuted
#'   differences, the settings and the group fits.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and one of:
#' * `lagdynamics_needs_events`: a fit was built from a transition matrix,
#'   so there are no sequences to permute.
#' * `lagdynamics_incompatible_fits`: the groups differ in their states,
#'   engine, lag or structural zeros.
#' * `lagdynamics_bad_input`: `adjust` is not a method of
#'   [stats::p.adjust()]; `x` is neither a grouped fit nor an `lsa` fit
#'   with a second fit in `y`; `y` is given together with a grouped fit;
#'   or the grouped fit has fewer than two groups.
#'
#' @examples
#' g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'          time = "Time", group = "Achiever")
#' cmp <- compare_lsa(g, R = 100, adjust = "BH")
#' cmp
#' transitions(cmp, significant = TRUE, sort = "strength") |> head()
#' as.data.frame(cmp, what = "global")
#' \donttest{
#' # A full analysis uses more permutations
#' compare_lsa(g, R = 1000, adjust = "BH")
#'
#' # Three courses: every pair, corrected across all pairs
#' courses <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'                time = "Time", group = "Course")
#' compare_lsa(courses, R = 100, adjust = "BH") |>
#'   as.data.frame(what = "global")
#' }
#'
#' @references
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never
#' be zero: Calculating exact p-values when permutations are randomly
#' drawn. \emph{Statistical Applications in Genetics and Molecular
#' Biology, 9}(1), Article 39.
#'
#' van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J.,
#' Tio, P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2023).
#' Comparing network structures on three aspects: A permutation test.
#' \emph{Psychological Methods, 28}(6), 1273-1285.
#'
#' @seealso [transitions()] and [as.data.frame.lsa_comparison()] to read
#'   the result; [plot.lsa_comparison()] to draw it; [bayes_compare_lsa()]
#'   for a Bayesian comparison; [permute_lsa()] to test one fit against
#'   no sequential structure; [lsa()] with `group =` to fit the groups.
#'
#' @export
compare_lsa <- function(x,
                        y = NULL,
                        R = 1000L,
                        measure = c("log_or", "adj_res", "yules_q",
                                    "prob", "count", "lift"),
                        adjust = "none",
                        min_count = 5L,
                        parallel = FALSE,
                        n_cores = NULL,
                        verbose = FALSE,
                        ...) {
  measure <- match.arg(measure)
  stopifnot(
    "`min_count` must be a single number >= 0" =
      is.numeric(min_count) && length(min_count) == 1L &&
      is.finite(min_count) && min_count >= 0
  )
  stopifnot(
    "`R` must be a single whole number >= 1" =
      is.numeric(R) && length(R) == 1L && is.finite(R) && R >= 1L &&
      R == floor(R)
  )
  R <- as.integer(R)
  if (!adjust %in% stats::p.adjust.methods) {
    .lsa_stop(sprintf("`adjust` must be one of: %s.",
                 paste(stats::p.adjust.methods, collapse = ", ")),
         call. = FALSE)
  }

  resolved <- .compare_resolve(x, y)
  fits <- resolved$fits
  nm <- resolved$names
  .compare_validate_all(fits)

  # Adjusted residuals scale with sqrt(N), so differencing them across
  # groups of unequal size detects sample size rather than behaviour.
  # Warn rather than forbid (the user may want it deliberately).
  if (measure == "adj_res") {
    sizes <- vapply(fits, function(f) sum(f$obs), numeric(1))
    if (is.finite(max(sizes)) && min(sizes) > 0 &&
        max(sizes) / min(sizes) > 1.5) {
      message(sprintf(
        paste0("compare_lsa(): adjusted residuals scale with sqrt(N); ",
               "comparing them across groups %.1fx apart in size mostly ",
               "reflects sample size, not behaviour. The default ",
               "measure = \"log_or\" is N-invariant."),
        max(sizes) / min(sizes)))
    }
  }

  if (length(fits) == 2L) {
    return(.compare_two(fits[[1L]], fits[[2L]], nm, R = R,
                        measure = measure, adjust = adjust,
                        min_count = min_count,
                        parallel = parallel, n_cores = n_cores,
                        verbose = verbose))
  }
  .compare_pairwise(fits, nm, R = R, measure = measure, adjust = adjust,
                    min_count = min_count,
                    parallel = parallel, n_cores = n_cores,
                    verbose = verbose)
}

# --- core two-group comparison ----------------------------------------

# Permutation comparison of exactly two (already validated) fits.
# `grp_names` is the length-two character vector labelling a and b.
.compare_two <- function(fit_a, fit_b, grp_names, R, measure, adjust,
                         min_count = 5L, parallel, n_cores, verbose) {
  labels <- fit_a$data$labels
  K <- length(labels)
  recipe <- fit_a$params

  # Pool the per-sequence event vectors.
  seqs_a <- .compare_per_seq(fit_a)
  seqs_b <- .compare_per_seq(fit_b)
  pooled <- c(seqs_a, seqs_b)
  n_a <- length(seqs_a)
  S <- length(pooled)
  seq_len_pooled <- vapply(pooled, length, integer(1))

  # Observed per-edge difference.
  meas_a <- .lsa_compare_matrix(fit_a, measure)
  meas_b <- .lsa_compare_matrix(fit_b, measure)
  obs_diff <- as.vector(meas_a) - as.vector(meas_b)

  # Refit one pseudo-group from a set of pooled-sequence indices.
  refit_idx <- function(idx) {
    seqs <- pooled[idx]
    new_events <- unlist(seqs, use.names = FALSE)
    new_seq_id <- rep.int(seq_along(idx), times = seq_len_pooled[idx])
    .refit_from_events(events = new_events, seq_id = new_seq_id,
                       labels = labels, recipe = recipe)
  }

  worker <- function(b) {
    perm <- sample.int(S)
    idx_a <- perm[seq_len(n_a)]
    idx_b <- perm[(n_a + 1L):S]
    ma <- .lsa_compare_matrix(refit_idx(idx_a), measure)
    mb <- .lsa_compare_matrix(refit_idx(idx_b), measure)
    as.vector(ma) - as.vector(mb)
  }

  if (verbose) message("Running ", R, " label permutations ...")
  results <- .run_parallel(worker, R = R, parallel = parallel,
                           n_cores = n_cores, verbose = verbose)

  perm_diff <- matrix(NA_real_, R, K * K)
  for (b in seq_len(R)) perm_diff[b, ] <- results[[b]]

  # Per-edge p-values (Phipson-Smyth, NA-disciplined).
  abs_obs <- abs(obs_diff)
  exceed <- colSums(abs(perm_diff) >= matrix(abs_obs, R, K * K,
                                             byrow = TRUE),
                    na.rm = TRUE)
  n_finite <- colSums(is.finite(perm_diff))
  p_perm <- (1 + exceed) / (1 + n_finite)
  p_perm[!is.finite(obs_diff) | n_finite == 0L] <- NA_real_

  # Minimum-support filter. An odds ratio (and most cell statistics) is
  # unstable for transitions that barely occur, and an all-but-empty cell
  # yields a spuriously small permutation p because its null is nearly
  # degenerate. Cells whose pooled observed count is below `min_count`
  # are not tested (p = NA) and are dropped from the family correction
  # and the omnibus, rather than flagged significant.
  pooled_n <- as.vector(fit_a$obs) + as.vector(fit_b$obs)
  p_perm[pooled_n < min_count] <- NA_real_
  # `tested` is the single source of truth for which cells enter the
  # family correction AND the omnibus, so observed and permuted statistics
  # are always computed over the identical cell set (a cell that is
  # non-estimable in the observed data must not contribute to the null).
  tested <- is.finite(p_perm)
  p_adj <- p_perm
  p_adj[tested] <- stats::p.adjust(p_perm[tested], method = adjust)

  alpha <- recipe$alpha
  grid <- expand.grid(from = labels, to = labels,
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  edges <- data.frame(
    from        = grid$from,
    to          = grid$to,
    a           = as.vector(meas_a),
    b           = as.vector(meas_b),
    diff        = obs_diff,
    p_perm      = p_perm,
    p_adj       = p_adj,
    significant = is.finite(p_adj) & p_adj < alpha,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  names(edges)[names(edges) == "a"] <- paste0(measure, "_a")
  names(edges)[names(edges) == "b"] <- paste0(measure, "_b")

  # Omnibus: Frobenius norm of the difference matrix, over the SAME tested
  # cells as the per-edge test (rare and non-estimable cells excluded, so
  # the global stat is not driven by unstable near-empty odds ratios and
  # observed/null share one cell set). NA when nothing is testable.
  if (any(tested)) {
    stat_obs <- sum(obs_diff[tested]^2, na.rm = TRUE)
    stat_null <- rowSums(perm_diff[, tested, drop = FALSE]^2, na.rm = TRUE)
    p_global <- (1 + sum(stat_null >= stat_obs)) / (1 + R)
  } else {
    stat_obs <- NA_real_
    p_global <- NA_real_
  }

  structure(
    list(
      edges     = edges,
      global    = list(statistic = stat_obs, p_value = p_global, R = R),
      perm_diff = perm_diff,
      measure   = measure,
      R         = R,
      adjust    = adjust,
      groups    = grp_names,
      fits      = stats::setNames(list(fit_a, fit_b), grp_names)
    ),
    class = c("lsa_comparison", "list")
  )
}

# --- all-pairwise comparison ------------------------------------------

# Run every C(k, 2) pairwise comparison, then apply the multiple-
# comparison correction ONCE across the pooled per-edge p-values of all
# pairs (family-wise), and separately across the per-pair omnibus tests.
.compare_pairwise <- function(fits, nm, R, measure, adjust,
                              min_count = 5L, parallel, n_cores, verbose) {
  pairs <- utils::combn(length(fits), 2L)
  npair <- ncol(pairs)
  alpha <- fits[[1L]]$params$alpha

  # Each pair fit with adjust = "none"; family correction applied below.
  comps <- vector("list", npair)
  pair_names <- character(npair)
  for (p in seq_len(npair)) {
    i <- pairs[1L, p]
    j <- pairs[2L, p]
    if (verbose) message("Pair ", p, "/", npair, ": ",
                         nm[i], " vs ", nm[j])
    comps[[p]] <- .compare_two(fits[[i]], fits[[j]], nm[c(i, j)], R = R,
                               measure = measure, adjust = "none",
                               min_count = min_count,
                               parallel = parallel, n_cores = n_cores,
                               verbose = verbose)
    pair_names[p] <- paste0(nm[i], "_vs_", nm[j])
  }
  names(comps) <- pair_names

  meas_a <- paste0(measure, "_a")
  meas_b <- paste0(measure, "_b")
  keep_cols <- c("from", "to", meas_a, meas_b, "diff", "p_perm")
  combined <- do.call(rbind, lapply(seq_len(npair), function(p) {
    e <- comps[[p]]$edges
    data.frame(
      group_a = nm[pairs[1L, p]],
      group_b = nm[pairs[2L, p]],
      e[, keep_cols],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))

  # Family-wide correction across all per-edge p-values.
  combined$p_adj <- NA_real_
  ok <- is.finite(combined$p_perm)
  combined$p_adj[ok] <- stats::p.adjust(combined$p_perm[ok],
                                        method = adjust)
  combined$significant <- is.finite(combined$p_adj) &
    combined$p_adj < alpha

  # Per-pair omnibus, corrected across pairs.
  global <- data.frame(
    group_a   = nm[pairs[1L, ]],
    group_b   = nm[pairs[2L, ]],
    statistic = vapply(comps, function(c) c$global$statistic, numeric(1)),
    p_value   = vapply(comps, function(c) c$global$p_value, numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  global$p_adj <- stats::p.adjust(global$p_value, method = adjust)

  structure(
    list(
      edges       = combined,
      global      = global,
      comparisons = comps,
      measure     = measure,
      R           = R,
      adjust      = adjust,
      groups      = nm
    ),
    class = c("lsa_comparison_pairwise", "list")
  )
}

# --- helpers ----------------------------------------------------------

# Resolve the (x, y) argument forms into a plain list of fits plus their
# group names. Accepts an lsa_group (>= 2 groups), or two bare lsa fits.
.compare_resolve <- function(x, y) {
  if (inherits(x, "lsa_group")) {
    if (!is.null(y)) {
      .lsa_stop("When `x` is an lsa_group, leave `y = NULL`: the groups are ",
           "taken from `x`.", call. = FALSE)
    }
    if (length(x) < 2L) {
      .lsa_stop("compare_lsa() needs at least two groups; the lsa_group has ",
           length(x), ".", call. = FALSE)
    }
    nm <- names(x)
    if (is.null(nm)) nm <- paste0("g", seq_along(x))
    fits <- lapply(seq_along(x), function(i) x[[i]])
    names(fits) <- nm
    return(list(fits = fits, names = nm))
  }
  if (!inherits(x, "lsa") || !inherits(y, "lsa")) {
    .lsa_stop("Supply either an `lsa_group` as `x`, or two `lsa` fits as `x` ",
         "and `y`.", call. = FALSE)
  }
  list(fits = list(x, y), names = c("a", "b"))
}

# Validate every fit: event-level, and sharing the first fit's state
# set / engine / lag / structural zeros so differences are meaningful.
.compare_validate_all <- function(fits) {
  ref <- fits[[1L]]
  for (f in fits) {
    if (identical(f$data$source, "transitions")) {
      .lsa_stop(class = "needs_events", "compare_lsa() requires event-level input. A fit built from ",
           "a pre-computed transition matrix cannot be permuted.",
           call. = FALSE)
    }
  }
  for (k in seq_along(fits)[-1L]) {
    f <- fits[[k]]
    if (!identical(ref$data$labels, f$data$labels)) {
      .lsa_stop(class = "incompatible_fits", "All groups must share the same state set (labels). Fit them ",
           "together with lsa(..., group = ) so labels are aligned.",
           call. = FALSE)
    }
    if (!identical(ref$params$engine, f$params$engine)) {
      .lsa_stop(class = "incompatible_fits", sprintf("Engine mismatch: '%s' vs '%s'.",
                   ref$params$engine, f$params$engine), call. = FALSE)
    }
    if (!identical(ref$params$lag, f$params$lag)) {
      .lsa_stop(class = "incompatible_fits", sprintf("Lag mismatch: %s vs %s.",
                   ref$params$lag, f$params$lag), call. = FALSE)
    }
    if (!identical(ref$params$structural_zeros, f$params$structural_zeros)) {
      .lsa_stop(class = "incompatible_fits", "Groups declare different structural zeros.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Split a fit's events into a list of per-sequence integer vectors, in
# original sequence order (mirrors .lsa_grouped()).
.compare_per_seq <- function(fit) {
  d <- fit$data
  split(d$events, factor(d$seq_id, levels = seq_len(d$n_sequences)))
}

# Raw measure matrix, NA/Inf preserved (NOT coerced to 0), so non-
# estimable cells stay out of the test. Works on a full lsa fit or on
# the bare engine-output list from .refit_from_events() (both carry
# $obs, $exp, $prob, $adj_res).
.lsa_compare_matrix <- function(obj, measure) {
  switch(measure,
    log_or  = .lsa_log_or(obj),
    adj_res = obj$adj_res,
    yules_q = obj$yules_q,
    prob    = obj$prob,
    count   = obj$obs,
    lift    = obj$obs / obj$exp
  )
}

# Per-cell log odds ratio of the 2x2 collapse for each transition i->j:
#   a = obs[i,j]              b = row_total[i] - a
#   c = col_total[j] - a      d = N - row_total[i] - col_total[j] + a
# log OR = log(a*d / (b*c)), with a Haldane-Anscombe +0.5 correction
# applied to all four cells whenever any is zero (standard for empty
# cells). Unbounded and additive on the log scale, so it is the
# N-invariant LSA effect size used to compare groups (Yule's Q saturates
# at +/-1 on any zero cell and compresses exactly the sparse tail LSA
# cares about). Cells the engine marks non-estimable (NA adj_res:
# structural zeros, zero-margin rows) stay NA. Works on a full lsa fit
# or the bare engine-output list (both carry $obs and $adj_res).
.lsa_log_or <- function(obj) {
  O <- obj$obs
  nr <- nrow(O); nc <- ncol(O)
  rt <- matrix(rowSums(O), nr, nc)            # row_total[i] in every col
  ct <- matrix(colSums(O), nr, nc, byrow = TRUE)  # col_total[j] in every row
  N <- sum(O)
  a <- O
  b <- rt - O
  cc <- ct - O
  d <- N - rt - ct + O
  zero <- a == 0 | b == 0 | cc == 0 | d == 0
  h <- 0.5 * zero
  lor <- log(((a + h) * (d + h)) / ((b + h) * (cc + h)))
  lor[is.na(obj$adj_res)] <- NA_real_         # honour engine NA discipline
  dimnames(lor) <- dimnames(O)
  lor
}

#' Read the tables of a group comparison
#'
#' @description
#' `as.data.frame()` returns the tables of a group comparison from
#' [compare_lsa()] or [bayes_compare_lsa()] as a plain data frame, and `print()`
#' shows a short summary of the comparison. A comparison holds two tables,
#' chosen with `what`: the comparison of every transition and the omnibus result
#' for the structure as a whole. The per-transition table is the one
#' [transitions()] reads; use [transitions()] when you want to keep only the
#' significant transitions or order them by the size of the difference.
#'
#' @param x A comparison from [compare_lsa()] or [bayes_compare_lsa()]
#'   (class `lsa_comparison` or `lsa_comparison_pairwise`).
#' @param what `"edges"` (default): one row per transition. `"global"`:
#'   the omnibus result, one row per pair of groups.
#' @param row.names,optional Not used; present for consistency with the
#'   generic.
#' @param ... Not used.
#'
#' @return `as.data.frame()` returns a data frame.
#'   * `what = "edges"`: one row per transition (and per pair of groups
#'     when more than two groups were compared, with leading `group_a`,
#'     `group_b` columns); the columns are listed under *Value* in
#'     [compare_lsa()] and [bayes_compare_lsa()].
#'   * `what = "global"`: one row per pair of groups with `group_a`,
#'     `group_b` and, for [compare_lsa()], the omnibus `statistic` (sum of
#'     squared differences), its permutation `p_value` and either the
#'     number of permutations `R` (two groups) or the across-pairs
#'     adjusted `p_adj` (more than two groups); for [bayes_compare_lsa()],
#'     `n_credible`, the number of credibly different transitions.
#'
#'   `print()` returns `x` invisibly.
#'
#' @examples
#' g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'          time = "Time", group = "Achiever")
#' cmp <- compare_lsa(g, R = 50, adjust = "BH")
#' cmp
#' as.data.frame(cmp, what = "global")
#' head(as.data.frame(cmp))
#'
#' @seealso [compare_lsa()], [bayes_compare_lsa()], [transitions()].
#'
#' @export
as.data.frame.lsa_comparison <- function(x, row.names = NULL,
                                         optional = FALSE,
                                         what = c("edges", "global"), ...) {
  what <- match.arg(what)
  if (identical(what, "edges")) return(x$edges)
  .comparison_global_df(x)
}

#' @rdname as.data.frame.lsa_comparison
#' @export
as.data.frame.lsa_comparison_pairwise <- function(x, row.names = NULL,
                                                  optional = FALSE,
                                                  what = c("edges", "global"),
                                                  ...) {
  what <- match.arg(what)
  if (identical(what, "edges")) return(x$edges)
  .comparison_global_df(x)
}

# The omnibus table of a comparison as one row per pair of groups. A
# two-group result stores it as a list of scalars; a pairwise result
# already stores a data frame.
.comparison_global_df <- function(x) {
  g <- x$global
  if (is.data.frame(g)) {
    out <- g
  } else {
    out <- data.frame(g, stringsAsFactors = FALSE)
    if (!"group_a" %in% names(out)) {
      out <- cbind(data.frame(group_a = x$groups[1L],
                              group_b = x$groups[2L],
                              stringsAsFactors = FALSE), out)
    }
  }
  rownames(out) <- NULL
  out
}

#' @rdname as.data.frame.lsa_comparison
#' @export
print.lsa_comparison <- function(x, ...) {
  cat("<lsa_comparison>\n")
  cat(sprintf("  groups:   %s vs %s\n", x$groups[1L], x$groups[2L]))
  cat(sprintf("  measure:  %s difference (%s - %s)\n",
              x$measure, x$groups[1L], x$groups[2L]))
  cat(sprintf("  R:        %d label permutations\n", x$R))
  tested <- sum(is.finite(x$edges$p_perm))
  nsig <- sum(x$edges$significant, na.rm = TRUE)
  cat(sprintf("  edges:    %d significant of %d tested (adjust = %s)\n",
              nsig, tested, x$adjust))
  cat(sprintf("  omnibus:  statistic = %.4g, p = %.4g\n",
              x$global$statistic, x$global$p_value))
  invisible(x)
}

#' @rdname as.data.frame.lsa_comparison
#' @export
print.lsa_comparison_pairwise <- function(x, ...) {
  cat("<lsa_comparison_pairwise>\n")
  cat(sprintf("  groups:   %d (%s)\n", length(x$groups),
              paste(x$groups, collapse = ", ")))
  cat(sprintf("  measure:  %s difference\n", x$measure))
  cat(sprintf("  R:        %d label permutations per pair\n", x$R))
  cat(sprintf("  pairs:    %d (adjust = %s, family-wide)\n",
              nrow(x$global), x$adjust))
  for (p in seq_len(nrow(x$global))) {
    ga <- x$global$group_a[p]
    gb <- x$global$group_b[p]
    rows <- x$edges$group_a == ga & x$edges$group_b == gb
    nsig <- sum(x$edges$significant[rows], na.rm = TRUE)
    tested <- sum(is.finite(x$edges$p_perm[rows]))
    cat(sprintf("    - %-20s %d/%d edges sig, omnibus p_adj = %.4g\n",
                paste0(ga, " vs ", gb), nsig, tested,
                x$global$p_adj[p]))
  }
  invisible(x)
}
