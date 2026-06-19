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

#' Compare Groups' Transition Structures
#'
#' Permutation test for whether groups produce different LSA transition
#' structures. For each pair of groups it pools their sequences,
#' repeatedly reassigns the group label of whole sequences (preserving
#' the original group sizes), refits each pseudo-group, and builds a
#' permutation distribution of the per-edge difference in the chosen
#' `measure`. The two-sided p-value is
#' `(1 + #{ |diff_perm| >= |diff_obs| }) / (1 + R)` (Phipson & Smyth
#' 2010). A single omnibus test of overall difference is reported from
#' the same permutations.
#'
#' With exactly two groups a single comparison is returned. With more
#' than two groups every pairwise comparison is run and the requested
#' `adjust` correction is applied **once across the whole family** of
#' per-edge p-values (and separately across the per-pair omnibus tests),
#' giving family-wise control rather than per-pair control.
#'
#' @param x Either an `lsa_group` object (from `lsa(..., group = )`)
#'   with two or more groups, or a single `lsa` fit for the first group.
#' @param y When `x` is a single `lsa` fit, the second group's `lsa`
#'   fit. Ignored (and must be `NULL`) when `x` is an `lsa_group`.
#' @param R Integer. Number of label permutations per comparison.
#'   Default `1000`.
#' @param measure Character. The per-edge quantity compared between
#'   groups. Default `"log_or"`: the per-cell log odds ratio of the 2x2
#'   transition collapse (Haldane-Anscombe corrected on empty cells) --
#'   an N-invariant LSA effect size, so the comparison reflects
#'   behaviour rather than sample size. Other options: `"yules_q"`
#'   (also N-invariant, but saturates at +/-1 on zero cells),
#'   `"adj_res"` (adjusted residuals -- the LSA *test statistic*, which
#'   scales with sqrt(N) and is therefore confounded by group size; a
#'   message warns when groups differ in size), `"prob"` (transition
#'   probabilities -- a raw rate, i.e. the TNA quantity, with no
#'   independence baseline), `"count"`, or `"lift"` (observed /
#'   expected).
#' @param adjust Multiple-comparison correction; any method accepted by
#'   [stats::p.adjust()] (e.g. `"holm"`, `"BH"`, `"bonferroni"`).
#'   Default `"none"`. For more than two groups it is applied across the
#'   pooled per-edge p-values of all pairs.
#' @param min_count Integer. Minimum pooled observed count (group a +
#'   group b) for a transition to be tested. Default `5`. Rarer cells
#'   carry an unstable odds ratio and a near-degenerate permutation null
#'   that produces spurious small p-values, so they get `p = NA` and are
#'   excluded from the multiple-comparison family and the omnibus rather
#'   than flagged significant. Set `0` to test every cell.
#' @param parallel Logical. Use multi-core resampling. Default `FALSE`.
#' @param n_cores Integer. Worker count when `parallel = TRUE`.
#' @param verbose Logical. Print progress every 100 permutations.
#' @param ... Reserved.
#'
#' **NA handling.** Non-estimable cells (structural zeros, zero-margin
#' rows in a permuted pseudo-group) carry `NA` in the measure matrix and
#' are never coerced to zero. Such cells get `p_perm = NA` rather than a
#' spurious significant flag, and the exceedance tally and omnibus
#' statistic are computed with `na.rm = TRUE`, matching [permute_lsa()].
#'
#' **Interpretation caveats.** The odds ratio is non-collapsible: the
#' per-group log odds ratios are group-specific departure-from-
#' independence measures and should not be pooled across groups that
#' have different marginal state distributions. As with any LSA, a
#' between-group difference can also be driven by subgroup composition
#' (Simpson's paradox); confirm with subgroup analysis when a confound
#' is plausible.
#'
#' @return For two groups, an object of class
#'   `c("lsa_comparison", "list")` with:
#' \describe{
#'   \item{edges}{Tidy per-edge data frame: `from`, `to`, the measure in
#'     each group (`<measure>_a`, `<measure>_b`), their difference
#'     `diff` (= a - b), the permutation p-value `p_perm`, the adjusted
#'     p-value `p_adj`, and a `significant` flag.}
#'   \item{global}{Omnibus test list: `statistic` (observed sum of
#'     squared edge differences), `p_value`, and `R`.}
#'   \item{perm_diff}{`R x K^2` matrix of permuted edge differences.}
#'   \item{measure, R, adjust, groups}{Call metadata; `groups` is the
#'     length-two character vector of group labels (a, b).}
#'   \item{fits}{The two original fits, named by group.}
#' }
#'   For more than two groups, an object of class
#'   `c("lsa_comparison_pairwise", "list")` with:
#' \describe{
#'   \item{edges}{Tidy per-edge data frame across all pairs, prefixed by
#'     `group_a`, `group_b`; `p_adj` and `significant` reflect the
#'     family-wide correction.}
#'   \item{global}{One row per pair: `group_a`, `group_b`, `statistic`,
#'     `p_value`, and the across-pairs adjusted `p_adj`.}
#'   \item{comparisons}{Named list of the underlying two-group
#'     `lsa_comparison` objects (each fit with `adjust = "none"`), for
#'     drill-down.}
#'   \item{measure, R, adjust, groups}{Call metadata; `groups` lists all
#'     group labels.}
#' }
#'
#' @examples
#' \donttest{
#' # group_regulation is wide sequences with no grouping column, so
#' # derive one: sessions whose first regulation act is planning vs not.
#' grp <- ifelse(group_regulation$T1 == "plan", "starts_plan", "other")
#' g <- lsa(group_regulation, group = grp)
#' cmp <- compare_lsa(g, R = 200)
#' head(cmp$edges)
#' cmp$global
#' }
#'
#' @references
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should
#' never be zero. \emph{Statistical Applications in Genetics and
#' Molecular Biology}, 9(1), Article 39.
#'
#' van Borkulo, C. D., et al. (2022). Comparing network structures on
#' three aspects: A permutation test. \emph{Psychological Methods}.
#'
#' @seealso [permute_lsa()], [bootstrap_lsa()]
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
  stopifnot(is.numeric(min_count), length(min_count) == 1L,
            is.finite(min_count), min_count >= 0)
  stopifnot(is.numeric(R), length(R) == 1L, is.finite(R),
            R >= 1L, R == floor(R))
  R <- as.integer(R)
  if (!adjust %in% stats::p.adjust.methods) {
    stop(sprintf("`adjust` must be one of: %s.",
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
      stop("When `x` is an lsa_group, leave `y = NULL`: the groups are ",
           "taken from `x`.", call. = FALSE)
    }
    if (length(x) < 2L) {
      stop("compare_lsa() needs at least two groups; the lsa_group has ",
           length(x), ".", call. = FALSE)
    }
    nm <- names(x)
    if (is.null(nm)) nm <- paste0("g", seq_along(x))
    fits <- lapply(seq_along(x), function(i) x[[i]])
    names(fits) <- nm
    return(list(fits = fits, names = nm))
  }
  if (!inherits(x, "lsa") || !inherits(y, "lsa")) {
    stop("Supply either an `lsa_group` as `x`, or two `lsa` fits as `x` ",
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
      stop("compare_lsa() requires event-level input. A fit built from ",
           "a pre-computed transition matrix cannot be permuted.",
           call. = FALSE)
    }
  }
  for (k in seq_along(fits)[-1L]) {
    f <- fits[[k]]
    if (!identical(ref$data$labels, f$data$labels)) {
      stop("All groups must share the same state set (labels). Fit them ",
           "together with lsa(..., group = ) so labels are aligned.",
           call. = FALSE)
    }
    if (!identical(ref$params$engine, f$params$engine)) {
      stop(sprintf("Engine mismatch: '%s' vs '%s'.",
                   ref$params$engine, f$params$engine), call. = FALSE)
    }
    if (!identical(ref$params$lag, f$params$lag)) {
      stop(sprintf("Lag mismatch: %s vs %s.",
                   ref$params$lag, f$params$lag), call. = FALSE)
    }
    if (!identical(ref$params$structural_zeros, f$params$structural_zeros)) {
      stop("Groups declare different structural zeros.", call. = FALSE)
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

#' Tidy a Group Comparison
#'
#' Returns the per-edge comparison table (the same data frame as
#' `x$edges`) so a comparison can be read with `as.data.frame()` like the
#' other result objects, without reaching into the object.
#'
#' @param x An `lsa_comparison` or `lsa_comparison_pairwise` object.
#' @param row.names,optional,... Standard [as.data.frame()] arguments
#'   (unused; present for method consistency).
#' @return The tidy per-edge data frame.
#' @export
as.data.frame.lsa_comparison <- function(x, row.names = NULL,
                                         optional = FALSE, ...) {
  x$edges
}

#' @rdname as.data.frame.lsa_comparison
#' @export
as.data.frame.lsa_comparison_pairwise <- function(x, row.names = NULL,
                                                  optional = FALSE, ...) {
  x$edges
}

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
