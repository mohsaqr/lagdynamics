# Bayesian Dirichlet-Multinomial comparison of two transition networks --
# the analytic Bayesian sibling of compare_lsa() (permutation), as
# certainty_lsa() is the sibling of bootstrap_lsa().
#
# Where the permutation test returns "is this difference more extreme than
# chance?", the Bayesian comparison returns "what is the plausible range of
# the true difference, and how precisely is it estimated from the counts?".
# Each state's outgoing transitions are modelled as Dirichlet-Multinomial
# (Jeffreys prior), so each edge probability is marginally Beta. The
# per-edge posterior mean difference (prob_a - prob_b) is closed form; a
# credible interval and the probability of direction come from a vectorised
# Monte Carlo draw on the Beta marginals.
#
# Ported/adapted from Nestimate::bayes_compare() (Johnston & Jendoubi 2026).
# Returns an lsa_comparison-compatible object so the barrel / heatmap plots
# and as.data.frame() apply unchanged.

#' Compare transition probabilities between groups with a Bayesian model
#'
#' @description
#' `bayes_compare_lsa()` compares the transition probabilities of groups with a
#' Bayesian model: given the counts, it estimates the plausible range of the
#' true difference between two groups' transition probabilities and how
#' precisely that difference is estimated. Under a Dirichlet-Multinomial model
#' the outgoing transitions of each state are a multinomial draw whose
#' probabilities have a Dirichlet posterior, so each transition probability has
#' a Beta posterior in each group. A transition whose source state is rarely
#' visited gets a wide interval, even when its estimated probabilities look far
#' apart. A permutation test ([compare_lsa()]) instead tests whether a group
#' difference is larger than chance.
#'
#' For every transition, `bayes_compare_lsa()` computes the posterior mean
#' difference of the transition probability between two groups (exact), its
#' credible interval and the probability of direction (the posterior probability
#' of the more probable sign of the difference), both from Monte Carlo draws of
#' the Beta posteriors, and a standardised effect size. A transition is credibly
#' different when its credible interval excludes zero, the absolute mean
#' difference exceeds `mean_threshold`, and the interval end nearest zero lies
#' further from zero than `bound_threshold`; the thresholds keep detectable but
#' negligible differences out. With more than two groups every pair is compared.
#' No resampling of sequences is involved, so the comparison is fast.
#'
#' The model treats transitions as independent draws and compares
#' transition probabilities only. For a test that respects the dependence
#' between events of one sequence, or that compares adjusted residuals or
#' odds ratios, use [compare_lsa()]. For the posterior of one group's
#' transition probabilities, use [certainty_lsa()].
#'
#' @param x A grouped fit from `lsa(data, group = )` with two or more
#'   groups, or the `lsa` fit of the first group. For the `print()`
#'   methods, the comparison result.
#' @param y The `lsa` fit of the second group when `x` is a single fit
#'   (the groups are then labelled `"a"` and `"b"`). Default `NULL`.
#' @param prior Dirichlet concentration added to every cell. Default
#'   `0.5` (the Jeffreys prior); `1` gives a uniform prior.
#' @param draws Number of Monte Carlo draws for the credible intervals.
#'   Default `10000`.
#' @param ci Credible-interval level. Default `0.95`.
#' @param mean_threshold Minimum absolute posterior mean difference for a
#'   transition to be credibly different. Default `0.01`.
#' @param bound_threshold Minimum distance from zero of the interval end
#'   nearest zero. Default `0.001`.
#' @param adjust Multiple-comparison correction applied to the two-sided
#'   Bayesian p-value `2 * (1 - pd)` across transitions (and across all
#'   pairs for more than two groups); any method of [stats::p.adjust()].
#'   Default `"none"`. It changes the `p_adj` column only: the
#'   `significant` flag is decided by the credible interval and the
#'   thresholds.
#' @param seed Optional seed for the Monte Carlo draws. Default `NULL`.
#'   The session's random number state is restored on exit.
#' @param ... Not used.
#'
#' @return With two groups, an object of class
#'   `c("lsa_bayes", "lsa_comparison", "list")`; with more than two, an
#'   object of class `c("lsa_bayes_pairwise", "lsa_comparison_pairwise",
#'   "list")`. Read the per-transition table with [transitions()]: one row
#'   per transition (`from -> to` pair), and for more than two groups one
#'   row per transition and pair of groups, with the columns
#'   * `group_a`, `group_b` (more than two groups only): the pair compared;
#'   * `from`, `to`: source and target state;
#'   * `prob_a`, `prob_b`: posterior mean transition probability in each
#'     group;
#'   * `diff`: `prob_a - prob_b`; `ci_low`, `ci_high`: its credible
#'     interval;
#'   * `pd`: probability of direction; `effect_size`: `diff` divided by
#'     its posterior standard deviation;
#'   * `p_value`: `2 * (1 - pd)` (`NA` for transitions that cannot be
#'     estimated); `p_adj`: after `adjust`;
#'   * `significant`: credibly different.
#'
#'   `transitions(x, significant = TRUE)` keeps the credibly different
#'   transitions; `sort = "strength"` orders by the absolute difference.
#'   `as.data.frame(x, what = "global")` gives the number of credibly
#'   different transitions per pair of groups (see
#'   [as.data.frame.lsa_comparison()]); [plot()][plot.lsa_comparison]
#'   draws the differences. `print()` shows the settings and the number
#'   of credibly different transitions and returns `x` invisibly.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and one of:
#' * `lagdynamics_needs_events`: a fit was built from a transition matrix.
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
#' bc <- bayes_compare_lsa(g, adjust = "BH", seed = 1)
#' bc
#' transitions(bc, significant = TRUE, sort = "strength") |> head()
#'
#' # Three courses: every pair
#' courses <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'                time = "Time", group = "Course")
#' bayes_compare_lsa(courses, draws = 2000, seed = 1) |>
#'   as.data.frame(what = "global")
#'
#' @references
#' Johnston, L., & Jendoubi, T. (2026). How delivery mode reshapes
#' resource engagement: A Bayesian differential network analysis.
#' TNA Workshop 2026.
#'
#' @seealso [transitions()] and [as.data.frame.lsa_comparison()] to read
#'   the result; [plot.lsa_comparison()] to draw it; [compare_lsa()] for
#'   the permutation test; [certainty_lsa()] for one group's posterior.
#'
#' @export
bayes_compare_lsa <- function(x,
                              y = NULL,
                              prior = 0.5,
                              draws = 10000L,
                              ci = 0.95,
                              mean_threshold = 0.01,
                              bound_threshold = 0.001,
                              adjust = "none",
                              seed = NULL) {
  stopifnot(
    "`prior` must be a single positive number" =
      is.numeric(prior) && length(prior) == 1L && is.finite(prior) &&
      prior > 0,
    "`draws` must be a single number >= 2" =
      is.numeric(draws) && length(draws) == 1L && is.finite(draws) &&
      draws >= 2,
    "`ci` must be a single number in (0, 1)" =
      is.numeric(ci) && length(ci) == 1L && is.finite(ci) && ci > 0 && ci < 1,
    "`mean_threshold` must be a single number" =
      is.numeric(mean_threshold) && length(mean_threshold) == 1L,
    "`bound_threshold` must be a single number" =
      is.numeric(bound_threshold) && length(bound_threshold) == 1L
  )
  draws <- as.integer(draws)
  if (!adjust %in% stats::p.adjust.methods) {
    .lsa_stop(sprintf("`adjust` must be one of: %s.",
                 paste(stats::p.adjust.methods, collapse = ", ")),
         call. = FALSE)
  }
  if (!is.null(seed)) {
    stopifnot("`seed` must be a single number" =
                is.numeric(seed) && length(seed) == 1L)
    old_seed <- .lsa_get_seed()
    on.exit(.lsa_restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  resolved <- .compare_resolve(x, y)        # reused from compare_lsa
  fits <- resolved$fits
  nm <- resolved$names
  .compare_validate_all(fits)

  if (length(fits) == 2L) {
    return(.bayes_compare_two(fits[[1L]], fits[[2L]], nm, prior = prior,
                              draws = draws, ci = ci,
                              mean_threshold = mean_threshold,
                              bound_threshold = bound_threshold,
                              adjust = adjust))
  }
  .bayes_compare_pairwise(fits, nm, prior = prior, draws = draws, ci = ci,
                          mean_threshold = mean_threshold,
                          bound_threshold = bound_threshold, adjust = adjust)
}

# --- core two-group Bayesian comparison -------------------------------

.bayes_compare_two <- function(fit_a, fit_b, grp_names, prior, draws, ci,
                               mean_threshold, bound_threshold, adjust) {
  labels <- rownames(fit_a$obs)
  K <- length(labels)
  dn <- list(labels, labels)

  # Dirichlet posterior per source row; each edge probability is Beta(a, b).
  aa <- fit_a$obs + prior
  ab <- fit_b$obs + prior
  rowA_a <- matrix(rowSums(aa), K, K, dimnames = dn)
  rowA_b <- matrix(rowSums(ab), K, K, dimnames = dn)
  prob_a <- aa / rowA_a
  prob_b <- ab / rowA_b
  diff_mat <- prob_a - prob_b

  # Credible interval / probability of direction via MC on the Beta
  # marginals (one vectorised rbeta call per group).
  n <- K * K
  ba <- as.vector(rowA_a - aa); bb <- as.vector(rowA_b - ab)
  da <- matrix(stats::rbeta(n * draws, rep(as.vector(aa), draws),
                            rep(ba, draws)), n, draws)
  db <- matrix(stats::rbeta(n * draws, rep(as.vector(ab), draws),
                            rep(bb, draws)), n, draws)
  dd <- da - db
  tail <- (1 - ci) / 2
  ci_low  <- apply(dd, 1L, stats::quantile, probs = tail, names = FALSE)
  ci_high <- apply(dd, 1L, stats::quantile, probs = 1 - tail, names = FALSE)
  prop_pos <- rowMeans(dd > 0)
  pd <- pmax(prop_pos, 1 - prop_pos)              # probability of direction
  p_bayes <- 2 * (1 - pd)                         # two-sided p-equivalent
  sd_d <- sqrt(pmax(rowMeans(dd^2) - rowMeans(dd)^2, 0))

  diff_v <- as.vector(diff_mat)
  eff <- diff_v / sd_d
  eff[!is.finite(eff)] <- 0

  # Non-estimable cells (zero-margin rows, structural zeros) stay out.
  valid <- is.finite(as.vector(fit_a$prob)) & is.finite(as.vector(fit_b$prob))
  sz <- fit_a$params$structural_zeros
  if (!is.null(sz)) valid <- valid & (as.vector(sz) != 0)

  ci_excl <- (ci_low > 0) | (ci_high < 0)
  nearest <- pmin(abs(ci_low), abs(ci_high))
  sig <- valid & ci_excl & abs(diff_v) > mean_threshold &
    nearest > bound_threshold

  p_bayes[!valid] <- NA_real_
  p_adj <- p_bayes
  ok <- is.finite(p_bayes)
  p_adj[ok] <- stats::p.adjust(p_bayes[ok], method = adjust)

  grid <- expand.grid(from = labels, to = labels,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  edges <- data.frame(
    from = grid$from, to = grid$to,
    a = as.vector(prob_a), b = as.vector(prob_b),
    diff = diff_v, ci_low = ci_low, ci_high = ci_high,
    pd = pd, effect_size = eff,
    p_value = p_bayes, p_adj = p_adj, significant = sig,
    stringsAsFactors = FALSE, row.names = NULL)
  names(edges)[names(edges) == "a"] <- "prob_a"
  names(edges)[names(edges) == "b"] <- "prob_b"

  structure(
    list(
      edges = edges,
      global = list(group_a = grp_names[1L], group_b = grp_names[2L],
                    n_credible = sum(sig, na.rm = TRUE)),
      measure = "prob", R = draws, adjust = adjust, groups = grp_names,
      fits = stats::setNames(list(fit_a, fit_b), grp_names),
      prior = prior, draws = draws, ci = ci,
      mean_threshold = mean_threshold, bound_threshold = bound_threshold
    ),
    class = c("lsa_bayes", "lsa_comparison", "list"))
}

# --- all-pairwise -----------------------------------------------------

.bayes_compare_pairwise <- function(fits, nm, prior, draws, ci,
                                    mean_threshold, bound_threshold, adjust) {
  pairs <- utils::combn(length(fits), 2L)
  npair <- ncol(pairs)
  comps <- vector("list", npair)
  pair_names <- character(npair)
  for (p in seq_len(npair)) {
    i <- pairs[1L, p]; j <- pairs[2L, p]
    comps[[p]] <- .bayes_compare_two(fits[[i]], fits[[j]], nm[c(i, j)],
                                     prior = prior, draws = draws, ci = ci,
                                     mean_threshold = mean_threshold,
                                     bound_threshold = bound_threshold,
                                     adjust = "none")
    pair_names[p] <- paste0(nm[i], "_vs_", nm[j])
  }
  names(comps) <- pair_names

  keep <- c("from", "to", "prob_a", "prob_b", "diff", "ci_low", "ci_high",
            "pd", "effect_size", "p_value")
  combined <- do.call(rbind, lapply(seq_len(npair), function(p) {
    data.frame(group_a = nm[pairs[1L, p]], group_b = nm[pairs[2L, p]],
               comps[[p]]$edges[, keep], stringsAsFactors = FALSE,
               row.names = NULL)
  }))
  combined$p_adj <- NA_real_
  ok <- is.finite(combined$p_value)
  combined$p_adj[ok] <- stats::p.adjust(combined$p_value[ok], method = adjust)
  # Recompute the credible-difference flag under the family-wide p, keeping
  # the same magnitude guards as the two-group test.
  nearest <- pmin(abs(combined$ci_low), abs(combined$ci_high))
  combined$significant <- ok & (combined$ci_low > 0 | combined$ci_high < 0) &
    abs(combined$diff) > mean_threshold & nearest > bound_threshold

  global <- data.frame(
    group_a = nm[pairs[1L, ]], group_b = nm[pairs[2L, ]],
    n_credible = vapply(comps, function(c) sum(c$edges$significant,
                                               na.rm = TRUE), integer(1)),
    stringsAsFactors = FALSE, row.names = NULL)

  structure(
    list(edges = combined, global = global, comparisons = comps,
         measure = "prob", R = draws, adjust = adjust, groups = nm,
         prior = prior, draws = draws, ci = ci,
         mean_threshold = mean_threshold, bound_threshold = bound_threshold),
    class = c("lsa_bayes_pairwise", "lsa_comparison_pairwise", "list"))
}

# --- print methods (Bayesian framing; dispatch before lsa_comparison) ---

#' @rdname bayes_compare_lsa
#' @export
print.lsa_bayes <- function(x, ...) {
  cat("<lsa_bayes>  (Bayesian Dirichlet-Multinomial comparison)\n")
  cat(sprintf("  groups:    %s vs %s\n", x$groups[1L], x$groups[2L]))
  cat(sprintf("  prior:     Dirichlet(%.2f)  |  draws: %d  |  CI: %.0f%%\n",
              x$prior, x$draws, 100 * x$ci))
  tested <- sum(is.finite(x$edges$p_value))
  ncred <- sum(x$edges$significant, na.rm = TRUE)
  cat(sprintf("  edges:     %d credibly different of %d compared\n",
              ncred, tested))
  invisible(x)
}

#' @rdname bayes_compare_lsa
#' @export
print.lsa_bayes_pairwise <- function(x, ...) {
  cat("<lsa_bayes_pairwise>\n")
  cat(sprintf("  groups:    %d (%s)\n", length(x$groups),
              paste(x$groups, collapse = ", ")))
  cat(sprintf("  prior:     Dirichlet(%.2f)  |  draws: %d  |  CI: %.0f%%\n",
              x$prior, x$draws, 100 * x$ci))
  for (p in seq_len(nrow(x$global))) {
    cat(sprintf("    - %-20s %d credibly-different edges\n",
                paste0(x$global$group_a[p], " vs ", x$global$group_b[p]),
                x$global$n_credible[p]))
  }
  invisible(x)
}
