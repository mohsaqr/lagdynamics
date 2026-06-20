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

#' Bayesian Comparison of Group Transition Structures (Dirichlet-Multinomial)
#'
#' Closed-form Bayesian alternative to [compare_lsa()] for comparing the
#' transition structures of two (or, pairwise, more) groups. Each state's
#' outgoing transitions are modelled as Dirichlet-Multinomial with a
#' Jeffreys prior, so each transition probability is marginally Beta. The
#' per-edge posterior mean difference `prob_a - prob_b` is exact; a
#' credible interval, the probability of direction `pd`, and a two-sided
#' Bayesian p-equivalent `2 * (1 - pd)` come from a Monte Carlo draw on the
#' Beta marginals.
#'
#' This complements [compare_lsa()]: the permutation test asks whether a
#' difference is more extreme than chance; the Bayesian comparison asks for
#' the plausible range of the true difference and how precisely it is
#' estimated. An edge whose source state is rarely visited gets a wide
#' credible interval even when its row-normalised probability looks
#' decisive.
#'
#' The result carries class `c("lsa_bayes", "lsa_comparison")` (and the
#' pairwise object `c("lsa_bayes_pairwise", "lsa_comparison_pairwise")`),
#' so [plot()][plot.lsa_comparison] (barrel / heatmap) and
#' [as.data.frame()] work as for a permutation comparison.
#'
#' @param x An `lsa_group` (two or more groups), or a single `lsa` fit for
#'   the first group.
#' @param y The second group's `lsa` fit when `x` is a single fit;
#'   otherwise `NULL`.
#' @param prior Numeric > 0. Dirichlet prior concentration added to every
#'   cell. Default `0.5` (Jeffreys). Use `1` for a uniform (Laplace) prior.
#' @param draws Integer. Monte Carlo posterior draws for the credible
#'   intervals. Default `10000`.
#' @param ci Numeric in (0, 1). Credible-interval mass. Default `0.95`.
#' @param mean_threshold,bound_threshold An edge is flagged credibly
#'   different only if its credible interval excludes zero, `|posterior
#'   mean diff|` exceeds `mean_threshold` (default `0.01`), and the
#'   credible bound nearest zero exceeds `bound_threshold` (default
#'   `0.001`). The thresholds guard against differences that are detectable
#'   but negligibly small.
#' @param adjust Multiple-comparison correction applied to the two-sided
#'   Bayesian p across edges (and family-wide across pairs); any method of
#'   [stats::p.adjust()]. Default `"none"`.
#' @param seed Optional integer for reproducible credible intervals.
#'
#' @return For two groups, class `c("lsa_bayes", "lsa_comparison", "list")`
#'   with an `edges` data frame (`from`, `to`, `prob_a`, `prob_b`, `diff`,
#'   `ci_low`, `ci_high`, `pd`, `effect_size`, `p_value`, `p_adj`,
#'   `significant`), the two `fits`, and the Bayesian settings. For more
#'   than two groups, an all-pairwise `c("lsa_bayes_pairwise",
#'   "lsa_comparison_pairwise", "list")`.
#'
#' @examples
#' \donttest{
#' g <- lsa(group_regulation,
#'          group = ifelse(group_regulation$T1 == "plan", "p", "o"))
#' bc <- bayes_compare_lsa(g, seed = 1)
#' head(as.data.frame(bc))
#' }
#'
#' @references
#' Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
#' Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.
#'
#' @seealso [compare_lsa()], [certainty_lsa()]
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
  stopifnot(is.numeric(prior), length(prior) == 1L, is.finite(prior),
            prior > 0,
            is.numeric(draws), length(draws) == 1L, draws >= 2,
            is.numeric(ci), length(ci) == 1L, ci > 0, ci < 1,
            is.numeric(mean_threshold), length(mean_threshold) == 1L,
            is.numeric(bound_threshold), length(bound_threshold) == 1L)
  draws <- as.integer(draws)
  if (!adjust %in% stats::p.adjust.methods) {
    stop(sprintf("`adjust` must be one of: %s.",
                 paste(stats::p.adjust.methods, collapse = ", ")),
         call. = FALSE)
  }
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L)
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
