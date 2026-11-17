# Analytic Dirichlet-Multinomial certainty -- the closed-form Bayesian
# sibling of bootstrap_lsa(). Where the bootstrap resamples sequences to
# quantify edge uncertainty, certainty_lsa() derives it in closed form:
# the outgoing transitions from each state follow a Dirichlet-Multinomial
# posterior (Jeffreys prior), so each transition probability is marginally
# Beta-distributed and its posterior mean, sd, credible interval and the
# stability decision all have exact formulas -- no resampling.
#
# Ported from Nestimate::certainty() (Johnston & Jendoubi 2026), adapted
# to lagdynamics's lsa fit (works off the fit's `obs` count matrix). Returns an
# lsa_bootstrap-compatible object (class c("lsa_certainty",
# "lsa_bootstrap")) so plot_forest(), as.data.frame() and the tidy
# workflow apply unchanged.

#' Compute the analytic certainty of transition probabilities
#'
#' @description
#' `certainty_lsa()` computes the analytic certainty of every transition
#' probability of a fit from [lsa()]: how precisely the data pin down each
#' transition probability, and which transitions are estimated precisely enough
#' to be called certain. Under a Dirichlet-Multinomial model the outgoing
#' transitions of each state are a multinomial draw whose probabilities have a
#' Dirichlet posterior, so each transition probability has a Beta posterior with
#' a closed-form mean, standard deviation and credible interval. A transition
#' probability estimated from counts is uncertain in proportion to how often its
#' source state was visited: a probability of 0.5 from 4 departures is far less
#' certain than the same value from 400.
#'
#' The posterior is computed with a Jeffreys prior by default and needs no
#' resampling. Each transition gets its posterior mean, standard deviation and
#' credible interval, and a posterior tail probability that decides whether it
#' is certain. With `inference = "stability"` a transition is certain when its
#' posterior stays within a band around the observed probability (by default 75%
#' to 125% of it); with `inference = "threshold"` it is certain when its
#' posterior lies above a minimum probability. A grouped fit gives one result
#' per group.
#'
#' The model treats the transitions as independent draws. When sequences
#' differ strongly from one another (for example, a mixture of latent
#' groups with long sequences) it reports more certainty than the data
#' support; use [bootstrap_lsa()], which resamples whole sequences, in
#' that case. To compare two groups' transition probabilities with the
#' same model, use [bayes_compare_lsa()].
#'
#' @param fit An `lsa` fit from [lsa()], or a grouped fit from
#'   `lsa(..., group = )`.
#' @param prior Dirichlet concentration added to every cell. Default
#'   `0.5` (the Jeffreys prior).
#' @param level_alpha Credible-interval level. Default `0.95`. The
#'   certainty decision uses `1 - level_alpha` as its threshold.
#' @param inference `"stability"` (default): a transition is certain when
#'   little posterior mass lies outside `consistency_range` times its
#'   observed probability. `"threshold"`: certain when little posterior
#'   mass lies below `edge_threshold`.
#' @param consistency_range Two multipliers of the observed probability
#'   that bound the band for `inference = "stability"`. Default
#'   `c(0.75, 1.25)`.
#' @param edge_threshold Minimum probability for
#'   `inference = "threshold"`. Default `NULL`, the 10% quantile of the
#'   non-zero transition probabilities.
#' @param x An `lsa_certainty` or `lsa_certainty_group` object, for the
#'   `print()` and `as.data.frame()` methods.
#' @param row.names,optional Not used; present for consistency with the
#'   generic.
#' @param ... Not used.
#'
#' @return An object of class `c("lsa_certainty", "lsa_bootstrap",
#'   "list")`. Read it with [transitions()], which returns a data frame
#'   with one row per transition (`from -> to` pair) and the columns
#'   `from`, `to`, `observed` (count), `prob_observed`, `prob_mean`
#'   (posterior mean), `prob_se` (posterior standard deviation),
#'   `prob_ci_low`, `prob_ci_high` (credible interval), `p_value`
#'   (posterior mass outside the band, or below the threshold), `stable`
#'   (certain), and `adj_res_observed`, `adj_res_stable` (the fit's
#'   adjusted residual and a copy of `stable`, used by [plot_forest()]).
#'   Transitions that cannot be estimated (a source state never left, or a
#'   structural zero) are `NA` and never certain.
#'   `transitions(x, significant = TRUE)` keeps the certain transitions.
#'
#'   For a grouped fit, an object of class `lsa_certainty_group`, a list
#'   with one `lsa_certainty` per group; `as.data.frame()` stacks their
#'   tables with a leading `group` column. `print()` shows the settings
#'   and the number of certain transitions and returns `x` invisibly.
#'
#' @examples
#' fit <- lsa(engagement)
#' cert <- certainty_lsa(fit)
#' cert
#' transitions(cert, significant = TRUE)
#'
#' # One table per group
#' g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'          time = "Time", group = "Achiever")
#' gcert <- certainty_lsa(g)
#' gcert
#' head(as.data.frame(gcert))
#'
#' @references
#' Johnston, L., & Jendoubi, T. (2026). How delivery mode reshapes
#' resource engagement: A Bayesian differential network analysis.
#' TNA Workshop 2026.
#'
#' @seealso [transitions()] to read the result; [plot.lsa_certainty()] and
#'   [plot_forest()] to draw it; [bootstrap_lsa()] for resampling
#'   intervals; [bayes_compare_lsa()] to compare groups under the same
#'   model.
#'
#' @export
certainty_lsa <- function(fit,
                          prior = 0.5,
                          level_alpha = 0.95,
                          inference = c("stability", "threshold"),
                          consistency_range = c(0.75, 1.25),
                          edge_threshold = NULL) {
  inference <- match.arg(inference)
  stopifnot(
    "`prior` must be a single positive number" =
      is.numeric(prior) && length(prior) == 1L && is.finite(prior) &&
      prior > 0,
    "`level_alpha` must be a single number in (0, 1)" =
      is.numeric(level_alpha) && length(level_alpha) == 1L &&
      level_alpha > 0 && level_alpha < 1,
    "`consistency_range` must be two positive numbers" =
      is.numeric(consistency_range) && length(consistency_range) == 2L &&
      all(consistency_range > 0)
  )

  if (inherits(fit, "lsa_group")) {
    out <- lapply(fit, certainty_lsa, prior = prior,
                  level_alpha = level_alpha, inference = inference,
                  consistency_range = consistency_range,
                  edge_threshold = edge_threshold)
    return(structure(out, class = c("lsa_certainty_group", "list")))
  }
  stopifnot("`fit` must be an `lsa` fit from lsa()" = inherits(fit, "lsa"))

  counts <- fit$obs
  labels <- rownames(counts)
  K <- length(labels)
  dn <- list(labels, labels)
  prob <- fit$prob

  # Beta posterior per edge (each Dirichlet row): a = count + prior,
  # b = (row concentration) - a. rowA is constant within a row.
  a <- counts + prior
  rowA <- matrix(rowSums(a), K, K, dimnames = dn)
  b <- rowA - a
  prob_mean <- a / rowA
  prob_se <- sqrt(a * b / (rowA^2 * (rowA + 1)))
  tail <- (1 - level_alpha) / 2
  ci_low  <- matrix(stats::qbeta(tail, a, b), K, K, dimnames = dn)
  ci_high <- matrix(stats::qbeta(1 - tail, a, b), K, K, dimnames = dn)

  # Decision via posterior mass (closed-form, no resampling).
  W <- prob
  if (inference == "stability") {
    cr_lo <- pmin(W * consistency_range[1], W * consistency_range[2])
    cr_hi <- pmax(W * consistency_range[1], W * consistency_range[2])
    p_val <- stats::pbeta(cr_lo, a, b) +
      (1 - stats::pbeta(cr_hi, a, b))            # mass outside the band
  } else {
    if (is.null(edge_threshold)) {
      nz <- W[is.finite(W) & W > 0]
      edge_threshold <- if (length(nz)) {
        unname(stats::quantile(nz, 0.10))
      } else 0
    }
    p_val <- stats::pbeta(edge_threshold, a, b)  # mass below the threshold
  }
  p_val <- matrix(p_val, K, K, dimnames = dn)
  alpha <- 1 - level_alpha

  # Non-estimable cells stay out: non-finite probability (zero-margin
  # rows) and declared structural zeros are NA and never "certain".
  valid <- is.finite(prob)
  sz <- fit$params$structural_zeros
  if (!is.null(sz)) valid <- valid & (sz != 0)
  prob_mean[!valid] <- NA; prob_se[!valid] <- NA
  ci_low[!valid] <- NA; ci_high[!valid] <- NA; p_val[!valid] <- NA
  stable <- valid & is.finite(W) & W > 0 & is.finite(p_val) & p_val < alpha

  grid <- expand.grid(from = labels, to = labels,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  edges <- data.frame(
    from             = grid$from,
    to               = grid$to,
    observed         = as.vector(counts),
    prob_observed    = as.vector(prob),
    prob_mean        = as.vector(prob_mean),
    prob_se          = as.vector(prob_se),
    prob_ci_low      = as.vector(ci_low),
    prob_ci_high     = as.vector(ci_high),
    p_value          = as.vector(p_val),
    stable           = as.vector(stable),
    # carried so plot_forest() can colour by direction / significance
    adj_res_observed = as.vector(fit$adj_res),
    adj_res_stable   = as.vector(stable),
    stringsAsFactors = FALSE, row.names = NULL
  )

  structure(
    list(
      edges             = edges,
      mean              = prob_mean,
      sd                = prob_se,
      ci_lower          = ci_low,
      ci_upper          = ci_high,
      p_values          = p_val,
      prior             = prior,
      level_alpha       = level_alpha,
      inference         = inference,
      consistency_range = consistency_range,
      edge_threshold    = edge_threshold,
      R                 = NA_integer_,
      level             = "analytic",
      fit               = fit
    ),
    class = c("lsa_certainty", "lsa_bootstrap", "list")
  )
}

#' @rdname certainty_lsa
#' @export
print.lsa_certainty <- function(x, ...) {
  cat("<lsa_certainty>  (analytic Dirichlet-Multinomial)\n")
  cat(sprintf("  engine:        %s\n", x$fit$method))
  cat(sprintf("  prior:         Dirichlet(%.2f)\n", x$prior))
  cat(sprintf("  CI level:      %.0f%%  |  inference: %s\n",
              100 * x$level_alpha, x$inference))
  n_stable <- sum(x$edges$stable, na.rm = TRUE)
  n_edges <- sum(is.finite(x$edges$prob_observed) &
                 x$edges$prob_observed > 0, na.rm = TRUE)
  cat(sprintf("  certain edges: %d of %d\n", n_stable, n_edges))
  invisible(x)
}

#' Plot credible intervals of transition probabilities
#'
#' @description
#' `plot()` draws the result of [certainty_lsa()] as a circular forest plot: one
#' spoke per transition, showing the observed transition probability and its
#' credible interval, so the precisely and imprecisely estimated transitions can
#' be seen side by side. The plot is drawn by [plot_forest()] with `metric =
#' "prob"`; see that page for the layout and the arguments that can be passed
#' through `...`.
#'
#' @param x An `lsa_certainty` object from [certainty_lsa()].
#' @param metric Interval to draw. Default `"prob"`, the only interval a
#'   certainty result has.
#' @param ... Passed to [plot_forest()].
#'
#' @return A `ggplot` object, drawn when printed.
#'
#' @examples
#' cert <- certainty_lsa(lsa(group_regulation))
#' plot(cert)
#'
#' @seealso [certainty_lsa()], [plot_forest()].
#'
#' @export
plot.lsa_certainty <- function(x, metric = "prob", ...) {
  plot_forest(x, metric = metric, ...)
}

#' @rdname certainty_lsa
#' @export
as.data.frame.lsa_certainty_group <- function(x, row.names = NULL,
                                              optional = FALSE, ...) {
  .grouped_df(x)
}

#' @rdname certainty_lsa
#' @export
print.lsa_certainty_group <- function(x, ...) {
  cat("<lsa_certainty_group>\n")
  cat(sprintf("  groups: %d (%s)\n", length(x),
              paste(names(x), collapse = ", ")))
  for (nm in names(x)) {
    n_stable <- sum(x[[nm]]$edges$stable, na.rm = TRUE)
    cat(sprintf("    - %-12s %d certain edges\n", paste0(nm, ":"),
                n_stable))
  }
  invisible(x)
}
