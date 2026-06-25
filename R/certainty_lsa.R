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

#' Analytic Certainty of Transition Edges (Dirichlet-Multinomial)
#'
#' Closed-form Bayesian alternative to [bootstrap_lsa()] for the
#' transition-probability edges of an `lsa` fit. Each state's outgoing
#' transitions are modelled as Dirichlet-Multinomial: with a Jeffreys
#' prior the posterior for a row is `Dirichlet(count + prior)`, so each
#' edge probability is marginally `Beta(a, b)` and its posterior mean,
#' standard deviation, credible interval and stability decision are
#' available analytically. No resampling, so it runs in microseconds.
#'
#' The result carries class `c("lsa_certainty", "lsa_bootstrap")` and an
#' `edges` table with the columns [plot_forest()] and
#' [as.data.frame()] expect, so it is a drop-in for a bootstrap result
#' (use `metric = "prob"`).
#'
#' **Certainty vs bootstrap.** Both answer "how precisely is this edge
#' pinned down?". They agree on homogeneous data. The Dirichlet posterior
#' treats transitions as independent, so on strongly heterogeneous data (a
#' mixture of latent classes with long sequences) it reports *more*
#' certainty than the sequence bootstrap -- prefer [bootstrap_lsa()] then.
#'
#' @param fit An `lsa` fit from [lsa()], or an `lsa_group`.
#' @param prior Numeric > 0. Dirichlet prior concentration added to every
#'   cell. Default `0.5` (the Jeffreys prior).
#' @param level_alpha Numeric in (0, 1). Credible-interval level. Default
#'   `0.95` (a 95% interval), matching [bootstrap_lsa()].
#' @param inference `"stability"` (default) flags an edge whose posterior
#'   keeps it within a multiplicative `consistency_range` of its observed
#'   probability; `"threshold"` flags an edge whose posterior mass lies
#'   above `edge_threshold`.
#' @param consistency_range Length-2 multiplicative bounds for stability
#'   inference. Default `c(0.75, 1.25)`.
#' @param edge_threshold Numeric or `NULL`. Fixed threshold for
#'   `inference = "threshold"`; `NULL` uses the 0.10 quantile of
#'   non-zero edge probabilities.
#'
#' @return An object of class `c("lsa_certainty", "lsa_bootstrap", "list")`
#'   with an `edges` data frame (`from`, `to`, `prob_observed`,
#'   `prob_mean`, `prob_se`, `prob_ci_low`, `prob_ci_high`, `p_value`,
#'   `stable`, plus `adj_res_observed`/`adj_res_stable` for plotting), the
#'   posterior matrices (`mean`, `sd`, `ci_lower`, `ci_upper`), and call
#'   metadata (`prior`, `level_alpha`, `inference`, ...). For an
#'   `lsa_group`, a named list of these (class `lsa_certainty_group`).
#'
#' @examples
#' \donttest{
#' fit  <- lsa(engagement)
#' cert <- certainty_lsa(fit)
#' cert
#' head(as.data.frame(cert))
#' }
#'
#' @references
#' Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
#' Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.
#'
#' @seealso [bootstrap_lsa()], [stability_lsa()], [plot_forest()]
#'
#' @export
certainty_lsa <- function(fit,
                          prior = 0.5,
                          level_alpha = 0.95,
                          inference = c("stability", "threshold"),
                          consistency_range = c(0.75, 1.25),
                          edge_threshold = NULL) {
  inference <- match.arg(inference)
  stopifnot(is.numeric(prior), length(prior) == 1L, is.finite(prior),
            prior > 0,
            is.numeric(level_alpha), length(level_alpha) == 1L,
            level_alpha > 0, level_alpha < 1,
            is.numeric(consistency_range), length(consistency_range) == 2L,
            all(consistency_range > 0))

  if (inherits(fit, "lsa_group")) {
    out <- lapply(fit, certainty_lsa, prior = prior,
                  level_alpha = level_alpha, inference = inference,
                  consistency_range = consistency_range,
                  edge_threshold = edge_threshold)
    return(structure(out, class = c("lsa_certainty_group", "list")))
  }
  stopifnot(inherits(fit, "lsa"))

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

#' Plot an Analytic-Certainty Result
#'
#' Circular forest of the per-edge transition-probability credible
#' intervals from [certainty_lsa()] (delegates to [plot_forest()] with
#' `metric = "prob"`).
#'
#' @param x An `lsa_certainty` object.
#' @param metric Which credible interval to draw. Default `"prob"`.
#' @param ... Passed to [plot_forest()].
#' @return A `ggplot` object (drawn when printed). Needs `ggplot2`.
#' @export
plot.lsa_certainty <- function(x, metric = "prob", ...) {
  plot_forest(x, metric = metric, ...)
}

#' @export
as.data.frame.lsa_certainty_group <- function(x, row.names = NULL,
                                              optional = FALSE, ...) {
  .grouped_df(x)
}

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
