# Transition-network plot via cograph::splot(). One call to draw the
# directed transition network with a chosen edge weight (counts,
# probabilities, residuals, or lift), white nodes and edge labels by
# default. cograph is a Suggests-level dependency; the function errors
# informatively if it is absent.

#' Plot a lag transition network
#'
#' @description
#' `plot_transitions()` draws a lag transition network with \pkg{cograph}: the
#' states are nodes and the transitions between them directed edges. Which
#' quantity weights the edges decides what the network shows. Weighted by
#' adjusted residuals it is a residual network: edges show where a transition is
#' more frequent (over-represented) or less frequent (avoided) than independence
#' predicts, the reading of lag sequential analysis. Weighted by transition
#' probabilities it is the network of Transition Network Analysis, which shows
#' where the process goes next, whether or not that is more than chance.
#'
#' With `weights = "residuals"` (default), over-represented transitions are
#' solid blue edges and avoided transitions dashed red edges; `"yules_q"` draws
#' the same signed network on the bounded Yule's Q scale, and `"lift"` the ratio
#' of observed to expected counts. With `"prob"` (or its aliases `"tna"` and
#' `"relative"`) or `"count"`, the network is drawn in the style of Transition
#' Network Analysis, each node ringed by its initial-state probability;
#' probability edges below 0.05 are hidden unless `edge_cutoff` is passed. Dense
#' networks are thinned with `significant` or, better at large sample sizes
#' where nearly every cell is significant, with `top`.
#'
#' The same display is `plot(fit, type = "network")`. For the full table
#' use the residual heatmap of [plot.lsa()]; for the per-source breakdown,
#' [plot_chords()] or [plot_polar()].
#'
#' @param fit An `lsa` fit from [lsa()].
#' @param weights Edge weight: `"residuals"` (default), `"tna"`,
#'   `"relative"`, `"prob"` (the three are the transition probability),
#'   `"count"`, `"lift"` or `"yules_q"`.
#' @param significant Keep only transitions whose residual p-value is below
#'   the fit's `alpha` (default `FALSE`).
#' @param top Keep only the strongest edges by absolute weight: a fraction
#'   below 1 keeps that share of the edges, a number of 1 or more keeps
#'   that many. Default `NULL` (all edges). Applied after `significant`.
#' @param decimals Decimal places of the edge labels (default `1`).
#' @param node_fill Node colour (default `"white"`). The probability and
#'   count networks colour nodes by state unless `node_fill` is set.
#' @param edge_labels Show the edge weights as labels (default `TRUE`).
#' @param ... Passed to [cograph::splot()], for example `layout`,
#'   `node_shape`, `edge_cutoff` or `curvature`.
#'
#' @return The network drawn by [cograph::splot()], invisibly; the plot is
#'   drawn on the current device.
#'
#' @examples
#' fit <- lsa(engagement)
#' plot_transitions(fit)
#' plot_transitions(fit, weights = "prob")
#'
#' fit9 <- lsa(group_regulation)
#' plot_transitions(fit9, top = 12)
#' plot_transitions(fit9, weights = "yules_q", significant = TRUE)
#'
#' @seealso [plot.lsa()], [plot_chords()], [transitions()],
#'   [transition_probabilities()]
#'
#' @export
plot_transitions <- function(fit,
                             weights = c("residuals", "tna", "relative",
                                          "count", "prob", "lift", "yules_q"),
                             significant = FALSE,
                             top = NULL,
                             decimals = 1,
                             node_fill = "white",
                             edge_labels = TRUE,
                             ...) {
  stopifnot("`fit` must be an `lsa` fit from lsa()" = inherits(fit, "lsa"))
  weights <- match.arg(weights)
  if (!is.null(top)) {
    stopifnot(
      "`top` must be a single positive number" =
        length(top) == 1L && is.numeric(top) && is.finite(top) && top > 0
    )
  }
  stopifnot(
    "`decimals` must be a single number >= 0" =
      length(decimals) == 1L && is.numeric(decimals) && is.finite(decimals) &&
      decimals >= 0
  )

  if (!requireNamespace("cograph", quietly = TRUE)) {
    .lsa_stop(class = "missing_package", "Package 'cograph' is required for plot_transitions(). ",
         "Install with install.packages('cograph').", call. = FALSE)
  }
  # Keep signed residuals (negatives are meaningful) when plotting them.
  weight_key <- if (weights %in% c("tna", "relative")) "prob" else weights
  wkey <- if (weight_key == "residuals") "adj_res" else weight_key
  m <- .lsa_weight_matrix(fit, wkey)

  if (isTRUE(significant)) {
    keep <- is.finite(fit$p) & fit$p < fit$params$alpha
    m[!keep] <- 0
  }
  # Effect-size pruning: keep the strongest edges by absolute weight, zero
  # the rest. Unlike `significant`, this stays a real filter at large N
  # (where every cell is significant) -- it declutters by magnitude. A
  # fractional `top` (0 < top < 1) keeps that proportion of the present
  # edges (top = 0.5 -> the strongest half); top >= 1 keeps that many edges.
  if (!is.null(top)) {
    mag <- abs(m)
    mag[!is.finite(mag)] <- 0
    n_edge <- sum(mag > 0)
    n_keep <- if (top < 1) ceiling(top * n_edge) else min(as.integer(top), n_edge)
    if (n_edge > n_keep) {
      kept <- order(mag, decreasing = TRUE)[seq_len(n_keep)]
      drop <- setdiff(which(mag > 0), kept)
      m[drop] <- 0
    }
  }
  # splot() colours edges BY SIGN. Residuals and Yule's Q are signed, so
  # map them to the shared TNA / dynalytics convention: BLUE = more
  # (over-represented, positive), RED = less (avoided, negative). Lift is
  # non-negative -- sign colouring would paint every edge one colour -- so
  # it is drawn in a single neutral blue with width carrying the magnitude.
  # Probabilities / counts branch off first, into the TNA-styled draw.
  # Probability / count networks are drawn in the Transition Network
  # Analysis (TNA) style using cograph's own `tna_styling` preset -- the
  # same styling tna applies, since tna renders through cograph. The
  # per-node donut ring carries each state's initial-state probability, and
  # edges are labelled with the transition probability / count.
  if (weight_key %in% c("prob", "count")) {
    defaults <- list(
      x = m, tna_styling = TRUE, directed = TRUE,
      edge_labels = edge_labels, weight_digits = decimals
    )
    # Respect an explicit node_fill; otherwise let tna_styling colour nodes.
    if (!identical(node_fill, "white")) defaults$node_fill <- node_fill
    # Initial-state probabilities become the node donut ring, matching the
    # tna plot. Omitted when the fit came from a bare transition matrix.
    if (!is.null(fit$inits)) {
      iv <- as.numeric(fit$inits[rownames(m)])
      iv[!is.finite(iv)] <- 0
      defaults$donut_values <- iv
    }
    # For probabilities, drop weak edges (< 0.05) by default so the network
    # stays legible, as tna does; the caller can override via edge_cutoff.
    if (weight_key == "prob") defaults$edge_cutoff <- 0.05
    return(do.call(cograph::splot, utils::modifyList(defaults, list(...))))
  }

  # splot() colours edges BY SIGN for the residual / Yule's Q networks.
  signed <- weight_key %in% c("residuals", "yules_q")
  # .cmp_high = blue (over-represented), .cmp_low = red (avoided): the same
  # convention as the comparison plots (see R/plot-comparison.R).
  epos <- if (signed) .cmp_high else "#4A6FA5"    # over-represented (blue)
  eneg <- if (signed) .cmp_low  else "#4A6FA5"    # avoided (red) / neutral
  # Plain node: filled circle with a single border (no donut ring). Any of
  # these (edge colours, node_*, node_shape, layout, ...) can be overridden
  # via `...`.
  defaults <- list(
    x = m, node_fill = node_fill, edge_labels = edge_labels,
    weight_digits = decimals,
    edge_positive_color = epos, edge_negative_color = eneg,
    node_border_color = "steelblue", node_border_width = 1.1
  )
  if (signed) {
    # Avoided edges (less than chance) are drawn DASHED and wrapped in a
    # soft red halo, so they read distinctly from the solid blue
    # over-represented edges -- the cograph bootstrap idiom (solid vs
    # dashed line + a translucent band). Over-represented edges stay solid
    # with no halo.
    neg <- is.finite(m) & m < 0
    mag <- abs(m); mag[!neg] <- 0
    denom <- max(mag, na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) denom <- 1
    defaults$edge_style <- ifelse(neg, 2L, 1L)        # dashed / solid (lty)
    # Halo only on avoided edges, normalised to [0, 1] so its size is
    # dataset-independent (residual magnitudes vary widely between fits).
    defaults$edge_ci <- mag / denom
    defaults$edge_ci_color <- eneg
    defaults$edge_ci_alpha <- 0.18
    defaults$edge_ci_scale <- 2
  }
  do.call(cograph::splot, utils::modifyList(defaults, list(...)))
}
