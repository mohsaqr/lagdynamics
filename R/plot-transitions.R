# Transition-network plot via cograph::splot(). One call to draw the
# directed transition network with a chosen edge weight (counts,
# probabilities, residuals, or lift), white nodes and edge labels by
# default. cograph is a Suggests-level dependency; the function errors
# informatively if it is absent.

#' Plot the Transition Network
#'
#' Draws the directed transition network of an `lsa` fit with
#' `cograph::splot()`. Pick the edge weight with `weights`; optionally
#' keep only significant edges. Nodes are white and edges are labelled
#' by default. Returns the `cograph` network invisibly.
#'
#' @param fit An `lsa` fit from [lsa()].
#' @param weights Which matrix becomes the edge weight, and how it is
#'   drawn:
#'   * `"residuals"` (default) -- a **residual network** (not a transition
#'     one): adjusted residuals coloured by sign on the TNA / Nestimate
#'     convention, **blue = more** (over-represented) solid and
#'     **red = less** (avoided) dashed with a soft halo.
#'   * `"prob"` / `"count"` -- the familiar **transition network** of
#'     Transition Network Analysis (TNA), drawn with
#'     `cograph::splot(tna_styling = TRUE)`: cograph's own TNA styling
#'     (coloured nodes, weighted directed edges) plus a donut ring per node
#'     carrying its initial-state probability, and edges labelled with the
#'     transition probability (`"prob"`) or observed count (`"count"`). For
#'     `"prob"`, edges below `0.05` are dropped by default so weak
#'     transitions do not clutter the plot (override with `edge_cutoff`).
#'   * `"lift"` -- observed / expected, drawn in a single neutral colour
#'     with magnitude carried by edge width.
#'   * `"yules_q"` -- a **signed association network**: Yule's Q on a fixed
#'     `[-1, 1]` scale, coloured by sign like the residual network (blue
#'     over-represented, red avoided) but bounded and not growing with
#'     sample size.
#' @param significant Logical. Keep only edges whose adjusted-residual
#'   p-value is below the fit's alpha; weaker cells are set to 0 (no
#'   edge). Default `FALSE`. Note that at large sample sizes almost every
#'   cell is significant, so this is a weak visual filter -- prefer `top`
#'   (effect-size pruning) to declutter a dense residual network.
#' @param top Numeric or `NULL`. Keep only the strongest edges by absolute
#'   weight (applied after `significant`); the rest are set to 0. A fraction
#'   `0 < top < 1` keeps that **proportion** of the present edges
#'   (`top = 0.5` -> the strongest half); a value `top >= 1` keeps that many
#'   edges (`top = 12` -> the 12 strongest). The legible way to thin a dense
#'   residual network: it prunes by effect size (`|adjusted residual|`)
#'   rather than by p-value. `NULL` (default) keeps every edge. Applies to
#'   every view; for the probability network it composes with the default
#'   `edge_cutoff = 0.05`.
#' @param decimals Number of decimal places for edge labels. Default `1`.
#' @param node_fill Node fill colour. Default `"white"`; the probability /
#'   count networks use a per-state palette instead unless `node_fill` is
#'   set explicitly.
#' @param edge_labels Logical (or a label vector). Show edge weights as
#'   labels. Default `TRUE`.
#' @param ... Passed to [cograph::splot()] (e.g. `node_shape`, `layout`,
#'   `edge_cutoff`, `curvature`).
#'
#' @return The `cograph_network` object, invisibly (drawn as a side
#'   effect).
#'
#' @examples
#' \dontrun{
#' fit <- lsa(group_regulation)
#' plot_transitions(fit)                                   # residual network
#' plot_transitions(fit, weights = "prob")                 # probabilities
#' plot_transitions(fit, weights = "residuals",            # residual network,
#'                  significant = TRUE)                     #   significant only
#' plot_transitions(fit, top = 12)                         # 12 strongest edges
#' plot_transitions(fit, top = 0.5)                        # strongest 50%
#' plot_transitions(fit, decimals = 2)                     # 2-dp edge labels
#' plot_transitions(fit, node_shape = "square")            # splot passthrough
#' }
#'
#' @seealso [plot.lsa()] (heatmap), [transitions()],
#'   [transition_probabilities()]
#'
#' @export
plot_transitions <- function(fit,
                             weights = c("residuals", "count", "prob",
                                          "lift", "yules_q"),
                             significant = FALSE,
                             top = NULL,
                             decimals = 1,
                             node_fill = "white",
                             edge_labels = TRUE,
                             ...) {
  stopifnot(inherits(fit, "lsa"))
  weights <- match.arg(weights)
  if (!is.null(top)) {
    stopifnot(length(top) == 1L, is.numeric(top), is.finite(top), top > 0)
  }
  stopifnot(length(decimals) == 1L, is.numeric(decimals),
            is.finite(decimals), decimals >= 0)

  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop("Package 'cograph' is required for plot_transitions(). ",
         "Install with install.packages('cograph').", call. = FALSE)
  }
  # Keep signed residuals (negatives are meaningful) when plotting them.
  wkey <- if (weights == "residuals") "adj_res" else weights
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
  if (weights %in% c("prob", "count")) {
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
    if (weights == "prob") defaults$edge_cutoff <- 0.05
    return(do.call(cograph::splot, utils::modifyList(defaults, list(...))))
  }

  # splot() colours edges BY SIGN for the residual / Yule's Q networks.
  signed <- weights %in% c("residuals", "yules_q")
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
