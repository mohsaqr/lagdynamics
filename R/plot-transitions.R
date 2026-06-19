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
#'   * `"prob"` / `"count"` -- the familiar **transition network**. This is
#'     a TNA model: the fit is converted to a `tna` object on the fly with
#'     [lsa_to_tna()] and rendered by tna's own plot method (coloured
#'     nodes, initial-probability arcs, weighted directed edges). Needs the
#'     `tna` package; `...` is forwarded to tna's plot.
#'   * `"lift"` -- observed / expected, drawn in a single neutral colour
#'     with magnitude carried by edge width.
#' @param significant Logical. Keep only edges whose adjusted-residual
#'   p-value is below the fit's alpha; weaker cells are set to 0 (no
#'   edge). Default `FALSE`.
#' @param node_fill Node fill colour. Default `"white"`.
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
#' plot_transitions(fit, node_shape = "square")            # splot passthrough
#' }
#'
#' @seealso [plot.lsa()] (heatmap), [transitions()], [lsa_to_tna()]
#'
#' @export
plot_transitions <- function(fit,
                             weights = c("residuals", "count", "prob",
                                          "lift"),
                             significant = FALSE,
                             node_fill = "white",
                             edge_labels = TRUE,
                             ...) {
  stopifnot(inherits(fit, "lsa"))
  weights <- match.arg(weights)

  # A probability- or count-weighted network IS a transition (TNA) model,
  # so build the `tna` object on the fly and let tna render it with its
  # own styling (coloured nodes, initial-probability arcs, weighted
  # directed edges) -- rather than lagseq's residual-network styling.
  if (weights %in% c("prob", "count")) {
    if (!requireNamespace("tna", quietly = TRUE)) {
      stop("Package 'tna' is required for the transition network ",
           "(weights = '", weights, "'). Install with ",
           "install.packages('tna'), or use weights = 'residuals' for ",
           "the residual network.", call. = FALSE)
    }
    return(plot(lsa_to_tna(fit, weights = weights), ...))
  }

  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop("Package 'cograph' is required for plot_transitions(). ",
         "Install with install.packages('cograph').", call. = FALSE)
  }
  # Keep signed residuals (negatives are meaningful) when plotting them.
  wkey <- if (weights == "residuals") "adj_res" else weights
  m <- .lsa_weight_matrix(fit, wkey, positive_residuals_only = FALSE)

  if (isTRUE(significant)) {
    keep <- is.finite(fit$p) & fit$p < fit$params$alpha
    m[!keep] <- 0
  }
  # splot() colours edges BY SIGN. Residuals are signed, so map them to
  # the shared TNA / Nestimate / dynalytics convention: BLUE = more
  # (over-represented, positive residual), RED = less (avoided, negative).
  # Counts/probs/lift are all non-negative -- sign colouring would paint
  # every edge one colour -- so draw them in a single neutral blue and let
  # width carry the magnitude.
  signed <- weights == "residuals"
  # .cmp_high = blue (over-represented), .cmp_low = red (avoided): the same
  # convention as the comparison plots (see R/plot-comparison.R).
  epos <- if (signed) .cmp_high else "#4A6FA5"    # over-represented (blue)
  eneg <- if (signed) .cmp_low  else "#4A6FA5"    # avoided (red) / neutral
  # Plain node: filled circle with a single border (no donut ring). Any of
  # these (edge colours, node_*, node_shape, layout, ...) can be overridden
  # via `...`.
  K <- nrow(m)
  defaults <- list(
    x = m, node_fill = node_fill, edge_labels = edge_labels,
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
