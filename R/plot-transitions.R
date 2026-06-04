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
#' @param weights Which matrix becomes the edge weight: `"residuals"`
#'   (default, adjusted residuals -- edges coloured warm = over-
#'   represented, cool = avoided, matching the heatmap), `"count"`
#'   (observed counts), `"prob"` (row-conditional probabilities), or
#'   `"lift"` (observed / expected). The non-residual weights are all
#'   non-negative, so they are drawn in a single neutral colour with the
#'   magnitude carried by edge width (rather than by sign).
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
  # the lagseq palette (over-represented = warm red, avoided = cool blue,
  # matching the heatmap). Counts/probs/lift are all non-negative -- sign
  # colouring would paint every edge one colour -- so draw them in a
  # single neutral blue and let width carry the magnitude.
  signed <- weights == "residuals"
  epos <- if (signed) .lsa_div_high else "#4A6FA5"   # over / neutral
  eneg <- if (signed) .lsa_div_low  else "#4A6FA5"   # avoided / neutral
  # Distinctive default node: white fill with a tight empty ring. Any of
  # these (edge colours, donut_*, node_shape, layout, ...) can be
  # overridden via `...`.
  K <- nrow(m)
  defaults <- list(
    x = m, node_fill = node_fill, edge_labels = edge_labels,
    edge_positive_color = epos, edge_negative_color = eneg,
    donut_fill = rep(NA, K), donut_empty = TRUE,
    donut_inner_ratio = 0.85, donut_border_color = "steelblue",
    donut_border_width = 1.1
  )
  do.call(cograph::splot, utils::modifyList(defaults, list(...)))
}
