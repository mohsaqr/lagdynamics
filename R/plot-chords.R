# Circular chord diagram for an lsa fit. The transition NETWORK is
# delegated to cograph::splot() (plot_transitions()); the CHORD diagram
# is delegated to cograph::plot_chord() here, exactly the same way -- so
# the heavy circular-ribbon drawing lives in cograph, and lagdynamics only
# decides what the ribbons mean.
#
# A chord ribbon's WIDTH is the transition frequency (or probability) and
# its FILL COLOUR is the adjusted residual: warm = over-represented, cool
# = avoided. A second fit (`compare`) fills each ribbon by the DIFFERENCE
# in the chosen colour metric between the two fits. Because an `lsa` fit
# carries the `cograph_network` class, cograph consumes it directly; we
# only build the per-ribbon colour vector, aligned to cograph's row-major
# nonzero edge order.

# Diverging (signed) and sequential (non-negative) anchors, identical to
# the ones plot.lsa() uses, so chords and heatmaps read as one family.
.lsa_div_low  <- "#2166AC"; .lsa_div_mid <- "#F7F7F7"; .lsa_div_high <- "#B2182B"
.lsa_seq_low  <- "#DEEBF7"; .lsa_seq_high <- "#08519C"

# Map a numeric vector to opaque ribbon colours (cograph applies the
# transparency via `chord_alpha`, so these stay fully opaque to avoid
# double-applying alpha). Diverging metrics (residuals, differences) are
# centred at 0 and symmetric; sequential metrics span their own range.
# Non-finite values render grey.
.lsa_chord_colors <- function(vals, diverging) {
  fin <- is.finite(vals)
  if (!any(fin)) return(rep("grey85", length(vals)))   # nothing to scale
  t01 <- rep(NA_real_, length(vals))
  if (diverging) {
    lim <- max(abs(vals[fin]), na.rm = TRUE)
    if (!is.finite(lim) || lim <= 0) lim <- 1
    anchors <- c(.lsa_div_low, .lsa_div_mid, .lsa_div_high)
    t01[fin] <- (pmin(pmax(vals[fin], -lim), lim) + lim) / (2 * lim)
  } else {
    rng <- range(vals[fin])
    anchors <- c(.lsa_seq_low, .lsa_seq_high)
    t01[fin] <- if (diff(rng) > 0) (vals[fin] - rng[1L]) / diff(rng)
                else rep(0.5, sum(fin))
  }
  ramp <- grDevices::colorRamp(anchors)
  out <- rep("grey85", length(vals))
  if (any(fin)) {
    m <- ramp(t01[fin])
    out[fin] <- grDevices::rgb(m[, 1L], m[, 2L], m[, 3L], maxColorValue = 255)
  }
  out
}

#' Plot a chord diagram of transitions
#'
#' @description
#' `plot_chords()` draws the transitions of a fit as a chord diagram with
#' \pkg{cograph}: the states are arcs around a circle and each transition is a
#' ribbon from its source arc to its target arc. The width of a ribbon shows how
#' much of the traffic the transition carries, and its colour can show whether
#' the transition is more or less frequent than independence predicts, so the
#' diagram shows where the transitions out of each state go, and which of the
#' frequent routes are also over-represented.
#'
#' Ribbon width is the count (default) or the transition probability; ribbon
#' colour is the adjusted residual on a diverging scale (red over-represented,
#' blue avoided), or the lift, probability or count on a sequential scale. With
#' `compare`, a second fit on the same states, each ribbon is coloured by the
#' difference between the two fits in the colour quantity, which shows where two
#' groups differ.
#'
#' The same display is `plot(fit, type = "chord")`. For the full table use
#' the residual heatmap of [plot.lsa()]; for a network layout,
#' [plot_transitions()]; for a test of group differences,
#' [compare_lsa()].
#'
#' @param fit An `lsa` fit from [lsa()].
#' @param compare Optional second `lsa` fit with the same states; ribbons
#'   are then coloured by `fit` minus `compare`. Default `NULL`.
#' @param width Ribbon width: `"count"` (default) or `"prob"`.
#' @param color Ribbon colour: `"residuals"` (default), `"lift"`, `"prob"`
#'   or `"count"`.
#' @param significant Draw only transitions whose residual p-value is
#'   below the fit's `alpha` (default `FALSE`); ignored with `compare`.
#' @param self_loops Draw self-transitions (default `TRUE`).
#' @param alpha Ribbon opacity (default `0.6`).
#' @param ... Passed to [cograph::plot_chord()], for example `ticks`,
#'   `segment_width`, `label_size` or `title`.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and
#' `lagdynamics_incompatible_fits` when `fit` and `compare` have different
#' states.
#'
#' @return The value of [cograph::plot_chord()] (the drawn segments and
#'   chords), invisibly; the plot is drawn on the current device.
#'
#' @examples
#' fit <- lsa(engagement)
#' plot_chords(fit)
#' plot_chords(fit, width = "prob", self_loops = FALSE)
#'
#' # Ribbons coloured by the difference between two groups
#' high <- lsa(subset(group_regulation_long, Achiever == "High"),
#'             actor = "Actor", action = "Action", time = "Time")
#' low <- lsa(subset(group_regulation_long, Achiever == "Low"),
#'            actor = "Actor", action = "Action", time = "Time")
#' plot_chords(high, compare = low)
#'
#' @seealso [plot.lsa()], [plot_transitions()], [plot_polar()],
#'   [compare_lsa()]
#'
#' @export
plot_chords <- function(fit, compare = NULL,
                        width = c("count", "prob"),
                        color = c("residuals", "lift", "prob", "count"),
                        significant = FALSE, self_loops = TRUE,
                        alpha = 0.6, ...) {
  stopifnot("`fit` must be an `lsa` fit from lsa()" = inherits(fit, "lsa"))
  width <- match.arg(width)
  color <- match.arg(color)
  if (!requireNamespace("cograph", quietly = TRUE)) {
    .lsa_stop(class = "missing_package", "Package 'cograph' is required for plot_chords(). ",
         "Install with install.packages('cograph').", call. = FALSE)
  }

  # Ribbon WIDTH: a non-negative weight matrix (frequency or probability).
  W <- .lsa_weight_matrix(fit, width, positive_residuals_only = FALSE)

  # Ribbon COLOUR: the chosen metric, or its between-fit difference.
  ckey <- if (color == "residuals") "adj_res" else color
  Cv <- .lsa_weight_matrix(fit, ckey, positive_residuals_only = FALSE)
  diverging <- color == "residuals"
  if (!is.null(compare)) {
    stopifnot("`compare` must be an `lsa` fit from lsa()" =
                inherits(compare, "lsa"))
    if (!setequal(rownames(fit$obs), rownames(compare$obs))) {
      .lsa_stop(class = "incompatible_fits", "`fit` and `compare` must have the same states.", call. = FALSE)
    }
    ord <- match(rownames(W), rownames(compare$obs))
    Cv <- Cv - .lsa_weight_matrix(compare, ckey,
                                  positive_residuals_only = FALSE)[ord, ord]
    diverging <- TRUE
    significant <- FALSE
  } else if (isTRUE(significant)) {
    keep <- is.finite(fit$p) & fit$p < fit$params$alpha
    W[!keep] <- 0                                  # drop non-significant ribbons
  }
  # cograph::plot_chord always preserves self-loops regardless of its
  # `self_loop` flag, so honour `self_loops = FALSE` by zeroing the
  # diagonal here (which also drops it from the nonzero colour set below).
  if (!isTRUE(self_loops)) diag(W) <- 0

  # cograph draws chords in row-major (i, j) order over nonzero widths;
  # build the colour vector in that exact order (transpose -> row-major).
  nz <- abs(as.vector(t(W))) > 0
  cols <- .lsa_chord_colors(as.vector(t(Cv))[nz], diverging)

  cograph::plot_chord(W, chord_color_by = cols, self_loop = self_loops,
                      chord_alpha = alpha, ...)
}
