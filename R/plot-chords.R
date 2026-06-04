# Circular chord diagram for an lsa fit. The transition NETWORK is
# delegated to cograph::splot() (plot_transitions()); the CHORD diagram
# is delegated to cograph::plot_chord() here, exactly the same way -- so
# the heavy circular-ribbon drawing lives in cograph, and lagseq only
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

#' Circular (Chord) Diagram of an LSA Fit
#'
#' Draws the transition structure as a chord diagram via
#' [cograph::plot_chord()]: states are arcs on an outer ring and each
#' transition is a curved ribbon whose **width** is its frequency (or
#' probability) and whose **fill colour** is its adjusted residual
#' (warm = over-represented, cool = avoided). Supply a second fit as
#' `compare` to fill each ribbon by the *difference* in the colour metric
#' between the two fits.
#'
#' This is the circular companion to the [plot.lsa()] heatmap and the
#' [plot_transitions()] network. Like them it delegates the drawing to
#' `cograph`; it needs the `cograph` package installed.
#'
#' @param fit An `lsa` fit from [lsa()].
#' @param compare Optional second `lsa` fit. When supplied, ribbon colour
#'   is `colour(fit) - colour(compare)` (a signed difference on the
#'   diverging scale). The two fits must share the same states. Default
#'   `NULL`.
#' @param width Which non-negative quantity sets ribbon width: `"count"`
#'   (default, transition frequency) or `"prob"` (row-conditional
#'   probability).
#' @param color Which quantity fills the ribbons: `"residuals"` (default,
#'   signed adjusted residual, diverging), `"lift"`, `"prob"`, or
#'   `"count"`. Non-residual metrics use a sequential scale unless
#'   `compare` makes them a signed difference.
#' @param significant Logical. Keep only significant transitions (drops
#'   the others' ribbons). Ignored when `compare` is set. Default `FALSE`.
#' @param self_loops Logical. Draw self-transition ribbons. Default
#'   `TRUE`.
#' @param alpha Ribbon fill opacity. Default `0.6`.
#' @param ... Passed to [cograph::plot_chord()] (e.g. `ticks`,
#'   `segment_width`, `label_size`, `title`).
#'
#' @return Invisibly, the list returned by [cograph::plot_chord()]
#'   (`segments` and `chords` data frames). Drawn as a side effect.
#'
#' @examples
#' \dontrun{
#' fit <- lsa(group_regulation)
#' plot_chords(fit)                          # ribbons filled by residual
#' plot_chords(fit, width = "prob")          # width = probability
#' plot_chords(fit, significant = TRUE, ticks = TRUE)
#'
#' # Compare two groups: ribbon colour = difference in residuals.
#' g <- lsa(group_regulation,
#'          group = rep(c("A", "B"), length.out = nrow(group_regulation)))
#' plot_chords(g$A, compare = g$B)
#' }
#'
#' @seealso [plot.lsa()] (heatmap), [plot_transitions()] (network),
#'   [transitions()]
#'
#' @export
plot_chords <- function(fit, compare = NULL,
                        width = c("count", "prob"),
                        color = c("residuals", "lift", "prob", "count"),
                        significant = FALSE, self_loops = TRUE,
                        alpha = 0.6, ...) {
  stopifnot(inherits(fit, "lsa"))
  width <- match.arg(width)
  color <- match.arg(color)
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop("Package 'cograph' is required for plot_chords(). ",
         "Install with install.packages('cograph').", call. = FALSE)
  }

  # Ribbon WIDTH: a non-negative weight matrix (frequency or probability).
  W <- .lsa_weight_matrix(fit, width, positive_residuals_only = FALSE)

  # Ribbon COLOUR: the chosen metric, or its between-fit difference.
  ckey <- if (color == "residuals") "adj_res" else color
  Cv <- .lsa_weight_matrix(fit, ckey, positive_residuals_only = FALSE)
  diverging <- color == "residuals"
  if (!is.null(compare)) {
    stopifnot(inherits(compare, "lsa"))
    if (!setequal(rownames(fit$obs), rownames(compare$obs))) {
      stop("`fit` and `compare` must have the same states.", call. = FALSE)
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
