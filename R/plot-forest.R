# Circular ("radial") forest of an lsa bootstrap. Each edge is a spoke
# placed clockwise around a ring; the spoke spans the metric's bootstrap
# confidence interval and a square marks the observed estimate. A dashed
# reference ring marks the null (0 for residuals). Edges whose residual
# is significant across resamples are coloured by direction (warm over-
# represented, cool avoided); the rest are grey. This is the circular
# read of bootstrap_lsa() -- "the bootstrap circles".
#
# Owned here in ggplot2 (sibling of plot.lsa()); it consumes the
# lsa_bootstrap object returned by bootstrap_lsa().

#' Circular Bootstrap Forest of an LSA Fit
#'
#' Draws a radial forest of an [bootstrap_lsa()] result: each transition
#' is a spoke around a ring, spanning its bootstrap confidence interval,
#' with a square at the observed estimate and a dashed reference ring at
#' the null. Spokes whose adjusted residual is significant across
#' resamples are coloured by direction (warm = over-represented, cool =
#' avoided); non-significant ones are grey. Needs `ggplot2`.
#'
#' @param boot An `lsa_bootstrap` object from [bootstrap_lsa()].
#' @param metric Which bootstrapped quantity to plot: `"residuals"`
#'   (default, adjusted residual), `"count"`, `"prob"`, or `"yules_q"`.
#' @param n_top Optional integer: keep only the `n_top` edges with the
#'   largest absolute estimate (the rest are dropped). Default `NULL`
#'   (all edges).
#' @param show_nonsig Logical. Draw non-significant edges (grey). Default
#'   `TRUE`; set `FALSE` to keep only significant transitions.
#' @param label_size Edge-label text size. Default `2.6`.
#'
#' @return A `ggplot` object (drawn when printed).
#'
#' @examples
#' \dontrun{
#' fit <- lsa(group_regulation)
#' b <- bootstrap_lsa(fit, R = 500)
#' plot_forest(b)                       # residual CIs, circular
#' plot_forest(b, metric = "prob")      # probability CIs
#' plot_forest(b, show_nonsig = FALSE)  # significant transitions only
#' }
#'
#' @seealso [bootstrap_lsa()], [plot.lsa()] (heatmap),
#'   [plot_polar()] (sunburst)
#'
#' @export
plot_forest <- function(boot, metric = c("residuals", "count", "prob",
                                         "yules_q"),
                        n_top = NULL, show_nonsig = TRUE, label_size = 2.6) {
  stopifnot(inherits(boot, "lsa_bootstrap"))
  metric <- match.arg(metric)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_forest(). ",
         "Install with install.packages('ggplot2').", call. = FALSE)
  }
  mkey <- switch(metric, residuals = "adj_res", metric)
  e <- boot$edges
  # The observed count lives in `observed`; every other metric uses
  # `<metric>_observed`.
  ecol <- paste0(mkey, "_observed")
  if (is.null(e[[ecol]])) ecol <- "observed"
  est <- e[[ecol]]
  lo  <- e[[paste0(mkey, "_ci_low")]]
  hi  <- e[[paste0(mkey, "_ci_high")]]
  if (is.null(est) || is.null(lo) || is.null(hi)) {
    stop("Bootstrap has no CI for '", metric, "'.", call. = FALSE)
  }

  # An edge is "significant" if its residual is stable across resamples;
  # its direction is the sign of the observed residual.
  sig  <- .lsa_true(e$adj_res_stable)
  over <- is.finite(e$adj_res_observed) & e$adj_res_observed > 0
  status <- ifelse(!sig, "n.s.",
            ifelse(over, "over-represented", "avoided"))

  df <- data.frame(edge = paste0(e$from, " \u2192 ", e$to),
                   est = est, lo = lo, hi = hi, sig = sig,
                   status = status, stringsAsFactors = FALSE)
  df <- df[is.finite(df$est) & is.finite(df$lo) & is.finite(df$hi), ,
           drop = FALSE]
  if (!show_nonsig) df <- df[df$sig, , drop = FALSE]
  if (!nrow(df)) stop("No edges to display.", call. = FALSE)
  df <- df[order(df$edge), , drop = FALSE]
  if (!is.null(n_top)) {
    stopifnot(is.numeric(n_top), length(n_top) == 1L, is.finite(n_top),
              n_top >= 1)
    n_top <- as.integer(n_top)
    keep <- order(abs(df$est), decreasing = TRUE)[seq_len(min(n_top, nrow(df)))]
    df <- df[sort(keep), , drop = FALSE]
  }
  df$status <- factor(df$status,
                      levels = c("over-represented", "avoided", "n.s."))

  n <- nrow(df)
  angles <- seq(pi / 2, pi / 2 - 2 * pi, length.out = n + 1)[seq_len(n)]
  df$angle <- angles

  # Radius scale: span the CI range, with the null (0) inside the ring.
  r_inner <- 0.55
  v_min <- min(c(df$lo, 0))
  v_max <- max(c(df$hi, 0)) * 1.02
  to_r <- function(v) {
    r_inner + pmin(pmax((v - v_min) / (v_max - v_min), 0), 1) * (1 - r_inner)
  }
  df$x_est <- to_r(df$est) * cos(angles); df$y_est <- to_r(df$est) * sin(angles)
  df$x_lo  <- to_r(df$lo)  * cos(angles); df$y_lo  <- to_r(df$lo)  * sin(angles)
  df$x_hi  <- to_r(df$hi)  * cos(angles); df$y_hi  <- to_r(df$hi)  * sin(angles)
  df$x_in  <- r_inner * cos(angles);      df$y_in  <- r_inner * sin(angles)
  df$x_out <- cos(angles);                df$y_out <- sin(angles)

  label_r <- 1.04
  flip <- cos(angles) < 0
  df$x_lab <- label_r * cos(angles); df$y_lab <- label_r * sin(angles)
  df$text_angle <- ifelse(flip, angles * 180 / pi + 180, angles * 180 / pi)
  df$hjust <- ifelse(flip, 1, 0)

  theta <- seq(0, 2 * pi, length.out = 300)
  ring_in  <- data.frame(x = r_inner * cos(theta), y = r_inner * sin(theta))
  ring_out <- data.frame(x = cos(theta), y = sin(theta))
  r_null <- to_r(0)
  ring_null <- data.frame(x = r_null * cos(theta), y = r_null * sin(theta))

  cols <- c("over-represented" = .lsa_div_high, "avoided" = .lsa_div_low,
            "n.s." = "grey70")

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = df,                      # faint guide spokes
      ggplot2::aes(x = .data$x_in, y = .data$y_in,
                   xend = .data$x_out, yend = .data$y_out),
      colour = "grey92", linewidth = 0.3) +
    ggplot2::geom_path(data = ring_in, ggplot2::aes(x = .data$x, y = .data$y),
      colour = "grey85", linewidth = 0.25) +
    ggplot2::geom_path(data = ring_out, ggplot2::aes(x = .data$x, y = .data$y),
      colour = "grey85", linewidth = 0.25) +
    ggplot2::geom_path(data = ring_null, ggplot2::aes(x = .data$x, y = .data$y),
      colour = "grey55", linewidth = 0.35, linetype = "dashed") +
    ggplot2::geom_segment(data = df,                      # CI bars
      ggplot2::aes(x = .data$x_lo, y = .data$y_lo,
                   xend = .data$x_hi, yend = .data$y_hi, colour = .data$status),
      linewidth = 0.8, lineend = "round") +
    ggplot2::geom_point(data = df,                        # estimate squares
      ggplot2::aes(x = .data$x_est, y = .data$y_est, colour = .data$status),
      shape = 15, size = 1.7) +
    ggplot2::geom_text(data = df,
      ggplot2::aes(x = .data$x_lab, y = .data$y_lab, label = .data$edge,
                   angle = .data$text_angle, hjust = .data$hjust,
                   colour = .data$status),
      size = label_size, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
    ggplot2::coord_equal(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4),
                         clip = "off") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40",
                                            size = 8.5),
      plot.margin = ggplot2::margin(8, 8, 8, 8)) +
    ggplot2::labs(
      title = sprintf("Bootstrap forest \u2014 %s",
        c(residuals = "adjusted residuals", count = "counts",
          prob = "probabilities", yules_q = "Yule's Q")[[metric]]),
      subtitle = sprintf(
        "%d resamples \u00b7 spoke = %.0f%% CI, square = estimate, dashed ring = null",
        boot$R, 100 * boot$level_alpha))
}

#' @rdname plot_forest
#' @param x An `lsa_bootstrap` object (for the `plot()` method).
#' @param ... Passed to [plot_forest()] (e.g. `metric`, `n_top`).
#' @export
plot.lsa_bootstrap <- function(x, ...) plot_forest(x, ...)

# Coerce a possibly-NA logical vector to a plain logical (NA -> FALSE).
.lsa_true <- function(x) !is.na(x) & x
