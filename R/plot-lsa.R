# Heatmap for an lsa fit, implemented directly in ggplot2 in the exact
# style of pathtree::plot_pathways(): white tile borders, the salient
# cell of each row annotated in bold white, lesser cells in grey,
# near-zero cells blank, the #deebf7 -> #08519c sequential ramp (fixed
# 0-1 for probabilities), monospace row labels, theme_minimal, no grid,
# a bold title and a descriptive subtitle. Adjusted residuals use a
# diverging blue-white-red scale centred at 0 (significant cells bold).
#
# Heatmaps are owned here; only the transition NETWORK is delegated to
# cograph::splot() (see plot_transitions()).

#' Plot an LSA Fit
#'
#' One entry point for every view of a fit; pick it with `type`:
#' `"heatmap"` (default, the `from x to` residual heatmap), `"network"`
#' (transition network via [cograph::splot()]), `"chord"` (chord diagram
#' via [cograph::plot_chord()]), or `"sunburst"` (polar rose). Extra
#' arguments are forwarded to the chosen view's worker
#' ([plot_transitions()], [plot_chords()], [plot_polar()]); see those for
#' view-specific options.
#'
#' @param x An `lsa` fit from [lsa()].
#' @param type Which view to draw: `"heatmap"` (default), `"network"`,
#'   `"chord"`, or `"sunburst"`.
#' @param ... Forwarded to the chosen view. For `"heatmap"`: `which`
#'   (`"residuals"` (default), `"prob"`, `"count"`, `"expected"`). For
#'   `"network"`/`"chord"`: `weights`. For `"sunburst"`: `style`, `fill`.
#'
#' @return A `ggplot` object for `"heatmap"` and `"sunburst"`; the
#'   (invisible) `cograph` object for `"network"` and `"chord"`. Drawn
#'   when printed.
#'
#' @examples
#' \dontrun{
#' fit <- lsa(group_regulation)
#' plot(fit)                     # residual heatmap (default)
#' plot(fit, which = "prob")     # heatmap of probabilities
#' plot(fit, type = "network")   # transition network
#' plot(fit, type = "chord")     # chord diagram
#' plot(fit, type = "sunburst")  # polar sunburst
#' }
#'
#' @seealso [plot_transitions()], [plot_chords()], [plot_polar()],
#'   [plot_forest()], [transitions()]
#'
#' @export
plot.lsa <- function(x, type = c("heatmap", "network", "chord", "sunburst"),
                     ...) {
  type <- match.arg(type)
  switch(type,
    heatmap  = .plot_lsa_heatmap(x, ...),
    network  = plot_transitions(x, ...),
    chord    = plot_chords(x, ...),
    sunburst = plot_polar(x, ...))
}

#' @rdname plot.lsa
#' @param combined Logical, for a grouped fit only. `FALSE` (default)
#'   draws each group as its own full-size figure; `TRUE` tiles all
#'   groups into a single figure (compact, but cramped for many groups).
#' @export
plot.lsa_group <- function(x, type = c("heatmap", "network", "chord",
                                       "sunburst"), combined = FALSE, ...) {
  type <- match.arg(type)
  groups <- names(x)
  is_gg <- type %in% c("heatmap", "sunburst")
  if (is_gg && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to plot this type.", call. = FALSE)
  }

  if (!isTRUE(combined)) {
    # One full-size figure per group (the default), each titled by group.
    if (is_gg) {
      plots <- lapply(groups, function(gname)
        plot(x[[gname]], type = type, ...) + ggplot2::ggtitle(gname))
      for (p in plots) print(p)
      return(invisible(plots))
    }
    res <- lapply(groups, function(gname) {
      out <- plot(x[[gname]], type = type, ...)
      graphics::title(main = gname)
      out
    })
    return(invisible(res))
  }

  # combined = TRUE: tile all groups into one figure.
  n <- length(groups)
  nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc)
  if (is_gg) {
    plots <- lapply(groups, function(gname)
      plot(x[[gname]], type = type, ...) + ggplot2::ggtitle(gname))
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nr, nc)))
    for (i in seq_along(plots)) {
      print(plots[[i]], vp = grid::viewport(
        layout.pos.row = ((i - 1L) %/% nc) + 1L,
        layout.pos.col = ((i - 1L) %% nc) + 1L))
    }
    return(invisible(plots))
  }
  op <- graphics::par(mfrow = c(nr, nc)); on.exit(graphics::par(op), add = TRUE)
  res <- lapply(groups, function(gname) {
    out <- plot(x[[gname]], type = type, ...)
    graphics::title(main = gname)
    out
  })
  invisible(res)
}

# Heatmap worker (type = "heatmap"), owned here in ggplot2: a from x to
# heatmap with rows ordered by outgoing volume, the salient cell of each
# row bold, lesser cells grey, 0-edge cells greyed and unlabelled.
.plot_lsa_heatmap <- function(x,
                     which = c("residuals", "prob", "count", "expected"),
                     ...) {
  which <- match.arg(which)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.lsa(). ",
         "Install with install.packages('ggplot2').", call. = FALSE)
  }
  col <- switch(which, residuals = "adj_res", prob = "prob",
                count = "count", expected = "expected")
  e <- x$edges
  if (!col %in% names(e) || all(!is.finite(e[[col]]))) {
    stop("This fit has no '", which, "' values to plot.", call. = FALSE)
  }
  e$.val <- e[[col]]

  # Rows (source states) ordered by outgoing volume, most active on top.
  flevels <- rownames(x$obs)[order(rowSums(x$obs))]   # asc -> top after rev
  e$.from <- factor(e$from_label, levels = flevels)
  e$.to   <- factor(e$to_label,   levels = rownames(x$obs))

  # 0-edge cells (no observed transition) render grey (the na.value) and
  # carry no label; every other cell is labelled.
  e$.val[!is.finite(e$count) | e$count == 0] <- NA_real_

  # Salient cell per row (modal next state, or largest |residual|) is
  # bold; all other non-zero cells are labelled in grey.
  rank_val <- if (which == "residuals") abs(e$.val) else e$.val
  salient_to <- tapply(seq_len(nrow(e)), e$from_label, function(i) {
    v <- rank_val[i]
    if (all(!is.finite(v))) NA_character_ else e$to_label[i][which.max(v)]
  })
  e$.bold <- is.finite(e$.val) &
    e$to_label == salient_to[as.character(e$from_label)]
  if (which == "residuals") e$.bold <- e$.bold &
    !is.na(e$significant) & e$significant
  e$.show <- is.finite(e$.val)            # label every non-zero cell
  d <- if (which == "count") 0L else 2L
  e$.lab <- formatC(e$.val, format = "f", digits = d)

  titles <- c(residuals = "Adjusted residuals",
              prob = "Transition probabilities",
              count = "Observed counts", expected = "Expected counts")
  caps <- c(residuals = "warm over-represented, cool avoided; significant cells bold",
            prob = "P(to | from); modal next state bold",
            count = "modal next state bold", expected = "modal next state bold")
  sub <- sprintf("%s engine, lag %d \u00b7 %d transitions \u00b7 %s",
                 x$method, x$params$lag, sum(x$obs), caps[[which]])

  gg <- ggplot2::ggplot(
    e, ggplot2::aes(x = .data$.to, y = .data$.from, fill = .data$.val)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::labs(title = titles[[which]], subtitle = sub,
                  x = "Next state", y = NULL,
                  fill = c(residuals = "z", prob = "P(to | from)",
                           count = "count", expected = "expected")[[which]]) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(family = "mono"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 8.5)
    )

  if (which == "residuals") {
    lim <- max(abs(e$.val[is.finite(e$.val)]))
    if (!is.finite(lim) || lim <= 0) lim <- 1   # degenerate (all 0/NA) fit
    gg <- gg + ggplot2::scale_fill_gradient2(
      low = .lsa_div_low, mid = .lsa_div_mid, high = .lsa_div_high,
      midpoint = 0, limits = c(-lim, lim), na.value = "grey85",
      name = "z")
  } else if (which == "prob") {
    gg <- gg + ggplot2::scale_fill_gradient(
      low = .lsa_seq_low, high = .lsa_seq_high, limits = c(0, 1),
      na.value = "grey85")
  } else {
    gg <- gg + ggplot2::scale_fill_gradient(
      low = .lsa_seq_low, high = .lsa_seq_high, na.value = "grey85")
  }

  soft <- e[e$.show & !e$.bold, , drop = FALSE]
  bold <- e[e$.show & e$.bold, , drop = FALSE]
  if (nrow(soft)) {
    gg <- gg + ggplot2::geom_text(data = soft,
      ggplot2::aes(label = .data$.lab), colour = "grey20", size = 3.0)
  }
  if (nrow(bold)) {
    gg <- gg + ggplot2::geom_text(data = bold,
      ggplot2::aes(label = .data$.lab), colour = "white",
      fontface = "bold", size = 3.2)
  }
  gg
}
