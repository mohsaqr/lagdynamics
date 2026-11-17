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

#' Plot a lag sequential analysis fit
#'
#' @description
#' `plot()` draws a lag sequential analysis fit as a residual heatmap, a lag
#' transition network, a chord diagram or a polar sunburst, chosen with `type`;
#' each display suits a different question. A residual heatmap lays out the full
#' table of antecedent by consequent states, coloured by the adjusted residual,
#' and shows at a glance which transitions are over-represented and which are
#' avoided. A lag transition network draws the states as nodes and the
#' transitions as directed edges, which shows the paths the process takes. A
#' chord diagram shows how the transitions out of each state are shared among
#' the targets, and a polar sunburst shows the same per-source breakdown on a
#' ring.
#'
#' The heatmap (default) shows the adjusted residuals, with significant cells in
#' bold; `which` switches it to transition probabilities, observed counts or
#' expected counts. The other types are drawn by [plot_transitions()],
#' [plot_chords()] and [plot_polar()], whose arguments can be passed through
#' `...`. A grouped fit is drawn once per group, or in one figure with `combined
#' = TRUE`.
#'
#' For the uncertainty of the residuals use [plot_forest()] on a
#' bootstrap; for differences between groups, `plot()` on the result of
#' [compare_lsa()] (see [plot.lsa_comparison()]).
#'
#' @param x An `lsa` fit from [lsa()], or a grouped fit from
#'   `lsa(..., group = )`.
#' @param type `"heatmap"` (default), `"network"`, `"chord"` or
#'   `"sunburst"`.
#' @param ... Passed to the chosen display. For `"heatmap"`: `which`, one
#'   of `"residuals"` (default), `"prob"`, `"count"` or `"expected"`. For
#'   the other types, the arguments of [plot_transitions()],
#'   [plot_chords()] or [plot_polar()].
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`
#' when the fit holds no finite values of the quantity requested for the
#' heatmap.
#'
#' @return For `"heatmap"` and `"sunburst"`, a `ggplot` object, drawn when
#'   printed. For `"network"` and `"chord"`, the value of
#'   [plot_transitions()] or [plot_chords()], invisibly, with the plot drawn
#'   on the current device. For a grouped fit, a list with one such value
#'   per group, invisibly, after drawing.
#'
#' @examples
#' fit <- lsa(engagement)
#' plot(fit)
#' plot(fit, which = "prob")
#' plot(fit, type = "network")
#' plot(fit, type = "chord")
#'
#' gfit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'             time = "Time", group = "Achiever")
#' plot(gfit, combined = TRUE)
#'
#' @seealso [plot_transitions()], [plot_chords()], [plot_polar()],
#'   [plot_forest()], [plot.lsa_comparison()], [transitions()]
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
#' @param combined For a grouped fit: `FALSE` (default) draws one figure
#'   per group; `TRUE` tiles the groups in one figure.
#' @export
plot.lsa_group <- function(x, type = c("heatmap", "network", "chord",
                                       "sunburst"), combined = FALSE, ...) {
  type <- match.arg(type)
  groups <- names(x)
  is_gg <- type %in% c("heatmap", "sunburst")
  if (is_gg && !requireNamespace("ggplot2", quietly = TRUE)) {
    .lsa_stop(class = "missing_package", "Package 'ggplot2' is required to plot this type.", call. = FALSE)
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
    .lsa_stop(class = "missing_package", "Package 'ggplot2' is required for plot.lsa(). ",
         "Install with install.packages('ggplot2').", call. = FALSE)
  }
  col <- switch(which, residuals = "adj_res", prob = "prob",
                count = "count", expected = "expected")
  e <- x$edges
  if (!col %in% names(e) || all(!is.finite(e[[col]]))) {
    .lsa_stop("This fit has no '", which, "' values to plot.", call. = FALSE)
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
  caps <- c(residuals = "warm over-represented, cool avoided; strongest significant cell of each row bold",
            prob = "P(to | from); modal next state bold",
            count = "modal next state bold", expected = "modal next state bold")
  sub <- .lsa_ascii(sprintf("%s engine, lag %d \u00b7 %d transitions \u00b7 %s",
                 x$method, x$params$lag, sum(x$obs), caps[[which]]))

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
