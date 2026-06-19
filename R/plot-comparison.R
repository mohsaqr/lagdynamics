# Plots for a group comparison from compare_lsa().
#
# The default is a back-to-back BARREL (pyramid): one row per
# transition, the first group's bar running left and the second's right
# from a central divider. Bar LENGTH is each group's transition
# probability (or count); bar FILL is each group's log odds ratio (the
# N-invariant LSA effect size); bar ends carry the observed count; the
# centre chip is the between-
# group difference p-value (bold + starred when significant). The
# pairwise object (more than two groups) draws one barrel per pair via
# facets, since every pairwise comparison is itself an a-vs-b.
#
# The palette follows the vcd::mosaic / Nestimate residual convention:
# BLUE = over-represented (positive residual, observed > expected),
# RED = avoided (negative residual). Same hues as Nestimate's mosaic.
#
# An alternate "heatmap" style draws the signed difference as a
# from x to grid (faceted by pair for the pairwise object).

# Comparison-specific diverging palette (vcd/Nestimate mosaic convention).
.cmp_low  <- "#D33F6A"   # rosy red        -> negative residual (avoided)
.cmp_mid  <- "#F7F7F7"   # near-white      -> at expectation
.cmp_high <- "#4A6FE3"   # periwinkle blue -> positive residual (over-represented)

#' Plot a Group Comparison
#'
#' Two views of a [compare_lsa()] result. The default `"barrel"` is a
#' back-to-back pyramid (one row per transition): the first group's bar
#' runs left and the second's right, bar length is each group's
#' transition probability, bar colour is each group's log odds ratio
#' (blue = over-represented, red = avoided, following the vcd / mosaic
#' convention), bar ends show the observed count, and the centre chip
#' shows the difference p-value (bold and starred when significant). The
#' bar of the group with the higher value gets a border that darkens with
#' the size of the difference (faint = small, dark = large). For more
#' than two groups, one
#' barrel is drawn per pair via facets. The `"heatmap"` style draws the
#' signed difference as a `from x to` grid on the same diverging scale.
#'
#' @param x An `lsa_comparison` or `lsa_comparison_pairwise` object.
#' @param style `"barrel"` (default) or `"heatmap"`.
#' @param value For `"barrel"`, the quantity mapped to bar length:
#'   `"prob"` (transition probability, default) or `"count"`.
#' @param rank For `"barrel"`, how to choose which transitions to show:
#'   `"frequency"` (default) ranks by pooled observed count -- the
#'   backbone transitions, which are mostly over-represented (blue);
#'   `"effect"` ranks by the strongest association in either group
#'   (`|log OR|`, among tested cells), surfacing both over- (blue) and
#'   under-represented / avoided (red) transitions.
#' @param top_n For `"barrel"`, how many transitions to show (highest
#'   `rank` on top; for the pairwise object the ranking is shared across
#'   facets so they line up). Default `12`.
#' @param ... Reserved.
#'
#' @return A `ggplot` object (drawn when printed). Needs `ggplot2`.
#'
#' @examples
#' \dontrun{
#' grp <- ifelse(group_regulation$T1 == "plan", "starts_plan", "other")
#' g <- lsa(group_regulation, group = grp)
#' cmp <- compare_lsa(g, R = 200)
#' plot(cmp)                    # back-to-back barrel
#' plot(cmp, style = "heatmap") # difference heatmap
#' }
#'
#' @seealso [compare_lsa()], [plot.lsa()]
#'
#' @export
plot.lsa_comparison <- function(x, style = c("barrel", "heatmap"),
                                value = c("prob", "count"),
                                rank = c("frequency", "effect"),
                                top_n = 12L, ...) {
  style <- match.arg(style)
  value <- match.arg(value)
  rank <- match.arg(rank)
  if (style == "heatmap") {
    labels <- rownames(x$fits[[1L]]$obs)
    e <- x$edges
    e$.diff <- e$diff
    e$.sig  <- !is.na(e$significant) & e$significant
    ga <- x$groups[1L]; gb <- x$groups[2L]
    nsig <- sum(e$.sig, na.rm = TRUE)
    sub <- sprintf(
      paste0("%s difference \u00b7 %d permutations \u00b7 %d significant edges ",
             "(adjust = %s)\nblue = %s higher, red = %s higher; ",
             "significant cells bold"),
      x$measure, x$R, nsig, x$adjust, ga, gb)
    return(.plot_compare_heatmap(
      e, labels, title = sprintf("Group difference: %s vs %s", ga, gb),
      subtitle = sub))
  }

  fit_a <- x$fits[[1L]]; fit_b <- x$fits[[2L]]
  ga <- x$groups[1L]; gb <- x$groups[2L]
  labels <- rownames(fit_a$obs)
  e <- x$edges
  ix <- cbind(match(e$from, labels), match(e$to, labels))
  if (rank == "frequency") {
    rank_val <- fit_a$obs[ix] + fit_b$obs[ix]
  } else {
    # Strongest association in either group (|log OR|), among tested cells
    # only, so the barrel surfaces both over- (blue) and under-represented
    # (red) transitions rather than just the frequent backbone.
    rank_val <- pmax(abs(.lsa_log_or(fit_a)[ix]),
                     abs(.lsa_log_or(fit_b)[ix]))
    rank_val[!is.finite(e$p_perm)] <- -Inf
  }
  sel <- .barrel_select(e, top_n, rank_val)
  one <- .barrel_one(fit_a, fit_b, sel$edges, sel$y, value)
  gutter <- .barrel_gutter(one$max_len)

  nsig <- sum(sel$edges$significant, na.rm = TRUE)
  sub <- sprintf(
    paste0("bar length = %s \u00b7 fill = log odds ratio \u00b7 centre = ",
           "difference p (adjust = %s)\n%s  \u2190\u2190    ",
           "\u2192\u2192  %s \u00b7 %d of %d shown significant \u00b7 rows by %s"),
    .barrel_value_name(value), x$adjust, ga, gb, nsig, nrow(sel$edges),
    rank)

  .barrel_plot(one$rect, one$p, sel$labels, sel$y_breaks, gutter,
               value, title = sprintf("Transition comparison: %s vs %s",
                                       ga, gb),
               subtitle = sub, facet = FALSE)
}

#' @rdname plot.lsa_comparison
#' @export
plot.lsa_comparison_pairwise <- function(x, style = c("barrel", "heatmap"),
                                         value = c("prob", "count"),
                                         rank = c("frequency", "effect"),
                                         top_n = 12L, ...) {
  style <- match.arg(style)
  value <- match.arg(value)
  rank <- match.arg(rank)
  if (style == "heatmap") {
    labels <- rownames(x$comparisons[[1L]]$fits[[1L]]$obs)
    e <- x$edges
    e$.diff <- e$diff
    e$.sig  <- !is.na(e$significant) & e$significant
    pair_levels <- paste0(x$global$group_a, " \u2212 ", x$global$group_b)
    e$.pair <- factor(paste0(e$group_a, " \u2212 ", e$group_b),
                      levels = pair_levels)
    nsig <- sum(e$.sig, na.rm = TRUE)
    sub <- sprintf(
      paste0("%s difference \u00b7 %d permutations/pair \u00b7 %d significant ",
             "edges (adjust = %s, family-wide)\nblue = first group ",
             "higher, red = second higher; significant cells bold"),
      x$measure, x$R, nsig, x$adjust)
    return(.plot_compare_heatmap(e, labels,
      title = "Pairwise group differences", subtitle = sub,
      facet = ".pair"))
  }

  # Shared row set: edges ranked by total transition frequency across all
  # groups, so every facet shows (and can be compared on) the same
  # backbone transitions rather than rare cells with noisy odds ratios.
  base <- x$comparisons[[1L]]$edges[, c("from", "to")]
  labels <- rownames(x$comparisons[[1L]]$fits[[1L]]$obs)
  gfits <- list()
  for (cp in x$comparisons) {
    gfits[[cp$groups[1L]]] <- cp$fits[[1L]]
    gfits[[cp$groups[2L]]] <- cp$fits[[2L]]
  }
  bix <- cbind(match(base$from, labels), match(base$to, labels))
  if (rank == "frequency") {
    rank_vec <- Reduce(`+`, lapply(gfits, function(f) f$obs))[bix]
  } else {
    # Strongest association across any group (|log OR|), among cells tested
    # in at least one pair, so red (avoided) transitions can surface.
    rank_vec <- Reduce(pmax, lapply(gfits,
                                    function(f) abs(.lsa_log_or(f))))[bix]
    tested_any <- tapply(is.finite(x$edges$p_perm),
                         paste(x$edges$from, x$edges$to, sep = "\r"), any)
    rank_vec[!tested_any[paste(base$from, base$to, sep = "\r")]] <- -Inf
  }
  rank_vec[!is.finite(rank_vec)] <- -Inf
  sel_idx <- utils::head(order(rank_vec, decreasing = TRUE),
                         min(as.integer(top_n), nrow(base)))
  shared <- base[sel_idx, , drop = FALSE]
  shared$y <- rev(seq_len(nrow(shared)))
  labels_y <- paste0(shared$from, " -> ", shared$to)
  pair_levels <- paste0(x$global$group_a, " \u2212 ", x$global$group_b)

  parts <- lapply(seq_along(x$comparisons), function(p) {
    cp <- x$comparisons[[p]]
    ce <- x$edges[x$edges$group_a == cp$groups[1L] &
                  x$edges$group_b == cp$groups[2L], ]
    m <- match(paste(shared$from, shared$to, sep = "\r"),
               paste(ce$from, ce$to, sep = "\r"))
    edges_p <- data.frame(
      from = shared$from, to = shared$to, diff = ce$diff[m],
      p_adj = ce$p_adj[m], significant = ce$significant[m],
      stringsAsFactors = FALSE)
    one <- .barrel_one(cp$fits[[1L]], cp$fits[[2L]], edges_p, shared$y,
                       value)
    facet_lab <- factor(pair_levels[p], levels = pair_levels)
    one$rect$.pair <- facet_lab
    one$p$.pair <- facet_lab
    one
  })
  rect <- do.call(rbind, lapply(parts, function(o) o$rect))
  pdf  <- do.call(rbind, lapply(parts, function(o) o$p))
  gutter <- .barrel_gutter(max(vapply(parts, function(o) o$max_len,
                                      numeric(1))))

  nsig <- sum(x$edges$significant, na.rm = TRUE)
  sub <- sprintf(
    paste0("bar length = %s \u00b7 fill = log odds ratio \u00b7 centre = ",
           "difference p (adjust = %s, family-wide)\nleft = first group, ",
           "right = second group of each pair \u00b7 %d significant"),
    .barrel_value_name(value), x$adjust, nsig)

  .barrel_plot(rect, pdf, labels_y, shared$y, gutter, value,
               title = "Pairwise transition comparison",
               subtitle = sub, facet = TRUE)
}

# --- barrel helpers ---------------------------------------------------

.barrel_value_name <- function(value) {
  if (value == "prob") "P(to | from)" else "count"
}

# Small offset that pushes both groups' bars away from the centre so the
# p-value chip has room. Scaled to the longest bar.
.barrel_gutter <- function(max_len) {
  g <- max_len * 0.06
  if (!is.finite(g) || g <= 0) 0.01 else g
}

# Pick the top_n edges by `rank_val` (highest on top) and assign row
# positions. Ranking by transition frequency keeps the barrel on the
# backbone transitions rather than rare cells with noisy odds ratios.
# Returns the edge subset, its y positions, axis labels/breaks.
.barrel_select <- function(edges, top_n, rank_val) {
  rank_val[!is.finite(rank_val)] <- -Inf
  keep <- utils::head(order(rank_val, decreasing = TRUE),
                      min(as.integer(top_n), nrow(edges)))
  sub <- edges[keep, , drop = FALSE]
  y <- rev(seq_len(nrow(sub)))
  list(edges = sub, y = y, y_breaks = y,
       labels = paste0(sub$from, " -> ", sub$to))
}

# Build the rectangle and p-value data frames for one a-vs-b pair over a
# given edge set. `edges` carries from, to, p_adj, significant; `y` is
# the shared row index per edge. Bar fill = log odds ratio (clamped
# to +/-3); length = prob or count, pulled from each fit.
.barrel_one <- function(fit_a, fit_b, edges, y, value) {
  labels <- rownames(fit_a$obs)
  ix <- cbind(match(edges$from, labels), match(edges$to, labels))
  len_a <- if (value == "prob") fit_a$prob[ix] else fit_a$obs[ix]
  len_b <- if (value == "prob") fit_b$prob[ix] else fit_b$obs[ix]
  len_a[!is.finite(len_a)] <- 0
  len_b[!is.finite(len_b)] <- 0
  # Fill = each group's log odds ratio (the N-invariant LSA effect size,
  # the default compare measure), clamped for the colour scale.
  clamp <- function(z) pmax(-3, pmin(3, z))
  res_a <- clamp(.lsa_log_or(fit_a)[ix]); res_b <- clamp(.lsa_log_or(fit_b)[ix])
  cnt_a <- fit_a$obs[ix]; cnt_b <- fit_b$obs[ix]
  bar_hw <- 0.40

  # Direction + magnitude mark: border the bar of the group with the
  # higher tested value (diff = measure_a - measure_b; group a is left),
  # darker the larger |diff| (faint for small differences, dark for
  # large), so a single mark shows which group is higher and by how much.
  dd <- edges$diff
  mag <- abs(dd); mag[!is.finite(mag)] <- 0
  denom <- max(mag); if (!is.finite(denom) || denom <= 0) denom <- 1
  dark <- mag / denom                              # 0 (faint) .. 1 (dark)
  win_col <- grDevices::grey(1 - 0.85 * dark)       # white .. near-black
  left_win  <- is.finite(dd) & dd > 0
  right_win <- is.finite(dd) & dd < 0
  border_l <- ifelse(left_win,  win_col, "white")
  border_r <- ifelse(right_win, win_col, "white")
  lw_l <- ifelse(left_win,  0.8, 0.3)
  lw_r <- ifelse(right_win, 0.8, 0.3)

  rect <- rbind(
    data.frame(y = y, ymin = y - bar_hw, ymax = y + bar_hw,
               len = len_a, resid = res_a, count = cnt_a,
               border = border_l, lwd = lw_l,
               side = "left", stringsAsFactors = FALSE),
    data.frame(y = y, ymin = y - bar_hw, ymax = y + bar_hw,
               len = len_b, resid = res_b, count = cnt_b,
               border = border_r, lwd = lw_r,
               side = "right", stringsAsFactors = FALSE)
  )
  sig <- !is.na(edges$significant) & edges$significant
  padj <- edges$p_adj
  p_txt <- ifelse(!is.finite(padj), "\u2014",
            ifelse(padj < 0.001, "<.001", sprintf("%.3f", padj)))
  p <- data.frame(
    y = y, label = ifelse(sig, paste0(p_txt, "*"), p_txt),
    face = ifelse(sig, "bold", "plain"),
    stringsAsFactors = FALSE)
  list(rect = rect, p = p, max_len = max(c(len_a, len_b), 0,
                                         na.rm = TRUE))
}

# Assemble the barrel ggplot from rect/p frames. Resolves left/right bar
# geometry from `side` and `gutter` here so callers stay declarative.
.barrel_plot <- function(rect, p, labels_y, y_breaks, gutter, value,
                         title, subtitle, facet) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to plot a comparison. ",
         "Install with install.packages('ggplot2').", call. = FALSE)
  }
  left <- rect$side == "left"
  rect$xmin <- ifelse(left, -rect$len - gutter, gutter)
  rect$xmax <- ifelse(left, -gutter, rect$len + gutter)
  rect$end_x <- ifelse(left, rect$xmin, rect$xmax)
  rect$hj <- ifelse(left, 1.2, -0.2)

  gg <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, colour = "grey75",
                        linewidth = 0.5) +
    ggplot2::geom_rect(data = rect,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax,
                   fill = .data$resid, colour = .data$border),
      linewidth = rect$lwd) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_text(data = rect,
      ggplot2::aes(x = .data$end_x, y = .data$y, label = .data$count),
      hjust = rect$hj, size = 2.6, colour = "grey35") +
    ggplot2::geom_label(data = p,
      ggplot2::aes(x = 0, y = .data$y, label = .data$label),
      fontface = p$face, size = 2.5, colour = "grey15",
      fill = "white", linewidth = 0.15, label.padding =
        ggplot2::unit(0.08, "lines")) +
    ggplot2::scale_fill_gradient2(
      low = .cmp_low, mid = .cmp_mid, high = .cmp_high,
      midpoint = 0, limits = c(-3, 3), name = "log OR") +
    ggplot2::scale_y_continuous(breaks = y_breaks, labels = labels_y,
      expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_x_continuous(
      labels = function(v) format(abs(v), scientific = FALSE)) +
    ggplot2::labs(x = .barrel_value_name(value), y = NULL,
                  title = title, subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(family = "mono", size = 8),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92",
                                                 linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(colour = "grey40",
                                            size = 8.5, hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold", size = 9))

  if (isTRUE(facet)) {
    gg <- gg + ggplot2::facet_wrap(~ .pair)
  }
  gg
}

# --- heatmap worker (alternate style) ---------------------------------

# `e` carries .diff (signed difference, NA where not estimable), .sig
# (logical significance), from, to, optional facet column. `labels` is
# the state order.
.plot_compare_heatmap <- function(e, labels, title, subtitle,
                                  facet = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to plot a comparison. ",
         "Install with install.packages('ggplot2').", call. = FALSE)
  }
  e$.from <- factor(e$from, levels = rev(labels))
  e$.to   <- factor(e$to,   levels = labels)
  e$.diff[!is.finite(e$.diff)] <- NA_real_
  e$.lab <- formatC(e$.diff, format = "f", digits = 2)
  e$.show <- is.finite(e$.diff)

  lim <- suppressWarnings(max(abs(e$.diff[is.finite(e$.diff)])))
  if (!is.finite(lim) || lim <= 0) lim <- 1

  gg <- ggplot2::ggplot(
    e, ggplot2::aes(x = .data$.to, y = .data$.from, fill = .data$.diff)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(
      low = .cmp_low, mid = .cmp_mid, high = .cmp_high,
      midpoint = 0, limits = c(-lim, lim), na.value = "grey85",
      name = "diff") +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = "Next state", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(family = "mono"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 8.5)
    )

  if (!is.null(facet)) {
    gg <- gg + ggplot2::facet_wrap(stats::as.formula(paste0("~ ", facet)))
  }

  soft <- e[e$.show & !e$.sig, , drop = FALSE]
  bold <- e[e$.show & e$.sig, , drop = FALSE]
  if (nrow(soft)) {
    gg <- gg + ggplot2::geom_text(data = soft,
      ggplot2::aes(label = .data$.lab), colour = "grey20", size = 3.0)
  }
  if (nrow(bold)) {
    gg <- gg + ggplot2::geom_text(data = bold,
      ggplot2::aes(label = .data$.lab), colour = "grey10",
      fontface = "bold", size = 3.2)
  }
  gg
}
