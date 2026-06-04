# Polar sunburst of an lsa fit, drawn directly in ggplot2 (a sibling of
# the plot.lsa() heatmap). Two styles share one large inner hole that
# carries the source-state names along its arc:
#
#   style = "rose"  (default): equal sectors, equal slots -- one per
#     target. A slot's radial BAR HEIGHT is the transition frequency
#     (sqrt-scaled) and its FILL is the adjusted residual. Equal slots
#     leave room to label every outer target, and nothing is crammed.
#   style = "wedge": sectors sized by a source's outgoing volume, each
#     target a wedge whose angular WIDTH is the transition's frequency
#     share, filled by the residual. Tiny wedges below `min_show` are
#     omitted and only wide wedges are labelled, so the dense small
#     transitions do not collapse into unreadable slivers.
#
# Owned here (like the heatmap); only the network and chord go to cograph.

# Readable tangential text angle for a label on the ring at angle `a`
# (radians): rotate along the arc, flipping the bottom half so it never
# reads upside-down.
.lsa_tangential <- function(a) {
  deg <- (a %% (2 * pi)) * 180 / pi
  rot <- deg - 90
  ifelse(deg > 180, rot + 180, rot)
}

# A filled polygon spanning [a0, a1] x [r0, r1] (a wedge / bar segment).
.lsa_wedge <- function(a0, a1, r0, r1, id, value) {
  arc <- seq(a0, a1, length.out = 14)
  data.frame(id = id, value = value,
             x = c(r0 * cos(arc), r1 * cos(rev(arc))),
             y = c(r0 * sin(arc), r1 * sin(rev(arc))),
             stringsAsFactors = FALSE)
}

# Larger inner circle (shared by both styles): empty centre + source band.
.lsa_polar_geom <- list(src_in = 0.34, src_out = 0.50, ring0 = 0.54, ring1 = 1)

#' Polar Sunburst of an LSA Fit
#'
#' Draws the transition structure as a polar sunburst with the source
#' states named along a large inner ring. Two styles: `"rose"` (default)
#' gives every target an equal angular slot and encodes frequency as the
#' radial **bar height** (so nothing crams); `"wedge"` sizes each
#' transition's angular **width** by its frequency share (the classic
#' look), omitting tiny wedges. Both fill by the adjusted residual (warm
#' = over-represented, cool = avoided), sharing the [plot.lsa()] heatmap
#' colour scale. Needs `ggplot2`.
#'
#' @param fit An `lsa` fit from [lsa()].
#' @param style `"rose"` (default, equal slots + bar height) or `"wedge"`
#'   (frequency-proportional wedge width).
#' @param fill Which quantity fills the bars/wedges: `"residuals"`
#'   (default, diverging), `"prob"`, or `"lift"`.
#' @param size For `style = "rose"`, which non-negative quantity sets bar
#'   height: `"count"` (default) or `"prob"`. Ignored for `"wedge"`.
#' @param significant Logical. Grey out non-significant cells (keeping
#'   their size). Default `FALSE`.
#' @param labels Which target cells to name: `"all"`, `"auto"`, or
#'   `"none"`. Default is `"all"` for `"rose"` (equal slots leave room)
#'   and `"auto"` for `"wedge"` (only wedges wide enough to fit a name).
#'   Source names are always shown.
#' @param min_show For `style = "wedge"`, drop wedges whose frequency
#'   share of their source's outflow is below this fraction. Default
#'   `0.01`; `0` keeps all.
#' @param label_size Label text size. Default `3`.
#' @param ... Ignored; accepted so `plot(fit, type = "sunburst", ...)`
#'   can forward arguments without error.
#'
#' @return A `ggplot` object (drawn when printed).
#'
#' @examples
#' \dontrun{
#' fit <- lsa(group_regulation)
#' plot_polar(fit)                          # rose: bars filled by residual
#' plot_polar(fit, style = "wedge")         # classic frequency wedges
#' plot_polar(fit, significant = TRUE)      # non-significant cells greyed
#' }
#'
#' @seealso [plot.lsa()] (heatmap), [plot_chords()] (chord),
#'   [plot_forest()] (bootstrap forest)
#'
#' @export
plot_polar <- function(fit, style = c("rose", "wedge"),
                       fill = c("residuals", "prob", "lift"),
                       size = c("count", "prob"), significant = FALSE,
                       labels = c("all", "auto", "none"),
                       min_show = 0.01, label_size = 3, ...) {
  stopifnot(inherits(fit, "lsa"))
  style <- match.arg(style)
  fill <- match.arg(fill); size <- match.arg(size)
  labels <- if (missing(labels)) (if (style == "rose") "all" else "auto")
            else match.arg(labels)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_polar(). ",
         "Install with install.packages('ggplot2').", call. = FALSE)
  }
  diverging <- fill == "residuals"

  states <- rownames(fit$obs)
  K <- length(states)
  fillmat <- switch(fill, residuals = fit$adj_res, prob = fit$prob,
                    lift = .lsa_weight_matrix(fit, "lift", FALSE))
  sig <- is.finite(fit$p) & fit$p < fit$params$alpha
  g <- .lsa_polar_geom

  built <- if (style == "rose") {
    .polar_rose(fit, states, K, fillmat, sig, size, significant, labels, g)
  } else {
    .polar_wedge(fit, states, K, fillmat, sig, significant, labels, min_show, g)
  }

  src_lab_df <- data.frame(
    label = states,
    x = ((g$src_in + g$src_out) / 2) * cos(built$sector_mid),
    y = ((g$src_in + g$src_out) / 2) * sin(built$sector_mid),
    angle = .lsa_tangential(built$sector_mid), stringsAsFactors = FALSE)

  theta <- seq(0, 2 * pi, length.out = 320)
  rings <- rbind(
    data.frame(x = g$ring0 * cos(theta), y = g$ring0 * sin(theta), grp = "r0"),
    data.frame(x = g$ring1 * cos(theta), y = g$ring1 * sin(theta), grp = "r1"))

  gg <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = built$src_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
      fill = "grey92", colour = "white", linewidth = 0.2) +
    ggplot2::geom_path(data = rings,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$grp),
      colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_polygon(data = built$poly_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                   fill = .data$value),
      colour = "white", linewidth = 0.12) +
    ggplot2::geom_text(data = src_lab_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                   angle = .data$angle),
      size = label_size, fontface = "bold", colour = "grey20") +
    ggplot2::coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25),
                         clip = "off") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40",
                                            size = 8.5),
      plot.margin = ggplot2::margin(8, 8, 8, 8))

  if (!is.null(built$lab_df)) {
    gg <- gg + ggplot2::geom_text(data = built$lab_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                   angle = .data$ang, hjust = .data$hjust),
      size = label_size - 0.5, colour = "grey30")
  }

  legname <- c(residuals = "z", prob = "P(to | from)", lift = "lift")[[fill]]
  if (diverging) {
    lim <- max(abs(built$poly_df$value[is.finite(built$poly_df$value)]),
               na.rm = TRUE)
    if (!is.finite(lim) || lim <= 0) lim <- 1
    gg <- gg + ggplot2::scale_fill_gradient2(
      low = .lsa_div_low, mid = .lsa_div_mid, high = .lsa_div_high,
      midpoint = 0, limits = c(-lim, lim), na.value = "grey85", name = legname)
  } else {
    gg <- gg + ggplot2::scale_fill_gradient(
      low = .lsa_seq_low, high = .lsa_seq_high, na.value = "grey85",
      name = legname)
  }

  fbase <- c(residuals = "adjusted residuals", prob = "transition probabilities",
             lift = "lift")[[fill]]
  sub <- if (style == "rose") {
    sprintf("%s engine, lag %d \u00b7 equal slots \u00b7 bar height = %s, fill = %s",
            fit$method, fit$params$lag, size, legname)
  } else {
    sprintf("%s engine, lag %d \u00b7 sector = source volume, wedge width = frequency, fill = %s",
            fit$method, fit$params$lag, legname)
  }
  gg + ggplot2::labs(title = sprintf("Transition sunburst \u2014 %s", fbase),
                     subtitle = sub)
}

# --- rose: equal sectors, equal slots, radial bar height = frequency ------
.polar_rose <- function(fit, states, K, fillmat, sig, size, significant,
                        labels, g) {
  sizemat <- switch(size, count = fit$obs, prob = fit$prob)
  gmax <- max(sizemat[is.finite(sizemat) & sizemat > 0], na.rm = TRUE)
  if (!is.finite(gmax) || gmax <= 0) {
    stop("Fit has no transitions to draw.", call. = FALSE)
  }
  gap <- 0.05
  sector_sz <- (2 * pi - gap * K) / K
  sector_start <- pi / 2 - (seq_len(K) - 1L) * (sector_sz + gap)
  sector_mid <- sector_start - sector_sz / 2
  to_r <- function(v) g$ring0 + sqrt(pmin(v / gmax, 1)) * (g$ring1 - g$ring0)

  pad <- sector_sz * 0.04
  slot_w <- (sector_sz - 2 * pad) / K
  bars <- list(); src <- list(); labs <- list(); id <- 0L
  for (i in seq_len(K)) {
    id <- id + 1L
    src[[i]] <- .lsa_wedge(sector_start[i] - sector_sz, sector_start[i],
                           g$src_in, g$src_out, id, NA_real_)
    modal <- which.max(sizemat[i, ])
    for (j in seq_len(K)) {
      v <- sizemat[i, j]
      if (!is.finite(v) || v <= 0) next
      c_ang <- sector_start[i] - pad - (j - 0.5) * slot_w
      id <- id + 1L
      val <- fillmat[i, j]
      if (isTRUE(significant) && !(is.finite(sig[i, j]) && sig[i, j])) {
        val <- NA_real_
      }
      bars[[id]] <- .lsa_wedge(c_ang - slot_w * 0.42, c_ang + slot_w * 0.42,
                               g$ring0, to_r(v), id, val)
      if (labels == "all" || (labels == "auto" && j == modal)) {
        labs[[id]] <- .polar_label(c_ang, states[j])
      }
    }
  }
  list(poly_df = do.call(rbind, bars), src_df = do.call(rbind, src),
       lab_df = if (length(labs)) do.call(rbind, labs) else NULL,
       sector_mid = sector_mid)
}

# --- wedge: sectors by volume, angular width = frequency share -------------
.polar_wedge <- function(fit, states, K, fillmat, sig, significant, labels,
                         min_show, g) {
  obs <- fit$obs
  out_vol <- rowSums(obs)
  out_vol[!is.finite(out_vol) | out_vol < 0] <- 0
  if (sum(out_vol) <= 0) stop("Fit has no outgoing transitions to draw.",
                              call. = FALSE)
  gap <- 0.05
  sector_sz <- out_vol / sum(out_vol) * (2 * pi - gap * K)
  sector_start <- numeric(K); sector_start[1L] <- pi / 2
  for (i in seq_len(K - 1L)) {
    sector_start[i + 1L] <- sector_start[i] - sector_sz[i] - gap
  }
  sector_mid <- sector_start - sector_sz / 2
  min_lab <- 0.05                          # ~2.9 deg: wedge must fit a name

  wedges <- list(); src <- list(); labs <- list(); id <- 0L
  for (i in seq_len(K)) {
    row_ct <- obs[i, ]; row_ct[!is.finite(row_ct) | row_ct < 0] <- 0
    tot <- sum(row_ct)
    if (tot <= 0) next
    id <- id + 1L
    src[[i]] <- .lsa_wedge(sector_start[i] - sector_sz[i], sector_start[i],
                           g$src_in, g$src_out, id, NA_real_)
    shares <- row_ct / tot
    keep <- which(shares >= min_show & shares > 0)
    if (!length(keep)) next
    usable <- sector_sz[i] - 2 * (sector_sz[i] * 0.02)
    widths <- shares[keep] / sum(shares[keep]) * usable
    cursor <- sector_start[i] - sector_sz[i] * 0.02
    for (k in seq_along(keep)) {
      j <- keep[k]; a_hi <- cursor; a_lo <- cursor - widths[k]; cursor <- a_lo
      id <- id + 1L
      val <- fillmat[i, j]
      if (isTRUE(significant) && !(is.finite(sig[i, j]) && sig[i, j])) {
        val <- NA_real_
      }
      wedges[[id]] <- .lsa_wedge(a_lo, a_hi, g$ring0, g$ring1, id, val)
      if (labels == "all" || (labels == "auto" && widths[k] >= min_lab)) {
        labs[[id]] <- .polar_label((a_lo + a_hi) / 2, states[j])
      }
    }
  }
  list(poly_df = do.call(rbind, wedges), src_df = do.call(rbind, src),
       lab_df = if (length(labs)) do.call(rbind, labs) else NULL,
       sector_mid = sector_mid)
}

# Radial outer label at angle `a` for target `lab`.
.polar_label <- function(a, lab) {
  data.frame(x = 1.03 * cos(a), y = 1.03 * sin(a), label = lab,
             ang = ifelse(cos(a) < 0, a * 180 / pi + 180, a * 180 / pi),
             hjust = ifelse(cos(a) < 0, 1, 0), stringsAsFactors = FALSE)
}
