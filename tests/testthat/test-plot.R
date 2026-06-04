# Heatmap (plot.lsa, ggplot2) and network (plot_transitions, cograph
# splot). We check they build/return without error (no visual assertion).

test_that("plot.lsa builds a heatmap for each matrix", {
  skip_if_not_installed("ggplot2")
  fit <- lsa(engagement, engine = "classical")
  for (w in c("residuals", "prob", "count", "expected")) {
    expect_s3_class(plot(fit, which = w), "ggplot")
  }
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 500, height = 500)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  expect_no_error(print(plot(fit)))
})

test_that("plot.lsa tolerates structural-zero (NA) cells", {
  skip_if_not_installed("ggplot2")
  fit <- lsa(engagement, structural_zeros = 1 - diag(3))
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 500, height = 500)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  expect_no_error(print(plot(fit)))
})

test_that("plot_transitions draws the network for each weight", {
  skip_if_not_installed("cograph")
  fit <- lsa(engagement, engine = "classical")
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 600, height = 500)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  for (w in c("count", "prob", "residuals", "lift")) {
    expect_no_error(plot_transitions(fit, weights = w))
  }
  expect_no_error(plot_transitions(fit, weights = "residuals",
                                   significant = TRUE))
  expect_error(plot_transitions(list(), weights = "count"), "lsa")
})

test_that("plot(fit, type=) dispatches to each view", {
  fit <- lsa(engagement, engine = "classical")
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 600, height = 600)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(fit), "ggplot")                       # heatmap default
  expect_s3_class(plot(fit, which = "prob"), "ggplot")       # heatmap arg
  expect_s3_class(plot(fit, type = "sunburst"), "ggplot")
  if (requireNamespace("cograph", quietly = TRUE)) {
    expect_no_error(plot(fit, type = "network"))
    expect_no_error(plot(fit, type = "chord"))
  }
  expect_error(plot(fit, type = "nope"))
})

test_that("plot() on a bootstrap draws the forest", {
  skip_if_not_installed("ggplot2")
  b <- bootstrap_lsa(lsa(engagement, engine = "classical"), R = 40)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 600, height = 600)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  expect_s3_class(plot(b), "ggplot")
})

test_that("plot() on a grouped fit draws every group", {
  skip_if_not_installed("ggplot2")
  gf <- lsa(engagement, group = rep(c("a", "b"), length.out = 136))
  expect_s3_class(gf, "lsa_group")
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 900, height = 500)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  # ggplot types return a list of per-group plots
  ph <- plot(gf, type = "heatmap")
  expect_length(ph, length(gf))
  expect_s3_class(ph[[1]], "ggplot")
  expect_s3_class(plot(gf, type = "sunburst")[[1]], "ggplot")
  if (requireNamespace("cograph", quietly = TRUE)) {
    expect_no_error(plot(gf, type = "network"))
    expect_no_error(plot(gf, type = "chord"))
  }
})
