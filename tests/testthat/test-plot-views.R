# Coverage for plot.lsa heatmap `which` variants and the probability-
# weighted (TNA-styled) network rendering, which other plot tests do not
# reach.

test_that("plot.lsa heatmap which = count / expected return ggplots", {
  skip_if_not_installed("ggplot2")
  fit <- lsa(engagement)
  expect_s3_class(plot(fit, which = "count"), "ggplot")
  expect_s3_class(plot(fit, which = "expected"), "ggplot")
})

test_that("probability/TNA/count-weighted network renders via cograph", {
  skip_if_not_installed("cograph")
  fit <- lsa(engagement)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # weights = "tna"/"relative"/"prob"/"count" are drawn in the TNA style.
  expect_no_error(plot(fit, type = "network", weights = "tna"))
  expect_no_error(plot(fit, type = "network", weights = "relative"))
  expect_no_error(plot(fit, type = "network", weights = "prob"))
  expect_no_error(plot(fit, type = "network", weights = "count"))
})
