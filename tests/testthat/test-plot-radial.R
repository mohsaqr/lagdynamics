# plot_polar() (residual sunburst) and plot_forest() (circular bootstrap
# CI forest). We check they build/return ggplots and print without error.

test_that("plot_polar builds a sunburst for each fill metric", {
  skip_if_not_installed("ggplot2")
  fit <- lsa(engagement, engine = "classical")
  for (f in c("residuals", "prob", "lift")) {
    expect_s3_class(plot_polar(fit, fill = f), "ggplot")
  }
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 600, height = 600)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  expect_no_error(print(plot_polar(fit)))
  expect_s3_class(plot_polar(fit, significant = TRUE), "ggplot")
})

test_that("plot_polar tolerates structural-zero cells and bad input", {
  skip_if_not_installed("ggplot2")
  fit <- lsa(engagement, structural_zeros = 1 - diag(3))
  expect_s3_class(plot_polar(fit), "ggplot")        # diagonal wedges absent
  expect_error(plot_polar(list()), "lsa")
})

test_that("plot_forest builds a radial bootstrap forest", {
  skip_if_not_installed("ggplot2")
  fit <- lsa(engagement, engine = "classical")
  b <- bootstrap_lsa(fit, R = 60)
  for (m in c("residuals", "count", "prob", "yules_q")) {
    expect_s3_class(plot_forest(b, metric = m), "ggplot")
  }
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 600, height = 600)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  expect_no_error(print(plot_forest(b)))
  expect_s3_class(plot_forest(b, n_top = 4), "ggplot")
})

test_that("plot_forest honours show_nonsig and validates input", {
  skip_if_not_installed("ggplot2")
  fit <- lsa(engagement, engine = "classical")
  b <- bootstrap_lsa(fit, R = 60)
  expect_s3_class(plot_forest(b, show_nonsig = FALSE), "ggplot")
  expect_error(plot_forest(fit), "lsa_bootstrap")    # needs a bootstrap
})

test_that(".lsa_true coerces NA logicals to FALSE", {
  expect_identical(lagseq:::.lsa_true(c(TRUE, FALSE, NA)),
                   c(TRUE, FALSE, FALSE))
})

test_that("plot_polar supports both styles and label modes", {
  skip_if_not_installed("ggplot2")
  fit <- lsa(engagement, engine = "classical")
  expect_s3_class(plot_polar(fit, style = "rose"), "ggplot")
  expect_s3_class(plot_polar(fit, style = "wedge"), "ggplot")
  expect_s3_class(plot_polar(fit, style = "wedge", labels = "all"), "ggplot")
  expect_s3_class(plot_polar(fit, labels = "none"), "ggplot")
  expect_s3_class(plot_polar(fit, style = "wedge", min_show = 0.1), "ggplot")
})

test_that("plot_forest validates n_top", {
  skip_if_not_installed("ggplot2")
  b <- bootstrap_lsa(lsa(engagement, engine = "classical"), R = 30)
  expect_error(plot_forest(b, n_top = 0))
  expect_error(plot_forest(b, n_top = NA_real_))
  expect_s3_class(plot_forest(b, n_top = 3), "ggplot")
})
