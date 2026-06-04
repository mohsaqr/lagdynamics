# plot_chords(): circular chord diagram, delegated to cograph::plot_chord().
# We check it draws without error and returns the chord/segment data
# invisibly (no visual assertion).

test_that("plot_chords draws a chord diagram for each colour metric", {
  skip_if_not_installed("cograph")
  fit <- lsa(engagement, engine = "classical")
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 600, height = 600)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  for (col in c("residuals", "prob", "count", "lift")) {
    out <- plot_chords(fit, color = col)
    expect_type(out, "list")
    expect_true(all(c("segments", "chords") %in% names(out)))
  }
})

test_that("plot_chords honours width, significant and self_loops", {
  skip_if_not_installed("cograph")
  fit <- lsa(engagement, engine = "classical")
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 600, height = 600)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  expect_no_error(plot_chords(fit, width = "prob"))
  expect_no_error(plot_chords(fit, significant = TRUE))
  expect_no_error(plot_chords(fit, self_loops = FALSE))
})

test_that("plot_chords compares two fits and validates inputs", {
  skip_if_not_installed("cograph")
  g <- lsa(engagement, group = rep(c("a", "b"), length.out = 136))
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 600, height = 600)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  expect_no_error(plot_chords(g$a, compare = g$b))
  expect_no_error(plot_chords(g$a, compare = g$b, color = "prob"))
  expect_error(plot_chords(list(), color = "count"), "lsa")
  other <- lsa(matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3,
                      dimnames = list(c("x", "y", "z"), c("x", "y", "z"))))
  expect_error(plot_chords(g$a, compare = other), "same states")
})

test_that(".lsa_chord_colors centres diverging metrics and greys NA", {
  cols <- lagseq:::.lsa_chord_colors(c(-2, 0, 2), diverging = TRUE)
  expect_length(cols, 3L)
  expect_equal(toupper(cols[2L]), "#F7F7F7")        # zero -> mid (white)
  expect_false(cols[1L] == cols[3L])                # opposite signs differ
  expect_identical(lagseq:::.lsa_chord_colors(NA_real_, TRUE), "grey85")
})
