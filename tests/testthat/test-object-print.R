# The lsa object's initial-state distribution slot and the Nestimate-
# style print.

test_that("lsa fit carries an inits slot that sums to 1", {
  fit <- lsa(engagement, engine = "classical")
  expect_type(fit$inits, "double")
  expect_named(fit$inits, fit$data$labels)
  expect_equal(sum(fit$inits), 1)
  # inits = proportion of sequences starting in each state.
  per <- split(fit$data$events,
               factor(fit$data$seq_id,
                      levels = seq_len(fit$data$n_sequences)))
  first <- vapply(per, function(s) s[1L], integer(1L))
  expect_equal(unname(fit$inits),
               unname(as.numeric(table(factor(first,
                 levels = seq_len(fit$data$n_states))) / length(first))))
})

test_that("transition-matrix fits have NULL inits", {
  fit <- lsa(matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3))
  expect_null(fit$inits)
})

test_that("print.lsa shows header, independence test, and inits bars", {
  fit <- lsa(engagement, engine = "classical")
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Lag Sequential Analysis")
  expect_match(out, "independence: G")
  expect_match(out, "Initial states:")
  expect_match(out, "█")   # the bar glyph
  expect_identical(withVisible(print(fit))$visible, FALSE)
})

test_that("print.lsa omits the initial-states block for matrix fits", {
  fit <- lsa(matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3))
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_false(grepl("Initial states", out))
})

test_that("summary.lsa adds node-activity bars and the full matrices", {
  fit <- lsa(engagement, engine = "classical")
  out <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(out, "Node activity")
  expect_match(out, "█")                    # bar glyph in the node section
  expect_match(out, "Adjusted residuals")   # full matrices still shown
  expect_identical(withVisible(summary(fit))$visible, FALSE)
})
