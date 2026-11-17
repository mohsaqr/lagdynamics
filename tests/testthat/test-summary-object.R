# summary.lsa() returns an object with two tidy tables and prints only
# when shown; as.data.frame(s, what = ) reads either table.

test_that("summary() returns a summary.lsa object without printing", {
  fit <- lsa(engagement)
  expect_silent(s <- summary(fit))
  expect_s3_class(s, "summary.lsa")
  expect_named(s, c("overview", "nodes"))
})

test_that("the overview table is one row that matches the fit", {
  fit <- lsa(engagement)
  o <- as.data.frame(summary(fit), what = "overview")
  expect_equal(nrow(o), 1L)
  expect_identical(o$engine, "classical")
  expect_equal(o$n_sequences, 136L)
  expect_equal(o$n_transitions, sum(fit$obs))
  expect_equal(o$n_significant, sum(transitions(fit)$significant))
  expect_equal(o$lrx2, tests(fit)$statistic[tests(fit)$test == "lrx2"])
})

test_that("the nodes table is one row per state with shares summing to 1", {
  fit <- lsa(engagement)
  n <- as.data.frame(summary(fit), what = "nodes")
  expect_identical(names(n), c("state", "outgoing", "incoming",
                               "out_share", "in_share", "init_prob"))
  expect_identical(n$state, fit$data$labels)
  expect_equal(sum(n$out_share), 1)
  expect_equal(sum(n$in_share), 1)
  expect_equal(sum(n$init_prob), 1)
  # Counts agree with nodes(), the reading verb for states.
  expect_equal(n$outgoing, nodes(fit)$outgoing)
  expect_equal(n$incoming, nodes(fit)$incoming)
})

test_that("a transition-matrix fit summarises with NA sizes and init", {
  fit <- lsa(matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3))
  s <- summary(fit)
  expect_true(is.na(as.data.frame(s)$n_sequences))
  expect_true(all(is.na(as.data.frame(s, what = "nodes")$init_prob)))
  expect_output(print(s), "from a transition matrix")
})

test_that("print.summary.lsa prints both tables and returns invisibly", {
  s <- summary(lsa(engagement))
  expect_output(print(s), "State activity")
  expect_output(print(s), "likelihood-ratio test")
  printed <- capture.output(visibility <- withVisible(print(s)))
  expect_false(visibility$visible)
  expect_gt(length(printed), 5L)
})

test_that("as.data.frame(summary) rejects an unknown table name", {
  s <- summary(lsa(engagement))
  expect_error(as.data.frame(s, what = "matrices"), "should be one of")
})
