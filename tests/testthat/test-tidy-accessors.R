# Reading verbs: transitions()/nodes()/tests()/initial() + summary().

test_that("nodes() is a tidy state frame keyed by state name", {
  fit <- lsa(engagement, engine = "classical")
  n <- nodes(fit)
  expect_s3_class(n, "data.frame")
  expect_setequal(names(n), c("state", "outgoing", "incoming"))
  expect_false(any(c("id", "node_id", "name", "label") %in% names(n)))
  # State names match the from/to endpoints of transitions().
  expect_setequal(n$state, unique(c(transitions(fit)$from, transitions(fit)$to)))
})

test_that("tests() gives one tidy row per tablewise test", {
  fit <- lsa(engagement, engine = "classical")
  t <- tests(fit)
  expect_setequal(names(t), c("test", "statistic", "df", "p"))
  expect_setequal(t$test, c("lrx2", "x2"))
  expect_true(all(t$p >= 0 & t$p <= 1))
})

test_that("initial() is a tidy state/init_prob frame; empty for matrix fits", {
  il <- initial(lsa(engagement, engine = "classical"))
  expect_setequal(names(il), c("state", "init_prob"))
  expect_equal(sum(il$init_prob), 1)
  mfit <- lsa(matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3))
  expect_equal(nrow(initial(mfit)), 0L)
})

test_that("summary(fit) returns a tidy one-row frame and still prints", {
  fit <- lsa(engagement, engine = "classical")
  expect_output(s <- summary(fit), "Lag Sequential Analysis")  # prints
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 1L)
  expect_true(all(c("engine", "n_transitions", "n_significant",
                    "lrx2", "x2") %in% names(s)))
  expect_identical(withVisible(summary(fit))$visible, FALSE)
})

test_that("summary() is one tidy row per fit, one row per group", {
  fit <- lsa(engagement, engine = "classical")
  expect_equal(nrow(summary(fit)), 1L)
  g <- lsa(engagement, group = rep(c("a", "b"), length.out = 136))
  fs <- summary(g)                       # grouped summary; no fit_summary()
  expect_s3_class(fs, "data.frame")
  expect_equal(nrow(fs), 2L)
  expect_identical(names(fs)[1L], "group")
})
