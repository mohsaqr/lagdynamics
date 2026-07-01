# transition_probabilities() exposes the row-stochastic P(to | from)
# matrix natively (the quantity a Transition Network Analysis reads),
# paired with initial() for the initial-state probabilities.

test_that("transition_probabilities() returns a square, labelled matrix", {
  fit <- lsa(engagement)
  tp <- transition_probabilities(fit)
  expect_true(is.matrix(tp))
  expect_equal(nrow(tp), ncol(tp))
  expect_identical(rownames(tp), colnames(tp))
  expect_setequal(rownames(tp), fit$data$labels)
})

test_that("each row of the transition matrix sums to 1 (no structural zeros)", {
  fit <- lsa(engagement)
  tp <- transition_probabilities(fit)
  expect_equal(unname(rowSums(tp)), rep(1, nrow(tp)), tolerance = 1e-8)
})

test_that("transition_probabilities() matches the fit's prob matrix", {
  fit <- lsa(engagement)
  expect_equal(transition_probabilities(fit), fit$prob)
})

test_that("initial() returns the initial-state probabilities (init P)", {
  fit <- lsa(engagement)
  ini <- initial(fit)
  expect_named(ini, c("state", "init_prob"))
  expect_equal(sum(ini$init_prob), 1, tolerance = 1e-8)
})
