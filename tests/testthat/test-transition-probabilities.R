# transition_probabilities() exposes the row-stochastic P(to | from) as a
# tidy one-row-per-pair frame, paired with initial() for the initial-state
# probabilities.

test_that("transition_probabilities() returns a tidy from/to/prob frame", {
  fit <- lsa(engagement)
  tp <- transition_probabilities(fit)
  expect_s3_class(tp, "data.frame")
  expect_false(is.matrix(tp))
  expect_identical(names(tp), c("from", "to", "prob"))
  K <- length(fit$data$labels)
  expect_equal(nrow(tp), K^2)
  expect_setequal(tp$from, fit$data$labels)
  expect_setequal(tp$to, fit$data$labels)
  expect_false(anyDuplicated(paste(tp$from, tp$to)) > 0)
})

test_that("probabilities out of each state sum to 1 (no structural zeros)", {
  fit <- lsa(engagement)
  tp <- transition_probabilities(fit)
  sums <- tapply(tp$prob, tp$from, sum)
  expect_equal(unname(as.numeric(sums)), rep(1, length(sums)),
               tolerance = 1e-8)
})

test_that("transition_probabilities() matches the fit's prob matrix cell by cell", {
  fit <- lsa(engagement)
  tp <- transition_probabilities(fit)
  ix <- cbind(match(tp$from, rownames(fit$prob)),
              match(tp$to, colnames(fit$prob)))
  expect_equal(tp$prob, fit$prob[ix])
})

test_that("hand-computed probabilities on a tiny sequence", {
  # a->b, b->a, a->b, b->b: from a: b = 2/2; from b: a = 1/2, b = 1/2.
  fit <- lsa(c("a", "b", "a", "b", "b"))
  tp <- transition_probabilities(fit)
  expect_equal(tp$prob, c(0, 1, 0.5, 0.5))
  expect_identical(tp$from, c("a", "a", "b", "b"))
  expect_identical(tp$to, c("a", "b", "a", "b"))
})

test_that("structural zeros are reported as probability 0", {
  fit <- lsa_forbidden_counts(engagement, loops = FALSE)
  tp <- transition_probabilities(fit)
  expect_true(all(tp$prob[tp$from == tp$to] == 0))
})

test_that("initial() returns the initial-state probabilities (init P)", {
  fit <- lsa(engagement)
  ini <- initial(fit)
  expect_named(ini, c("state", "init_prob"))
  expect_equal(sum(ini$init_prob), 1, tolerance = 1e-8)
})
