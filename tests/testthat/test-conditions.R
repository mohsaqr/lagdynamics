# User-facing errors carry a specific lagdynamics_* class plus the umbrella
# class lagdynamics_error. Tests assert the class, not the message text.

test_that("every classed error also inherits lagdynamics_error", {
  expect_error(lsa(engagement, engine = "no_such_engine"),
               class = "lagdynamics_error")
})

test_that("an unknown engine raises lagdynamics_unknown_engine", {
  expect_error(lsa(engagement, engine = "no_such_engine"),
               class = "lagdynamics_unknown_engine")
  expect_error(unregister_lsa_engine("no_such_engine"),
               class = "lagdynamics_unknown_engine")
})

test_that("resampling a transition-matrix fit raises lagdynamics_needs_events", {
  mfit <- lsa(matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3))
  expect_error(bootstrap_lsa(mfit, R = 2), class = "lagdynamics_needs_events")
  expect_error(permute_lsa(mfit, R = 2), class = "lagdynamics_needs_events")
  expect_error(reliability_lsa(mfit, R = 2),
               class = "lagdynamics_needs_events")
  expect_error(stability_lsa(mfit, R = 2), class = "lagdynamics_needs_events")
})

test_that("data with no transitions raises lagdynamics_no_transitions", {
  expect_error(lsa(list("a", "b", "c")),
               class = "lagdynamics_no_transitions")
  expect_error(lsa(list("a", "b", "c"), engine = "two_cell"),
               class = "lagdynamics_no_transitions")
})

test_that("structural zeros on a non-classical engine are unsupported", {
  expect_error(lsa(engagement, engine = "two_cell", loops = FALSE),
               class = "lagdynamics_unsupported_option")
  expect_error(lsa(engagement, engine = "bidirectional", loops = FALSE),
               class = "lagdynamics_unsupported_option")
})

test_that("comparing fits at different lags raises lagdynamics_incompatible_fits", {
  expect_error(compare_lsa(lsa(engagement, lag = 1), lsa(engagement, lag = 2),
                           R = 2),
               class = "lagdynamics_incompatible_fits")
})

test_that("invalid input values raise lagdynamics_bad_input", {
  expect_error(lsa(c(1.5, 2, 1)), class = "lagdynamics_bad_input")
  expect_error(lsa(engagement, structural_zeros = matrix(1, 2, 2)),
               class = "lagdynamics_bad_input")
})

test_that("argument contracts use named stopifnot messages", {
  expect_error(lsa(engagement, alpha = 2), "`alpha` must be a single number")
  expect_error(bootstrap_lsa(lsa(engagement), R = 0),
               "`R` must be a single whole number")
  expect_error(register_lsa_engine("x", fn = 1, description = "d"),
               "`fn` must be a function")
})

test_that("counts in structural-zero cells raise a classed warning", {
  # engagement has many self-transitions; loops = FALSE forbids them.
  expect_warning(lsa(engagement, loops = FALSE),
                 class = "lagdynamics_structural_zero_counts")
  # A table whose forbidden cells are empty fits without the warning.
  m <- matrix(c(0, 30, 12, 25, 0, 18, 9, 22, 0), 3, 3,
              dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  expect_no_warning(fit <- lsa(m, loops = FALSE))
  expect_gte(tests(fit)$statistic[1], 0)
})
