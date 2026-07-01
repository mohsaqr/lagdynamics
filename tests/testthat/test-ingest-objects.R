# lsa() can duck-type sequence-bearing objects from sibling packages
# (tna / group_tna, Nestimate nestimate_data, TraMineR stslist) via
# inherits(), so a user who has those packages can fit directly from
# their objects. Those packages are not lagdynamics dependencies, so the
# round-trips are exercised by the sibling projects, not here. This file
# guards the one rule that must hold with no sibling package installed:
# object dispatch must not hijack ordinary list / vector input.

test_that("a plain list is treated as sequences, not misread as an object", {
  fit <- lsa(list(c("a", "b", "a"), c("b", "c", "a", "b")))
  expect_equal(fit$data$n_sequences, 2L)
  expect_setequal(fit$data$labels, c("a", "b", "c"))
})

test_that("a plain character vector is a single sequence", {
  fit <- lsa(c("a", "b", "a", "c", "b"))
  expect_s3_class(fit, "lsa")
  expect_equal(fit$data$n_sequences, 1L)
})
