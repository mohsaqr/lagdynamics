# Fit with structural zeros on data that contain transitions in the
# forbidden cells. Such a fit warns with lagdynamics_structural_zero_counts
# (the forbidden counts are left out of the quasi-independence fit). The tests that use this helper
# pin the behaviour of exactly such fits, so the warning is asserted here
# rather than hidden.
lsa_forbidden_counts <- function(...) {
  fit <- NULL
  expect_warning(fit <- lsa(...),
                 class = "lagdynamics_structural_zero_counts")
  fit
}
