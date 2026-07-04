test_that("the Dynalytics evidence surface is available from one fit", {
  fit <- lsa(engagement)

  expect_s3_class(fit, "lsa")
  expect_s3_class(transitions(fit), "data.frame")
  expect_s3_class(nodes(fit), "data.frame")
  expect_s3_class(tests(fit), "data.frame")
  expect_s3_class(initial(fit), "data.frame")
  fit_summary <- NULL
  capture.output(fit_summary <- summary(fit))
  expect_s3_class(fit_summary, "data.frame")
  expect_true(any(transitions(fit)$significant))

  expect_s3_class(certainty_lsa(fit), "lsa_certainty")
  expect_s3_class(bootstrap_lsa(fit, R = 4, seed = 1), "lsa_bootstrap")
  expect_s3_class(reliability_lsa(fit, R = 4, seed = 1), "lsa_reliability")
  expect_s3_class(stability_lsa(fit, R = 4, seed = 1), "lsa_stability")
  expect_s3_class(permute_lsa(fit, R = 4, seed = 1), "lsa_permutation")
})

test_that("the Dynalytics group-comparison evidence surface is available", {
  seqs <- list(
    c("A", "B", "A", "B", "C"),
    c("A", "B", "C", "C", "B"),
    c("C", "B", "C", "A", "B"),
    c("C", "A", "C", "B", "A")
  )
  gfit <- lsa(seqs, group = c("high", "high", "low", "low"))

  expect_s3_class(gfit, "lsa_group")
  expect_s3_class(transitions(gfit), "data.frame")
  expect_s3_class(compare_lsa(gfit, R = 4, seed = 1), "lsa_comparison")
  expect_s3_class(bayes_compare_lsa(gfit, draws = 40, seed = 1),
                  "lsa_bayes")
})
