test_that("significant_transitions returns rows with p < alpha", {
  fit <- lsa(engagement, engine = "classical")
  sig <- significant_transitions(fit, alpha = 0.05)
  expect_s3_class(sig, "data.frame")
  expect_identical(colnames(sig), colnames(fit$edges))
  expect_true(all(sig$p < 0.05))
  expect_true(all(is.finite(sig$p)))
})

test_that("significant_transitions defaults to the fit's alpha", {
  fit <- lsa(engagement, engine = "classical", alpha = 0.01)
  sig_default <- significant_transitions(fit)
  sig_explicit <- significant_transitions(fit, alpha = 0.01)
  expect_identical(sig_default, sig_explicit)
  expect_true(all(sig_default$p < 0.01))
})

test_that("significant_transitions excludes NA p-values", {
  # Structural-zero fit produces NA p-values on forbidden cells.
  # Those rows must not appear in the significant set.
  K <- 3
  sz <- 1 - diag(K)            # forbid self-transitions
  fit <- lsa(engagement, engine = "classical", structural_zeros = sz)
  sig <- significant_transitions(fit)
  expect_false(any(is.na(sig$p)))
})

test_that("over/under partition significant transitions by sign", {
  fit <- lsa(engagement, engine = "classical")
  sig  <- significant_transitions(fit)
  over <- overrepresented_transitions(fit)
  under <- underrepresented_transitions(fit)
  # Every over row has adj_res > 0; every under row has adj_res < 0.
  expect_true(all(over$adj_res > 0))
  expect_true(all(under$adj_res < 0))
  # Disjoint and exhaustive (modulo zero residuals, which are rare on
  # real data but should still be handled).
  expect_equal(nrow(over) + nrow(under),
               sum(is.finite(sig$adj_res) & sig$adj_res != 0))
})

test_that("common_transitions filters by minimum observed count", {
  fit <- lsa(engagement, engine = "classical")
  c1 <- common_transitions(fit, min_count = 1L)
  c5 <- common_transitions(fit, min_count = 5L)
  expect_true(all(c1$count >= 1L))
  expect_true(all(c5$count >= 5L))
  expect_true(nrow(c5) <= nrow(c1))
})

test_that("filter helpers reject non-lsa input", {
  expect_error(significant_transitions(list(edges = data.frame())),
               class = "simpleError")
  expect_error(common_transitions("not a fit"),
               class = "simpleError")
})

test_that("significant_transitions rejects out-of-range alpha", {
  fit <- lsa(engagement, engine = "classical")
  expect_error(significant_transitions(fit, alpha = 0))
  expect_error(significant_transitions(fit, alpha = 1))
  expect_error(significant_transitions(fit, alpha = -0.1))
  expect_error(significant_transitions(fit, alpha = c(0.05, 0.01)))
})

test_that("common_transitions rejects min_count < 1", {
  fit <- lsa(engagement, engine = "classical")
  expect_error(common_transitions(fit, min_count = 0))
  expect_error(common_transitions(fit, min_count = -1))
})
