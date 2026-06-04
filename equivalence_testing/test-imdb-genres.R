# Validation on the IMDB primary-genre sequence dataset
# (1,000 highly-rated films, 1970-2024, K = 16). No published LSA
# result exists for this slice; the test cross-validates lagseq's
# classical engine against stats::chisq.test()$stdres (the
# authoritative Haberman residual oracle) on the same input.

test_that("imdb_genres data object is well-formed", {
  expect_type(imdb_genres$sequence, "character")
  expect_equal(length(imdb_genres$sequence), 1000L)
  expect_equal(imdb_genres$n_events, 1000L)
  expect_equal(imdb_genres$k_states, 16L)
  expect_equal(length(imdb_genres$alphabet), 16L)
  expect_true(all(imdb_genres$sequence %in% imdb_genres$alphabet))
})

test_that("imdb: lagseq classical engine fits without error", {
  fit <- lsa(imdb_genres$sequence, engine = "classical")
  expect_s3_class(fit, "lsa")
  expect_equal(nrow(fit$obs), 16L)
  expect_equal(sum(fit$obs), 999L)   # T - 1 transitions at lag 1
  expect_true(is.finite(fit$lrx2$statistic))
  expect_true(is.finite(fit$lrx2$p))
})

test_that("imdb: adj_res == chisq.test()$stdres at 1e-10", {
  fit <- lsa(imdb_genres$sequence, engine = "classical")
  suppressWarnings({
    ct <- stats::chisq.test(fit$obs, correct = FALSE)
  })
  ok <- !is.na(fit$adj_res) & !is.na(ct$stdres)
  expect_equal(unname(fit$adj_res[ok]), unname(ct$stdres[ok]),
               tolerance = 1e-10)
})

test_that("imdb: expected freqs == outer(R, C) / N at 1e-12", {
  fit <- lsa(imdb_genres$sequence, engine = "classical")
  R <- rowSums(fit$obs); C <- colSums(fit$obs); N <- sum(fit$obs)
  expect_equal(unname(fit$exp), unname(outer(R, C) / N),
               tolerance = 1e-12)
})

test_that("imdb: LR p-value == pchisq() at 1e-12", {
  fit <- lsa(imdb_genres$sequence, engine = "classical")
  expect_equal(fit$lrx2$p,
               stats::pchisq(fit$lrx2$statistic,
                             df = fit$lrx2$df,
                             lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("imdb: structural-zero variant zeros exp and NAs residuals on diagonal", {
  S <- 1 - diag(16)
  fit <- lsa(imdb_genres$sequence, engine = "classical",
             structural_zeros = S)
  expect_true(all(diag(fit$exp) == 0))
  # Forbidden cells are non-estimable: residuals are NA.
  expect_true(all(is.na(diag(fit$adj_res))))
})

test_that("imdb: large dataset doesn't break Yule's Q computation", {
  fit <- lsa(imdb_genres$sequence, engine = "classical")
  # All cells with a > 0 and finite marginals should yield finite Q in [-1, 1].
  ok <- is.finite(fit$yules_q)
  expect_true(all(fit$yules_q[ok] >= -1 - 1e-10 &
                  fit$yules_q[ok] <=  1 + 1e-10))
})
