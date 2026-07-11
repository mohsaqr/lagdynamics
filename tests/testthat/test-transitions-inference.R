# transitions() must read every inference result, not only a fit. These
# guard the tidy contract the vignettes rely on: one verb, named arguments,
# a tidy one-row-per-transition data.frame.

make_fit <- function() lsa(engagement)
make_gfit <- function() {
  lsa(engagement, group = rep(c("a", "b"), length.out = 136))
}

inference_results <- function() {
  set.seed(1)
  fit <- make_fit()
  gfit <- make_gfit()
  list(
    bootstrap   = bootstrap_lsa(fit, R = 30),
    certainty   = certainty_lsa(fit),
    stability   = stability_lsa(fit, R = 30),
    permutation = permute_lsa(fit, R = 30),
    comparison  = compare_lsa(gfit),
    bayes       = bayes_compare_lsa(gfit),
    lags        = lsa_lags(engagement, lags = 1:2)
  )
}

test_that("transitions() returns a tidy data.frame for every inference result", {
  results <- inference_results()
  lapply(names(results), function(nm) {
    out <- transitions(results[[nm]])
    expect_s3_class(out, "data.frame")
    expect_true(all(c("from", "to") %in% names(out)), info = nm)
    expect_gt(nrow(out), 0L)
    expect_identical(rownames(out), as.character(seq_len(nrow(out))), info = nm)
  })
})

test_that("transitions() agrees with the as.data.frame() view", {
  results <- inference_results()
  lapply(names(results), function(nm) {
    expect_equal(transitions(results[[nm]]), as.data.frame(results[[nm]]),
                 ignore_attr = TRUE, info = nm)
  })
})

test_that("significant = TRUE keeps only the flagged transitions", {
  results <- inference_results()
  flags <- c(bootstrap = "adj_res_stable", certainty = "stable",
             stability = "stable", permutation = "significant",
             comparison = "significant", bayes = "significant",
             lags = "significant")
  lapply(names(flags), function(nm) {
    kept <- transitions(results[[nm]], significant = TRUE)
    all_tx <- transitions(results[[nm]])
    expect_lte(nrow(kept), nrow(all_tx))
    if (nrow(kept) > 0L) {
      expect_true(all(kept[[flags[[nm]]]]), info = nm)
    }
    expect_equal(nrow(kept), sum(all_tx[[flags[[nm]]]], na.rm = TRUE),
                 info = nm)
  })
})

test_that("sort = 'strength' orders by descending absolute effect", {
  results <- inference_results()
  strengths <- c(bootstrap = "adj_res_observed", certainty = "adj_res_observed",
                 stability = "stability", permutation = "observed_adj_res",
                 comparison = "diff", bayes = "diff")
  lapply(names(strengths), function(nm) {
    out <- transitions(results[[nm]], sort = "strength")
    key <- abs(out[[strengths[[nm]]]])
    expect_false(is.unsorted(rev(key), na.rm = TRUE), info = nm)
  })
})

test_that("transitions() on a fit is unchanged by the widened generic", {
  fit <- make_fit()
  expect_s3_class(transitions(fit), "data.frame")
  expect_lte(nrow(transitions(fit, significant = TRUE)), nrow(transitions(fit)))
  expect_true(all(transitions(fit, direction = "over")$adj_res > 0))
  expect_true(all(transitions(fit, direction = "under")$adj_res < 0))
  expect_true(all(transitions(fit, min_count = 100)$count >= 100))
})

test_that("transitions() on lsa_lags stacks the per-lag fits", {
  lg <- lsa_lags(engagement, lags = 1:2)
  out <- transitions(lg)
  expect_setequal(unique(out$lag), c(1L, 2L))
  expect_equal(nrow(out), nrow(transitions(lsa(engagement))) * 2L)
})

test_that("significant must be a single non-missing logical", {
  cert <- certainty_lsa(make_fit())
  expect_error(transitions(cert, significant = NA))
  expect_error(transitions(cert, significant = c(TRUE, FALSE)))
})
