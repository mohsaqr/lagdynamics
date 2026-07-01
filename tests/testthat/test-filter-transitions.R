# The transitions() verb and its filter arguments.

test_that("transitions(fit) returns the tidy name-keyed edge frame", {
  fit <- lsa(engagement, engine = "classical")
  tr <- transitions(fit)
  expect_s3_class(tr, "data.frame")
  # Name-keyed endpoints, distinct per-cell stats, no protocol cruft.
  expect_identical(colnames(tr),
    c("from", "to", "lag", "count", "expected", "prob", "prob_col",
      "adj_res", "p", "yules_q", "kappa", "kappa_z", "kappa_p", "lift",
      "sign", "significant"))
  expect_false(any(c("from_label", "to_label", "weight", "edge") %in%
                     colnames(tr)))
  # from/to carry state names, matching nodes()$state.
  expect_type(tr$from, "character")
  expect_setequal(unique(c(tr$from, tr$to)), nodes(fit)$state)
  expect_equal(nrow(tr), nrow(fit$obs)^2)
  expect_identical(rownames(tr), as.character(seq_len(nrow(tr))))
})

test_that("transitions(sort=) orders rows strongest-first, columns intact", {
  fit <- lsa(engagement, engine = "classical")
  base <- transitions(fit)
  strg <- transitions(fit, sort = "strength")
  cnt  <- transitions(fit, sort = "count")
  # same rows and columns, only the order changes
  expect_identical(colnames(strg), colnames(base))
  expect_setequal(paste(strg$from, strg$to), paste(base$from, base$to))
  # monotone non-increasing on the sort key
  expect_false(is.unsorted(rev(abs(strg$adj_res[is.finite(strg$adj_res)]))))
  expect_false(is.unsorted(rev(cnt$count)))
  expect_error(transitions(fit, sort = "nonsense"))
})

test_that("transitions(significant = TRUE) keeps p < alpha", {
  fit <- lsa(engagement, engine = "classical")
  sig <- transitions(fit, significant = TRUE, alpha = 0.05)
  expect_true(all(sig$p < 0.05))
  expect_true(all(is.finite(sig$p)))
})

test_that("transitions defaults to the fit's alpha", {
  fit <- lsa(engagement, engine = "classical", alpha = 0.01)
  d <- transitions(fit, significant = TRUE)
  e <- transitions(fit, significant = TRUE, alpha = 0.01)
  expect_identical(d, e)
  expect_true(all(d$p < 0.01))
})

test_that("transitions(significant=) excludes NA p-values", {
  fit <- lsa(engagement, engine = "classical",
             structural_zeros = 1 - diag(3))
  expect_false(any(is.na(transitions(fit, significant = TRUE)$p)))
})

test_that("direction over/under partition significant transitions by sign", {
  fit <- lsa(engagement, engine = "classical")
  sig  <- transitions(fit, significant = TRUE)
  over <- transitions(fit, direction = "over")
  under <- transitions(fit, direction = "under")
  expect_true(all(over$adj_res > 0))
  expect_true(all(under$adj_res < 0))
  expect_equal(nrow(over) + nrow(under),
               sum(is.finite(sig$adj_res) & sig$adj_res != 0))
})

test_that("min_count filters by minimum observed count", {
  fit <- lsa(engagement, engine = "classical")
  c1 <- transitions(fit, min_count = 1L)
  c5 <- transitions(fit, min_count = 5L)
  expect_true(all(c1$count >= 1L))
  expect_true(all(c5$count >= 5L))
  expect_true(nrow(c5) <= nrow(c1))
})

test_that("transitions rejects out-of-range alpha", {
  fit <- lsa(engagement, engine = "classical")
  expect_error(transitions(fit, significant = TRUE, alpha = 0))
  expect_error(transitions(fit, significant = TRUE, alpha = 1))
  expect_error(transitions(fit, significant = TRUE, alpha = c(0.05, 0.01)))
})

test_that("transitions(min_count=) rejects NA / non-positive thresholds", {
  fit <- lsa(engagement, engine = "classical")
  expect_error(transitions(fit, min_count = NA_real_))   # was silent all-NA rows
  expect_error(transitions(fit, min_count = 0))
  expect_error(transitions(fit, min_count = -5))
  expect_s3_class(transitions(fit, min_count = 100), "data.frame")  # valid
})
