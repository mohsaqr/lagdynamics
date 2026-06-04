# GSEQ-parity additions: arbitrary integer lags, multi-lag profiles,
# column-conditional probabilities, and the Pearson tablewise chi-square.

test_that("negative lag yields the predecessor (transpose) table", {
  seqs <- list(c("a", "b", "a", "c", "b", "a"), c("b", "a", "c", "a", "b"))
  fp <- lsa(seqs, lag = 1)
  fn <- lsa(seqs, lag = -1)
  expect_equal(unname(fn$obs), unname(t(fp$obs)))
})

test_that("lag 0 is the degenerate self-pairing diagonal", {
  seqs <- list(c("a", "b", "a", "c", "b", "a"), c("b", "a", "c", "a", "b"))
  f0 <- lsa(seqs, lag = 0)
  expect_equal(sum(f0$obs[lower.tri(f0$obs)]) +
                 sum(f0$obs[upper.tri(f0$obs)]), 0)            # off-diag all 0
  expect_equal(unname(diag(f0$obs)),
               unname(as.integer(table(factor(unlist(seqs),
                 levels = c("a", "b", "c"))))))               # diag = counts
})

test_that("lag accepts any single integer; matrix input still lag 1 only", {
  expect_s3_class(lsa(c("a", "b", "a", "c", "b"), lag = 2), "lsa")
  tm <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3)
  expect_error(lsa(tm, lag = 2), "lag = 1")
})

test_that("prob_col is the column-conditional P(from | to)", {
  fit <- lsa(engagement, engine = "classical")
  # Each estimable target column sums to 1.
  cs <- colSums(fit$prob_col, na.rm = TRUE)
  expect_equal(unname(cs), rep(1, ncol(fit$obs)))
  # prob_col[i,j] = obs[i,j] / colSums(obs)[j].
  ref <- sweep(fit$obs, 2L, colSums(fit$obs), "/")
  expect_equal(unname(fit$prob_col), unname(ref))
  expect_true("prob_col" %in% names(fit$edges))
})

test_that("Pearson x2 accompanies G2 with the same df", {
  fit <- lsa(engagement, engine = "classical")
  expect_type(fit$x2, "list")
  expect_equal(fit$x2$df, fit$lrx2$df)
  # Pearson statistic matches a direct computation on estimable cells.
  ok <- is.finite(fit$exp) & fit$exp > 0
  expect_equal(fit$x2$statistic,
               sum((fit$obs[ok] - fit$exp[ok])^2 / fit$exp[ok]))
  expect_true(fit$x2$p >= 0 && fit$x2$p <= 1)
})

test_that("lsa_lags builds a per-lag profile and stacks edges", {
  prof <- lsa_lags(engagement, lags = 1:3)
  expect_s3_class(prof, "lsa_lags")
  expect_length(prof, 3L)
  expect_identical(attr(prof, "lags"), 1:3)
  expect_true(all(vapply(prof, inherits, logical(1L), "lsa")))
  d <- as.data.frame(prof)
  expect_equal(nrow(d), 9L * 3L)                  # 3x3 cells x 3 lags
  expect_setequal(unique(d$lag), 1:3)
  expect_output(print(prof), "lsa_lags")
})

test_that("lag_profile returns a tidy one-row-per-lag frame for a transition", {
  p <- lag_profile(engagement, "Active", "Average", lags = 1:3)
  expect_s3_class(p, "data.frame")
  expect_equal(nrow(p), 3L)
  expect_identical(p$lag, 1:3)
  expect_setequal(names(p),
    c("lag", "from", "to", "count", "prob", "adj_res", "p", "significant"))
  expect_true(all(p$from == "Active" & p$to == "Average"))
  expect_identical(rownames(p), as.character(1:3))   # clean row names
  # Accepts an existing lsa_lags object too.
  prof <- lsa_lags(engagement, lags = 1:2)
  expect_equal(nrow(lag_profile(prof, "Active", "Average")), 2L)
  # Unknown transition errors clearly.
  expect_error(lag_profile(engagement, "Active", "Nope"), "not found")
})
