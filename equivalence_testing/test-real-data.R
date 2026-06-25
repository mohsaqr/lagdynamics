# Multi-oracle validation on published, real-world categorical sequence
# datasets. The strategy:
#   1. Use lagdynamics::engagement (shipped) — 138 students x 15 weeks, K = 3.
#   2. If Nestimate is installed, also validate against
#      group_regulation_long (K = 9, ~27,500 events). This sweeps the
#      multi-K / long-format / long-sequence path.
#
# All oracles are base-R primitives (chisq.test, loglin, pchisq, pnorm).
# No prior LSA package is used as oracle.

test_that("engagement: wide-matrix input produces correct counts", {
  fit <- lsa(engagement, engine = "classical")
  # Hand check: total transitions should equal the number of
  # (week_i, week_{i+1}) cells where neither is NA.
  expected_n <- 0L
  for (r in seq_len(nrow(engagement))) {
    row_vals <- engagement[r, ]
    valid_idx <- which(!is.na(row_vals))
    if (length(valid_idx) < 2L) next
    # Consecutive valid positions
    pairs_n <- sum(diff(valid_idx) == 1L)
    expected_n <- expected_n + pairs_n
  }
  expect_equal(sum(fit$obs), expected_n)
})

test_that("engagement: adjusted residuals match chisq.test()$stdres", {
  fit <- lsa(engagement, engine = "classical")
  suppressWarnings({
    ct <- stats::chisq.test(fit$obs, correct = FALSE)
  })
  expect_equal(unname(fit$adj_res), unname(ct$stdres),
               tolerance = 1e-10)
})

test_that("engagement: expected freqs match outer(R, C) / N", {
  fit <- lsa(engagement, engine = "classical")
  R <- rowSums(fit$obs)
  C <- colSums(fit$obs)
  N <- sum(fit$obs)
  expect_equal(unname(fit$exp), unname(outer(R, C) / N),
               tolerance = 1e-12)
})

test_that("engagement: LR p-value matches pchisq()", {
  fit <- lsa(engagement, engine = "classical")
  expect_equal(
    fit$lrx2$p,
    stats::pchisq(fit$lrx2$statistic, df = fit$lrx2$df,
                  lower.tail = FALSE),
    tolerance = 1e-12
  )
})

test_that("engagement: structural-zero variant zeros exp and NAs residuals on diagonal", {
  S <- 1 - diag(3)
  fit <- lsa(engagement, engine = "classical", structural_zeros = S)
  expect_true(all(diag(fit$exp) == 0))
  # Forbidden cells are non-estimable: residuals and p-values are NA.
  expect_true(all(is.na(diag(fit$adj_res))))
  expect_true(all(is.na(diag(fit$p))))
})

test_that("engagement: edges frame nrows = K * K and signs are coherent", {
  fit <- lsa(engagement, engine = "classical")
  expect_equal(nrow(fit$edges), 3L * 3L)
  # The "sign" column should be `over` iff count > expected (when both
  # finite).
  for (i in seq_len(nrow(fit$edges))) {
    cnt <- fit$edges$count[i]
    exp_v <- fit$edges$expected[i]
    if (is.finite(cnt) && is.finite(exp_v) && cnt != exp_v) {
      expected_sign <- if (cnt > exp_v) "over" else "under"
      expect_equal(fit$edges$sign[i], expected_sign)
    }
  }
})

test_that("group_regulation_long: long-sequence path matches chisq.test", {
  skip_if_not_installed("Nestimate")
  e <- new.env()
  utils::data("group_regulation_long", package = "Nestimate", envir = e)
  df <- get("group_regulation_long", envir = e)
  # Treat the whole stream as one long sequence (ignore actor breaks for
  # this oracle test — same code path as a single learner).
  fit <- lsa(df$Action, engine = "classical")
  suppressWarnings({
    ct <- stats::chisq.test(fit$obs, correct = FALSE)
  })
  expect_equal(unname(fit$adj_res), unname(ct$stdres),
               tolerance = 1e-10)
})

test_that("group_regulation_long: per-actor list-of-sequences input", {
  skip_if_not_installed("Nestimate")
  e <- new.env()
  utils::data("group_regulation_long", package = "Nestimate", envir = e)
  df <- get("group_regulation_long", envir = e)
  # Split into one sequence per actor (canonical multi-sequence input).
  seqs <- split(df$Action, df$Actor)
  fit <- lsa(seqs, engine = "classical")
  # Lag-1 within-sequence: total transitions = events - sequences.
  expect_equal(sum(fit$obs),
               fit$data$n_events - fit$data$n_sequences)
  # Marginals are finite, totals consistent.
  expect_true(all(is.finite(fit$obs)))
  expect_gt(sum(fit$obs), 0)
})
