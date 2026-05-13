test_that("bootstrap_lsa returns the right class and shape", {
  set.seed(1L)
  fit <- lsa(engagement, engine = "classical")
  bs <- bootstrap_lsa(fit, R = 50)
  expect_s3_class(bs, "lsa_bootstrap")
  expect_equal(nrow(bs$edges), 9L)         # K^2 = 9 for K=3
  expect_equal(dim(bs$boot_obs), c(50L, 9L))
  expect_equal(dim(bs$boot_adj_res), c(50L, 9L))
  expect_true("adj_res_stable" %in% names(bs$edges))
})

test_that("bootstrap_lsa reads recipe from fit$params", {
  set.seed(2L)
  fit <- lsa(engagement, engine = "classical", alternative = "greater")
  bs <- bootstrap_lsa(fit, R = 30)
  expect_equal(bs$fit$params$engine, "classical")
  expect_equal(bs$fit$params$alternative, "greater")
})

test_that("bootstrap_lsa CIs match stats::quantile() on the boot matrix", {
  set.seed(3L)
  fit <- lsa(engagement, engine = "classical")
  bs <- bootstrap_lsa(fit, R = 100, level_alpha = 0.95)
  # Pick cell 1: hand-verify the CI matches quantile()
  cell1 <- bs$boot_adj_res[, 1]
  expect_equal(unname(bs$edges$adj_res_ci_low[1]),
               unname(stats::quantile(cell1, 0.025, na.rm = TRUE,
                                      names = FALSE)),
               tolerance = 1e-10)
  expect_equal(unname(bs$edges$adj_res_ci_high[1]),
               unname(stats::quantile(cell1, 0.975, na.rm = TRUE,
                                      names = FALSE)),
               tolerance = 1e-10)
})

test_that("bootstrap_lsa indices= reproducibility hook", {
  fit <- lsa(engagement, engine = "classical")
  set.seed(4L)
  bs1 <- bootstrap_lsa(fit, R = 20)
  # Replay the exact same indices in a fresh session-state.
  set.seed(999L)   # different seed
  bs2 <- bootstrap_lsa(fit, R = 20, indices = bs1$indices_used)
  expect_equal(bs1$boot_obs, bs2$boot_obs)
  expect_equal(bs1$boot_adj_res, bs2$boot_adj_res, tolerance = 1e-12)
})

test_that("bootstrap_lsa falls back to event-level for single sequence", {
  set.seed(5L)
  seq1 <- sample(c("a","b","c"), 200, replace = TRUE)
  fit <- lsa(seq1, engine = "classical")
  expect_warning(
    bs <- bootstrap_lsa(fit, R = 30, level = "sequence"),
    "Only one sequence"
  )
  expect_equal(bs$level, "event")
})

test_that("bootstrap_lsa errors on transition-matrix input", {
  obs <- matrix(c(0,3,1, 2,0,4, 5,1,0), 3, 3,
                dimnames = list(c("a","b","c"), c("a","b","c")))
  fit <- lsa(obs, engine = "classical")
  expect_error(bootstrap_lsa(fit), "requires event-level input")
})

test_that("bootstrap_lsa: stable edges match observed sign", {
  set.seed(6L)
  fit <- lsa(engagement, engine = "classical")
  bs <- bootstrap_lsa(fit, R = 100)
  # For stable edges, the observed residual must lie within the CI.
  stable_rows <- which(bs$edges$adj_res_stable)
  for (i in stable_rows) {
    expect_true(
      bs$edges$adj_res_observed[i] >= bs$edges$adj_res_ci_low[i] - 1e-6 &&
      bs$edges$adj_res_observed[i] <= bs$edges$adj_res_ci_high[i] + 1e-6,
      label = sprintf("row %d observed within CI", i)
    )
  }
})

test_that("bootstrap_lsa: as.data.frame returns the edges frame", {
  set.seed(7L)
  fit <- lsa(engagement, engine = "classical")
  bs <- bootstrap_lsa(fit, R = 20)
  df <- as.data.frame(bs)
  expect_identical(df, bs$edges)
})
