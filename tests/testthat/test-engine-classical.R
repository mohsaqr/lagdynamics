# Unit tests for the classical engine against hand-computed values.

test_that("classical engine: 2-state hand example", {
  # Sequence: a b a b a b a b a   (length 9, lag 1 → 8 transitions)
  # All pairs alternate (a,b) or (b,a):
  #   a -> b : 4
  #   b -> a : 4
  d  <- lsa_data(c("a", "b", "a", "b", "a", "b", "a", "b", "a"))
  tx <- lsa_transitions(d)
  expect_equal(unname(tx$obs["a", "b"]), 4L)
  expect_equal(unname(tx$obs["b", "a"]), 4L)
  expect_equal(sum(tx$obs), 8L)

  fit <- lsa(d, engine = "classical")

  # Row totals: a = 4, b = 4. Column totals: a = 4, b = 4. N = 8.
  # E[a,a] = 4 * 4 / 8 = 2.
  expect_equal(unname(fit$exp["a", "a"]), 2)
  expect_equal(unname(fit$exp["a", "b"]), 2)
  expect_equal(unname(fit$exp["b", "a"]), 2)
  expect_equal(unname(fit$exp["b", "b"]), 2)

  # Transitional probabilities: a always goes to b => 1, never to a => 0.
  expect_equal(unname(fit$prob["a", "b"]), 1)
  expect_equal(unname(fit$prob["a", "a"]), 0)
  expect_equal(unname(fit$prob["b", "a"]), 1)
  expect_equal(unname(fit$prob["b", "b"]), 0)
})

test_that("classical engine: adjusted residuals match Haberman formula", {
  d  <- lsa_data(c("a", "b", "a", "b", "a", "b", "a", "b", "a"))
  fit <- lsa(d, engine = "classical")

  # Hand: z[a,a] = (0 - 2) / sqrt(2 * (1 - 4/8) * (1 - 4/8))
  #              = -2 / sqrt(2 * 0.5 * 0.5)
  #              = -2 / sqrt(0.5)
  #              = -2.828427
  expect_equal(unname(fit$adj_res["a", "a"]), -2 / sqrt(0.5),
               tolerance = 1e-10)
  expect_equal(unname(fit$adj_res["a", "b"]),  2 / sqrt(0.5),
               tolerance = 1e-10)
})

test_that("classical engine: Yule's Q hand check", {
  d  <- lsa_data(c("a", "b", "a", "b", "a", "b", "a", "b", "a"))
  fit <- lsa(d, engine = "classical")
  # 2x2 collapse for [a, a]:
  #   a = O[a,a] = 0; b = R[a] - 0 = 4; c = C[a] - 0 = 4;
  #   d = N - R[a] - C[a] + 0 = 8 - 4 - 4 + 0 = 0
  # Q = (a*d - b*c)/(a*d + b*c) = (0 - 16)/(0 + 16) = -1
  expect_equal(unname(fit$yules_q["a", "a"]), -1)
  expect_equal(unname(fit$yules_q["a", "b"]),  1)
})

test_that("classical engine: LR p-value matches pchisq", {
  d  <- lsa_data(c("a", "b", "a", "b", "a", "b", "a", "b", "a"))
  fit <- lsa(d, engine = "classical")
  # df for K=2 with no structural zeros: (K-1)^2 = 1.
  expect_equal(fit$lrx2$df, 1)
  expect_equal(
    fit$lrx2$p,
    stats::pchisq(fit$lrx2$statistic, df = 1, lower.tail = FALSE),
    tolerance = 1e-12
  )
})

test_that("classical engine: edges frame is consistent with matrices", {
  d   <- lsa_data(c("a", "b", "a", "c", "b", "a", "c", "a", "b"))
  fit <- lsa(d, engine = "classical")
  for (i in seq_len(nrow(fit$edges))) {
    from_i <- fit$edges$from[i]
    to_i   <- fit$edges$to[i]
    expect_identical(fit$edges$count[i], unname(fit$obs[from_i, to_i]))
    expect_equal(fit$edges$expected[i], unname(fit$exp[from_i, to_i]),
                 tolerance = 1e-12)
    expect_equal(fit$edges$adj_res[i], unname(fit$adj_res[from_i, to_i]),
                 tolerance = 1e-12)
  }
})

test_that("classical engine: structural zeros via IPF stay zero", {
  d <- lsa_data(c("a", "b", "a", "c", "b", "a", "c", "b", "a", "c"))
  S <- 1 - diag(3)
  fit <- lsa(d, engine = "classical", structural_zeros = S)
  expect_true(all(fit$exp[S == 0] == 0))
  # Marginals of E equal marginals of O.
  expect_equal(rowSums(fit$exp), rowSums(fit$obs), tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_equal(colSums(fit$exp), colSums(fit$obs), tolerance = 1e-6,
               ignore_attr = TRUE)
  # IPF metadata is captured.
  expect_true(fit$meta$ipf$converged)
})
