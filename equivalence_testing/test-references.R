# Equivalence with base-R statistical primitives. These tests demonstrate
# that lagdynamics's classical engine produces the same numerical values as
# long-established R primitives (stats::chisq.test, stats::loglin,
# stats::pchisq), with no dependence on any prior LSA package.

set.seed(20260512)

# A reproducible random sequence with K = 4 codes, length 200.
seq200 <- sample(c("a", "b", "c", "d"), 200, replace = TRUE)

test_that("expected frequencies equal outer(R, C) / N", {
  fit <- lsa(seq200, engine = "classical")
  R   <- rowSums(fit$obs)
  C   <- colSums(fit$obs)
  N   <- sum(fit$obs)
  expected_hand <- outer(R, C) / N
  expect_equal(unname(fit$exp), unname(expected_hand),
               tolerance = 1e-12)
})

test_that("adjusted residuals equal chisq.test()$stdres (no structural zeros)", {
  fit <- lsa(seq200, engine = "classical")
  suppressWarnings({
    ct <- stats::chisq.test(fit$obs, correct = FALSE)
  })
  # chisq.test()$stdres uses the same Haberman formula as §4.1 of
  # inst/REFERENCES.md. Tolerance is tight because this is an analytic
  # identity, not an approximation.
  expect_equal(unname(fit$adj_res), unname(ct$stdres),
               tolerance = 1e-10)
})

test_that("LR p-value matches pchisq()", {
  fit <- lsa(seq200, engine = "classical")
  expect_equal(
    fit$lrx2$p,
    stats::pchisq(fit$lrx2$statistic, df = fit$lrx2$df,
                  lower.tail = FALSE),
    tolerance = 1e-12
  )
})

test_that("IPF matches stats::loglin() with structural zeros", {
  # Hand-built observed counts on 4 codes with a forbidden-diagonal
  # pattern. We compare lagdynamics's IPF against base-R's loglin() which
  # implements the same iterative-proportional-fitting algorithm via a
  # different code path.
  set.seed(99)
  obs <- matrix(sample.int(20, 16, replace = TRUE), 4, 4)
  diag(obs) <- 0L

  S <- 1 - diag(4)
  lagdynamics_ipf <- lsa_ipf(obs, structure = S)

  # loglin() with `start` initialized at the structural-zero pattern and
  # the same tight tolerance lagdynamics_ipf used (default loglin eps = 0.1
  # is far too loose for a numerical-equivalence test).
  base_fit <- suppressWarnings(
    stats::loglin(obs, margin = list(1, 2), start = S,
                  fit = TRUE, print = FALSE, eps = 1e-10)
  )

  # Both algorithms converge to the same MLE under the same model. With
  # matched tolerances they agree to the convergence precision.
  expect_equal(lagdynamics_ipf$fit, base_fit$fit, tolerance = 1e-6)

  # Marginals match observed exactly.
  expect_equal(rowSums(lagdynamics_ipf$fit), rowSums(obs), tolerance = 1e-6)
  expect_equal(colSums(lagdynamics_ipf$fit), colSums(obs), tolerance = 1e-6)
  # Structural zeros stay zero.
  expect_true(all(lagdynamics_ipf$fit[S == 0] == 0))
})

test_that("Yule's Q matches hand 2x2 collapse for every cell", {
  fit <- lsa(seq200, engine = "classical")
  obs <- fit$obs
  R   <- unname(rowSums(obs))
  C   <- unname(colSums(obs))
  N   <- sum(obs)
  K   <- nrow(obs)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      a <- unname(obs[i, j])
      b <- R[i] - a
      c <- C[j] - a
      d <- N - R[i] - C[j] + a
      num <- a * d - b * c
      den <- a * d + b * c
      q   <- if (den > 0) num / den else NA_real_
      expect_equal(unname(fit$yules_q[i, j]), q, tolerance = 1e-12,
                   label = sprintf("cell (%s, %s)",
                                   rownames(obs)[i], colnames(obs)[j]))
    }
  }
})

test_that("p-values match pnorm() under each alternative", {
  fit_two  <- lsa(seq200, engine = "classical", alternative = "two.sided")
  fit_grt  <- lsa(seq200, engine = "classical", alternative = "greater")
  fit_less <- lsa(seq200, engine = "classical", alternative = "less")

  expect_equal(unname(fit_two$p),
               unname(2 * stats::pnorm(-abs(fit_two$adj_res))),
               tolerance = 1e-12)
  expect_equal(unname(fit_grt$p),
               unname(stats::pnorm(fit_grt$adj_res, lower.tail = FALSE)),
               tolerance = 1e-12)
  expect_equal(unname(fit_less$p),
               unname(stats::pnorm(fit_less$adj_res)),
               tolerance = 1e-12)
})

test_that("transition matrix input + sequence input agree on identical counts", {
  fit_seq <- lsa(seq200, engine = "classical")
  fit_mat <- lsa(fit_seq$obs, engine = "classical")
  expect_equal(unname(fit_seq$obs),     unname(fit_mat$obs))
  expect_equal(unname(fit_seq$exp),     unname(fit_mat$exp),
               tolerance = 1e-12)
  expect_equal(unname(fit_seq$adj_res), unname(fit_mat$adj_res),
               tolerance = 1e-10)
  expect_equal(unname(fit_seq$yules_q), unname(fit_mat$yules_q),
               tolerance = 1e-12)
})
