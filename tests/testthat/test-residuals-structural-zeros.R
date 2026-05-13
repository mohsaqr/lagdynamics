# Christensen (1997, p. 357) design-matrix residual tests.
#
# The most important test in this file is the *reduction property*: when
# structural_zeros = matrix(1, K, K) (i.e. no actual zeros), Christensen's
# hat-matrix formula must collapse to Haberman's closed-form residual
# from §4.1. This is an analytic identity, not an approximation, so the
# tolerance is tight (1e-10).

set.seed(20260512)

test_that("residuals with all-ones structure match no-structure residuals", {
  # Build a real LSA dataset, then fit with and without a trivial
  # structural_zeros = all-ones matrix. The hat-matrix engine path
  # should produce numerically identical residuals to the Haberman path
  # because the underlying math is the same.
  set.seed(11L)
  seq200 <- sample(c("a", "b", "c", "d"), 200, replace = TRUE)
  fit_haberman   <- lsa(seq200, engine = "classical")
  fit_christensen <- lsa(seq200, engine = "classical",
                          structural_zeros = matrix(1, 4, 4))

  expect_equal(
    unname(fit_christensen$adj_res),
    unname(fit_haberman$adj_res),
    tolerance = 1e-10
  )
})

test_that("structural-zero cells have residual exactly zero", {
  set.seed(12L)
  seq300 <- sample(c("a", "b", "c", "d"), 300, replace = TRUE)
  S <- 1 - diag(4)              # forbid self-transitions
  fit <- lsa(seq300, engine = "classical", structural_zeros = S)
  expect_true(all(diag(fit$adj_res) == 0))
  expect_true(all(diag(fit$exp) == 0))
})

test_that("expected freqs under structural zeros match IPF directly", {
  # The engine's expected matrix should be the IPF fit. This pins the
  # contract: engine calls lsa_ipf() (not its own ad-hoc routine).
  set.seed(13L)
  seq300 <- sample(c("a", "b", "c", "d"), 300, replace = TRUE)
  S <- 1 - diag(4)
  fit <- lsa(seq300, engine = "classical", structural_zeros = S)
  ipf <- lsa_ipf(fit$obs, structure = S)
  expect_equal(unname(fit$exp), unname(ipf$fit), tolerance = 1e-8)
})

test_that("structural-zero marginals match observed marginals", {
  set.seed(14L)
  seq300 <- sample(c("a", "b", "c"), 300, replace = TRUE)
  S <- 1 - diag(3)
  fit <- lsa(seq300, engine = "classical", structural_zeros = S)
  expect_equal(rowSums(fit$exp), rowSums(fit$obs),
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(colSums(fit$exp), colSums(fit$obs),
               tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("structural-zero residuals have correct sign", {
  # Where O > E, the residual must be > 0; where O < E, residual < 0.
  set.seed(15L)
  seq500 <- sample(c("a", "b", "c", "d"), 500, replace = TRUE)
  S <- 1 - diag(4)
  fit <- lsa(seq500, engine = "classical", structural_zeros = S)
  # Compare signs in off-diagonal cells only.
  off <- S == 1
  signs_obs_minus_exp <- sign(fit$obs[off] - fit$exp[off])
  signs_residual      <- sign(fit$adj_res[off])
  # When residual is finite and nonzero, signs must agree.
  ok <- is.finite(signs_residual) & signs_residual != 0
  expect_true(all(signs_obs_minus_exp[ok] == signs_residual[ok]))
})

test_that("structural-zero LR test uses correct df", {
  set.seed(16L)
  seq200 <- sample(c("a", "b", "c", "d"), 200, replace = TRUE)
  S <- 1 - diag(4)               # 4 structural zeros
  fit <- lsa(seq200, engine = "classical", structural_zeros = S)
  K <- 4L
  s_zeros <- sum(S == 0)
  expect_equal(fit$lrx2$df, (K - 1L)^2 - s_zeros)
})
