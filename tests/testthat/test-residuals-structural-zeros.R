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

test_that("structural-zero cells have non-estimable (NA) residuals", {
  # Cells forbidden by structural_zeros are NOT "exactly expected" --
  # they are non-estimable. The engine must encode that as NA in
  # adj_res (and the derived p-value), not 0 (which falsely yields
  # p = 1 from 2*pnorm(0) and lights up as "not significant" rather
  # than "untestable").
  set.seed(12L)
  seq300 <- sample(c("a", "b", "c", "d"), 300, replace = TRUE)
  S <- 1 - diag(4)              # forbid self-transitions
  fit <- lsa(seq300, engine = "classical", structural_zeros = S)
  expect_true(all(is.na(diag(fit$adj_res))))
  expect_true(all(diag(fit$exp) == 0))
  # p-values on those cells must also be NA, not 1.
  expect_true(all(is.na(diag(fit$p))))
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

test_that("structural-zero LR test uses correct df (non-degenerate)", {
  set.seed(16L)
  seq200 <- sample(c("a", "b", "c", "d"), 200, replace = TRUE)
  S <- 1 - diag(4)               # 4 structural zeros, no degeneracy
  fit <- lsa(seq200, engine = "classical", structural_zeros = S)
  K <- 4L
  s_zeros <- sum(S == 0)
  # For a non-degenerate pattern the rank-based formula reduces to
  # the textbook (K-1)^2 - s_zeros result.
  expect_equal(fit$lrx2$df, (K - 1L)^2 - s_zeros)
})

test_that(".lrx2 uses rank-based df under degenerate structural-zero patterns", {
  # Forbid the entire first row: state "a" is never a "from" state.
  # The row-"a" effect drops out of the design, so rank(X) loses one
  # column. The naive formula (K-1)^2 - s_zeros = 9 - 4 = 5 is wrong;
  # the correct df = n_estimable - rank(X) = 12 - 6 = 6.
  obs <- matrix(c(
    0, 0, 0, 0,
    5, 0, 4, 3,
    2, 6, 0, 4,
    3, 2, 4, 0
  ), 4L, 4L, byrow = TRUE,
  dimnames = list(c("a", "b", "c", "d"), c("a", "b", "c", "d")))
  S <- matrix(c(
    0, 0, 0, 0,
    1, 1, 1, 1,
    1, 1, 1, 1,
    1, 1, 1, 1
  ), 4L, 4L, byrow = TRUE)
  # The residual design matrix is singular under this pattern (row-a
  # effect is non-identifiable), so .adjusted_residuals warns and
  # falls back. The LR df is computed independently via .lrx2 and
  # uses the correct rank-based path.
  expect_warning(
    fit <- lsa(obs, engine = "classical", structural_zeros = S),
    "Design matrix is singular"
  )
  expect_equal(fit$lrx2$df, 6L)
  # Naive formula would give 5; the test confirms we are NOT using it.
  expect_false(identical(fit$lrx2$df, (4L - 1L)^2 - sum(S == 0)))
})

test_that(".lrx2 returns NA p when no cells are estimable", {
  # All-zero structural pattern: every cell forbidden. df is
  # undefined; p must be NA, not forced to df = 1.
  obs <- matrix(c(0, 0, 0,
                  0, 0, 0,
                  0, 0, 0), 3L, 3L, byrow = TRUE)
  # We can't construct this via lsa() because obs would have N == 0
  # and the engine errors before reaching .lrx2. Instead exercise
  # .lrx2 directly.
  S <- matrix(0L, 3L, 3L)
  res <- lagseq:::.lrx2(obs = matrix(1, 3, 3), exp_mat = matrix(1, 3, 3),
                       K = 3L, structural_zeros = S)
  expect_true(is.na(res$p))
})

test_that("non-classical engines reject structural_zeros explicitly", {
  S <- 1 - diag(3)
  for (eng in c("two_cell", "bidirectional",
                "parallel_dominance", "nonparallel_dominance")) {
    expect_error(
      lsa(engagement, engine = eng, structural_zeros = S),
      "does not support structural_zeros",
      info = eng
    )
  }
})
