# Property tests for iterative proportional fitting. These hold for
# *every* feasible input, not just hand-picked oracle pairs, so they
# pin algorithmic invariants that a future refactor could break.

set.seed(20260512)

# Small random table with diagonal zeros — the canonical lagseq IPF case.
make_random_obs <- function(K = 4L, max_count = 25L) {
  m <- matrix(sample.int(max_count, K * K, replace = TRUE), K, K)
  diag(m) <- 0L
  m
}

test_that("IPF marginals match observed marginals exactly", {
  for (trial in seq_len(20)) {
    obs <- make_random_obs(K = sample(3:6, 1))
    out <- lsa_ipf(obs)
    expect_equal(rowSums(out$fit), rowSums(obs), tolerance = 1e-7)
    expect_equal(colSums(out$fit), colSums(obs), tolerance = 1e-7)
  }
})

test_that("IPF preserves structural zeros exactly", {
  for (trial in seq_len(20)) {
    K <- sample(3:6, 1)
    obs <- make_random_obs(K = K)
    S <- 1 - diag(K)
    out <- lsa_ipf(obs, structure = S)
    expect_true(all(out$fit[S == 0] == 0))
  }
})

test_that("IPF converges in fewer iterations on small tables", {
  obs <- make_random_obs(K = 4L)
  out <- lsa_ipf(obs, tol = 1e-10)
  expect_true(out$converged)
  expect_lt(out$iterations, 200L)
})

test_that("IPF on a no-zero pattern reduces to outer(R, C) / N", {
  set.seed(7L)
  obs <- matrix(sample.int(30, 9), 3, 3)   # no zeros forced
  S <- matrix(1, 3, 3)                     # no structural zeros
  out <- lsa_ipf(obs, structure = S)
  R <- rowSums(obs)
  C <- colSums(obs)
  N <- sum(obs)
  expect_equal(out$fit, outer(R, C) / N, tolerance = 1e-7)
})

test_that("IPF errors when a row has positive total but no estimable cells", {
  obs <- matrix(c(0, 5, 0,   # row 1 has total 5, but if all S[1,] = 0...
                  3, 0, 2,
                  1, 4, 0), 3, 3, byrow = TRUE)
  S <- matrix(c(0, 0, 0,
                1, 0, 1,
                1, 1, 0), 3, 3, byrow = TRUE)
  expect_error(lsa_ipf(obs, structure = S),
               "row has positive marginal but no estimable cells")
})

test_that("IPF total mass equals N when feasible", {
  for (trial in seq_len(10)) {
    obs <- make_random_obs(K = sample(3:5, 1))
    out <- lsa_ipf(obs)
    expect_equal(sum(out$fit), sum(obs), tolerance = 1e-7)
  }
})
