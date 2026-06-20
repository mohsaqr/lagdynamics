# Tests for the four non-classical engines: two_cell, bidirectional,
# parallel_dominance, nonparallel_dominance.
#
# Each engine is tested for:
#   1. Existence in the registry and callable through lsa()
#   2. Internal consistency (symmetry / antisymmetry / sign agreement)
#   3. Hand-computed values on a small example
#   4. Cross-check against base-R primitives where applicable

# Shared small example: 2x2 with known cell values
LITTLE <- matrix(c(
  10, 20,
  30, 40
), 2, 2, byrow = TRUE, dimnames = list(c("a","b"), c("a","b")))

# A 3x3 example with no zero cells
THREE <- matrix(c(
  12,  8,  5,
   4, 15,  6,
   3,  7, 20
), 3, 3, byrow = TRUE, dimnames = list(c("x","y","z"), c("x","y","z")))

# --- Registry ----------------------------------------------------------

test_that("five engines are registered after package load", {
  reg <- list_lsa_engines()
  expect_setequal(reg$name,
                  c("classical", "two_cell", "bidirectional",
                    "parallel_dominance", "nonparallel_dominance"))
})

test_that("convenience wrappers dispatch to the right engine", {
  expect_equal(lsa_two_cell(LITTLE)$method, "two_cell")
  expect_equal(lsa_bidirectional(LITTLE)$method, "bidirectional")
  expect_equal(lsa_parallel_dominance(LITTLE)$method, "parallel_dominance")
  expect_equal(lsa_nonparallel_dominance(LITTLE)$method,
               "nonparallel_dominance")
})

# --- two_cell ---------------------------------------------------------

test_that("two_cell: hand-computed odds ratio on a 2x2 input", {
  # Feed a 2x2 transition matrix directly; the 2x2 collapse for cell
  # [1, 1] reproduces the original cells exactly:
  #   a = obs[1,1] = 10
  #   b = R[1] - a = 30 - 10 = 20
  #   c = C[1] - a = 40 - 10 = 30
  #   d = N - R[1] - C[1] + a = 100 - 30 - 40 + 10 = 40
  # OR = (a*d) / (b*c) = (10*40) / (20*30) = 400/600 = 0.6667
  fit <- lsa(LITTLE, engine = "two_cell")
  expect_equal(unname(fit$meta$extra$odds_ratio["a", "a"]),
               (10 * 40) / (20 * 30), tolerance = 1e-12)
  # log(OR) Wald SE = sqrt(1/10 + 1/20 + 1/30 + 1/40) = sqrt(0.1 + 0.05 + 0.0333 + 0.025) = sqrt(0.2083)
  expected_se <- sqrt(1/10 + 1/20 + 1/30 + 1/40)
  expect_equal(unname(fit$meta$extra$log_or_se["a", "a"]),
               expected_se, tolerance = 1e-12)
})

test_that("two_cell: zero-cell continuity correction kicks in", {
  m <- matrix(c(0, 5, 8, 12), 2, 2,
              dimnames = list(c("a","b"), c("a","b")))
  fit <- lsa(m, engine = "two_cell")
  # With a=0, default continuity=0.5 should produce a finite log_or
  expect_true(is.finite(fit$meta$extra$odds_ratio["a", "a"]))
  expect_true(is.finite(fit$meta$extra$log_or["a", "a"]))
})

# Base-R equivalence checks (pnorm / chisq.test / binom.test agreement)
# and the oconnor_couple oracle live in
# equivalence_testing/test-base-r-engines.R, out of the shipped suite.

# --- bidirectional ----------------------------------------------------

test_that("bidirectional: residual matrix is symmetric", {
  fit <- lsa(THREE, engine = "bidirectional")
  expect_equal(unname(fit$adj_res),
               unname(t(fit$adj_res)),
               tolerance = 1e-12)
  expect_equal(unname(fit$p),
               unname(t(fit$p)),
               tolerance = 1e-12)
})

test_that("bidirectional: symmetric_obs = O + t(O)", {
  fit <- lsa(THREE, engine = "bidirectional")
  W <- fit$meta$extra$symmetric_obs
  expect_equal(unname(W), unname(THREE + t(THREE)))
})

# --- parallel_dominance ----------------------------------------------

test_that("parallel_dominance: residual matrix is antisymmetric", {
  fit <- lsa(THREE, engine = "parallel_dominance")
  expect_equal(unname(fit$adj_res),
               unname(-t(fit$adj_res)),
               tolerance = 1e-12)
})

test_that("parallel_dominance: diagonal residuals are zero", {
  fit <- lsa(THREE, engine = "parallel_dominance")
  expect_true(all(diag(fit$adj_res) == 0))
})

test_that("parallel_dominance: hand-computed z[1, 2] on THREE", {
  fit <- lsa(THREE, engine = "parallel_dominance")
  # D[x, y] = O[x, y] - O[y, x] = 8 - 4 = 4
  # E[x, y] + E[y, x]: need full table marginals
  R <- rowSums(THREE); C <- colSums(THREE); N <- sum(THREE)
  E <- outer(R, C) / N
  expected_z <- (THREE["x","y"] - THREE["y","x"]) / sqrt(E["x","y"] + E["y","x"])
  expect_equal(unname(fit$adj_res["x", "y"]), unname(expected_z),
               tolerance = 1e-12)
})

# --- nonparallel_dominance --------------------------------------------

test_that("nonparallel_dominance: residual matrix is antisymmetric", {
  fit <- lsa(THREE, engine = "nonparallel_dominance")
  expect_equal(unname(fit$adj_res),
               unname(-t(fit$adj_res)),
               tolerance = 1e-12)
})

test_that("nonparallel_dominance: hand-computed z[1, 2] on THREE", {
  fit <- lsa(THREE, engine = "nonparallel_dominance")
  # D[x, y] = 8 - 4 = 4; SE = sqrt(8 + 4) = sqrt(12) = 3.464
  # z = 4 / 3.464 = 1.155
  expected_z <- (THREE["x","y"] - THREE["y","x"]) /
                sqrt(THREE["x","y"] + THREE["y","x"])
  expect_equal(unname(fit$adj_res["x", "y"]), unname(expected_z),
               tolerance = 1e-12)
})

# nonparallel_dominance binom.test agreement and the oconnor_couple
# cross-engine oracle live in equivalence_testing/test-base-r-engines.R.
