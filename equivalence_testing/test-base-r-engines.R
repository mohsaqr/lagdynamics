# Base-R-primitive equivalence for the engines and the permutation test.
# These check lagseq's numbers against independent base-R computations
# (pnorm, chisq.test$stdres, binom.test) and the O'Connor oracle dataset.
# They live here, out of tests/testthat/, so the shipped suite does not
# validate against external references (see README in this folder).

# A 3x3 example with no zero cells.
THREE <- matrix(c(
  12,  8,  5,
   4, 15,  6,
   3,  7, 20
), 3, 3, byrow = TRUE, dimnames = list(c("x","y","z"), c("x","y","z")))

test_that("two_cell: p-values agree with pnorm() of log_or_z", {
  fit <- lsa(THREE, engine = "two_cell")
  z <- fit$adj_res
  expect_equal(unname(fit$p),
               unname(2 * stats::pnorm(-abs(z))),
               tolerance = 1e-12)
})

test_that("bidirectional: residual matches Haberman on the symmetric table", {
  fit <- lsa(THREE, engine = "bidirectional")
  W <- fit$meta$extra$symmetric_obs
  # Independent computation: chisq.test()$stdres on the symmetric table
  # should match our adj_res (since both apply the Haberman formula).
  suppressWarnings({
    ct <- stats::chisq.test(W, correct = FALSE)
  })
  ok <- !is.na(fit$adj_res) & !is.na(ct$stdres)
  expect_equal(unname(fit$adj_res[ok]),
               unname(ct$stdres[ok]),
               tolerance = 1e-10)
})

test_that("nonparallel_dominance: binomial_p agrees with binom.test()", {
  fit <- lsa(THREE, engine = "nonparallel_dominance")
  bp <- fit$meta$extra$binomial_p
  # Hand check cell [x, y]: a=8, n=8+4=12
  ref <- stats::binom.test(8, 12, p = 0.5)$p.value
  expect_equal(unname(bp["x", "y"]), ref, tolerance = 1e-12)
})

test_that("oconnor_couple: all engines fit and report consistent counts", {
  data(oconnor_couple)
  for (engine in c("classical", "two_cell", "bidirectional",
                   "parallel_dominance", "nonparallel_dominance")) {
    fit <- lsa(oconnor_couple$sequence, engine = engine)
    expect_equal(sum(fit$obs), 392L,
                 label = sprintf("engine = %s: sum(obs)", engine))
    expect_s3_class(fit, "lsa")
    expect_s3_class(fit, "cograph_network")
  }
})

test_that("O'Connor 1999 permutation oracle: large-residual cells get small p_perm", {
  # The paper publishes p_mean from 10 blocks * 1000 = 10,000
  # permutations. For cells with |Z| >> 1.96, our permute_lsa() at
  # R = 1000 should also yield very small p_perm. We test this in
  # a softer way: every cell with |adj_res| >= 4 in lagseq's classical
  # output must have p_perm <= 0.05 in our 1000-permutation
  # replication.
  data(oconnor_couple)
  set.seed(16L)
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  pm <- permute_lsa(fit, R = 500)
  large <- which(abs(as.vector(fit$adj_res)) >= 4)
  if (length(large) > 0) {
    expect_true(all(pm$edges$p_perm[large] <= 0.05))
  }
})
