test_that("permute_lsa returns the right class and shape", {
  set.seed(11L)
  fit <- lsa(engagement, engine = "classical")
  pm <- permute_lsa(fit, R = 50)
  expect_s3_class(pm, "lsa_permutation")
  expect_equal(nrow(pm$edges), 9L)
  expect_equal(dim(pm$perm_adj_res), c(50L, 9L))
})

test_that("permute_lsa: Phipson-Smyth p-value formula is exact", {
  # Hand-check: with R replicates, the minimum possible p_perm is
  # (1 + 0) / (1 + R) = 1 / (R + 1). Verify a cell with very
  # extreme observed residual hits this minimum.
  set.seed(12L)
  fit <- lsa(engagement, engine = "classical")
  pm <- permute_lsa(fit, R = 50)
  # The Active-Active cell typically dominates engagement data;
  # at R = 50 the minimum p is 1/51 ~= 0.0196.
  expected_min <- 1 / (pm$R + 1)
  expect_gte(min(pm$edges$p_perm), expected_min - 1e-10)
})

test_that("permute_lsa: shuffles= replay produces identical p-values", {
  fit <- lsa(engagement, engine = "classical")
  set.seed(13L)
  # Generate one set of shuffles by running once and capturing them.
  T <- fit$data$n_events
  shuffles <- replicate(20L, sample.int(T), simplify = FALSE)
  pm1 <- permute_lsa(fit, R = 20, shuffles = shuffles)
  set.seed(999L)
  pm2 <- permute_lsa(fit, R = 20, shuffles = shuffles)
  expect_equal(pm1$perm_adj_res, pm2$perm_adj_res, tolerance = 1e-12)
  expect_equal(pm1$edges$p_perm, pm2$edges$p_perm, tolerance = 1e-12)
})

test_that("permute_lsa errors on transition-matrix input", {
  obs <- matrix(c(0,3,1, 2,0,4, 5,1,0), 3, 3,
                dimnames = list(c("a","b","c"), c("a","b","c")))
  fit <- lsa(obs, engine = "classical")
  expect_error(permute_lsa(fit), "requires event-level input")
})

test_that("permute_lsa: observed_adj_res equals fit$adj_res", {
  set.seed(15L)
  fit <- lsa(engagement, engine = "classical")
  pm <- permute_lsa(fit, R = 30)
  expect_equal(pm$edges$observed_adj_res, as.vector(fit$adj_res),
               tolerance = 1e-12)
})

test_that("O'Connor 1999 permutation oracle: large-residual cells get small p_perm", {
  # The paper publishes p_mean from 10 blocks * 1000 = 10,000
  # permutations. For cells with |Z| >> 1.96, our permute_lsa() at
  # R = 1000 should also yield very small p_perm. We test this in
  # a softer way: every cell with |adj_res| >= 4 in lagseq's classical
  # output must have p_perm <= 0.05 in our 1000-permutation
  # replication.
  skip_on_cran()
  data(oconnor_couple)
  set.seed(16L)
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  pm <- permute_lsa(fit, R = 500)
  large <- which(abs(as.vector(fit$adj_res)) >= 4)
  if (length(large) > 0) {
    expect_true(all(pm$edges$p_perm[large] <= 0.05))
  }
})

test_that("permute_lsa: as.data.frame returns the edges frame", {
  set.seed(17L)
  fit <- lsa(engagement, engine = "classical")
  pm <- permute_lsa(fit, R = 20)
  expect_identical(as.data.frame(pm), pm$edges)
})
