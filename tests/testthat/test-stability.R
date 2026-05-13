test_that("stability_lsa returns the right class and shape", {
  set.seed(21L)
  fit <- lsa(engagement, engine = "classical")
  st <- stability_lsa(fit, R = 30)
  expect_s3_class(st, "lsa_stability")
  expect_equal(nrow(st$edges), 9L)
  expect_equal(dim(st$stability_matrix), c(30L, 9L))
})

test_that("stability_lsa: stability is in [0, 1]", {
  set.seed(22L)
  fit <- lsa(engagement, engine = "classical")
  st <- stability_lsa(fit, R = 30)
  expect_true(all(st$edges$stability >= 0))
  expect_true(all(st$edges$stability <= 1))
})

test_that("stability_lsa: stable flag matches min_stable threshold", {
  set.seed(23L)
  fit <- lsa(engagement, engine = "classical")
  st <- stability_lsa(fit, R = 30, min_stable = 0.8)
  expect_equal(st$edges$stable, st$edges$stability >= 0.8)
})

test_that("stability_lsa: edges with strongly significant observed stay stable", {
  set.seed(24L)
  fit <- lsa(engagement, engine = "classical")
  st <- stability_lsa(fit, R = 50)
  # Cells with extremely small observed p (< 1e-6) should be stable
  # under 80% subsampling.
  obs_p <- as.vector(fit$p)
  very_sig <- which(is.finite(obs_p) & obs_p < 1e-6)
  for (i in very_sig) {
    expect_gte(st$edges$stability[i], 0.7,
               label = sprintf("cell %d very-sig stability", i))
  }
})

test_that("stability_lsa errors on transition-matrix input", {
  obs <- matrix(c(0,3,1, 2,0,4, 5,1,0), 3, 3,
                dimnames = list(c("a","b","c"), c("a","b","c")))
  fit <- lsa(obs, engine = "classical")
  expect_error(stability_lsa(fit), "requires event-level input")
})

test_that("stability_lsa: as.data.frame returns the edges frame", {
  set.seed(25L)
  fit <- lsa(engagement, engine = "classical")
  st <- stability_lsa(fit, R = 20)
  expect_identical(as.data.frame(st), st$edges)
})
