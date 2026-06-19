# A homogeneous multi-sequence fixture: the analytic Dirichlet certainty
# and the sequence bootstrap should agree closely here.
.cert_fit <- function(seed = 5L) {
  set.seed(seed)
  events <- sample(c("a", "b", "c"), 900L, replace = TRUE,
                   prob = c(0.5, 0.3, 0.2))
  seqs <- split(events, rep(seq_len(45L), each = 20L))
  lsa(seqs, engine = "classical")
}

test_that("certainty_lsa returns the documented shape and class", {
  cert <- certainty_lsa(.cert_fit())
  expect_s3_class(cert, "lsa_certainty")
  expect_s3_class(cert, "lsa_bootstrap")          # drop-in compatible
  expect_identical(nrow(cert$edges), 9L)          # K = 3 -> 9 edges
  expect_true(all(c("from", "to", "prob_observed", "prob_mean", "prob_se",
                    "prob_ci_low", "prob_ci_high", "p_value", "stable",
                    "adj_res_observed", "adj_res_stable") %in%
                  names(cert$edges)))
  # as.data.frame inherits the bootstrap method.
  expect_identical(as.data.frame(cert), cert$edges)
  expect_identical(cert$prior, 0.5)
  expect_true(is.na(cert$R))                       # no iterations
})

test_that("posterior is a valid Beta: ci_low <= mean <= ci_high, near observed", {
  e <- certainty_lsa(.cert_fit())$edges
  ok <- is.finite(e$prob_mean)
  expect_true(all(e$prob_ci_low[ok] <= e$prob_mean[ok] + 1e-9))
  expect_true(all(e$prob_mean[ok] <= e$prob_ci_high[ok] + 1e-9))
  # With a weak prior the posterior mean sits close to the observed prob.
  expect_lt(max(abs(e$prob_mean[ok] - e$prob_observed[ok])), 0.05)
})

test_that("analytic CIs agree with the sequence bootstrap on homogeneous data", {
  fit <- .cert_fit()
  cert <- certainty_lsa(fit)
  set.seed(1); bs <- bootstrap_lsa(fit, R = 1000L)
  m <- merge(cert$edges[, c("from", "to", "prob_ci_low", "prob_ci_high")],
             bs$edges[, c("from", "to", "prob_ci_low", "prob_ci_high")],
             by = c("from", "to"), suffixes = c("_c", "_b"))
  expect_lt(mean(abs(m$prob_ci_low_c  - m$prob_ci_low_b),  na.rm = TRUE), 0.03)
  expect_lt(mean(abs(m$prob_ci_high_c - m$prob_ci_high_b), na.rm = TRUE), 0.03)
})

test_that("loops = FALSE makes the diagonal non-estimable and never certain", {
  set.seed(9L)
  seqs <- split(sample(c("a", "b", "c"), 600L, replace = TRUE),
                rep(seq_len(30L), each = 20L))
  cert <- certainty_lsa(lsa(seqs, loops = FALSE))
  diag_rows <- cert$edges$from == cert$edges$to
  expect_true(all(is.na(cert$edges$prob_mean[diag_rows])))
  expect_false(any(cert$edges$stable[diag_rows]))
})

test_that("a stronger prior shrinks edges toward uniform (more regularised)", {
  fit <- .cert_fit()
  weak   <- certainty_lsa(fit, prior = 0.5)$edges
  strong <- certainty_lsa(fit, prior = 50)$edges
  # The dominant self-ish edge is pulled down toward 1/K by the strong prior.
  i <- which.max(weak$prob_observed)
  expect_lt(strong$prob_mean[i], weak$prob_mean[i])
})

test_that("threshold inference and grouped fits work", {
  fit <- .cert_fit()
  ct <- certainty_lsa(fit, inference = "threshold")
  expect_identical(ct$inference, "threshold")
  expect_true(is.finite(ct$edge_threshold))

  g <- lsa(split(sample(c("a", "b"), 400L, replace = TRUE),
                 rep(seq_len(20L), each = 20L)),
           group = rep(c("x", "y"), each = 10L), engine = "classical")
  cg <- certainty_lsa(g)
  expect_s3_class(cg, "lsa_certainty_group")
  expect_identical(names(cg), c("x", "y"))
  expect_s3_class(cg[["x"]], "lsa_certainty")
})

test_that("plot.lsa_certainty returns a ggplot", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(certainty_lsa(.cert_fit())), "ggplot")
})
