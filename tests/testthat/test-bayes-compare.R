# Two groups with a genuinely different A-> routing (shared state set).
.bayes_groups <- function(n_a, n_b, bias_a, bias_b, seed = 1L) {
  set.seed(seed)
  states <- c("A", "B", "C")
  gen <- function(n, bias) replicate(n, {
    L <- sample(8:14, 1L); s <- character(L); s[1L] <- sample(states, 1L)
    for (i in 2:L) { p <- rep(1, 3); if (s[i - 1L] == "A") p <- p + bias
                     s[i] <- sample(states, 1L, prob = p) }
    s
  }, simplify = FALSE)
  list(seqs = c(gen(n_a, bias_a), gen(n_b, bias_b)),
       group = rep(c("g1", "g2"), c(n_a, n_b)))
}

test_that("bayes_compare_lsa returns the documented shape and class", {
  d <- .bayes_groups(40, 40, c(0, 5, 0), c(0, 0, 5), seed = 1L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  bc <- bayes_compare_lsa(fit, draws = 2000, seed = 1)

  expect_s3_class(bc, "lsa_bayes")
  expect_s3_class(bc, "lsa_comparison")          # drop-in for the plots
  expect_identical(nrow(bc$edges), 9L)           # K = 3 -> 9 edges
  expect_true(all(c("from", "to", "prob_a", "prob_b", "diff", "ci_low",
                    "ci_high", "pd", "effect_size", "p_value", "p_adj",
                    "significant") %in% names(bc$edges)))
  expect_identical(as.data.frame(bc), bc$edges)  # inherits as.data.frame
  expect_identical(bc$measure, "prob")
})

test_that("posterior is well-formed: diff in CI, pd in [.5,1], p in [0,1]", {
  d <- .bayes_groups(40, 40, c(0, 5, 0), c(0, 0, 5), seed = 2L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  e <- bayes_compare_lsa(fit, draws = 3000, seed = 2)$edges
  ok <- is.finite(e$diff)
  # The closed-form posterior mean difference lies inside the MC interval.
  expect_true(all(e$ci_low[ok] <= e$diff[ok] + 1e-2))
  expect_true(all(e$diff[ok] <= e$ci_high[ok] + 1e-2))
  expect_true(all(e$pd[ok] >= 0.5 - 1e-9 & e$pd[ok] <= 1 + 1e-9))
  pp <- e$p_value[ok]
  expect_true(all(pp >= 0 & pp <= 1))
  # diff is exactly prob_a - prob_b (posterior means).
  expect_equal(e$diff, e$prob_a - e$prob_b, tolerance = 1e-12)
})

test_that("credibly-different edges have a CI that excludes zero", {
  d <- .bayes_groups(45, 45, c(0, 6, 0), c(0, 0, 6), seed = 3L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  e <- bayes_compare_lsa(fit, draws = 3000, seed = 3)$edges
  sig <- e$significant & !is.na(e$significant)
  expect_true(all(e$ci_low[sig] > 0 | e$ci_high[sig] < 0))
  # The engineered A->B / A->C edges are detected as credibly different.
  ab <- e[e$from == "A" & e$to == "B", ]
  expect_true(ab$significant)
})

test_that("direction agrees with the permutation compare_lsa", {
  d <- .bayes_groups(45, 45, c(0, 6, 0), c(0, 0, 6), seed = 4L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  set.seed(4); bc <- bayes_compare_lsa(fit, draws = 2000)
  set.seed(4); cp <- compare_lsa(fit, R = 300, measure = "prob")
  m <- merge(bc$edges[, c("from", "to", "diff")],
             cp$edges[, c("from", "to", "diff")],
             by = c("from", "to"), suffixes = c("_b", "_p"))
  ok <- is.finite(m$diff_b) & is.finite(m$diff_p) & m$diff_p != 0
  expect_gt(mean(sign(m$diff_b[ok]) == sign(m$diff_p[ok])), 0.8)
})

test_that("seed makes the credible intervals reproducible", {
  d <- .bayes_groups(30, 30, c(0, 4, 0), c(0, 0, 4), seed = 5L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  a <- bayes_compare_lsa(fit, draws = 1500, seed = 7)
  b <- bayes_compare_lsa(fit, draws = 1500, seed = 7)
  expect_equal(a$edges$ci_low, b$edges$ci_low, tolerance = 1e-12)
})

test_that("more than two groups runs all pairwise comparisons", {
  set.seed(6L)
  d3 <- list(
    seqs = c(.bayes_groups(30, 30, c(0, 6, 0), c(0, 0, 6), seed = 6L)$seqs,
             .bayes_groups(30, 0, c(0, 0, 0), c(0, 0, 0), seed = 7L)$seqs),
    group = rep(c("AB", "AC", "flat"), each = 30))
  fit <- lsa(d3$seqs, group = d3$group, engine = "classical")
  bc <- bayes_compare_lsa(fit, draws = 1500, seed = 1)
  expect_s3_class(bc, "lsa_bayes_pairwise")
  expect_identical(nrow(bc$edges), 27L)          # 3 pairs x 9 edges
  expect_identical(nrow(bc$global), 3L)
  expect_identical(names(bc$comparisons),
                   c("AB_vs_AC", "AB_vs_flat", "AC_vs_flat"))
})

test_that("plot reuses the comparison barrel", {
  skip_if_not_installed("ggplot2")
  d <- .bayes_groups(30, 30, c(0, 5, 0), c(0, 0, 5), seed = 8L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  expect_s3_class(plot(bayes_compare_lsa(fit, draws = 1500, seed = 8)),
                  "ggplot")
})
