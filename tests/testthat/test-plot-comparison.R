# Coverage for the comparison plot paths not exercised by test-compare.R:
# the "heatmap" style (two-group and pairwise -> .plot_compare_heatmap),
# the count-valued barrel, the pairwise effect ranking, and the Bayesian
# object routed through the same plot methods.

# Three groups of {A,B,C} sequences with distinct A-> biases, so a
# pairwise comparison has real structure. Sequences are built with Reduce
# (Markov accumulate) rather than a loop.
.cmpplot_groups <- function(seed = 1L) {
  set.seed(seed)
  states <- c("A", "B", "C")
  gen <- function(n, bias) replicate(n, {
    L <- sample(8:14, 1L)
    Reduce(function(prev, .) {
      p <- rep(1, 3)
      if (prev == "A") p <- p + bias
      sample(states, 1L, prob = p)
    }, seq_len(L - 1L), accumulate = TRUE, init = sample(states, 1L))
  }, simplify = FALSE)
  list(seqs = c(gen(20, c(0, 4, 0)), gen(20, c(0, 0, 4)), gen(20, c(4, 0, 0))),
       group = rep(c("g1", "g2", "g3"), each = 20))
}

.two_group_cmp <- function(seed) {
  d <- .cmpplot_groups(seed)
  fit <- lsa(d$seqs[1:40], group = d$group[1:40], engine = "classical")
  compare_lsa(fit, R = 100L)
}

test_that("plot.lsa_comparison style = 'heatmap' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  cmp <- .two_group_cmp(2)
  expect_s3_class(cmp, "lsa_comparison")
  expect_s3_class(plot(cmp, style = "heatmap"), "ggplot")
})

test_that("barrel value = 'count' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  cmp <- .two_group_cmp(3)
  expect_s3_class(plot(cmp, value = "count"), "ggplot")
})

test_that("pairwise plot covers heatmap, effect rank and count", {
  skip_if_not_installed("ggplot2")
  d <- .cmpplot_groups(4)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 100L)
  expect_s3_class(cmp, "lsa_comparison_pairwise")
  expect_s3_class(plot(cmp, style = "heatmap"), "ggplot")  # faceted heatmap
  expect_s3_class(plot(cmp, rank = "effect"), "ggplot")    # pairwise effect rank
  expect_s3_class(plot(cmp, value = "count"), "ggplot")
})

test_that("bayes comparison routes through both plot styles", {
  skip_if_not_installed("ggplot2")
  d <- .cmpplot_groups(5)
  fit <- lsa(d$seqs[1:40], group = d$group[1:40], engine = "classical")
  bc <- bayes_compare_lsa(fit, draws = 1500, seed = 3)
  expect_s3_class(plot(bc), "ggplot")
  expect_s3_class(plot(bc, style = "heatmap"), "ggplot")
})
