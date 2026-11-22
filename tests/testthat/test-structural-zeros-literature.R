# Regressions for two literature fixes (2026-11-22).
#
# 1. Quasi-independence with structural zeros (Bakeman & Quera, 1995,
#    Psychological Bulletin 118, 272-284; Bishop, Fienberg & Holland, 1975,
#    ch. 5) is a model of the possible cells: its margins are the totals of
#    those cells. Counts the data put in a forbidden cell are left out. Before
#    the fix the margins kept them, and lsa(engagement, loops = FALSE)
#    reported G2 = -1184.633 with p = 1.
# 2. certainty_lsa() computes the upper tail of the Beta posterior directly
#    (pbeta(lower.tail = FALSE); DiDonato & Morris, 1992), not 1 - pbeta().

# A 3 x 3 table with self-transitions; loops = FALSE forbids the diagonal.
.qi_table <- function() {
  matrix(c(40, 12, 7,
           9, 55, 14,
           6, 16, 61), 3, 3, byrow = TRUE,
         dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
}

# The same model as a Poisson log-linear fit on the possible cells only.
.qi_glm <- function(obs) {
  cells <- as.data.frame(as.table(obs), stringsAsFactors = FALSE)
  names(cells) <- c("from", "to", "count")
  cells <- subset(cells, from != to)
  stats::glm(count ~ from + to, family = stats::poisson(), data = cells,
             control = stats::glm.control(epsilon = 1e-14, maxit = 100L))
}

test_that("lsa_ipf() reproduces the margins of the possible cells", {
  obs <- .qi_table()
  possible <- 1 - diag(3L)
  out <- lsa_ipf(obs, structure = possible, tol = 1e-12, max_iter = 1000L)
  expect_true(out$converged)
  expect_equal(unname(rowSums(out$fit)), unname(rowSums(obs * possible)),
               tolerance = 1e-10)
  expect_equal(unname(colSums(out$fit)), unname(colSums(obs * possible)),
               tolerance = 1e-10)
  expect_equal(unname(diag(out$fit)), c(0, 0, 0))
  # Invariant: counts in a forbidden cell do not change the fit.
  obs_without <- obs * possible
  expect_equal(out$fit, lsa_ipf(obs_without, structure = possible,
                                tol = 1e-12, max_iter = 1000L)$fit,
               tolerance = 1e-12)
})

test_that("G2 with loops = FALSE and self-transitions equals the glm quasi-independence deviance", {
  obs <- .qi_table()
  model <- .qi_glm(obs)
  fit <- lsa_forbidden_counts(obs, loops = FALSE)
  g2 <- subset(tests(fit), test == "lrx2")
  expect_equal(g2$statistic, stats::deviance(model), tolerance = 1e-6)
  expect_identical(as.integer(g2$df), 1L)
  expect_gte(g2$statistic, 0)
  expect_equal(g2$p, stats::pchisq(stats::deviance(model), df = 1,
                                   lower.tail = FALSE), tolerance = 1e-6)
})

test_that("engagement with loops = FALSE: G2 is 5.9937 on 1 df, not -1184.633", {
  fit <- lsa_forbidden_counts(engagement, loops = FALSE)
  g2 <- subset(tests(fit), test == "lrx2")
  # glm quasi-independence deviance on the six possible cells: 5.993710.
  expect_equal(g2$statistic, 5.993710, tolerance = 1e-6)
  expect_identical(as.integer(g2$df), 1L)
  expect_lt(g2$p, 0.05)
})

test_that("adjusted residuals use the possible-cell fit (equal glm standardized residuals)", {
  obs <- .qi_table()
  model <- .qi_glm(obs)
  fit <- lsa_forbidden_counts(obs, loops = FALSE)
  ours <- subset(transitions(fit), from != to)
  theirs <- data.frame(from = model$data$from, to = model$data$to,
                       expected = unname(stats::fitted(model)),
                       adj_res = unname(stats::rstandard(model, type = "pearson")))
  both <- merge(ours, theirs, by = c("from", "to"), suffixes = c("", "_glm"))
  expect_identical(nrow(both), 6L)
  expect_equal(both$expected, both$expected_glm, tolerance = 1e-6)
  expect_equal(both$adj_res, both$adj_res_glm, tolerance = 1e-6)
})

test_that("certainty_lsa() upper tail is pbeta(lower.tail = FALSE) to full relative accuracy", {
  fit <- lsa(engagement)
  prior <- 0.5
  cert <- transitions(certainty_lsa(fit, prior = prior, inference = "stability"))
  counts <- fit$obs
  a <- counts + prior
  b <- matrix(rowSums(a), nrow(a), ncol(a)) - a
  share <- counts / rowSums(counts)
  expected <- stats::pbeta(0.75 * share, a, b) +
    stats::pbeta(1.25 * share, a, b, lower.tail = FALSE)
  lookup <- data.frame(from = rep(rownames(counts), times = ncol(counts)),
                       to = rep(colnames(counts), each = nrow(counts)),
                       reference = as.vector(expected))
  both <- merge(cert, lookup, by = c("from", "to"))
  expect_identical(nrow(both), length(expected))
  positive <- subset(both, reference > 0)
  expect_gt(nrow(positive), 0L)
  relative <- abs(positive$p_value - positive$reference) / positive$reference
  # Before the fix Average -> Average was 3.14102e-17 against 3.14168e-17
  # (relative 2.1e-4).
  expect_lt(max(relative), 1e-12)
})
