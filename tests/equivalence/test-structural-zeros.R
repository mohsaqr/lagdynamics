# Structural zeros. lsa_ipf() fits quasi-independence by iterative proportional fitting; the references are
# stats::loglin() (the same model by another routine) and stats::glm() (the model as a Poisson regression on
# the cells that are possible). Then the two findings recorded in equivalence/FINDINGS.md for loops = FALSE:
# G2 with self-transitions present is a known difference (owner, 2026-11-22); the print header is a declared
# defect (`defect =`): recorded as a disagreement, and the test fails once it no longer shows.

test_that("lsa_ipf() equals stats::loglin() started at the pattern of structural zeros", {
  set.seed(99)
  observed <- matrix(sample.int(20L, 16L, replace = TRUE), 4L, 4L, dimnames = list(letters[1:4], letters[1:4]))
  diag(observed) <- 0L
  possible <- 1 - diag(4L)
  reference <- suppressWarnings(stats::loglin(observed, margin = list(1, 2), start = possible, fit = TRUE, print = FALSE, eps = 1e-10))
  # Two iterative fits of one model: lsa_ipf() stops at a margin difference of 1e-8, loglin() at 1e-10.
  expect_reference(lsa_ipf(observed, structure = possible)$fit, reference$fit, reference = "stats", covers = "lsa_ipf()",
                   quantity = "expected count", tolerance = 1e-6)
})

test_that("lsa(loops = FALSE) on a table without self-transitions equals quasi-independence fitted by stats::glm()", {
  observed <- ref_counts_without_repeats(group_regulation)
  cells <- subset(ref_long(observed, "count"), from != to)
  model <- stats::glm(count ~ from + to, family = stats::poisson(), data = cells,
                      control = stats::glm.control(epsilon = 1e-14, maxit = 100L))
  expected <- data.frame(from = cells$from, to = cells$to, expected = unname(stats::fitted(model)),
                         adj_res = unname(stats::rstandard(model, type = "pearson")))
  fit <- lsa(observed, lag = 1, loops = FALSE)
  # IPF to a margin difference of 1e-8 against IRLS: the fits agree to the convergence of the first.
  expect_reference(subset(transitions(fit), from != to), expected, reference = "stats",
                   covers = c("lsa()", "lsa_ipf()", "transitions()"), quantity = "quasi-independence", tolerance = 1e-6,
                   key = c("from", "to"))
  # df = 72 possible cells - (1 + 8 + 8) parameters; G2 is the model's residual deviance
  expect_reference(subset(tests(fit), test == "lrx2", select = c(test, statistic, df)),
                   data.frame(test = "lrx2", statistic = stats::deviance(model), df = stats::df.residual(model)),
                   reference = "stats", covers = c("tests()", "lsa()"), quantity = "quasi-independence", tolerance = 1e-6, key = "test")
})

test_that("with loops = FALSE on data that has self-transitions, G2 is reported negative (known difference)", {
  # 1074 of the 1734 engagement transitions are self-transitions. Declaring them impossible while they are in
  # the margins gives expected counts that do not sum to the observed ones, and G2 = 2 sum O log(O / E) < 0.
  # The package warns (class lagdynamics_structural_zero_counts) and still reports the number.
  fit <- suppressWarnings(lsa(engagement, loops = FALSE), classes = "lagdynamics_structural_zero_counts")
  cells <- subset(ref_long(ref_pairs(ref_sequences(engagement)), "count"), from != to)
  model <- stats::glm(count ~ from + to, family = stats::poisson(), data = cells)
  # A known difference (owner, 2026-11-22): documented and warned as not interpretable; the glm reference fits
  # another model (the self-transitions left out of the margins).
  expect_known_difference(
    subset(tests(fit), test == "lrx2", select = c(test, statistic)), data.frame(test = "lrx2", statistic = stats::deviance(model)),
    reference = "stats", covers = c("tests()", "lsa()"), quantity = "G2 with loops = FALSE and self-transitions present",
    tolerance = 1e-6, key = "test",
    reason = "With loops = FALSE and self-transitions present, lsa() keeps them in the margins and reports G2 = -1184.633 with a documented, classed warning (lagdynamics_structural_zero_counts) that the tablewise test is not interpretable. The glm reference fits a different model: quasi-independence on the six possible cells, with the self-transitions left out of the margins (deviance 5.99).",
    agrees_with = "lsa(loops = FALSE) on a table without self-transitions (stats::glm quasi-independence, compared above)")
})

test_that("the print header counts significant transitions out of all cells when every cell is possible", {
  fit <- lsa(group_regulation_long, actor = "Actor", action = "Action", time = "Time")
  header <- grep("Significant transitions", utils::capture.output(print(fit)), value = TRUE)
  states <- length(unique(group_regulation_long$Action))
  expect_reference(as.numeric(sub(".* of ([0-9]+).*", "\\1", header)), states * states, reference = "independent: cells of the table, K * K",
                   covers = "print.lsa()", quantity = "transitions the header counts significance out of, all cells possible",
                   tolerance = 0)
})

test_that("with loops = FALSE the print header counts significant transitions out of the cells tested", {
  fit <- suppressWarnings(lsa(group_regulation_long, actor = "Actor", action = "Action", time = "Time", loops = FALSE),
                          classes = "lagdynamics_structural_zero_counts")
  header <- grep("Significant transitions", utils::capture.output(print(fit)), value = TRUE)
  printed_total <- as.numeric(sub(".* of ([0-9]+).*", "\\1", header))
  states <- length(unique(group_regulation_long$Action))
  expect_reference(
    printed_total, states * states - states, reference = "independent: cells that are possible, K * K - K",
    covers = "print.lsa()", quantity = "transitions the header counts significance out of", tolerance = 0)
})
