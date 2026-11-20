# lsa(), transitions(): the classical lag sequential analysis of the bundled engagement sequences (138
# students, 15 weeks, 3 states). Counts and probabilities against a plain table of consecutive pairs and
# against tna; expected counts, adjusted residuals (Haberman 1973) and p-values against stats::chisq.test()
# and stats::pnorm(); the association indices against their definitions.

counts <- ref_pairs(ref_sequences(engagement), lag = 1L)
independence <- suppressWarnings(stats::chisq.test(counts, correct = FALSE))

test_that("lag-1 counts and transition probabilities equal a plain table of consecutive pairs", {
  fit <- lsa(engagement, lag = 1)
  expected <- ref_cells(count = counts, prob = prop.table(counts, margin = 1L))
  # Counts are integers and a probability is one division of two integers: rounding of the last bit only.
  expect_reference(transitions(fit), expected, reference = "independent: table of consecutive pairs, base R",
                   covers = c("lsa()", "transitions()", "transitions.lsa()"), quantity = "lag 1", tolerance = 1e-12,
                   key = c("from", "to"))
  expect_reference(as.data.frame(fit), ref_labelled(expected), reference = "independent: table of consecutive pairs, base R",
                   covers = "as.data.frame.lsa()", quantity = "lag 1, as.data.frame()", tolerance = 1e-12,
                   key = c("from_label", "to_label"))
})

test_that("lag-1 counts and transition probabilities equal tna's frequency and relative models", {
  skip_without_reference("tna")
  expected <- ref_cells(count = tna::ftna(engagement)$weights, prob = tna::tna(engagement)$weights)
  expect_reference(transitions(lsa(engagement, lag = 1)), expected, reference = "tna",
                   covers = c("lsa()", "transitions()", "transitions.lsa()"), quantity = "lag 1", tolerance = 1e-12,
                   key = c("from", "to"))
})

test_that("expected counts, adjusted residuals and two-sided p-values equal stats::chisq.test()", {
  fit <- lsa(engagement, lag = 1, engine = "classical", alternative = "two.sided")
  expected <- ref_cells(expected = independence$expected, adj_res = independence$stdres,
                        p = 2 * stats::pnorm(-abs(independence$stdres)))
  # An analytic identity. Relative, because the p-values run down to 1e-104, where an absolute 1e-10
  # would compare nothing.
  expect_reference(transitions(fit), expected, reference = "stats",
                   covers = c("lsa()", "transitions()", "transitions.lsa()"), quantity = "independence", tolerance = 1e-10,
                   measure = "relative", key = c("from", "to"),
                   note = "stats::chisq.test(correct = FALSE): $expected, $stdres, and 2 * pnorm(-|stdres|).")
  expect_reference(transitions(lsa_classical(engagement)), expected, reference = "stats", covers = "lsa_classical()",
                   quantity = "independence, lsa_classical()", tolerance = 1e-10, measure = "relative", key = c("from", "to"))
})

test_that("one-sided p-values equal stats::pnorm() of the standardised residual", {
  greater <- lsa(engagement, alternative = "greater")
  less <- lsa(engagement, alternative = "less")
  expect_reference(transitions(greater), ref_cells(p = stats::pnorm(independence$stdres, lower.tail = FALSE)),
                   reference = "stats", covers = c("lsa()", "transitions()"), quantity = "alternative greater",
                   tolerance = 1e-10, measure = "relative", key = c("from", "to"))
  expect_reference(transitions(less), ref_cells(p = stats::pnorm(independence$stdres)),
                   reference = "stats", covers = c("lsa()", "transitions()"), quantity = "alternative less",
                   tolerance = 1e-10, measure = "relative", key = c("from", "to"))
})

test_that("Yule's Q, lift and the column-conditional probability equal their definitions", {
  cell <- ref_collapse(counts)
  observed <- unclass(counts)
  expected <- ref_cells(
    # Q = (ad - bc) / (ad + bc) on the 2 x 2 collapse of each cell
    yules_q = (cell$a * cell$d - cell$b * cell$c) / (cell$a * cell$d + cell$b * cell$c),
    # lift = observed / expected under independence
    lift = observed / (outer(rowSums(observed), colSums(observed)) / sum(observed)),
    # P(antecedent | consequent)
    prob_col = prop.table(observed, margin = 2L))
  expect_reference(transitions(lsa(engagement)), expected, reference = "independent: the definitions, base R",
                   covers = c("lsa()", "transitions()", "transitions.lsa()"), quantity = "association", tolerance = 1e-12,
                   key = c("from", "to"))
})

test_that("a random sequence of four codes: adjusted residuals equal stats::chisq.test()", {
  set.seed(20260512)
  events <- sample(c("a", "b", "c", "d"), 200L, replace = TRUE)
  pairs <- ref_pairs(list(events))
  test <- suppressWarnings(stats::chisq.test(pairs, correct = FALSE))
  expect_reference(transitions(lsa(events)), ref_cells(count = pairs, expected = test$expected, adj_res = test$stdres),
                   reference = "stats", covers = c("lsa()", "transitions()"), quantity = "random sequence", tolerance = 1e-10,
                   key = c("from", "to"))
})

test_that("a matrix of counts gives the same residuals as stats::chisq.test() on that matrix", {
  three <- ref_three()
  test <- stats::chisq.test(three, correct = FALSE)
  expect_reference(transitions(lsa(three)), ref_cells(count = three, expected = test$expected, adj_res = test$stdres),
                   reference = "stats", covers = c("lsa()", "transitions()"), quantity = "3 x 3 table", tolerance = 1e-10,
                   key = c("from", "to"))
})
