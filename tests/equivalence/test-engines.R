# The four engines beside the classical one, on the 3 x 3 table of the engine tests (no zero cell).
# Two-cell: odds ratio, its logarithm, Wald standard error and z of each cell's 2 x 2 collapse (Bakeman and
# Gottman 1997). Bidirectional: Haberman residuals of the symmetrised table O + t(O). Parallel dominance:
# (O - t(O)) / sqrt(E + t(E)). Non-parallel dominance: (O - t(O)) / sqrt(O + t(O)), with an exact binomial test.

three <- ref_three()

test_that("two-cell: odds ratio, log odds ratio, standard error, z and p equal the definitions", {
  cell <- ref_collapse(three)
  log_or <- log((cell$a * cell$d) / (cell$b * cell$c))
  log_or_se <- sqrt(1 / cell$a + 1 / cell$b + 1 / cell$c + 1 / cell$d)
  z <- log_or / log_or_se
  expected <- ref_cells(odds_ratio = exp(log_or), log_or = log_or, log_or_se = log_or_se, adj_res = z)
  expect_reference(transitions(lsa_two_cell(three)), expected, reference = "independent: the definitions, base R",
                   covers = c("lsa_two_cell()", "lsa()", "transitions()"), quantity = "two-cell", tolerance = 1e-12,
                   key = c("from", "to"))
  # The p-value of the z written out above: a formula of this test, not a separate program (until 0.35.4 it was
  # labelled "stats" because it calls stats::pnorm()).
  expect_reference(transitions(lsa(three, engine = "two_cell")), ref_cells(p = 2 * stats::pnorm(-abs(z))),
                   reference = "independent: 2 * pnorm(-|z|) of the z from the definition, base R",
                   covers = c("lsa_two_cell()", "lsa()"), quantity = "two-cell p", tolerance = 1e-12,
                   key = c("from", "to"))
})

# Outside reference for the two-cell engine: vcd::loddsratio() on each cell's 2 x 2 collapse gives the log odds
# ratio and its Wald standard error, and its summary() the z test (lmtest's coeftest). No collapse of this table
# holds a zero, so vcd adds nothing to the cells (its default `correct = any(x == 0)`) and neither does
# lsa_two_cell().
test_that("two-cell: log odds ratio, standard error, z and p equal vcd::loddsratio()", {
  skip_without_reference("vcd")
  cell <- ref_collapse(three)
  tests <- Map(\(a, b, c, d) summary(vcd::loddsratio(matrix(c(a, c, b, d), 2L), log = TRUE)),
               as.vector(cell$a), as.vector(cell$b), as.vector(cell$c), as.vector(cell$d))
  column <- \(name) matrix(vapply(tests, \(test) unname(test[1L, name]), numeric(1)), nrow(three),
                           dimnames = dimnames(three))
  expected <- ref_cells(log_or = column("Estimate"), log_or_se = column("Std. Error"), adj_res = column("z value"),
                        p = column("Pr(>|z|)"))
  expect_reference(transitions(lsa_two_cell(three)), expected, reference = "vcd",
                   covers = c("lsa_two_cell()", "lsa()", "transitions()"), quantity = "two-cell, vcd",
                   tolerance = 1e-12, key = c("from", "to"))
})

test_that("bidirectional: residuals equal stats::chisq.test() on the symmetrised table", {
  symmetric <- stats::chisq.test(three + t(three), correct = FALSE)
  expected <- ref_cells(adj_res = symmetric$stdres, symmetric_obs = three + t(three), symmetric_exp = symmetric$expected)
  expect_reference(transitions(lsa_bidirectional(three)), expected, reference = "stats",
                   covers = c("lsa_bidirectional()", "lsa()", "transitions()"), quantity = "bidirectional", tolerance = 1e-10,
                   key = c("from", "to"))
})

test_that("parallel dominance: the difference of the two directions over sqrt(E + t(E))", {
  expected_counts <- outer(rowSums(three), colSums(three)) / sum(three)
  difference <- three - t(three)
  expected <- ref_cells(dominance = difference, dominance_se = sqrt(expected_counts + t(expected_counts)),
                        adj_res = difference / sqrt(expected_counts + t(expected_counts)))
  expect_reference(transitions(lsa_parallel_dominance(three)), expected, reference = "independent: the definitions, base R",
                   covers = c("lsa_parallel_dominance()", "lsa()", "transitions()"), quantity = "parallel dominance",
                   tolerance = 1e-12, key = c("from", "to"))
})

test_that("non-parallel dominance: the difference over sqrt(O + t(O)), and stats::binom.test() for each pair", {
  difference <- three - t(three)
  totals <- three + t(three)
  expected <- subset(ref_cells(dominance = difference, pair_totals = totals, adj_res = difference / sqrt(totals)), from != to)
  fit <- lsa_nonparallel_dominance(three)
  expect_reference(subset(transitions(fit), from != to), expected, reference = "independent: the definitions, base R",
                   covers = c("lsa_nonparallel_dominance()", "lsa()", "transitions()"), quantity = "non-parallel dominance",
                   tolerance = 1e-12, key = c("from", "to"))
  exact <- mapply(\(successes, trials) stats::binom.test(successes, trials, p = 0.5)$p.value, as.vector(three), as.vector(totals))
  binomial <- subset(ref_cells(binomial_p = matrix(exact, 3L, dimnames = dimnames(three))), from != to)
  expect_reference(subset(transitions(fit), from != to), binomial, reference = "stats", covers = c("lsa_nonparallel_dominance()", "lsa()"),
                   quantity = "non-parallel dominance", tolerance = 1e-12, key = c("from", "to"))
})
