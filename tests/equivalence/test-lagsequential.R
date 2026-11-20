# lsa() against LagSequential (O'Connor; CRAN 0.1.1), a second implementation of classical lag sequential analysis
# (Bakeman and Quera 2011). Its sequential() takes one sequence of numeric codes and returns, per cell of the lagged
# table, the frequency, expected frequency, transitional probability, adjusted residual and its p-value, Yule's Q,
# and Wampold's transformed kappa with its z and p. lsa() on the same single sequence (one row of a wide table)
# names numeric codes as LagSequential does ("Code 1", ...). The sequences are LagSequential's own example data.

#' LagSequential's cell tables as rows keyed by from and to (its dimnames carry a leading space)
lagsequential_cells <- function(codes, lag, adjacent = TRUE) {
  printed <- utils::capture.output(result <- LagSequential::sequential(codes, lag = lag, adjacent = adjacent))
  tables <- list(count = result$freqs, expected = result$expfreqs, prob = result$probs, adj_res = result$adjres,
                 p = result$p, yules_q = result$YulesQ, kappa = result$kappas, kappa_z = result$z,
                 kappa_p = result$pk)
  states <- trimws(rownames(result$expfreqs))
  named <- lapply(tables, \(cells) matrix(cells, length(states), dimnames = list(states, states)))
  do.call(ref_cells, named)
}

our_cells <- function(codes, lag, loops = TRUE) {
  transitions(lsa(as.data.frame(t(codes)), lag = lag, loops = loops))
}

test_that("lag 1 table of LagSequential's example sequence equals LagSequential::sequential()", {
  skip_without_reference("LagSequential")
  codes <- as.vector(LagSequential::data_sequential)
  # Ratios, square roots and normal tails of counts below 125.
  expect_reference(our_cells(codes, lag = 1L), lagsequential_cells(codes, lag = 1L), reference = "LagSequential",
                   covers = c("lsa()", "transitions()"), quantity = "data_sequential, lag 1", tolerance = 1e-10,
                   key = c("from", "to"))
})

test_that("lag 2: the table, residuals and Yule's Q equal LagSequential; Wampold's kappa follows another event count", {
  skip_without_reference("LagSequential")
  codes <- as.vector(LagSequential::data_sequential)
  ours <- our_cells(codes, lag = 2L)
  theirs <- lagsequential_cells(codes, lag = 2L)
  expect_reference(ours, subset(theirs, select = c(from, to, count, expected, prob, adj_res, p, yules_q)),
                   reference = "LagSequential", covers = c("lsa()", "transitions()"),
                   quantity = "data_sequential, lag 2", tolerance = 1e-10, key = c("from", "to"))
  expect_known_difference(
    ours, subset(theirs, select = c(from, to, kappa, kappa_z, kappa_p)), reference = "LagSequential",
    covers = c("lsa()", "transitions()"), quantity = "data_sequential, lag 2, Wampold's kappa", tolerance = 1e-10,
    key = c("from", "to"),
    reason = paste("Wampold's transformed kappa needs each code's number of events. lsa() counts every event of",
                   "the sequence (n = 122); LagSequential's sequential() takes the senders of the lag-2 pairs plus",
                   "the last event (n = number of pairs + 1 = 121), which equals the event count at lag 1 only. The",
                   "two agree at lag 1 (the test above)."))
})

test_that("Wampold's (1984) sequence equals LagSequential::sequential()", {
  skip_without_reference("LagSequential")
  codes <- as.vector(as.matrix(LagSequential::data_Wampold_1984))
  expect_reference(our_cells(codes, lag = 1L), lagsequential_cells(codes, lag = 1L), reference = "LagSequential",
                   covers = c("lsa()", "transitions()"), quantity = "data_Wampold_1984, lag 1", tolerance = 1e-10,
                   key = c("from", "to"))
})

test_that("codes that cannot follow themselves: loops = FALSE equals LagSequential's adjacent = FALSE", {
  skip_without_reference("LagSequential")
  # Repeats collapsed (rle()), so that no code follows itself on either side.
  codes <- rle(as.vector(LagSequential::data_sequential))$values
  ours <- subset(our_cells(codes, lag = 1L, loops = FALSE), from != to)
  theirs <- subset(lagsequential_cells(codes, lag = 1L, adjacent = FALSE), from != to)
  expect_reference(ours, subset(theirs, select = c(from, to, count, prob)), reference = "LagSequential",
                   covers = c("lsa()", "transitions()"), quantity = "data_sequential without repeats, lag 1",
                   tolerance = 1e-10, key = c("from", "to"))
  expect_known_difference(
    ours, subset(theirs, select = c(from, to, expected, adj_res, p)), reference = "LagSequential",
    covers = c("lsa()", "transitions()"), quantity = "data_sequential without repeats, lag 1, expected and residuals",
    tolerance = 1e-8, measure = "relative", key = c("from", "to"),
    reason = paste("With repeats impossible, LagSequential's sequential(adjacent = FALSE) takes the one-step",
                   "formula E[i, j] = R[i] C[j] / (N - R[j]) and residuals from it; lsa(loops = FALSE) fits",
                   "quasi-independence by iterative proportional fitting, the maximum likelihood fit (equal to",
                   "stats::glm() and stats::loglin() in test-structural-zeros.R)."))
})
