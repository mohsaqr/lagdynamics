# lsa_lags(), lag_profile(): the engagement sequences at lags 1 to 3, each event paired with the event one,
# two and three steps later in the same sequence. Against a plain table of lagged pairs and
# stats::chisq.test() on each lag's table.

sequences <- ref_sequences(engagement)
by_lag <- do.call(rbind, lapply(1:3, \(lag) {
  counts <- ref_pairs(sequences, lag = lag)
  test <- stats::chisq.test(counts, correct = FALSE)
  cbind(lag = lag, ref_cells(count = counts, prob = prop.table(counts, margin = 1L), adj_res = test$stdres,
                             p = 2 * stats::pnorm(-abs(test$stdres))))
}))

test_that("counts, probabilities and adjusted residuals at lags 1 to 3", {
  fits <- lsa_lags(engagement, lags = 1:3)
  # Counts are exact; residuals and p-values are an analytic identity, relative because p runs to 1e-104.
  expect_reference(transitions(fits), by_lag, reference = "stats",
                   covers = c("lsa_lags()", "transitions()", "transitions.lsa_lags()"), quantity = "lags 1 to 3",
                   tolerance = 1e-10, measure = "relative", key = c("lag", "from", "to"),
                   note = "A plain table of lagged pairs in base R, with stats::chisq.test(correct = FALSE)$stdres on each lag's table.")
  expect_reference(as.data.frame(fits), subset(by_lag, select = c(lag, from, to, count, adj_res)), reference = "stats",
                   covers = "as.data.frame.lsa_lags()", quantity = "lags 1 to 3, as.data.frame()", tolerance = 1e-10,
                   measure = "relative", key = c("lag", "from", "to"))
})

test_that("the lag profile of one transition is that transition's row at each lag", {
  expected <- subset(by_lag, from == "Active" & to == "Average")
  expect_reference(lag_profile(lsa_lags(engagement, lags = 1:3), from = "Active", to = "Average"), expected, reference = "stats",
                   covers = "lag_profile()", quantity = "Active -> Average", tolerance = 1e-10, measure = "relative",
                   key = c("lag", "from", "to"))
})
