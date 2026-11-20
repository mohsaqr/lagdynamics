# tests(): the tablewise tests of independence of a fit (likelihood-ratio G2 and Pearson X2), against
# stats::chisq.test() with G2 from its definition, and against stats::loglin().

counts <- ref_pairs(ref_sequences(engagement), lag = 1L)

test_that("G2 and X2 of the engagement table equal stats::chisq.test() and the definition of G2", {
  pearson <- stats::chisq.test(counts, correct = FALSE)
  g2 <- 2 * sum(counts * log(counts / pearson$expected))
  df <- unname(pearson$parameter)
  expected <- data.frame(test = c("lrx2", "x2"), statistic = c(g2, unname(pearson$statistic)), df = c(df, df),
                         p = c(stats::pchisq(g2, df, lower.tail = FALSE), pearson$p.value))
  # Closed-form statistics on the same table; relative because the p-values are of the order 1e-132.
  expect_reference(tests(lsa(engagement)), expected, reference = "stats", covers = c("tests()", "tests.lsa()", "lsa()"),
                   quantity = "independence of the table", tolerance = 1e-10, measure = "relative", key = "test",
                   note = "X2 from stats::chisq.test(correct = FALSE); G2 = 2 sum O log(O / E) with its expected counts; p by stats::pchisq().")
})

test_that("G2 and X2 of the engagement table equal the independence model of stats::loglin()", {
  model <- stats::loglin(counts, margin = list(1, 2), fit = TRUE, eps = 1e-10, print = FALSE)
  expected <- data.frame(test = c("lrx2", "x2"), statistic = c(model$lrt, model$pearson), df = c(model$df, model$df),
                         p = stats::pchisq(c(model$lrt, model$pearson), model$df, lower.tail = FALSE))
  expect_reference(tests(lsa(engagement)), expected, reference = "stats", covers = c("tests()", "tests.lsa()"),
                   quantity = "independence of the table, loglin()", tolerance = 1e-10, measure = "relative", key = "test")
})

test_that("the summary of a fit carries the same tablewise tests and the size of the data", {
  pearson <- stats::chisq.test(counts, correct = FALSE)
  g2 <- 2 * sum(counts * log(counts / pearson$expected))
  sequences <- ref_sequences(engagement)
  expected <- data.frame(engine = "classical", lrx2 = g2, x2 = unname(pearson$statistic), lrx2_df = unname(pearson$parameter),
                         n_transitions = sum(counts), n_events = length(unlist(sequences)), n_sequences = length(sequences))
  expect_reference(as.data.frame(summary(lsa(engagement))), expected, reference = "stats",
                   covers = c("summary.lsa()", "as.data.frame.summary.lsa()"), quantity = "summary", tolerance = 1e-10,
                   measure = "relative", key = "engine")
})

test_that("a grouped fit tests each group's own table", {
  sequences <- ref_log_sequences(group_regulation_long, actor = "Actor", action = "Action", time = "Time")
  achiever <- tapply(group_regulation_long$Achiever, group_regulation_long$Actor, \(level) level[[1L]])
  states <- sort(unique(group_regulation_long$Action))
  by_group <- lapply(split(sequences, achiever[names(sequences)]), ref_pairs, lag = 1L, states = states)
  expected <- do.call(rbind, lapply(names(by_group), \(level) {
    pearson <- suppressWarnings(stats::chisq.test(by_group[[level]], correct = FALSE))
    data.frame(group = level, test = "x2", statistic = unname(pearson$statistic), df = unname(pearson$parameter))
  }))
  fit <- lsa(group_regulation_long, actor = "Actor", action = "Action", time = "Time", group = "Achiever")
  expect_reference(subset(tests(fit), test == "x2"), expected, reference = "stats", covers = c("tests()", "tests.lsa_group()", "lsa()"),
                   quantity = "Pearson X2 by group", tolerance = 1e-10, measure = "relative", key = c("group", "test"))
  expect_reference(subset(summary(fit), select = c(group, x2, x2_df)),
                   data.frame(group = expected$group, x2 = expected$statistic, x2_df = expected$df), reference = "stats",
                   covers = "summary.lsa_group()", quantity = "summary by group", tolerance = 1e-10, measure = "relative", key = "group")
})
