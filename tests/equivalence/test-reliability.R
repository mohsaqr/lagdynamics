# reliability_lsa(): split-half reliability (Epskamp, Borsboom and Fried 2018). The halves are drawn inside the
# function with no argument to supply them, so the evidence is of two kinds.
#  - The summary of a run is recomputed from its own replicates (as.data.frame()): mean, standard deviation and the
#    2.5% and 97.5% quantiles of the finite correlations, and the number of failed replicates. Exact.
#  - The distribution is compared statistically with the procedure written out in base R (random halves of the
#    sequences, lag-1 tables of each half, the correlation of the two halves' values over all transitions): a Monte
#    Carlo summary. The tolerance for the mean correlation is four standard errors of the difference of two
#    independent means of R draws, 4 * sqrt(2) * sd / sqrt(R), with sd the reference's own spread; for the standard
#    deviation, four standard errors of the difference of two sample standard deviations, with the standard error
#    of one sd sqrt((kurtosis - 1) / (4 R)) * sd from the reference's own sample kurtosis (the large-sample variance
#    of a standard deviation; it is sd / sqrt(2 R) for normal draws). A first rule that assumed normal draws failed
#    2 of 48 comparisons in a 24-seed sweep (the split-half correlations are skewed towards one), so the kurtosis was
#    put in before these tests were run as declared comparisons.

replicates <- 300L
seeds <- 1:4
split_half <- "independent: split-half resampling of whole sequences (Epskamp et al. 2018), base R"

#' Values of one half: transition probabilities, counts or adjusted residuals of its lag-1 table
half_values <- function(members, states, weights) {
  counts <- unclass(ref_pairs(members, states = states))
  switch(weights,
         prob = as.vector(prop.table(counts, margin = 1L)),
         count = as.vector(counts),
         adj_res = as.vector(suppressWarnings(stats::chisq.test(counts, correct = FALSE))$stdres))
}

#' The split-half correlations of R random splits, written out
reference_split_half <- function(sequences, replicates, weights, method) {
  states <- sort(unique(unlist(sequences)))
  n <- length(sequences)
  vapply(seq_len(replicates), \(r) {
    first <- sample.int(n, n %/% 2L)
    one <- half_values(sequences[first], states, weights)
    two <- half_values(sequences[-first], states, weights)
    ok <- is.finite(one) & is.finite(two)
    stats::cor(one[ok], two[ok], method = method)
  }, numeric(1))
}

#' The large-sample standard error of a sample standard deviation, from the draws' own kurtosis
sd_standard_error <- function(draws) {
  centred <- draws - mean(draws)
  kurtosis <- mean(centred^4) / mean(centred^2)^2
  sqrt((kurtosis - 1) / (4 * length(draws))) * stats::sd(draws)
}

monte_carlo_rows <- function(fit, sequences, weights, method, label) {
  lapply(seeds, \(s) {
    set.seed(s)
    ours <- as.data.frame(reliability_lsa(fit, R = replicates, weights = weights, method = method),
                          what = "summary")
    set.seed(1000 + s)
    theirs <- reference_split_half(sequences, replicates, weights, method)
    list(ours = data.frame(case = sprintf("%s, seed %d", label, s), mean = ours$mean, sd = ours$sd),
         theirs = data.frame(case = sprintf("%s, seed %d", label, s), mean = mean(theirs), sd = stats::sd(theirs)),
         spread = stats::sd(theirs), sd_error = sd_standard_error(theirs))
  })
}

declare_monte_carlo <- function(rows, covers, label) {
  ours <- do.call(rbind, lapply(rows, `[[`, "ours"))
  theirs <- do.call(rbind, lapply(rows, `[[`, "theirs"))
  spread <- max(vapply(rows, `[[`, numeric(1), "spread"))
  sd_error <- max(vapply(rows, `[[`, numeric(1), "sd_error"))
  expect_reference(subset(ours, select = c(case, mean)), subset(theirs, select = c(case, mean)),
                   reference = split_half, covers = covers,
                   quantity = sprintf("Monte Carlo summary: mean split-half correlation, %s", label),
                   tolerance = 4 * sqrt(2) * spread / sqrt(replicates), key = "case")
  expect_reference(subset(ours, select = c(case, sd)), subset(theirs, select = c(case, sd)),
                   reference = split_half, covers = covers,
                   quantity = sprintf("Monte Carlo summary: sd of split-half correlations, %s", label),
                   tolerance = 4 * sqrt(2) * sd_error, key = "case")
}

test_that("the summary of a run is the mean, sd and 2.5% / 97.5% quantiles of its own replicates", {
  rows <- lapply(seeds, \(s) {
    set.seed(s)
    result <- reliability_lsa(lsa(engagement), R = 60L)
    correlations <- as.data.frame(result)$correlation
    finite <- correlations[is.finite(correlations)]
    list(ours = transform(as.data.frame(result, what = "summary"), seed = s),
         theirs = data.frame(seed = s, n_failed = sum(!is.finite(correlations)), mean = mean(finite),
                             sd = stats::sd(finite), ci_low = unname(stats::quantile(finite, 0.025)),
                             ci_high = unname(stats::quantile(finite, 0.975))))
  })
  expect_reference(do.call(rbind, lapply(rows, `[[`, "ours")), do.call(rbind, lapply(rows, `[[`, "theirs")),
                   reference = "independent: mean, sd and quantiles of the replicates, base R",
                   covers = c("reliability_lsa()", "reliability_lsa.lsa()", "as.data.frame.lsa_reliability()"),
                   quantity = "summary of the replicates", tolerance = 1e-12, key = "seed")
})

test_that("split-half correlations of transition probabilities match the procedure written out (Monte Carlo)", {
  rows <- monte_carlo_rows(lsa(engagement), ref_sequences(engagement), "prob", "pearson", "engagement, prob")
  declare_monte_carlo(rows, c("reliability_lsa()", "reliability_lsa.lsa()", "as.data.frame.lsa_reliability()"),
                      "engagement, prob, pearson")
})

test_that("split-half Spearman correlations of adjusted residuals match the procedure written out (Monte Carlo)", {
  rows <- monte_carlo_rows(lsa(group_regulation), ref_sequences(group_regulation), "adj_res", "spearman",
                           "group_regulation, adj_res")
  declare_monte_carlo(rows, c("reliability_lsa()", "reliability_lsa.lsa()"),
                      "group_regulation, adj_res, spearman")
})

test_that("a grouped fit: each group's split-half correlations match the procedure on its own sequences", {
  log_sequences <- ref_log_sequences(group_regulation_long, actor = "Actor", action = "Action", time = "Time")
  achiever <- tapply(group_regulation_long$Achiever, group_regulation_long$Actor, \(level) level[[1L]])
  groups <- split(log_sequences, achiever[names(log_sequences)])
  grouped <- lsa(group_regulation_long, actor = "Actor", action = "Action", time = "Time", group = "Achiever")
  rows <- lapply(seeds, \(s) {
    set.seed(s)
    ours <- as.data.frame(reliability_lsa(grouped, R = replicates), what = "summary")
    set.seed(1000 + s)
    per_group <- lapply(names(groups), \(g) reference_split_half(groups[[g]], replicates, "prob", "pearson"))
    list(ours = data.frame(case = sprintf("%s, seed %d", ours$group, s), mean = ours$mean, sd = ours$sd),
         theirs = data.frame(case = sprintf("%s, seed %d", names(groups), s),
                             mean = vapply(per_group, mean, numeric(1)), sd = vapply(per_group, stats::sd, numeric(1))),
         spread = max(vapply(per_group, stats::sd, numeric(1))),
         sd_error = max(vapply(per_group, sd_standard_error, numeric(1))))
  })
  declare_monte_carlo(rows, c("reliability_lsa.lsa_group()", "as.data.frame.lsa_reliability_group()"),
                      "group_regulation_long by achievement")
})
