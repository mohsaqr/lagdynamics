# bootstrap_lsa() and permute_lsa(), with the resamples supplied, so that an independent computation can
# follow the same draws: the bootstrap of whole sequences (percentile intervals by stats::quantile()) and the
# permutation test of the adjusted residuals (Phipson and Smyth 2010: p = (1 + exceedances) / (1 + R)).

sequences <- ref_sequences(engagement)
states <- sort(unique(unlist(sequences)))
cell_names <- ref_long(unclass(ref_pairs(sequences, states = states)), "count")

test_that("the sequence bootstrap: means, standard errors and percentile intervals over the same resamples", {
  resamples <- 200L
  set.seed(7)
  indices <- matrix(sample.int(length(sequences), resamples * length(sequences), replace = TRUE), resamples, length(sequences))
  # each sequence's own table of pairs, one row per sequence; a resample's table is the sum of its sequences'
  per_sequence <- t(vapply(sequences, \(events) as.vector(ref_pairs(list(events), states = states)), numeric(length(states)^2)))
  boot_counts <- t(apply(indices, 1L, \(drawn) colSums(per_sequence[drawn, , drop = FALSE])))
  as_table <- \(cells) matrix(cells, length(states), dimnames = list(states, states))
  boot_residuals <- t(apply(boot_counts, 1L, \(cells) as.vector(stats::chisq.test(as_table(cells), correct = FALSE)$stdres)))
  boot_probabilities <- t(apply(boot_counts, 1L, \(cells) as.vector(prop.table(as_table(cells), margin = 1L))))
  low <- \(draws) apply(draws, 2L, stats::quantile, probs = 0.025, names = FALSE)
  high <- \(draws) apply(draws, 2L, stats::quantile, probs = 0.975, names = FALSE)
  expected <- data.frame(
    from = cell_names$from, to = cell_names$to,
    count_mean = colMeans(boot_counts), count_se = apply(boot_counts, 2L, stats::sd),
    count_ci_low = low(boot_counts), count_ci_high = high(boot_counts),
    adj_res_mean = colMeans(boot_residuals), adj_res_se = apply(boot_residuals, 2L, stats::sd),
    adj_res_ci_low = low(boot_residuals), adj_res_ci_high = high(boot_residuals),
    prob_mean = colMeans(boot_probabilities), prob_ci_low = low(boot_probabilities), prob_ci_high = high(boot_probabilities))
  ours <- bootstrap_lsa(lsa(engagement), R = resamples, level = "sequence", level_alpha = 0.95, indices = indices)
  by_resampling <- "independent: bootstrap of whole sequences with the same resamples, stats::chisq.test() and stats::quantile()"
  expect_reference(transitions(ours), expected, reference = by_resampling,
                   covers = c("bootstrap_lsa()", "transitions()", "transitions.lsa_bootstrap()"), quantity = "bootstrap", tolerance = 1e-10,
                   key = c("from", "to"))
  expect_reference(as.data.frame(ours), expected, reference = by_resampling, covers = "as.data.frame.lsa_bootstrap()",
                   quantity = "bootstrap, as.data.frame()", tolerance = 1e-10, key = c("from", "to"))
})

test_that("the permutation test: p = (1 + exceedances) / (1 + R) over the same shuffles", {
  shuffles_count <- 99L
  events <- unlist(sequences)
  owner <- rep(seq_along(sequences), lengths(sequences))
  set.seed(11)
  shuffles <- replicate(shuffles_count, sample.int(length(events)), simplify = FALSE)
  residuals_of <- \(shuffled) {
    pairs <- ref_pairs(split(shuffled, owner), states = states)
    as.vector(stats::chisq.test(pairs, correct = FALSE)$stdres)
  }
  observed <- residuals_of(events)
  shuffled <- vapply(shuffles, \(order) residuals_of(events[order]), numeric(length(observed)))
  exceed <- rowSums(abs(shuffled) >= abs(observed))
  expected <- data.frame(from = cell_names$from, to = cell_names$to, observed_count = cell_names$count,
                         observed_adj_res = observed, p_perm = (1 + exceed) / (1 + shuffles_count))
  ours <- permute_lsa(lsa(engagement), R = shuffles_count, shuffles = shuffles)
  by_shuffling <- "independent: permutation test with the same shuffles, stats::chisq.test()"
  expect_reference(transitions(ours), expected, reference = by_shuffling,
                   covers = c("permute_lsa()", "transitions()", "transitions.lsa_permutation()"), quantity = "permutation",
                   tolerance = 1e-10, key = c("from", "to"))
  expect_reference(as.data.frame(ours), expected, reference = by_shuffling, covers = "as.data.frame.lsa_permutation()",
                   quantity = "permutation, as.data.frame()", tolerance = 1e-10, key = c("from", "to"))
})
