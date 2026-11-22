# stability_lsa(): case-drop stability of significant transitions (Epskamp, Borsboom and Fried 2018). The dropped
# cases are drawn inside the function with no argument to supply them, so the evidence is:
#  - exact: whether each transition is significant in the full fit, from stats::chisq.test()'s adjusted residuals
#    and 2 * pnorm(-|z|) at the fit's alpha;
#  - statistical: each transition's stability (the share of subsets of 80% of the sequences in which it is
#    significant) against the same case-drop written out in base R, a Monte Carlo summary. Each share is a mean of R
#    Bernoulli draws, so the difference of two independent shares has a standard error of at most
#    sqrt(2 * 0.25 / R) (the variance of a Bernoulli draw is at most 1/4); the tolerance is four of those, fixed
#    before any run: 4 * sqrt(0.5 / 400) = 0.141.

subsets <- 400L
seeds <- 1:3
alpha <- 0.05
sequences <- ref_sequences(group_regulation)
states <- sort(unique(unlist(sequences)))

significant_cells <- function(members) {
  counts <- unclass(ref_pairs(members, states = states))
  residuals <- suppressWarnings(stats::chisq.test(counts, correct = FALSE))$stdres
  p <- 2 * stats::pnorm(-abs(residuals))
  is.finite(p) & p < alpha
}

test_that("the full fit's significance decision of each transition equals stats::chisq.test()", {
  set.seed(1)
  result <- stability_lsa(lsa(group_regulation), R = 20L)
  expected <- ref_cells(observed_sig = significant_cells(sequences))
  expect_reference(transitions(result), expected, reference = "stats",
                   covers = c("stability_lsa()", "transitions.lsa_stability()"),
                   quantity = "significant in the full fit", tolerance = 0, key = c("from", "to"))
  expect_reference(as.data.frame(result), expected, reference = "stats", covers = "as.data.frame.lsa_stability()",
                   quantity = "significant in the full fit, as.data.frame()", tolerance = 0, key = c("from", "to"))
})

test_that("stability of each transition matches the case-drop written out (Monte Carlo summary)", {
  n_keep <- round(length(sequences) * 0.8)
  rows <- lapply(seeds, \(s) {
    set.seed(s)
    ours <- as.data.frame(stability_lsa(lsa(group_regulation), R = subsets, proportion = 0.8))
    set.seed(1000 + s)
    decisions <- vapply(seq_len(subsets), \(r) as.vector(significant_cells(sequences[sample.int(length(sequences),
                                                                                                  n_keep)])),
                        logical(length(states)^2))
    theirs <- ref_cells(stability = matrix(rowMeans(decisions), length(states), dimnames = list(states, states)))
    list(ours = transform(subset(ours, select = c(from, to, stability)), seed = s),
         theirs = transform(theirs, seed = s))
  })
  expect_reference(do.call(rbind, lapply(rows, `[[`, "ours")), do.call(rbind, lapply(rows, `[[`, "theirs")),
                   reference = "independent: case-drop resampling of whole sequences (Epskamp et al. 2018), base R",
                   covers = c("stability_lsa()", "as.data.frame.lsa_stability()", "transitions.lsa_stability()"),
                   quantity = "Monte Carlo summary: share of subsets significant",
                   tolerance = 4 * sqrt(0.5 / subsets), key = c("seed", "from", "to"))
})
