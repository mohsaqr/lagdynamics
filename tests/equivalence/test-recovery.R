# Parameter recovery (simulab, the owner's simulation package; the truth is the parameters it stores with the data,
# not its code): Markov chains drawn by simulab::simulate_markov() from a known transition matrix and known initial
# probabilities, fitted by lsa(). The tolerances are fixed before any run from the sampling error of the design.
# 3000 chains of length 30 give 87000 transitions; under this matrix each state is left at least 20000 times
# (asserted below), so the standard error of an estimated transition probability is at most sqrt(0.25 / 20000) and
# the tolerance is four of those, 0.0141. The initial state of 3000 chains: at most sqrt(0.25 / 3000) per state,
# tolerance four of those, 0.0365. Three seeds.

truth <- matrix(c(0.60, 0.30, 0.10,
                  0.20, 0.50, 0.30,
                  0.25, 0.25, 0.50), 3L, 3L, byrow = TRUE, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
start <- c(A = 0.5, B = 0.3, C = 0.2)
seeds <- 1:3
recovery <- "independent: recovery of the parameters of simulab simulate_markov()"

simulated <- lapply(seeds, \(s) simulab::simulate_markov(n = 3000L, transition = truth, chain_length = 30L,
                                                         initial = start, states = rownames(truth), seed = s))

test_that("the transition probabilities of lsa() recover simulab's transition matrix", {
  rows <- Map(\(data, s) {
    fit <- lsa(data, actor = "id", action = "state", time = "period")
    stored <- as.data.frame(data, what = "transitions")
    list(ours = cbind(seed = s, transition_probabilities(fit)),
         truth = data.frame(seed = s, from = stored$from, to = stored$to, prob = stored$probability),
         left = nodes(fit)$outgoing)
  }, simulated, seeds)
  expect_true(all(unlist(lapply(rows, `[[`, "left")) >= 20000))
  expect_reference(do.call(rbind, lapply(rows, `[[`, "ours")), do.call(rbind, lapply(rows, `[[`, "truth")),
                   reference = recovery, covers = c("lsa()", "transition_probabilities()", "transition_probabilities.lsa()"),
                   quantity = "transition probability, three seeds", tolerance = 4 * sqrt(0.25 / 20000),
                   key = c("seed", "from", "to"))
})

test_that("the initial probabilities of lsa() recover simulab's initial distribution", {
  rows <- Map(\(data, s) {
    stored <- as.data.frame(data, what = "initial_probabilities")
    list(ours = cbind(seed = s, initial(lsa(data, actor = "id", action = "state", time = "period"))),
         truth = data.frame(seed = s, state = stored$state, init_prob = stored$probability))
  }, simulated, seeds)
  expect_reference(do.call(rbind, lapply(rows, `[[`, "ours")), do.call(rbind, lapply(rows, `[[`, "truth")),
                   reference = recovery, covers = c("initial()", "initial.lsa()"),
                   quantity = "initial probability, three seeds", tolerance = 4 * sqrt(0.25 / 3000),
                   key = c("seed", "state"))
})

test_that("the sequence bootstrap interval of each transition probability covers the true value in most cells", {
  # The bootstrap's own standard error of a transition probability is close to its sampling error: across the 27
  # cells of three seeds, the 95% interval misses the truth in few cells. Recorded per cell as the share covered:
  # the expected share is 0.95; with 27 cells, fewer than 23 covered would be below 0.95 by more than 4 binomial SE
  # (sqrt(0.95 * 0.05 / 27) = 0.042). The comparison is the coverage share against 0.95 with tolerance 4 * 0.042.
  covered <- unlist(Map(\(data, s) {
    set.seed(100 + s)
    boot <- transitions(bootstrap_lsa(lsa(data, actor = "id", action = "state", time = "period"), R = 200L,
                                      level = "sequence", level_alpha = 0.95))
    stored <- as.data.frame(data, what = "transitions")
    merged <- merge(boot, data.frame(from = stored$from, to = stored$to, true_prob = stored$probability),
                    by = c("from", "to"))
    merged$prob_ci_low <= merged$true_prob & merged$true_prob <= merged$prob_ci_high
  }, simulated, seeds))
  expect_reference(mean(covered), 0.95, reference = recovery, covers = c("bootstrap_lsa()", "transitions.lsa_bootstrap()"),
                   quantity = "coverage of the 95% percentile interval, 27 cells", tolerance = 4 * sqrt(0.95 * 0.05 / 27))
})
