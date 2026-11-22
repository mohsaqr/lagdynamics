# certainty_lsa(): the analytic uncertainty of each transition probability. With a Dirichlet(prior) prior on
# each row of the table, a transition probability has a Beta(a, b) posterior, a = count + prior and
# b = (row total + K * prior) - a. Mean and standard deviation from the moments of the Beta distribution;
# interval and the mass outside 0.75 to 1.25 times the observed probability from stats::qbeta() and pbeta().

counts <- unclass(ref_pairs(ref_sequences(engagement), lag = 1L))
prior <- 0.5
a <- counts + prior
row_total <- matrix(rowSums(a), nrow(a), ncol(a), dimnames = dimnames(a))
b <- row_total - a

test_that("posterior mean and standard deviation equal the moments of the Beta posterior", {
  expected <- ref_cells(observed = counts, prob_mean = a / row_total, prob_se = sqrt(a * b / (row_total^2 * (row_total + 1))))
  expect_reference(transitions(certainty_lsa(lsa(engagement), prior = prior)), expected,
                   reference = "independent: moments of the Beta posterior, base R",
                   covers = c("certainty_lsa()", "transitions()", "transitions.lsa_certainty()"), quantity = "posterior", tolerance = 1e-12,
                   key = c("from", "to"))
})

test_that("the 95% interval and the mass outside the consistency range equal stats::qbeta() and stats::pbeta()", {
  observed <- prop.table(counts, margin = 1L)
  outside <- stats::pbeta(0.75 * observed, a, b) + stats::pbeta(1.25 * observed, a, b, lower.tail = FALSE)
  ours <- transitions(certainty_lsa(lsa(engagement), prior = prior, level_alpha = 0.95, inference = "stability"))
  expect_reference(ours, ref_cells(prob_ci_low = stats::qbeta(0.025, a, b), prob_ci_high = stats::qbeta(0.975, a, b)),
                   reference = "stats", covers = c("certainty_lsa()", "transitions.lsa_certainty()"), quantity = "posterior",
                   tolerance = 1e-9, measure = "relative", key = c("from", "to"))
  # The mass is a probability: as one, it agrees to 1e-12 in every cell.
  expect_reference(ours, ref_cells(p_value = outside), reference = "stats", covers = c("certainty_lsa()", "transitions.lsa_certainty()"),
                   quantity = "posterior mass", tolerance = 1e-12, key = c("from", "to"))
  # As a relative quantity it does not, where the mass is tiny (first run of this test, 2026-11-19: Average ->
  # Average, 3.14102e-17 against 3.14168e-17). A known difference (owner, 2026-11-22): floating point, not a defect.
  expect_known_difference(
    ours, ref_cells(p_value = outside), reference = "stats", covers = c("certainty_lsa()", "transitions.lsa_certainty()"),
    quantity = "posterior mass, relative", tolerance = 1e-9, measure = "relative", key = c("from", "to"),
    reason = "certainty_lsa() computes the upper tail as 1 - pbeta(), which loses relative accuracy near 1e-16 (3.14102e-17 against 3.14168e-17). As probabilities the two agree to 1e-12 and no certainty decision changes.",
    agrees_with = "stats::pbeta(lower.tail = FALSE)")
})

test_that("a grouped fit: each group's Beta posterior moments, interval and mass outside the range", {
  log_sequences <- ref_log_sequences(group_regulation_long, actor = "Actor", action = "Action", time = "Time")
  achiever <- tapply(group_regulation_long$Achiever, group_regulation_long$Actor, \(level) level[[1L]])
  group_states <- sort(unique(group_regulation_long$Action))
  groups <- split(log_sequences, achiever[names(log_sequences)])
  expected <- do.call(rbind, lapply(names(groups), \(level) {
    observed <- unclass(ref_pairs(groups[[level]], states = group_states))
    a <- observed + prior
    total <- matrix(rowSums(a), nrow(a), ncol(a), dimnames = dimnames(a))
    b <- total - a
    share <- prop.table(observed, margin = 1L)
    outside <- stats::pbeta(0.75 * share, a, b) + stats::pbeta(1.25 * share, a, b, lower.tail = FALSE)
    cbind(group = level, ref_cells(observed = observed, prob_observed = share, prob_mean = a / total,
              prob_se = sqrt(a * b / (total^2 * (total + 1))), prob_ci_low = stats::qbeta(0.025, a, b),
              prob_ci_high = stats::qbeta(0.975, a, b), p_value = outside))
  }))
  ours <- as.data.frame(certainty_lsa(lsa(group_regulation_long, actor = "Actor", action = "Action", time = "Time",
                                          group = "Achiever"), prior = prior))
  # Quantiles of the Beta distribution by stats::qbeta() on both sides; the masses as probabilities (absolute).
  expect_reference(ours, expected, reference = "independent: moments and quantiles of the Beta posterior per group",
                   covers = "as.data.frame.lsa_certainty_group()", quantity = "posterior per group",
                   tolerance = 1e-12, key = c("group", "from", "to"))
})

# markovchain::markovchainFit(method = "map") is a Bayesian fit of a Markov chain with a Dirichlet prior on each row
# (hyperparameter matrix): its expectedValue, standardError and confidenceInterval are the posterior mean, standard
# deviation and equal-tailed interval of each transition probability. markovchain takes hyperparameters of at
# least 1, so the comparison is at prior = 1 (uniform) and prior = 2.

test_that("posterior mean, sd and 95% interval equal markovchain::markovchainFit(method = 'map')", {
  skip_without_reference("markovchain")
  sequences <- ref_sequences(engagement)
  rows <- lapply(c(1, 2), \(pseudo) {
    states <- sort(unique(unlist(sequences)))
    fit <- markovchain::markovchainFit(data = sequences, method = "map", confidencelevel = 0.95,
                                       hyperparam = matrix(pseudo, length(states), length(states),
                                                           dimnames = list(states, states)))
    interval <- fit$confidenceInterval
    theirs <- cbind(prior = pseudo, ref_cells(prob_mean = fit$expectedValue, prob_se = fit$standardError,
                                              prob_ci_low = interval$lowerEndpointMatrix,
                                              prob_ci_high = interval$upperEndpointMatrix))
    ours <- cbind(prior = pseudo, transitions(certainty_lsa(lsa(engagement), prior = pseudo, level_alpha = 0.95)))
    list(ours = ours, theirs = theirs)
  })
  # Closed-form moments and Beta quantiles on both sides (markovchain by stats::qbeta()); relative for the quantiles.
  expect_reference(do.call(rbind, lapply(rows, `[[`, "ours")), do.call(rbind, lapply(rows, `[[`, "theirs")),
                   reference = "markovchain", covers = c("certainty_lsa()", "transitions.lsa_certainty()"),
                   quantity = "Dirichlet posterior, priors 1 and 2", tolerance = 1e-9, measure = "relative",
                   key = c("prior", "from", "to"))
})

test_that("a grouped fit: each group's posterior mean, sd and interval equal markovchain on its sequences", {
  skip_without_reference("markovchain")
  log_sequences <- ref_log_sequences(group_regulation_long, actor = "Actor", action = "Action", time = "Time")
  achiever <- tapply(group_regulation_long$Achiever, group_regulation_long$Actor, \(level) level[[1L]])
  groups <- split(log_sequences, achiever[names(log_sequences)])
  states <- sort(unique(group_regulation_long$Action))
  expected <- do.call(rbind, lapply(names(groups), \(level) {
    fit <- markovchain::markovchainFit(data = unname(groups[[level]]), method = "map", confidencelevel = 0.95,
                                       possibleStates = states,
                                       hyperparam = matrix(1, length(states), length(states),
                                                           dimnames = list(states, states)))
    interval <- fit$confidenceInterval
    cbind(group = level, ref_cells(prob_mean = fit$expectedValue, prob_se = fit$standardError,
                                   prob_ci_low = interval$lowerEndpointMatrix,
                                   prob_ci_high = interval$upperEndpointMatrix))
  }))
  ours <- as.data.frame(certainty_lsa(lsa(group_regulation_long, actor = "Actor", action = "Action", time = "Time",
                                          group = "Achiever"), prior = 1, level_alpha = 0.95))
  expect_reference(ours, expected, reference = "markovchain", covers = "as.data.frame.lsa_certainty_group()",
                   quantity = "Dirichlet posterior per group, prior 1", tolerance = 1e-9, measure = "relative",
                   key = c("group", "from", "to"))
})
