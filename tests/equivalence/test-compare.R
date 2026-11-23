# compare_lsa() and bayes_compare_lsa(): two groups of the bundled regulation log (high and low achievers).
# The observed part of each is closed form: the log odds ratio of every transition in each group (with 0.5
# added to the four cells of a 2 x 2 collapse that holds a zero, Haldane 1956, Anscombe 1956) and the
# posterior mean of every transition probability under a Dirichlet prior. The drawn parts (compare_lsa()'s
# permutation p-values, bayes_compare_lsa()'s credible intervals, probability of direction and effect size) are
# compared at the end of the file as Monte Carlo summaries, each number at four of its own standard errors.

log_sequences <- ref_log_sequences(group_regulation_long, actor = "Actor", action = "Action", time = "Time")
achiever <- tapply(group_regulation_long$Achiever, group_regulation_long$Actor, \(level) level[[1L]])
states <- sort(unique(group_regulation_long$Action))
tables <- lapply(split(log_sequences, achiever[names(log_sequences)]), \(members) unclass(ref_pairs(members, states = states)))
grouped <- lsa(group_regulation_long, actor = "Actor", action = "Action", time = "Time", group = "Achiever")

test_that("the log odds ratio of each transition in each group, and their difference", {
  log_odds <- function(observed) {
    cell <- ref_collapse(observed)
    add <- 0.5 * (cell$a == 0 | cell$b == 0 | cell$c == 0 | cell$d == 0)
    log(((cell$a + add) * (cell$d + add)) / ((cell$b + add) * (cell$c + add)))
  }
  expected <- ref_cells(log_or_a = log_odds(tables$High), log_or_b = log_odds(tables$Low), diff = log_odds(tables$High) - log_odds(tables$Low))
  set.seed(3)
  ours <- compare_lsa(grouped, R = 19L, measure = "log_or")
  by_definition <- "independent: log odds ratio of the 2 x 2 collapse, base R"
  expect_reference(transitions(ours), expected, reference = by_definition,
                   covers = c("compare_lsa()", "transitions()", "transitions.lsa_comparison()"), quantity = "High against Low",
                   tolerance = 1e-10, key = c("from", "to"))
  expect_reference(as.data.frame(ours), expected, reference = by_definition, covers = "as.data.frame.lsa_comparison()",
                   quantity = "High against Low, as.data.frame()", tolerance = 1e-10, key = c("from", "to"))
})

test_that("the posterior mean of each transition probability in each group, and their difference", {
  prior <- 0.5
  posterior_mean <- \(observed) (observed + prior) / (rowSums(observed) + prior * ncol(observed))
  expected <- ref_cells(prob_a = posterior_mean(tables$High), prob_b = posterior_mean(tables$Low),
                        diff = posterior_mean(tables$High) - posterior_mean(tables$Low))
  ours <- bayes_compare_lsa(grouped, prior = prior, draws = 500L, seed = 1)
  expect_reference(transitions(ours), expected, reference = "independent: mean of the Dirichlet posterior, base R",
                   covers = c("bayes_compare_lsa()", "transitions.lsa_comparison()"), quantity = "High against Low, posterior",
                   tolerance = 1e-12, key = c("from", "to"))
})

test_that("three groups: every pair's log odds ratios and their difference, as the two-group computation per pair", {
  three <- transform(group_regulation_long, Third = sprintf("g%d", as.integer(factor(Actor)) %% 3L))
  third <- tapply(three$Third, three$Actor, \(level) level[[1L]])
  per_third <- lapply(split(log_sequences, third[names(log_sequences)]),
                      \(members) unclass(ref_pairs(members, states = states)))
  log_odds <- function(observed) {
    cell <- ref_collapse(observed)
    add <- 0.5 * (cell$a == 0 | cell$b == 0 | cell$c == 0 | cell$d == 0)
    log(((cell$a + add) * (cell$d + add)) / ((cell$b + add) * (cell$c + add)))
  }
  pairs <- utils::combn(names(per_third), 2L)
  expected <- do.call(rbind, lapply(seq_len(ncol(pairs)), \(i) {
    a <- per_third[[pairs[1L, i]]]
    b <- per_third[[pairs[2L, i]]]
    cbind(group_a = pairs[1L, i], group_b = pairs[2L, i],
          ref_cells(log_or_a = log_odds(a), log_or_b = log_odds(b), diff = log_odds(a) - log_odds(b)))
  }))
  set.seed(5)
  ours <- compare_lsa(lsa(three, actor = "Actor", action = "Action", time = "Time", group = "Third"), R = 19L,
                      measure = "log_or")
  by_definition <- "independent: log odds ratio of the 2 x 2 collapse, base R"
  expect_reference(transitions(ours), expected, reference = by_definition,
                   covers = c("compare_lsa()", "transitions.lsa_comparison_pairwise()"), quantity = "every pair of three groups",
                   tolerance = 1e-10, key = c("group_a", "group_b", "from", "to"))
  expect_reference(as.data.frame(ours), expected, reference = by_definition, covers = "as.data.frame.lsa_comparison_pairwise()",
                   quantity = "every pair of three groups, as.data.frame()", tolerance = 1e-10,
                   key = c("group_a", "group_b", "from", "to"))
})

# The same observed parts against outside software: vcd::loddsratio() computes the log odds ratio of a 2 x 2 table,
# adding 0.5 to every cell when one is zero (its default `correct = any(x == 0)`, Haldane 1956); markovchain's
# Bayesian fit gives the posterior mean of every transition probability under a Dirichlet prior (hyperparameters of
# at least 1, so prior = 1).

#' The log odds ratio of each cell's 2 x 2 collapse by vcd::loddsratio()
vcd_log_odds <- function(observed) {
  cell <- ref_collapse(observed)
  values <- mapply(\(a, b, c, d) unname(stats::coef(vcd::loddsratio(matrix(c(a, c, b, d), 2L), log = TRUE))),
                   cell$a, cell$b, cell$c, cell$d)
  matrix(values, nrow(observed), dimnames = dimnames(observed))
}

test_that("the log odds ratio of each transition in each group equals vcd::loddsratio()", {
  skip_without_reference("vcd")
  high <- vcd_log_odds(tables$High)
  low <- vcd_log_odds(tables$Low)
  expected <- ref_cells(log_or_a = high, log_or_b = low, diff = high - low)
  set.seed(3)
  ours <- compare_lsa(grouped, R = 19L, measure = "log_or")
  expect_reference(transitions(ours), expected, reference = "vcd",
                   covers = c("compare_lsa()", "transitions.lsa_comparison()"), quantity = "High against Low, vcd",
                   tolerance = 1e-10, key = c("from", "to"))
  expect_reference(as.data.frame(ours), expected, reference = "vcd", covers = "as.data.frame.lsa_comparison()",
                   quantity = "High against Low, as.data.frame(), vcd", tolerance = 1e-10, key = c("from", "to"))
})

test_that("three groups: every pair's log odds ratios equal vcd::loddsratio()", {
  skip_without_reference("vcd")
  three <- transform(group_regulation_long, Third = sprintf("g%d", as.integer(factor(Actor)) %% 3L))
  third <- tapply(three$Third, three$Actor, \(level) level[[1L]])
  per_third <- lapply(split(log_sequences, third[names(log_sequences)]),
                      \(members) vcd_log_odds(unclass(ref_pairs(members, states = states))))
  pairs <- utils::combn(names(per_third), 2L)
  expected <- do.call(rbind, lapply(seq_len(ncol(pairs)), \(i) {
    a <- per_third[[pairs[1L, i]]]
    b <- per_third[[pairs[2L, i]]]
    cbind(group_a = pairs[1L, i], group_b = pairs[2L, i], ref_cells(log_or_a = a, log_or_b = b, diff = a - b))
  }))
  set.seed(5)
  ours <- compare_lsa(lsa(three, actor = "Actor", action = "Action", time = "Time", group = "Third"), R = 19L,
                      measure = "log_or")
  expect_reference(transitions(ours), expected, reference = "vcd",
                   covers = c("compare_lsa()", "transitions.lsa_comparison_pairwise()"),
                   quantity = "every pair of three groups, vcd", tolerance = 1e-10,
                   key = c("group_a", "group_b", "from", "to"))
  expect_reference(as.data.frame(ours), expected, reference = "vcd", covers = "as.data.frame.lsa_comparison_pairwise()",
                   quantity = "every pair of three groups, as.data.frame(), vcd", tolerance = 1e-10,
                   key = c("group_a", "group_b", "from", "to"))
})

test_that("the posterior mean of each transition probability in each group equals markovchain's Bayesian fit", {
  skip_without_reference("markovchain")
  groups <- split(log_sequences, achiever[names(log_sequences)])
  posterior <- lapply(groups, \(members) markovchain::markovchainFit(
    data = unname(members), method = "map", possibleStates = states,
    hyperparam = matrix(1, length(states), length(states), dimnames = list(states, states)))$expectedValue)
  expected <- ref_cells(prob_a = posterior$High, prob_b = posterior$Low, diff = posterior$High - posterior$Low)
  ours <- bayes_compare_lsa(grouped, prior = 1, draws = 500L, seed = 1)
  expect_reference(transitions(ours), expected, reference = "markovchain",
                   covers = c("bayes_compare_lsa()", "transitions.lsa_comparison()"),
                   quantity = "High against Low, posterior mean, prior 1", tolerance = 1e-12, key = c("from", "to"))
})

# ---- the drawn parts: Monte Carlo summaries ------------------------------------------------------------------
# compare_lsa(): the group labels of whole sequences are permuted, the log odds ratio of every transition is
# recomputed in both pseudo-groups, and each transition's p-value is (1 + the number of permutations whose
# |difference| reaches the observed one) / (1 + R) (Phipson and Smyth 2010); transitions with fewer than five
# events in the two groups together are not tested; the omnibus statistic is the sum of squared differences over
# the tested transitions. The reference below writes that out in base R on per-sequence counts, with its own seeds.
# A p-value of R permutations is a share of R draws: tolerance 4 standard errors of the difference of two shares,
# per transition, at the reference's share (ref_share_se()), fixed before the run.

test_that("permutation p-values of each transition and the omnibus match label permutation written out (Monte Carlo)", {
  permutations <- 1000L
  per_sequence <- t(vapply(log_sequences, \(events) as.vector(unclass(ref_pairs(list(events), states = states))),
                           numeric(length(states)^2)))
  group <- unname(achiever[names(log_sequences)])
  log_odds_of <- function(counts) {
    observed <- matrix(counts, length(states), dimnames = list(states, states))
    cell <- ref_collapse(observed)
    add <- 0.5 * (cell$a == 0 | cell$b == 0 | cell$c == 0 | cell$d == 0)
    as.vector(log(((cell$a + add) * (cell$d + add)) / ((cell$b + add) * (cell$c + add))))
  }
  difference_of <- \(labels) log_odds_of(colSums(per_sequence[labels == "High", , drop = FALSE])) -
    log_odds_of(colSums(per_sequence[labels == "Low", , drop = FALSE]))
  observed <- difference_of(group)
  tested <- colSums(per_sequence) >= 5
  cells <- \(values) ref_cells(value = matrix(values, length(states), dimnames = list(states, states)))
  rows <- lapply(1:3, \(s) {
    set.seed(s)
    ours <- compare_lsa(grouped, R = permutations, measure = "log_or")
    set.seed(1000 + s)
    null <- t(replicate(permutations, difference_of(sample(group))))
    p_cell <- (1 + colSums(abs(null) >= matrix(abs(observed), permutations, length(observed), byrow = TRUE))) /
      (1 + permutations)
    p_cell[!tested] <- NA_real_
    statistic <- sum(observed[tested]^2)
    p_global <- (1 + sum(rowSums(null[, tested, drop = FALSE]^2) >= statistic)) / (1 + permutations)
    list(ours = transform(subset(transitions(ours), select = c(from, to, p_perm)), seed = s),
         theirs = transform(setNames(cells(p_cell), c("from", "to", "p_perm")), seed = s),
         ours_global = data.frame(seed = s, p_value = as.data.frame(ours, what = "global")$p_value),
         theirs_global = data.frame(seed = s, p_value = p_global))
  })
  bind <- \(part) do.call(rbind, lapply(rows, `[[`, part))
  by_permutation <- "independent: label permutation of whole sequences, Phipson & Smyth (2010) p-value, base R"
  theirs <- bind("theirs")
  expect_each_cell(bind("ours"), theirs, tolerance = 4 * ref_share_se(theirs$p_perm, permutations, permutations),
                   key = c("seed", "from", "to"), value = "p_perm", reference = by_permutation,
                   covers = c("compare_lsa()", "transitions.lsa_comparison()"),
                   quantity = "Monte Carlo summary: permutation p-value, High against Low")
  theirs_global <- bind("theirs_global")
  expect_each_cell(bind("ours_global"), theirs_global,
                   tolerance = 4 * ref_share_se(theirs_global$p_value, permutations, permutations),
                   key = "seed", value = "p_value", reference = by_permutation,
                   covers = c("compare_lsa()", "as.data.frame.lsa_comparison()"),
                   quantity = "Monte Carlo summary: omnibus permutation p-value, High against Low")
})

# bayes_compare_lsa(): each source state's outgoing transitions have the posterior Dirichlet(counts + prior) in each
# group; the difference of a transition's probabilities between the groups is summarised by its 2.5% and 97.5%
# quantiles, the probability of direction (the posterior mass of its more probable sign) and the effect size (the
# exact mean difference over the posterior standard deviation). The reference draws whole rows from the Dirichlet
# posterior with gtools::rdirichlet() (seeded, three seeds). Tolerances per cell, four standard errors of the
# difference of two independent Monte Carlo estimates of `draws` draws each, from the reference's draws, fixed before
# the run: a share (probability of direction) by ref_share_se(); a quantile q of level t by the distance between the
# quantiles at t -/+ sqrt(t (1 - t) / draws), halved (the order-statistic standard error); the effect size by its
# delta-method standard error |effect| sqrt((kurtosis - 1) / (4 draws)), from the standard error of a standard
# deviation.

#' Draws of the Dirichlet posterior of every row, one column per cell in the order of as.vector() of the table
dirichlet_cells <- function(observed, prior, draws) {
  rows <- lapply(seq_len(nrow(observed)), \(i) gtools::rdirichlet(draws, observed[i, ] + prior))
  matrix(aperm(array(unlist(rows), c(draws, ncol(observed), nrow(observed))), c(1L, 3L, 2L)), draws)
}

test_that("bayes_compare_lsa() intervals, probability of direction and effect size match a gtools::rdirichlet() simulation", {
  skip_without_reference("gtools")
  draws <- 20000L
  prior <- 0.5
  level <- 0.025
  quantile_se <- \(values, t) {
    step <- sqrt(t * (1 - t) / draws)
    diff(stats::quantile(values, c(t - step, t + step), names = FALSE)) / 2
  }
  posterior_mean <- \(observed) as.vector((observed + prior) / (rowSums(observed) + prior * ncol(observed)))
  exact_difference <- posterior_mean(tables$High) - posterior_mean(tables$Low)
  cells <- \(values) ref_cells(value = matrix(values, length(states), dimnames = list(states, states)))
  rows <- lapply(1:3, \(s) {
    ours <- transitions(bayes_compare_lsa(grouped, prior = prior, draws = draws, seed = s))
    set.seed(1000 + s)
    difference <- dirichlet_cells(tables$High, prior, draws) - dirichlet_cells(tables$Low, prior, draws)
    positive <- colMeans(difference > 0)
    spread <- apply(difference, 2L, \(values) sqrt(mean((values - mean(values))^2)))
    kurtosis <- apply(difference, 2L, \(values) mean((values - mean(values))^4) / mean((values - mean(values))^2)^2)
    effect <- exact_difference / spread
    reference <- list(
      ci_low = list(value = apply(difference, 2L, stats::quantile, probs = level, names = FALSE),
                    se = sqrt(2) * apply(difference, 2L, quantile_se, t = level)),
      ci_high = list(value = apply(difference, 2L, stats::quantile, probs = 1 - level, names = FALSE),
                     se = sqrt(2) * apply(difference, 2L, quantile_se, t = 1 - level)),
      pd = list(value = pmax(positive, 1 - positive), se = ref_share_se(pmax(positive, 1 - positive), draws, draws)),
      effect_size = list(value = effect, se = sqrt(2) * abs(effect) * sqrt((kurtosis - 1) / (4 * draws))))
    list(ours = transform(ours, seed = s),
         theirs = lapply(reference, \(part) {
           table <- cells(part$value)
           transform(table, se = cells(part$se)$value, seed = s)
         }))
  })
  ours <- do.call(rbind, lapply(rows, `[[`, "ours"))
  invisible(lapply(c("ci_low", "ci_high", "pd", "effect_size"), \(column) {
    theirs <- do.call(rbind, lapply(rows, \(row) row$theirs[[column]]))
    names(theirs)[names(theirs) == "value"] <- column
    expect_each_cell(ours, theirs, tolerance = 4 * theirs$se, key = c("seed", "from", "to"), value = column,
                     reference = "gtools", covers = c("bayes_compare_lsa()", "transitions.lsa_comparison()"),
                     quantity = sprintf("Monte Carlo summary: %s, High against Low, prior 0.5", column))
  }))
})
