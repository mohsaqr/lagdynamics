# compare_lsa() and bayes_compare_lsa(): two groups of the bundled regulation log (high and low achievers).
# The observed part of each is closed form: the log odds ratio of every transition in each group (with 0.5
# added to the four cells of a 2 x 2 collapse that holds a zero, Haldane 1956, Anscombe 1956) and the
# posterior mean of every transition probability under a Dirichlet prior. The permutation p-values and the
# simulated credible intervals are not compared here: they are drawn inside the functions.

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
