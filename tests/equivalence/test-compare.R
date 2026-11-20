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
