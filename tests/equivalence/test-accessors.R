# The reading verbs of a fit and the two data verbs: nodes(), initial(), transition_probabilities(), nobs(),
# lsa_transitions(), lsa_data(). Each against the same quantity taken from a plain table in base R.

sequences <- ref_sequences(engagement)
counts <- unclass(ref_pairs(sequences, lag = 1L))
by_hand <- "independent: table of consecutive pairs, base R"

test_that("nodes(): outgoing and incoming totals are the margins of the table", {
  expected <- data.frame(state = rownames(counts), outgoing = unname(rowSums(counts)), incoming = unname(colSums(counts)))
  expect_reference(nodes(lsa(engagement)), expected, reference = by_hand, covers = c("nodes()", "nodes.lsa()"),
                   quantity = "state totals", tolerance = 0, key = "state")
})

test_that("initial(): the share of sequences that start in each state", {
  first <- table(factor(vapply(sequences, \(events) events[[1L]], character(1)), rownames(counts)))
  expected <- data.frame(state = names(first), init_prob = as.vector(first) / length(sequences))
  expect_reference(initial(lsa(engagement)), expected, reference = by_hand, covers = c("initial()", "initial.lsa()"),
                   quantity = "initial state", tolerance = 1e-12, key = "state")
})

test_that("initial(): equal to the initial probabilities tna estimates", {
  skip_without_reference("tna")
  inits <- tna::tna(engagement)$inits
  expect_reference(initial(lsa(engagement)), data.frame(state = names(inits), init_prob = unname(inits)), reference = "tna",
                   covers = c("initial()", "initial.lsa()"), quantity = "initial state", tolerance = 1e-12, key = "state")
})

test_that("transition_probabilities(): each row of the table over its total", {
  expect_reference(transition_probabilities(lsa(engagement)), ref_cells(prob = prop.table(counts, margin = 1L)), reference = by_hand,
                   covers = c("transition_probabilities()", "transition_probabilities.lsa()"), quantity = "transition probability",
                   tolerance = 1e-12, key = c("from", "to"))
})

test_that("nobs(): the number of transitions", {
  expect_reference(stats::nobs(lsa(engagement)), sum(counts), reference = by_hand, covers = "nobs.lsa()",
                   quantity = "transitions", tolerance = 0)
})

test_that("lsa_transitions(): counts and margins at lag 1 and lag 2", {
  expected <- ref_cells(count = counts, row_total = matrix(rowSums(counts), 3L, 3L, dimnames = dimnames(counts)),
                        col_total = matrix(colSums(counts), 3L, 3L, byrow = TRUE, dimnames = dimnames(counts)))
  expect_reference(as.data.frame(lsa_transitions(engagement)), expected, reference = by_hand,
                   covers = c("lsa_transitions()", "as.data.frame.lsa_transitions()"), quantity = "lag 1", tolerance = 0,
                   key = c("from", "to"))
  expect_reference(as.data.frame(lsa_transitions(engagement, lag = 2)), ref_cells(count = ref_pairs(sequences, lag = 2L)),
                   reference = by_hand, covers = c("lsa_transitions()", "as.data.frame.lsa_transitions()"), quantity = "lag 2",
                   tolerance = 0, key = c("from", "to"))
})

test_that("lsa_data(): one row per event, in order within its sequence", {
  expected <- data.frame(seq_id = rep(seq_along(sequences), lengths(sequences)), index = sequence(lengths(sequences)),
                         state = unlist(sequences))
  expect_reference(as.data.frame(lsa_data(engagement)), expected, reference = "independent: the rows of the data, base R",
                   covers = c("lsa_data()", "as.data.frame.lsa_data()"), quantity = "event", tolerance = 0, key = c("seq_id", "index"))
})

test_that("a grouped fit reads each group's own table", {
  log_sequences <- ref_log_sequences(group_regulation_long, actor = "Actor", action = "Action", time = "Time")
  achiever <- tapply(group_regulation_long$Achiever, group_regulation_long$Actor, \(level) level[[1L]])
  states <- sort(unique(group_regulation_long$Action))
  groups <- split(log_sequences, achiever[names(log_sequences)])
  per_group <- function(table_of) do.call(rbind, lapply(names(groups), \(level) cbind(group = level, table_of(groups[[level]]))))
  cells <- per_group(\(members) {
    pairs <- unclass(ref_pairs(members, lag = 1L, states = states))
    ref_cells(count = pairs, prob = prop.table(pairs, margin = 1L))
  })
  totals <- per_group(\(members) {
    pairs <- unclass(ref_pairs(members, lag = 1L, states = states))
    data.frame(state = states, outgoing = unname(rowSums(pairs)), incoming = unname(colSums(pairs)))
  })
  starts <- per_group(\(members) {
    first <- table(factor(vapply(members, \(events) events[[1L]], character(1)), states))
    data.frame(state = states, init_prob = as.vector(first) / length(members))
  })
  fit <- lsa(group_regulation_long, actor = "Actor", action = "Action", time = "Time", group = "Achiever")
  by_log <- "independent: table of consecutive pairs per group, base R"
  expect_reference(transitions(fit), cells, reference = by_log, covers = c("transitions()", "transitions.lsa_group()", "lsa()"),
                   quantity = "by group", tolerance = 1e-12, key = c("group", "from", "to"))
  expect_reference(as.data.frame(fit), ref_labelled(cells), reference = by_log, covers = "as.data.frame.lsa_group()",
                   quantity = "by group, as.data.frame()", tolerance = 1e-12, key = c("group", "from_label", "to_label"))
  expect_reference(transition_probabilities(fit), subset(cells, select = -count), reference = by_log,
                   covers = c("transition_probabilities()", "transition_probabilities.lsa_group()"), quantity = "by group, probability",
                   tolerance = 1e-12, key = c("group", "from", "to"))
  expect_reference(nodes(fit), totals, reference = by_log, covers = c("nodes()", "nodes.lsa_group()"), quantity = "by group",
                   tolerance = 0, key = c("group", "state"))
  expect_reference(initial(fit), starts, reference = by_log, covers = c("initial()", "initial.lsa_group()"), quantity = "by group",
                   tolerance = 1e-12, key = c("group", "state"))
})
