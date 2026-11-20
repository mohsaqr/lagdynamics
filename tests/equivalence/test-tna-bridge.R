# lsa_to_tna(): a fit handed to tna as a transition network. With probability weights it must be the model
# tna itself estimates from the same sequences; with count weights, tna's frequency model.
# (Moved from tests/testthat/test-tna-bridge.R, where the class and the warnings are still tested.)

test_that("probability weights and initial probabilities equal tna::tna()", {
  skip_without_reference("tna")
  ours <- lsa_to_tna(lsa(engagement))
  theirs <- tna::tna(engagement)
  expect_reference(ours$weights, theirs$weights, reference = "tna", covers = c("lsa_to_tna()", "lsa_to_tna.lsa()"),
                   quantity = "transition probability", tolerance = 1e-12)
  expect_reference(ours$inits, theirs$inits, reference = "tna", covers = c("lsa_to_tna()", "lsa_to_tna.lsa()"),
                   quantity = "initial probability", tolerance = 1e-12)
})

test_that("count weights equal tna::ftna()", {
  skip_without_reference("tna")
  expect_reference(lsa_to_tna(lsa(engagement), weights = "count")$weights, tna::ftna(engagement)$weights, reference = "tna",
                   covers = c("lsa_to_tna()", "lsa_to_tna.lsa()"), quantity = "transition count", tolerance = 0)
})

test_that("a grouped fit gives each group's probability weights and initial probabilities as tna::group_tna()", {
  skip_without_reference("tna")
  log_sequences <- ref_log_sequences(group_regulation_long, actor = "Actor", action = "Action", time = "Time")
  achiever <- tapply(group_regulation_long$Achiever, group_regulation_long$Actor, \(level) level[[1L]])
  longest <- max(lengths(log_sequences))
  wide <- as.data.frame(do.call(rbind, lapply(log_sequences, \(events) c(events, rep(NA, longest - length(events))))))
  theirs <- tna::group_tna(wide, group = achiever[names(log_sequences)])
  ours <- lsa_to_tna(lsa(group_regulation_long, actor = "Actor", action = "Action", time = "Time", group = "Achiever"))
  per_group <- \(object, part) do.call(rbind, lapply(c("High", "Low"), \(g) {
    values <- object[[g]][[part]]
    if (is.matrix(values)) cbind(group = g, ref_long(values, "value")) else
      data.frame(group = g, state = object[[g]]$labels, value = unname(values))
  }))
  expect_reference(per_group(ours, "weights"), per_group(theirs, "weights"), reference = "tna",
                   covers = "lsa_to_tna.lsa_group()", quantity = "transition probability per group", tolerance = 1e-12,
                   key = c("group", "from", "to"))
  expect_reference(per_group(ours, "inits"), per_group(theirs, "inits"), reference = "tna",
                   covers = "lsa_to_tna.lsa_group()", quantity = "initial probability per group", tolerance = 1e-12,
                   key = c("group", "state"))
})
