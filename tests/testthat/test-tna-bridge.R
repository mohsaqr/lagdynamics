# lsa_to_tna() hands a fit to tna: with prob weights it must reproduce the
# model tna itself estimates from the same sequences.

# The comparison of the weights and initial probabilities with tna::tna() is a declared equivalence test:
# tests/equivalence/test-tna-bridge.R.

test_that("lsa_to_tna() returns a tna object with the states as labels", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement)
  net <- lsa_to_tna(fit)
  expect_s3_class(net, "tna")
  expect_identical(net$labels, nodes(fit)$state)
})

test_that("lsa_to_tna() count weights are the observed counts", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement)
  net <- lsa_to_tna(fit, weights = "count")
  expect_identical(attr(net, "type"), "frequency")
  expect_equal(sum(net$weights), nobs(fit))
})

test_that("non-resampleable weights warn with a class and carry no data", {
  skip_if_not_installed("tna")
  expect_warning(net <- lsa_to_tna(lsa(engagement), weights = "adj_res"),
                 class = "lagdynamics_tna_no_data")
  expect_null(net$data)
  expect_true(all(net$weights >= 0))
})

test_that("a grouped fit becomes a group_tna with one network per group", {
  skip_if_not_installed("tna")
  g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
           time = "Time", group = "Achiever")
  nets <- lsa_to_tna(g)
  expect_s3_class(nets, "group_tna")
  expect_setequal(names(nets), c("High", "Low"))
})

test_that("sequence-data colours come from the Okabe-Ito palette", {
  skip_if_not_installed("tna")
  net <- lsa_to_tna(lsa(engagement))
  expect_true(all(attr(net$data, "colors") %in% .lsa_okabe_ito))
})
