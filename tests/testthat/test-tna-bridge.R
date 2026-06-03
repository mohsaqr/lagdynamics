test_that("lsa_to_tna generic dispatches to lsa_to_tna.lsa", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement, engine = "classical")
  net <- lsa_to_tna(fit)
  expect_s3_class(net, "tna")
})

test_that("lsa_to_tna.lsa exposes prob weights by default", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement, engine = "classical")
  net <- lsa_to_tna(fit)
  expect_identical(dim(net$weights), dim(fit$prob))
  # The prob matrix carried over with non-finite cells coerced to 0.
  fit_prob <- fit$prob
  fit_prob[!is.finite(fit_prob)] <- 0
  expect_equal(unname(net$weights), unname(fit_prob))
  expect_identical(attr(net, "type"), "relative")
})

test_that("lsa_to_tna.lsa weights = 'count' carries observed counts and frequency type", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement, engine = "classical")
  net <- lsa_to_tna(fit, weights = "count")
  # tna::build_model(type = "frequency") demands integer storage, so
  # we compare values without comparing storage mode.
  expect_equal(unname(net$weights), unname(fit$obs),
               ignore_attr = TRUE)
  expect_identical(attr(net, "type"), "frequency")
})

test_that("lsa_to_tna.lsa weights = 'adj_res' clips negatives (tna requires non-negative)", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement, engine = "classical")
  # adj_res omits $data and warns (see dedicated test below).
  net <- suppressWarnings(lsa_to_tna(fit, weights = "adj_res"))
  expect_true(all(net$weights >= 0))
  # Cells where original adj_res was positive should appear unchanged.
  ref <- fit$adj_res
  pos_cells <- is.finite(ref) & ref > 0
  expect_equal(net$weights[pos_cells], ref[pos_cells])
  # Cells where original was negative should be exactly 0.
  neg_cells <- is.finite(ref) & ref < 0
  expect_true(all(net$weights[neg_cells] == 0))
})

test_that("lsa_to_tna.lsa zeros out non-finite cells from structural-zero fits", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement, engine = "classical",
             structural_zeros = 1 - diag(3))
  # Diagonal of adj_res is NA after the priority-2 fix.
  expect_true(all(is.na(diag(fit$adj_res))))
  net <- lsa_to_tna(fit, weights = "prob")
  expect_true(all(diag(net$weights) == 0))
})

test_that("lsa_to_tna.lsa errors with helpful message when tna is absent", {
  if (requireNamespace("tna", quietly = TRUE)) {
    skip("tna is installed; cannot test absence path here")
  }
  fit <- lsa(engagement, engine = "classical")
  expect_error(lsa_to_tna(fit), "Package 'tna' is required")
})

test_that("as.igraph.lsa builds a directed weighted igraph", {
  skip_if_not_installed("igraph")
  fit <- lsa(engagement, engine = "classical")
  g <- igraph::as.igraph(fit)
  expect_true(inherits(g, "igraph"))
  expect_true(igraph::is_directed(g))
  expect_true(igraph::is_weighted(g))
  expect_equal(as.integer(igraph::vcount(g)), nrow(fit$obs))
})

test_that("lsa_to_tna output is consumable by tna::centralities", {
  skip_if_not_installed("tna")
  # lagseq's contract is conversion, not centrality computation: the
  # object from lsa_to_tna() must be accepted by tna's own verbs. Computing
  # centralities is tna's job, so we only check the bridge is fit for it.
  fit <- lsa(engagement, engine = "classical")
  cents <- tna::centralities(lsa_to_tna(fit))
  expect_true(nrow(cents) == nrow(fit$obs))
})

test_that("lsa_to_tna attaches $data + inits so tna's sequence verbs run", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement, engine = "classical")
  net <- lsa_to_tna(fit, weights = "prob")
  # $data is the reconstructed tna_seq_data; inits sum to 1.
  expect_s3_class(net$data, "tna_seq_data")
  expect_equal(nrow(net$data), fit$data$n_sequences)
  expect_equal(sum(net$inits), 1)
  # The resampling verbs that need sequences now work end to end.
  expect_no_error(tna::bootstrap(net, iter = 20))
  expect_no_error(tna::permutation_test(net, net, iter = 20))
})

test_that("lsa_to_tna leaves $data NULL for transition-matrix fits", {
  skip_if_not_installed("tna")
  tm <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3)
  expect_null(lsa_to_tna(lsa(tm))$data)
})

test_that("lsa_to_tna omits $data + warns for residual/lift scales", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement, engine = "classical")
  # prob / count are resampleable -> $data attached, no warning.
  expect_silent(p <- lsa_to_tna(fit, weights = "prob"))
  expect_false(is.null(p$data))
  expect_silent(c <- lsa_to_tna(fit, weights = "count"))
  expect_false(is.null(c$data))
  # adj_res / lift would resample on the wrong scale -> $data omitted + warn.
  for (w in c("adj_res", "lift")) {
    expect_warning(net <- lsa_to_tna(fit, weights = w), "omits the .data slot")
    expect_null(net$data)
    expect_null(net$inits)
  }
})
