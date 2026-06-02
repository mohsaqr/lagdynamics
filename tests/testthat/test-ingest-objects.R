# Recovering event sequences from sibling-package objects so lsa() can
# analyse them directly (tna / group_tna, Nestimate nestimate_data,
# TraMineR stslist). These guard the input-side interoperability and
# the safety rule for sequence-less tna objects.

test_that("lsa() ingests a tna object via its $data sequences", {
  skip_if_not_installed("tna")
  m <- tna::tna(tna::group_regulation)
  fit <- lsa(m)
  expect_s3_class(fit, "lsa")
  expect_identical(fit$data$source, "events")
  # Recovered labels are exactly the tna alphabet.
  expect_setequal(fit$data$labels, attr(m$data, "alphabet"))
  # One sequence per row of the tna_seq_data matrix.
  expect_equal(fit$data$n_sequences, nrow(m$data))
})

test_that("lsa() refuses a tna built from a bare matrix (no sequences)", {
  skip_if_not_installed("tna")
  # A transition-matrix fit has no sequences, so lsa_to_tna() leaves
  # data = NULL -- and feeding that back to lsa() must error rather than
  # fabricate counts from probability weights. (An event-sourced fit, by
  # contrast, round-trips fine because its sequences are recovered.)
  net <- lsa_to_tna(lsa(matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3)))
  expect_null(net$data)
  expect_error(lsa(net), "carries no sequence data")
})

test_that("an event-sourced fit round-trips lsa -> tna -> lsa", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement, engine = "classical")
  back <- lsa(lsa_to_tna(fit))
  expect_s3_class(back, "lsa")
  expect_equal(back$data$n_sequences, fit$data$n_sequences)
  expect_setequal(back$data$labels, fit$data$labels)
})

test_that("lsa() pools sequences from a group_tna", {
  skip_if_not_installed("tna")
  gm <- tryCatch(
    tna::group_model(tna::group_regulation,
                     group = rep(c("a", "b"),
                                 length.out = nrow(tna::group_regulation))),
    error = function(e) skip("group_model unavailable in this tna build")
  )
  skip_if_not(inherits(gm, "group_tna"))
  fit <- lsa(gm)
  expect_s3_class(fit, "lsa")
  # Pooled count equals the sum of per-group sequence counts.
  expect_equal(fit$data$n_sequences,
               sum(vapply(unclass(gm), function(g) nrow(g$data), integer(1L))))
})

test_that("lsa() ingests a Nestimate nestimate_data object", {
  skip_if_not_installed("Nestimate")
  data("group_regulation_long", package = "Nestimate")
  nd <- Nestimate::prepare(group_regulation_long,
                           actor = "Actor", action = "Action")
  expect_s3_class(nd, "nestimate_data")
  fit <- lsa(nd)
  expect_s3_class(fit, "lsa")
  expect_identical(fit$data$source, "events")
  expect_equal(fit$data$n_sequences, nrow(nd$sequence_data))
})

test_that("lsa() ingests a TraMineR stslist", {
  skip_if_not_installed("TraMineR")
  data("mvad", package = "TraMineR")
  ss <- suppressMessages(TraMineR::seqdef(mvad, var = 17:20))
  fit <- lsa(ss)
  expect_s3_class(fit, "lsa")
  expect_equal(fit$data$n_sequences, nrow(ss))
  expect_setequal(fit$data$labels, attr(ss, "alphabet"))
})

test_that("a plain list is still treated as sequences, not misread", {
  # Regression: object dispatch must not hijack ordinary list input.
  fit <- lsa(list(c("a", "b", "a"), c("b", "c", "a", "b")))
  expect_equal(fit$data$n_sequences, 2L)
  expect_setequal(fit$data$labels, c("a", "b", "c"))
})
