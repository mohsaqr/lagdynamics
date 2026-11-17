# Console glyphs fall back to ASCII outside UTF-8; resampling verbs catch
# only the documented degenerate-split condition and count failures; a
# `seed` argument never changes the caller's random stream.

# --- ASCII fallback ------------------------------------------------------

test_that("print.lsa uses ASCII glyphs in a non-UTF-8 session", {
  local_mocked_bindings(.lsa_utf8 = function() FALSE)
  out <- paste(capture.output(print(lsa(engagement))), collapse = "\n")
  expect_match(out, "G^2", fixed = TRUE)
  expect_match(out, "####", fixed = TRUE)
  expect_false(grepl("[^\\x01-\\x7F]", out, perl = TRUE))
})

test_that("print.summary.lsa uses ASCII glyphs in a non-UTF-8 session", {
  local_mocked_bindings(.lsa_utf8 = function() FALSE)
  out <- paste(capture.output(print(summary(lsa(engagement)))),
               collapse = "\n")
  expect_match(out, "X^2", fixed = TRUE)
  expect_false(grepl("[^\\x01-\\x7F]", out, perl = TRUE))
})

test_that("print.lsa uses the UTF-8 glyphs in a UTF-8 session", {
  # Captured output is re-encoded by the session locale, so this half of
  # the contract can only be observed in a UTF-8 session.
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "session is not UTF-8")
  local_mocked_bindings(.lsa_utf8 = function() TRUE)
  out <- paste(capture.output(print(lsa(engagement))), collapse = "\n")
  expect_match(out, "G²", fixed = TRUE)
  expect_match(out, "█", fixed = TRUE)
})

test_that("plot subtitles drop the middle dot outside UTF-8", {
  local_mocked_bindings(.lsa_utf8 = function() FALSE)
  expect_identical(.lsa_ascii("a · b"), "a | b")
  gg <- plot(lsa(engagement), type = "heatmap")
  expect_false(grepl("·", gg$labels$subtitle, fixed = TRUE))
})

# --- failed replicates ---------------------------------------------------

# One sequence with transitions and nine singletons: every split leaves one
# half made of singletons only, so every split-half replicate fails.
singleton_fit <- function() {
  lsa(c(list(c("a", "b", "a", "b", "a")), as.list(rep(c("a", "b"), 4)),
        list("a")))
}

test_that("reliability_lsa() counts failed replicates and reports them", {
  rel <- reliability_lsa(singleton_fit(), R = 12)
  expect_equal(rel$n_failed, 12L)
  expect_equal(as.data.frame(rel, what = "summary")$n_failed, 12L)
  expect_output(print(rel), "failed:        12 of 12 replicates")
})

test_that("reliability_lsa() reports zero failures on ordinary data", {
  set.seed(11)
  rel <- reliability_lsa(lsa(engagement), R = 10)
  expect_equal(rel$n_failed, 0L)
  expect_output(print(rel), "failed:        0 of 10 replicates")
})

test_that("reliability_lsa() lets errors other than no-transitions through", {
  local_mocked_bindings(
    .refit_from_events = function(...) stop("unexpected engine failure")
  )
  expect_error(reliability_lsa(lsa(engagement), R = 2),
               "unexpected engine failure")
})

test_that("stability_lsa() counts failed replicates and reports them", {
  set.seed(5)
  st <- stability_lsa(singleton_fit(), R = 40, proportion = 0.1)
  expect_gt(st$n_failed, 0L)
  expect_lte(st$n_failed, 40L)
  expect_output(print(st), sprintf("failed:        %d of 40", st$n_failed))
})

test_that("stability_lsa() lets errors other than no-transitions through", {
  local_mocked_bindings(
    .refit_from_events = function(...) stop("unexpected engine failure")
  )
  expect_error(stability_lsa(lsa(engagement), R = 2),
               "unexpected engine failure")
})

# --- seed scoping --------------------------------------------------------

test_that("transfer_entropy(seed = ) leaves the caller's RNG stream intact", {
  set.seed(1)
  expected <- stats::runif(1)
  set.seed(1)
  transfer_entropy(engagement, R = 5, seed = 99)
  expect_identical(stats::runif(1), expected)
})

test_that("bayes_compare_lsa(seed = ) leaves the caller's RNG stream intact", {
  g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
           time = "Time", group = "Achiever")
  set.seed(2)
  expected <- stats::runif(1)
  set.seed(2)
  bayes_compare_lsa(g, draws = 20, seed = 7)
  expect_identical(stats::runif(1), expected)
})

test_that("a seeded call is reproducible", {
  a <- transfer_entropy(engagement, R = 5, seed = 3)
  b <- transfer_entropy(engagement, R = 5, seed = 3)
  expect_identical(a, b)
})

# --- comparison and reliability accessors ---------------------------------

test_that("as.data.frame(comparison, what = 'global') is one tidy row", {
  g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
           time = "Time", group = "Achiever")
  cmp <- compare_lsa(g, R = 20)
  gl <- as.data.frame(cmp, what = "global")
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("group_a", "group_b", "statistic", "p_value") %in%
                    names(gl)))
  expect_identical(as.data.frame(cmp), as.data.frame(cmp, what = "edges"))
})

test_that("the engine-registration example engine fits", {
  pearson_engine <- function(transitions, ...) {
    obs <- transitions$obs
    expected <- outer(rowSums(obs), colSums(obs)) / sum(obs)
    residual <- (obs - expected) / sqrt(expected)
    list(obs = obs, exp = expected, prob = obs / rowSums(obs),
         adj_res = residual, p = 2 * stats::pnorm(-abs(residual)))
  }
  register_lsa_engine("pearson_test", pearson_engine,
                      description = "Pearson residuals")
  on.exit(unregister_lsa_engine("pearson_test"), add = TRUE)
  fit <- lsa(engagement, engine = "pearson_test")
  tx <- transitions(fit)
  expected <- outer(rowSums(fit$obs), colSums(fit$obs)) / sum(fit$obs)
  ix <- cbind(match(tx$from, rownames(fit$obs)),
              match(tx$to, colnames(fit$obs)))
  expect_equal(tx$adj_res,
               ((fit$obs - expected) / sqrt(expected))[ix])
})
