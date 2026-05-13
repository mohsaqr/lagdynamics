test_that("lsa_data() accepts a character vector", {
  d <- lsa_data(c("a", "b", "a", "c", "b"))
  expect_s3_class(d, "lsa_data")
  expect_identical(d$labels, c("a", "b", "c"))
  expect_identical(d$events, c(1L, 2L, 1L, 3L, 2L))
  expect_identical(d$seq_id, rep(1L, 5L))
  expect_identical(d$n_sequences, 1L)
  expect_identical(d$n_events, 5L)
  expect_identical(d$source, "events")
})

test_that("lsa_data() accepts an integer vector", {
  d <- lsa_data(c(1L, 2L, 1L, 3L, 2L))
  expect_identical(d$labels, c("Code 1", "Code 2", "Code 3"))
  expect_identical(d$events, c(1L, 2L, 1L, 3L, 2L))
})

test_that("lsa_data() accepts a list of sequences", {
  d <- lsa_data(list(c("a", "b", "a"), c("b", "c", "a")))
  expect_identical(d$n_sequences, 2L)
  expect_identical(d$n_events, 6L)
  expect_identical(d$seq_id, c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_identical(d$transitions_per_seq, c(2L, 2L))
})

test_that("lsa_data() accepts a wide matrix", {
  m <- matrix(c("a", "b", "a", NA,
                "b", "a", NA,  NA,
                "c", "b", "a", "b"), nrow = 3, byrow = TRUE)
  d <- lsa_data(m)
  expect_identical(d$n_sequences, 3L)
  expect_identical(d$n_events, 3L + 2L + 4L)
})

test_that("lsa_data() accepts a wide data.frame", {
  df <- data.frame(t1 = c("a", "b"), t2 = c("b", "a"),
                   t3 = c("a", "c"), stringsAsFactors = FALSE)
  d <- lsa_data(df)
  expect_identical(d$n_sequences, 2L)
  expect_identical(d$n_events, 6L)
})

test_that("lsa_data() preserves a pre-computed transition matrix", {
  tm <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3,
               dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  d <- lsa_data(tm)
  expect_identical(d$source, "transitions")
  expect_identical(d$n_states, 3L)
  expect_identical(d$labels, c("a", "b", "c"))
  expect_identical(d$obs_input, tm)
})

test_that("lsa_data() supplied labels override defaults", {
  d <- lsa_data(c(1L, 2L, 1L), labels = c("foo", "bar"))
  expect_identical(d$labels, c("foo", "bar"))
})

test_that("lsa_data() with extra labels treats them as a fixed alphabet", {
  d <- lsa_data(c("a", "b"), labels = c("a", "b", "c"))
  expect_identical(d$labels, c("a", "b", "c"))
  expect_identical(d$n_states, 3L)
})

test_that("lsa_data() rejects duplicated labels", {
  expect_error(lsa_data(c("a", "b"), labels = c("a", "a")),
               "must be unique")
})

test_that("lsa_data() rejects characters not in supplied labels", {
  expect_error(lsa_data(c("a", "z"), labels = c("a", "b")),
               "Some events do not match supplied labels")
})

test_that("lsa_data() rejects integer codes exceeding supplied labels", {
  expect_error(lsa_data(c(1L, 5L), labels = c("a", "b")),
               "exceeds `labels` length")
})

test_that("lsa_data() rejects empty input", {
  expect_error(lsa_data(character(0)), "zero events")
  expect_error(lsa_data(list(character(0), character(0))),
               "All empty after NA removal|zero events|No usable sequences")
})

test_that("lsa_data() is idempotent", {
  d <- lsa_data(c("a", "b", "a"))
  expect_identical(lsa_data(d), d)
})
