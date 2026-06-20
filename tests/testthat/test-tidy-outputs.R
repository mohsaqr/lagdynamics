# Every result object should yield a tidy one-row-per-observation data
# frame via as.data.frame(); grouped objects carry a `group` column.

.tidy_seqs <- function(seed = 1L) {
  set.seed(seed)
  split(sample(c("a", "b", "c"), 600L, replace = TRUE),
        rep(seq_len(30L), each = 20L))
}

test_that("as.data.frame returns a tidy data frame for every result object", {
  seqs <- .tidy_seqs()
  grp  <- rep(c("x", "y"), each = 15L)
  m <- matrix(c(0, 5, 3, 2, 0, 4, 1, 6, 0), 3, 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  objs <- list(
    lsa(seqs), lsa_data(seqs), lsa_data(m), lsa_transitions(seqs),
    lsa_lags(seqs, lags = 1:2), bootstrap_lsa(lsa(seqs), R = 30L),
    permute_lsa(lsa(seqs), R = 30L), stability_lsa(lsa(seqs), R = 30L),
    reliability_lsa(lsa(seqs), R = 15L), certainty_lsa(lsa(seqs)),
    compare_lsa(lsa(seqs, group = grp), R = 30L),
    bayes_compare_lsa(lsa(seqs, group = grp), draws = 300L, seed = 1))
  for (o in objs) {
    df <- as.data.frame(o)
    expect_s3_class(df, "data.frame")
    expect_gt(nrow(df), 0L)
  }
})

test_that("grouped objects tidy to long format with a group column", {
  seqs <- .tidy_seqs()
  grp  <- rep(c("x", "y"), each = 15L)

  gfit <- as.data.frame(lsa(seqs, group = grp))
  expect_true("group" %in% names(gfit))
  expect_setequal(unique(gfit$group), c("x", "y"))
  # One block of edges per group (K^2 = 9 each).
  expect_identical(nrow(gfit), 18L)

  gcert <- as.data.frame(certainty_lsa(lsa(seqs, group = grp)))
  expect_true("group" %in% names(gcert))
  expect_identical(nrow(gcert), 18L)

  grel <- as.data.frame(reliability_lsa(lsa(seqs, group = grp), R = 15L))
  expect_true("group" %in% names(grel))
})

test_that("lsa_data tidies to one row per event (or per cell for a matrix)", {
  d <- lsa_data(list(c("a", "b", "a"), c("b", "c")))
  df <- as.data.frame(d)
  expect_identical(names(df), c("seq_id", "index", "state"))
  expect_identical(nrow(df), 5L)                       # 3 + 2 events
  expect_identical(df$state, c("a", "b", "a", "b", "c"))
  expect_identical(df$index, c(1L, 2L, 3L, 1L, 2L))    # within-sequence

  m <- matrix(c(0, 5, 3, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  dm <- as.data.frame(lsa_data(m))
  expect_identical(names(dm), c("from", "to", "count"))
  expect_identical(nrow(dm), 4L)
})
