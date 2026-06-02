# Grouped (multi-group) lsa fits and the verbs that dispatch on them.
# A grouped fit is a named list of single-group `lsa` fits sharing one
# global label set; the downstream verbs map over it.

make_group <- function() rep(c("low", "high"), length.out = 136L)

test_that("lsa(group=) returns an lsa_group of per-group lsa fits", {
  fit <- lsa(engagement, engine = "classical", group = make_group())
  expect_s3_class(fit, "lsa_group")
  expect_length(fit, 2L)
  expect_setequal(names(fit), c("high", "low"))
  expect_true(all(vapply(fit, inherits, logical(1L), "lsa")))
})

test_that("group fits share one global label set even if a group lacks a state", {
  # Force a group whose sequences only ever visit a subset of states by
  # building a tiny dataset where group "b" never sees state "c".
  seqs <- list(c("a", "b", "c", "a"), c("a", "b", "a", "b"),
               c("b", "a", "b", "a"))
  fit <- lsa(seqs, group = c("a", "b", "b"))
  labs_a <- rownames(fit$a$obs)
  labs_b <- rownames(fit$b$obs)
  expect_identical(labs_a, labs_b)
  expect_identical(labs_a, c("a", "b", "c"))
  # Group b never visits "c", so its row/col sums there are zero, but
  # the matrix is still full 3 x 3.
  expect_equal(dim(fit$b$obs), c(3L, 3L))
  expect_equal(unname(rowSums(fit$b$obs)["c"]), 0)
})

test_that("group attributes record levels and sizes", {
  fit <- lsa(engagement, group = make_group())
  expect_identical(attr(fit, "levels"), c("high", "low"))
  expect_identical(attr(fit, "group_sizes"), c(68L, 68L))
  expect_identical(attr(fit, "engine"), "classical")
})

test_that("print.lsa_group is informative and returns invisibly", {
  fit <- lsa(engagement, group = make_group())
  expect_output(print(fit), "lsa_group")
  expect_output(print(fit), "groups:")
  expect_identical(withVisible(print(fit))$visible, FALSE)
})

test_that("lsa(group=) validates the grouping vector", {
  expect_error(lsa(engagement, group = c("a", "b")),
               "one entry per sequence")
  expect_error(lsa(engagement, group = c(NA, rep("a", 135L))),
               "must not contain NA")
})

test_that("lsa(group=) rejects transition-matrix input", {
  tm <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3)
  expect_error(lsa(tm, group = c("a", "b", "c")),
               "requires event-level sequence data")
})

test_that("single-group lsa() is unchanged when group is NULL", {
  fit <- lsa(engagement, engine = "classical")
  expect_s3_class(fit, "lsa")
  expect_false(inherits(fit, "lsa_group"))
})

# --- grouped filter helpers --------------------------------------------

test_that("significant_transitions.lsa_group binds a long frame with a group col", {
  fit <- lsa(engagement, group = make_group())
  st <- significant_transitions(fit)
  expect_s3_class(st, "data.frame")
  expect_false(inherits(st, "tbl_df"))   # base-R house style, not tibble
  expect_true("group" %in% names(st))
  # group must be the leading column.
  expect_identical(names(st)[1L], "group")
  # Every group label seen is a real level.
  expect_true(all(st$group %in% names(fit)))
  # Same downstream columns as the single-group edge frame, plus group.
  single_cols <- names(significant_transitions(fit$high))
  expect_setequal(names(st), c("group", single_cols))
})

test_that("grouped over/under filters keep only the right residual sign", {
  fit <- lsa(engagement, group = make_group())
  ov <- overrepresented_transitions(fit)
  un <- underrepresented_transitions(fit)
  expect_true(all(ov$adj_res > 0))
  expect_true(all(un$adj_res < 0))
  expect_true("group" %in% names(ov))
  expect_true("group" %in% names(un))
})

test_that("grouped common_transitions respects min_count per group", {
  fit <- lsa(engagement, group = make_group())
  cm <- common_transitions(fit, min_count = 3L)
  expect_true(all(cm$count >= 3L))
  expect_true("group" %in% names(cm))
})

test_that("grouped filter returns a stable zero-row frame when nothing passes", {
  fit <- lsa(engagement, group = make_group())
  # No transition is observed a billion times, so every group is empty.
  cm <- common_transitions(fit, min_count = 1e9)
  expect_equal(nrow(cm), 0L)
  expect_identical(names(cm)[1L], "group")
  # Columns still match the single-group shape (+ group).
  expect_setequal(names(cm),
                  c("group", names(common_transitions(fit$high))))
})

# --- grouped tna / igraph bridge ---------------------------------------

test_that("lsa_to_tna.lsa_group builds a group_tna with named groups", {
  skip_if_not_installed("tna")
  fit <- lsa(engagement, group = make_group())
  gt <- lsa_to_tna(fit)
  expect_s3_class(gt, "group_tna")
  expect_setequal(names(gt), c("high", "low"))
  expect_true(all(vapply(gt, inherits, logical(1L), "tna")))
})

test_that("as.igraph.lsa_group returns a named list of igraph graphs", {
  skip_if_not_installed("igraph")
  fit <- lsa(engagement, group = make_group())
  gl <- igraph::as.igraph(fit)
  expect_type(gl, "list")
  expect_setequal(names(gl), c("high", "low"))
  expect_true(all(vapply(gl, igraph::is_igraph, logical(1L))))
})

# --- grouped reliability ------------------------------------------------

test_that("reliability_lsa.lsa_group returns a grouped reliability container", {
  set.seed(42)
  fit <- lsa(engagement, group = make_group())
  rel <- reliability_lsa(fit, R = 10L)
  expect_s3_class(rel, "lsa_reliability_group")
  expect_setequal(names(rel), c("high", "low"))
  expect_true(all(vapply(rel, inherits, logical(1L), "lsa_reliability")))
  expect_output(print(rel), "lsa_reliability_group")
})
