# The lsa object's initial-state distribution slot and the Nestimate-
# style print.

test_that("lsa fit carries an inits slot that sums to 1", {
  fit <- lsa(engagement, engine = "classical")
  expect_type(fit$inits, "double")
  expect_named(fit$inits, fit$data$labels)
  expect_equal(sum(fit$inits), 1)
  # inits = proportion of sequences starting in each state.
  per <- split(fit$data$events,
               factor(fit$data$seq_id,
                      levels = seq_len(fit$data$n_sequences)))
  first <- vapply(per, function(s) s[1L], integer(1L))
  expect_equal(unname(fit$inits),
               unname(as.numeric(table(factor(first,
                 levels = seq_len(fit$data$n_states))) / length(first))))
})

test_that("transition-matrix fits have NULL inits", {
  fit <- lsa(matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3))
  expect_null(fit$inits)
})

test_that("print.lsa shows header, independence test, and inits bars", {
  fit <- lsa(engagement, engine = "classical")
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Lag Sequential Analysis")
  expect_match(out, "independence: G")
  expect_match(out, "Initial states:")
  expect_match(out, .lsa_glyph("\u2588", "#"))   # the bar glyph
  expect_identical(withVisible(print(fit))$visible, FALSE)
})

test_that("print.lsa omits the initial-states block for matrix fits", {
  fit <- lsa(matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3))
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_false(grepl("Initial states", out))
})

test_that("print.summary.lsa shows the overview and the state activity", {
  fit <- lsa(engagement, engine = "classical")
  out <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(out, "Lag sequential analysis: classical engine, lag 1")
  expect_match(out, "State activity")
  expect_match(out, "136 sequences")
  # summary() itself prints nothing: the object prints only when shown.
  expect_silent(s <- summary(fit))
  expect_identical(withVisible(summary(fit))$visible, TRUE)
})

test_that("print() counts significant transitions out of the cells tested", {
  # Three states and no self-transition in the data. With loops = FALSE the
  # three self-transitions are structural zeros: 9 - 3 = 6 cells are tested.
  seqs <- data.frame(a = c("A", "B", "C", "A", "C", "B", "A"),
                     b = c("B", "C", "A", "C", "B", "A", "B"),
                     c = c("C", "A", "B", "B", "A", "C", "C"),
                     d = c("A", "B", "C", "A", "B", "A", "A"))
  header <- function(fit) {
    grep("Significant transitions", utils::capture.output(print(fit)), value = TRUE)
  }
  all_cells <- lsa(seqs)
  expect_match(header(all_cells), " of 9$")
  no_loops <- lsa(seqs, loops = FALSE)
  expect_match(header(no_loops), " of 6$")
})
