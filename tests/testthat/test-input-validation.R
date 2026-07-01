# Regression tests for input-validation hardening. Each guard rejects a
# shape-valid but value-invalid input that previously slipped through to
# the math layer and produced a cryptic error or a plausible-looking but
# wrong result. See CODE_REVIEW_TRANSITION_NETWORKS.md / audit.md.

# --- lsa_ipf: finite, non-negative count table -------------------------

test_that("lsa_ipf() rejects negative observed counts", {
  expect_error(lsa_ipf(matrix(c(-1, 2, 3, 4), 2, 2)),
               "non-negative")
})

test_that("lsa_ipf() rejects non-finite observed counts", {
  expect_error(lsa_ipf(matrix(c(NA, 2, 3, 4), 2, 2)), "finite")
  expect_error(lsa_ipf(matrix(c(Inf, 2, 3, 4), 2, 2)), "finite")
  expect_error(lsa_ipf(matrix(c(NaN, 2, 3, 4), 2, 2)), "finite")
})

test_that("lsa_ipf() still accepts a valid count table", {
  obs <- matrix(c(0, 4, 6, 3, 0, 5, 7, 2, 0), 3, 3, byrow = TRUE)
  out <- lsa_ipf(obs)
  expect_true(out$converged)
  expect_equal(rowSums(out$fit), rowSums(obs), tolerance = 1e-7)
})

# --- lsa() / lsa_data(): finite transition-matrix input ----------------

test_that("transition-matrix input rejects non-finite cells", {
  expect_error(lsa(matrix(c(1, Inf, 2, 3), 2, 2)), "finite")
  expect_error(lsa(matrix(c(1, NA, 2, 3), 2, 2)), "finite")
})

test_that("transition-matrix input still rejects negatives", {
  expect_error(lsa(matrix(c(1, -2, 3, 4), 2, 2)), "non-negative")
})

# --- lsa_data(): whole-number event codes ------------------------------

test_that("fractional numeric event codes are rejected, not truncated", {
  expect_error(lsa_data(c(1.2, 2.8, 1.9), labels = c("A", "B")),
               "whole numbers")
  expect_error(lsa_data(c(1.5, 2, 1)), "whole numbers")
})

test_that("whole-number numeric codes are still accepted", {
  d <- lsa_data(c(1, 2, 1, 3))
  expect_identical(d$events, c(1L, 2L, 1L, 3L))
})

test_that("non-finite numeric event codes are rejected", {
  expect_error(lsa_data(c(1, Inf, 2)), "finite")
})

# --- lsa_transitions(): finite lag -------------------------------------

test_that("non-finite lag is rejected", {
  d <- lsa_data(c("a", "b", "a", "b"))
  expect_error(lsa_transitions(d, lag = Inf))
  expect_error(lsa_transitions(d, lag = NA_real_))
})

test_that("negative and zero lags are still accepted", {
  d <- lsa_data(c("a", "b", "a", "b"))
  expect_s3_class(lsa_transitions(d, lag = -1), "lsa_transitions")
  expect_s3_class(lsa_transitions(d, lag = 0), "lsa_transitions")
})

# --- permute_lsa(): shuffles must be genuine permutations ---------------

test_that("permute_lsa() rejects non-permutation shuffles", {
  fit <- lsa(list(c("a", "b", "a", "b", "a"), c("b", "a", "b", "a")),
             engine = "classical")
  T <- fit$data$n_events
  # duplicated index (drops/duplicates events)
  expect_error(permute_lsa(fit, R = 1, shuffles = list(rep(1L, T))),
               "permutation")
  # wrong length
  expect_error(permute_lsa(fit, R = 1, shuffles = list(seq_len(T - 1L))),
               "permutation")
  # out of range
  bad <- seq_len(T); bad[1] <- T + 1L
  expect_error(permute_lsa(fit, R = 1, shuffles = list(bad)),
               "permutation")
  # NA
  na_perm <- seq_len(T); na_perm[1] <- NA_integer_
  expect_error(permute_lsa(fit, R = 1, shuffles = list(na_perm)),
               "permutation")
})

test_that("permute_lsa() accepts valid permutation shuffles", {
  set.seed(11)
  fit <- lsa(list(c("a", "b", "a", "b", "a"), c("b", "a", "b", "a")),
             engine = "classical")
  T <- fit$data$n_events
  pm <- permute_lsa(fit, R = 2,
                    shuffles = list(sample(T), sample(T)))
  expect_s3_class(pm, "lsa_permutation")
})

# --- bootstrap_lsa(): indices and block_length validation --------------

test_that("bootstrap_lsa() validates replay-index dimensions and range", {
  set.seed(7)
  fit <- lsa(list(c("a", "b", "a", "b", "a"), c("b", "a", "b", "a")),
             engine = "classical")
  S <- fit$data$n_sequences
  # wrong number of columns for sequence-level replay
  expect_error(bootstrap_lsa(fit, R = 1, indices = matrix(1L, 1, 1)),
               "columns")
  # out-of-range sequence index
  expect_error(
    bootstrap_lsa(fit, R = 1, indices = matrix(c(1L, 99L), 1, S)),
    "1\\.\\.")
  # NA index
  expect_error(
    bootstrap_lsa(fit, R = 1, indices = matrix(c(1L, NA), 1, S)),
    "finite")
})

test_that("bootstrap_lsa() accepts a well-formed index matrix", {
  fit <- lsa(list(c("a", "b", "a", "b", "a"), c("b", "a", "b", "a")),
             engine = "classical")
  S <- fit$data$n_sequences
  idx <- matrix(c(1L, 2L, 2L, 1L), nrow = 2, ncol = S, byrow = TRUE)
  bs <- bootstrap_lsa(fit, R = 2, indices = idx)
  expect_equal(dim(bs$indices_used), c(2L, S))
})

test_that("bootstrap_lsa() validates event-level replay width", {
  fit <- lsa(c("a", "b", "a", "b", "a", "b"), engine = "classical")
  T <- fit$data$n_events
  expect_error(
    bootstrap_lsa(fit, R = 1, level = "event",
                  indices = matrix(seq_len(T - 1L), nrow = 1)),
    "event-level replay")
  idx <- matrix(seq_len(T), nrow = 1)
  bs <- bootstrap_lsa(fit, R = 1, level = "event", indices = idx)
  expect_equal(dim(bs$indices_used), c(1L, T))
})

test_that("bootstrap_lsa() rejects an invalid block_length", {
  fit <- lsa(c("a", "b", "a", "b", "a", "b"), engine = "classical")
  expect_error(bootstrap_lsa(fit, R = 1, block_length = -2), "block_length")
  expect_error(bootstrap_lsa(fit, R = 1, block_length = NA), "block_length")
})

test_that("resampling helpers reject non-scalar and non-integer R values", {
  fit <- lsa(list(c("a", "b", "a", "b", "a"), c("b", "a", "b", "a")),
             engine = "classical")
  expect_error(bootstrap_lsa(fit, R = c(1, 2)))
  expect_error(bootstrap_lsa(fit, R = 1.5))
  expect_error(permute_lsa(fit, R = c(1, 2)))
  expect_error(permute_lsa(fit, R = 1.5))
  expect_error(stability_lsa(fit, R = c(1, 2)))
  expect_error(stability_lsa(fit, R = 1.5))
  expect_error(reliability_lsa(fit, R = c(1, 2)))
  expect_error(reliability_lsa(fit, R = 1.5))
})

test_that("bootstrap_lsa() rejects non-scalar confidence levels", {
  fit <- lsa(list(c("a", "b", "a", "b", "a"), c("b", "a", "b", "a")),
             engine = "classical")
  expect_error(bootstrap_lsa(fit, R = 1, level_alpha = c(0.9, 0.95)))
  expect_error(bootstrap_lsa(fit, R = 1, level_alpha = NA_real_))
})

# --- stability_lsa(): degenerate subsamples do not abort the run -------

test_that("stability_lsa() tolerates all-singleton subsamples", {
  set.seed(3)
  fit <- lsa(list("A", c("A", "B", "A", "B")), engine = "classical")
  # Half the sequences are singletons; a subsample can be transition-free.
  # Previously this aborted with "No transitions in input."
  st <- expect_no_error(stability_lsa(fit, R = 30, proportion = 0.5))
  expect_s3_class(st, "lsa_stability")
  expect_equal(nrow(st$edges), fit$data$n_states^2)
})

test_that("loops = FALSE matches 1 - diag(K) and forbids the diagonal", {
  set.seed(7L)
  seqs <- split(sample(c("a", "b", "c"), 400L, replace = TRUE),
                rep(seq_len(20L), each = 20L))
  noloop <- lsa(seqs, loops = FALSE)
  mat    <- lsa(seqs, structural_zeros = 1 - diag(3))
  expect_equal(noloop$exp, mat$exp)
  expect_true(all(diag(noloop$exp) == 0))        # self-transitions forbidden
  expect_true(all(is.na(diag(noloop$adj_res))))  # and not tested
  # loops = TRUE (default) keeps the diagonal estimable.
  expect_false(any(is.na(diag(lsa(seqs)$adj_res))))
  # loops = FALSE also zeros the diagonal of an explicit matrix.
  combo <- lsa(seqs, loops = FALSE, structural_zeros = matrix(1, 3, 3))
  expect_true(all(diag(combo$exp) == 0))
})

test_that("transfer_entropy() validates scalar integer arguments and alignment", {
  expect_error(transfer_entropy(c("a", "b"), lag = c(1, 2), test = "none"),
               "`lag`")
  expect_error(transfer_entropy(c("a", "b"), history = 1.5, test = "none"),
               "`history`")
  expect_error(transfer_entropy(c("a", "b"), R = Inf),
               "`R`")
  expect_error(transfer_entropy(c("a", "b", "c"), c("a", "b"),
                                test = "none"),
               "same length")
})
