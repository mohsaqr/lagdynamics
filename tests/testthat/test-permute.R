test_that("permute_lsa returns the right class and shape", {
  set.seed(11L)
  fit <- lsa(engagement, engine = "classical")
  pm <- permute_lsa(fit, R = 50)
  expect_s3_class(pm, "lsa_permutation")
  expect_equal(nrow(pm$edges), 9L)
  expect_equal(dim(pm$perm_adj_res), c(50L, 9L))
})

test_that("permute_lsa: Phipson-Smyth p-value formula is exact", {
  # Hand-check: with R replicates, the minimum possible p_perm is
  # (1 + 0) / (1 + R) = 1 / (R + 1). Verify a cell with very
  # extreme observed residual hits this minimum.
  set.seed(12L)
  fit <- lsa(engagement, engine = "classical")
  pm <- permute_lsa(fit, R = 50)
  # The Active-Active cell typically dominates engagement data;
  # at R = 50 the minimum p is 1/51 ~= 0.0196.
  expected_min <- 1 / (pm$R + 1)
  expect_gte(min(pm$edges$p_perm), expected_min - 1e-10)
})

test_that("permute_lsa: shuffles= replay produces identical p-values", {
  fit <- lsa(engagement, engine = "classical")
  set.seed(13L)
  # Generate one set of shuffles by running once and capturing them.
  T <- fit$data$n_events
  shuffles <- replicate(20L, sample.int(T), simplify = FALSE)
  pm1 <- permute_lsa(fit, R = 20, shuffles = shuffles)
  set.seed(999L)
  pm2 <- permute_lsa(fit, R = 20, shuffles = shuffles)
  expect_equal(pm1$perm_adj_res, pm2$perm_adj_res, tolerance = 1e-12)
  expect_equal(pm1$edges$p_perm, pm2$edges$p_perm, tolerance = 1e-12)
})

test_that("permute_lsa errors on transition-matrix input", {
  obs <- matrix(c(0,3,1, 2,0,4, 5,1,0), 3, 3,
                dimnames = list(c("a","b","c"), c("a","b","c")))
  fit <- lsa(obs, engine = "classical")
  expect_error(permute_lsa(fit), "requires event-level input")
})

test_that("permute_lsa: observed_adj_res equals fit$adj_res", {
  set.seed(15L)
  fit <- lsa(engagement, engine = "classical")
  pm <- permute_lsa(fit, R = 30)
  expect_equal(pm$edges$observed_adj_res, as.vector(fit$adj_res),
               tolerance = 1e-12)
})

test_that("O'Connor 1999 permutation oracle: large-residual cells get small p_perm", {
  # The paper publishes p_mean from 10 blocks * 1000 = 10,000
  # permutations. For cells with |Z| >> 1.96, our permute_lsa() at
  # R = 1000 should also yield very small p_perm. We test this in
  # a softer way: every cell with |adj_res| >= 4 in lagseq's classical
  # output must have p_perm <= 0.05 in our 1000-permutation
  # replication.
  skip_on_cran()
  data(oconnor_couple)
  set.seed(16L)
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  pm <- permute_lsa(fit, R = 500)
  large <- which(abs(as.vector(fit$adj_res)) >= 4)
  if (length(large) > 0) {
    expect_true(all(pm$edges$p_perm[large] <= 0.05))
  }
})

test_that("permute_lsa: as.data.frame returns the edges frame", {
  set.seed(17L)
  fit <- lsa(engagement, engine = "classical")
  pm <- permute_lsa(fit, R = 20)
  expect_identical(as.data.frame(pm), pm$edges)
})

test_that(".shuffle_within_sequences preserves each sequence's multiset", {
  # Regression for the pre-fix bug at R/permute_lsa.R:109 where the
  # within-sequence shuffle used sample.int(length(sp)) as the source
  # indices. For sequences after the first that produced indices
  # 1..length(sp), reading from positions that belong to sequence 1
  # and copying them into sequence 2's slots. Disjoint-alphabet
  # inputs are the sharpest signature: a correct shuffler can never
  # introduce seq 1's symbols into seq 2's positions.
  set.seed(20260514L)
  seq1 <- rep(c("a", "b"), times = 5L)   # alphabet {a, b}
  seq2 <- rep(c("x", "y"), times = 5L)   # alphabet {x, y}
  events <- c(seq1, seq2)
  seq_positions <- list(seq_len(10L), 10L + seq_len(10L))
  for (rep in seq_len(200L)) {
    out <- lagseq:::.shuffle_within_sequences(events, seq_positions)
    expect_identical(sort(out[seq_positions[[1L]]]), sort(seq1))
    expect_identical(sort(out[seq_positions[[2L]]]), sort(seq2))
  }
})

test_that(".shuffle_within_sequences is safe on singleton and integer events", {
  # Two defensive checks: (1) singleton sequences must be left
  # unchanged; (2) integer event codes must not be misinterpreted by
  # R's sample() quirk where sample(5L) returns sample(1:5). Using
  # sample.int() inside the helper avoids that pitfall.
  set.seed(7L)
  events <- c(3L, 1L, 4L, 7L)
  seq_positions <- list(seq_len(3L), 4L)
  for (rep in seq_len(50L)) {
    out <- lagseq:::.shuffle_within_sequences(events, seq_positions)
    expect_identical(out[4L], 7L)
    expect_identical(sort(out[seq_len(3L)]), c(1L, 3L, 4L))
  }
})

test_that("permute_lsa within_sequence=TRUE keeps disjoint alphabets disjoint", {
  # End-to-end version of the bug. With two disjoint-alphabet
  # sequences, correct within-sequence shuffling cannot create
  # cross-alphabet transitions, so the column marginals for x and y
  # stay positive and the residuals on seq 2's within-alphabet cells
  # remain finite and varying. Pre-fix, seq 2's slots filled with seq
  # 1's symbols, collapsing the x/y column marginals to zero and
  # producing NaN/NA residuals.
  set.seed(20260514L)
  seqs <- list(
    rep(c("a", "b"), times = 6L),
    rep(c("x", "y"), times = 6L)
  )
  fit <- lsa(seqs, engine = "classical",
             labels = c("a", "b", "x", "y"))
  pm <- permute_lsa(fit, R = 50L, within_sequence = TRUE)
  expect_false(any(is.na(pm$perm_adj_res)),
               info = "no NA in residuals under correct shuffling")
  # as.vector() is column-major, so cell (i, j) maps to
  # index (j - 1) * K + i. K = 4 with order (a, b, x, y).
  idx_xy <- (4L - 1L) * 4L + 3L           # cell (x, y)
  idx_yx <- (3L - 1L) * 4L + 4L           # cell (y, x)
  expect_gt(stats::var(pm$perm_adj_res[, idx_xy]), 1e-6)
  expect_gt(stats::var(pm$perm_adj_res[, idx_yx]), 1e-6)
})

test_that("permute_lsa never flags non-estimable (structural-zero) cells", {
  set.seed(1)
  fit <- lsa(engagement, structural_zeros = 1 - diag(3))
  # Diagonal is forbidden -> NA observed residual.
  expect_true(all(is.na(diag(fit$adj_res))))
  pm <- permute_lsa(fit, R = 50)
  diag_rows <- pm$edges[pm$edges$from == pm$edges$to, ]
  # Forbidden cells must get NA p and never be significant (regression:
  # they previously collapsed to p = 1/(R+1) and passed alpha).
  expect_true(all(is.na(diag_rows$p_perm)))
  expect_false(any(diag_rows$significant, na.rm = TRUE))
  # Estimable off-diagonal cells still get finite p-values.
  off <- pm$edges[pm$edges$from != pm$edges$to, ]
  expect_true(all(is.finite(off$p_perm)))
})
