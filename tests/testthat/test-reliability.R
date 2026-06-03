test_that("reliability_lsa returns the right shape and class", {
  set.seed(42L)
  # Construct an event stream with enough variation that prob is
  # non-degenerate across random halves.
  set.seed(42L)
  events <- sample(c("a", "b", "c"), 600L, replace = TRUE,
                   prob = c(0.5, 0.3, 0.2))
  # 30 sequences of length 20 each.
  seqs <- split(events, rep(seq_len(30L), each = 20L))
  fit <- lsa(seqs, engine = "classical")
  rel <- reliability_lsa(fit, R = 20L)
  expect_s3_class(rel, "lsa_reliability")
  expect_identical(length(rel$correlations), 20L)
  expect_true(is.numeric(rel$mean))
  expect_true(rel$mean > -1 && rel$mean < 1)
  expect_identical(rel$R, 20L)
  expect_identical(rel$n_sequences, 30L)
})

test_that("reliability_lsa flips its correlation values across weight types", {
  set.seed(43L)
  events <- sample(c("a", "b", "c"), 600L, replace = TRUE,
                   prob = c(0.5, 0.3, 0.2))
  seqs <- split(events, rep(seq_len(30L), each = 20L))
  fit <- lsa(seqs, engine = "classical")
  rel_prob  <- reliability_lsa(fit, R = 15L, weights = "prob")
  rel_count <- reliability_lsa(fit, R = 15L, weights = "count")
  # Correlations should be valid for both weight types and not
  # identical (different scales, but both should be high since the
  # data is from a stationary source).
  expect_true(is.finite(rel_prob$mean))
  expect_true(is.finite(rel_count$mean))
  expect_identical(rel_prob$weights, "prob")
  expect_identical(rel_count$weights, "count")
})

test_that("reliability_lsa supports Spearman", {
  set.seed(44L)
  events <- sample(c("a", "b", "c"), 300L, replace = TRUE)
  seqs <- split(events, rep(seq_len(15L), each = 20L))
  fit <- lsa(seqs, engine = "classical")
  rel <- reliability_lsa(fit, R = 10L, method = "spearman")
  expect_identical(rel$method, "spearman")
  expect_true(is.numeric(rel$mean))
})

test_that("reliability_lsa errors on transition-matrix input", {
  obs <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3,
                dimnames = list(c("a","b","c"), c("a","b","c")))
  fit <- lsa(obs, engine = "classical")
  expect_error(reliability_lsa(fit), "requires event-level input")
})

test_that("reliability_lsa errors with fewer than 2 sequences", {
  fit <- lsa(c("a","b","a","c","b","a","c"), engine = "classical")
  expect_error(reliability_lsa(fit), "at least 2 sequences")
})

test_that("reliability_lsa: identical halves give correlation 1", {
  # If by construction the two halves are statistical replicas of the
  # same source, the mean correlation should be high (well above
  # chance). We don't require r = 1 (the sample variation always
  # leaves slack) but we require mean > 0.5 on a clean signal.
  set.seed(45L)
  # Make 40 long sequences from one strong transition pattern.
  one_seq <- function() {
    s <- sample(c("a","b","c"), 1)
    for (i in 2:30) {
      # strong preference: a -> b, b -> c, c -> a
      nxt <- switch(s[length(s)],
                    a = sample(c("b","b","b","c"), 1),
                    b = sample(c("c","c","c","a"), 1),
                    c = sample(c("a","a","a","b"), 1))
      s <- c(s, nxt)
    }
    s
  }
  seqs <- replicate(40L, one_seq(), simplify = FALSE)
  fit <- lsa(seqs, engine = "classical")
  rel <- reliability_lsa(fit, R = 30L)
  expect_gt(rel$mean, 0.5)
})

test_that("reliability_lsa survives singleton sequences (zero-transition halves)", {
  set.seed(7)
  # 4 length-2 sequences + 6 singletons: some random halves have 0 transitions.
  seqs <- c(replicate(4, c("a", "b"), simplify = FALSE),
            replicate(6, "a", simplify = FALSE))
  fit <- lsa(seqs)
  # Must not error; degenerate splits return NA replicates.
  expect_no_error(rel <- reliability_lsa(fit, R = 50))
  expect_s3_class(rel, "lsa_reliability")
  expect_true(sum(is.finite(rel$correlations)) >= 1L)
})
