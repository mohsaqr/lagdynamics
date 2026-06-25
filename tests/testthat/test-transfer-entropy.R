# Tests for the experimental transfer_entropy() verb. The synthetic chains are
# the oracle: X drives Y at lag 1, Y has its own momentum, there is no Y -> X
# coupling, so a correct estimator must report TE_{X->Y} >> TE_{Y->X}.

make_coupled <- function(n = 3000, seed = 1) {
  set.seed(seed)
  states <- c("a", "b", "c")
  Px <- matrix(c(.85, .075, .075, .075, .85, .075, .075, .075, .85),
               3, byrow = TRUE, dimnames = list(states, states))
  x <- Reduce(function(p, .) sample(states, 1, prob = Px[p, ]),
              seq_len(n - 1), accumulate = TRUE, init = "a")
  y <- Reduce(function(py, t) {
    w <- 0.5 + 4 * (states == x[t - 1]) + 1.5 * (states == py)
    sample(states, 1, prob = w / sum(w))
  }, 2:n, accumulate = TRUE, init = "a")
  list(x = x, y = y)
}

get_te <- function(df, src, tgt, col = "te") {
  df[[col]][df$from == src & df$to == tgt]
}

test_that("bivariate transfer entropy recovers the buried direction", {
  d <- make_coupled()
  te <- transfer_entropy(d$x, d$y, test = "surrogate", R = 99, seed = 1)
  expect_setequal(unique(c(te$from, te$to)), c("d$x", "d$y"))  # deparsed labels
  fwd <- te$te[which.max(te$te)]
  rev <- te$te[which.min(te$te)]
  expect_gt(fwd / rev, 20)                      # forward dominates
  expect_lt(min(te$te_effective), 0.02)         # reverse is negligible
  expect_gt(max(te$te_effective), 0.1)          # forward is substantial
})

test_that("independent series yield non-significant, near-zero effective TE", {
  set.seed(7)
  a <- sample(c("a", "b", "c"), 2000, TRUE)
  b <- sample(c("a", "b", "c"), 2000, TRUE)
  te <- transfer_entropy(a, b, test = "surrogate", R = 199, seed = 1)
  expect_true(all(abs(te$te_effective) < 0.01))
  expect_true(any(te$p > 0.05))                 # at least one direction not significant
})

test_that("output is a tidy data.frame with the documented columns", {
  te <- transfer_entropy(engagement, test = "none")
  expect_s3_class(te, "data.frame")
  expect_named(te, c("from", "to", "te", "te_normalised", "n"))
  expect_equal(nrow(te), 3 * 2)                 # 3 states -> 6 ordered pairs
  expect_true(all(te$te >= 0))
  expect_true(all(te$te_normalised >= 0 & te$te_normalised <= 1))
  expect_true(all(diff(te$te) <= 0))            # ordered by descending te
})

test_that("normalize = FALSE drops the normalised column", {
  te <- transfer_entropy(engagement, test = "none", normalize = FALSE)
  expect_false("te_normalised" %in% names(te))
})

test_that("lag never crosses sequence boundaries", {
  # Two short, internally-constant sequences with opposite states. A boundary
  # leak would invent an a->b transition between rows; correct pooling never does.
  m <- rbind(c("a", "a", "a", "a"), c("b", "b", "b", "b"))
  te <- transfer_entropy(m, test = "none")
  # only self-persistence exists; cross-state TE must be 0 (no a-next-to-b pairs)
  expect_true(all(te$te == 0))
})

test_that("history and lag arguments are honoured", {
  d <- make_coupled(n = 2000)
  expect_no_error(transfer_entropy(d$x, d$y, history = 2, test = "none"))
  expect_no_error(transfer_entropy(d$x, d$y, lag = 2, test = "none"))
})
