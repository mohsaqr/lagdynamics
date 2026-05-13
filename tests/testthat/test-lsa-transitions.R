test_that("transitions count correctly on a hand-computed example", {
  # Sequence: a b a c b a   (length 6, lag 1 → 5 transitions)
  # Pairs: (a,b), (b,a), (a,c), (c,b), (b,a)
  # Counts:
  #   a -> b : 1
  #   a -> c : 1
  #   b -> a : 2
  #   c -> b : 1
  d <- lsa_data(c("a", "b", "a", "c", "b", "a"))
  tx <- lsa_transitions(d, lag = 1)
  expect_equal(unname(tx$obs["a", "b"]), 1L)
  expect_equal(unname(tx$obs["a", "c"]), 1L)
  expect_equal(unname(tx$obs["b", "a"]), 2L)
  expect_equal(unname(tx$obs["c", "b"]), 1L)
  expect_equal(sum(tx$obs), 5L)
  expect_equal(tx$n_transitions, 5L)
})

test_that("transitions never span sequence boundaries", {
  d <- lsa_data(list(c("a", "b"), c("c", "a")))
  tx <- lsa_transitions(d, lag = 1)
  # Within-sequence: a->b (seq1), c->a (seq2). No b->c boundary.
  expect_equal(unname(tx$obs["a", "b"]), 1L)
  expect_equal(unname(tx$obs["c", "a"]), 1L)
  expect_equal(sum(tx$obs), 2L)
})

test_that("transitions accept lag = 2", {
  # Sequence: a b c a b c     (length 6, lag 2 → 4 transitions)
  # Pairs at lag 2: (a,c), (b,a), (c,b), (a,c)
  d <- lsa_data(c("a", "b", "c", "a", "b", "c"))
  tx <- lsa_transitions(d, lag = 2)
  expect_equal(unname(tx$obs["a", "c"]), 2L)
  expect_equal(unname(tx$obs["b", "a"]), 1L)
  expect_equal(unname(tx$obs["c", "b"]), 1L)
  expect_equal(sum(tx$obs), 4L)
})

test_that("transitions edge frame layout matches as.vector(obs)", {
  d <- lsa_data(c("a", "b", "a", "c", "b", "a"))
  tx <- lsa_transitions(d)
  # The tidy frame should have count == as.vector(obs)
  expect_identical(tx$edges$count, as.vector(tx$obs))
  expect_identical(nrow(tx$edges), length(tx$obs))
  # And the (from, to) pairs should index correctly
  for (i in seq_len(nrow(tx$edges))) {
    expect_identical(
      tx$edges$count[i],
      tx$obs[tx$edges$from[i], tx$edges$to[i]]
    )
  }
})

test_that("pre-computed transition matrix passes through", {
  tm <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3,
               dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  tx <- lsa_transitions(tm)
  expect_identical(tx$obs, tm)
})

test_that("pre-computed transition matrix rejects non-unit lag", {
  tm <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3)
  expect_error(lsa_transitions(tm, lag = 2),
               "only supports lag = 1")
})
