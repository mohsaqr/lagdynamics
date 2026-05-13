# Third-party validation: Qi An, Wanli Xing, Yuzhe Wang, Xiuyu Li (2026)
# "Behavioural Trajectories and Spatial Responses", Sustainability,
# 18(5), 2326. doi:10.3390/su18052326.
#
# This test runs lagseq's classical engine on the EXACT published
# input matrix (Table 4) and compares its outputs to the EXACT
# published Z-scores and Yule's Q (Table 5), via the shipped
# `qi2026_grandmother` data object.
#
# Result: lagseq agrees with `stats::chisq.test()$stdres` on the
# same input at floating-point precision (< 1e-12). lagseq agrees
# with the paper on 86.7% of cells within 0.05; the 4 cells where it
# disagrees are documented paper typos catalogued in
# `qi2026_grandmother$known_typos`.

LABS <- c("CO","DH","WO","CK","SM","EA","CM","ED","RE","UO")

test_that("qi2026_grandmother data object is well-formed", {
  expect_named(qi2026_grandmother,
               c("obs", "adj_res", "yules_q", "known_typos",
                 "code_descriptions", "source", "license",
                 "n_transitions", "k_states", "notes"),
               ignore.order = TRUE)
  expect_equal(dim(qi2026_grandmother$obs), c(10L, 10L))
  expect_equal(rownames(qi2026_grandmother$obs), LABS)
  expect_equal(sum(qi2026_grandmother$obs), 1531L)
  expect_equal(qi2026_grandmother$n_transitions, 1531L)
  expect_equal(qi2026_grandmother$k_states, 10L)
})

test_that("qi2026: lagseq adjusted residuals == chisq.test()$stdres at 1e-12", {
  # The actual oracle. lagseq's math must equal the standardized
  # Pearson residual formula in stats::chisq.test() on this input.
  fit <- lsa(qi2026_grandmother$obs, engine = "classical")
  suppressWarnings({
    ct <- stats::chisq.test(qi2026_grandmother$obs, correct = FALSE)
  })
  ours <- fit$adj_res[LABS, LABS]
  ok <- !is.na(ours) & !is.na(ct$stdres)
  expect_equal(unname(ours[ok]), unname(ct$stdres[ok]),
               tolerance = 1e-12)
})

test_that("qi2026: lagseq agrees with paper Z-scores on >85% of cells within 0.05", {
  fit <- lsa(qi2026_grandmother$obs, engine = "classical")
  ours <- fit$adj_res[LABS, LABS]
  pap  <- qi2026_grandmother$adj_res[LABS, LABS]
  ok <- !is.na(ours) & !is.na(pap)
  within <- abs(ours[ok] - pap[ok]) < 0.05
  expect_gt(mean(within), 0.85)
})

test_that("qi2026: known paper typos diverge from paper but match the math", {
  fit <- lsa(qi2026_grandmother$obs, engine = "classical")
  typos <- qi2026_grandmother$known_typos
  for (i in seq_len(nrow(typos))) {
    rr <- typos$from[i]; cc <- typos$to[i]
    expect_equal(unname(fit$adj_res[rr, cc]),
                 typos$math_computed[i],
                 tolerance = 0.005,
                 label = sprintf("typo cell %s -> %s: lagseq vs math",
                                 rr, cc))
    expect_gt(abs(unname(fit$adj_res[rr, cc]) -
                  typos$paper_printed[i]),
              0.1,
              label = sprintf("typo cell %s -> %s: lagseq diverges from paper",
                              rr, cc))
  }
})

test_that("qi2026: sign agreement on cells the paper got right", {
  fit <- lsa(qi2026_grandmother$obs, engine = "classical")
  ours <- fit$adj_res[LABS, LABS]
  pap  <- qi2026_grandmother$adj_res[LABS, LABS]
  typo_mask <- matrix(FALSE, 10, 10, dimnames = list(LABS, LABS))
  typos <- qi2026_grandmother$known_typos
  for (i in seq_len(nrow(typos))) {
    typo_mask[typos$from[i], typos$to[i]] <- TRUE
  }
  sig_paper <- abs(pap) >= 1.96
  ok <- !is.na(ours) & sig_paper & !typo_mask
  if (sum(ok) > 0) {
    expect_true(all(sign(ours[ok]) == sign(pap[ok])))
  }
})

test_that("qi2026: Yule's Q matches hand 2x2 collapse on all 100 cells", {
  fit <- lsa(qi2026_grandmother$obs, engine = "classical")
  ours <- fit$yules_q[LABS, LABS]
  obs  <- qi2026_grandmother$obs
  R <- unname(rowSums(obs)); C <- unname(colSums(obs)); N <- sum(obs)
  for (i in seq_along(LABS)) for (j in seq_along(LABS)) {
    a <- unname(obs[i, j])
    b <- R[i] - a; c <- C[j] - a
    d <- N - R[i] - C[j] + a
    num <- a * d - b * c; den <- a * d + b * c
    expected_q <- if (den > 0) num / den else NA_real_
    expect_equal(unname(ours[i, j]), expected_q, tolerance = 1e-12)
  }
})
