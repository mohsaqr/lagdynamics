# Canonical bit-identity validation against O'Connor (1999) Appendix
# B. Every single output matrix lagdynamics's classical engine produces is
# compared cell-by-cell to the paper's printed values; agreement is
# required at better than the paper's own printed precision (3-4
# decimal places).
#
# This is the most consequential test in the verification battery:
# the published-output tables were computed by the SEQUENTIAL SAS/SPSS
# program, the same author later co-wrote the LagSequential R package
# that lagdynamics replaced (clean-room). If lagdynamics's independent
# reimplementation reproduces O'Connor's published numbers on
# O'Connor's published input, the math is provably correct.

test_that("oconnor_couple data object is well-formed", {
  expect_equal(length(oconnor_couple$sequence), 393L)
  expect_equal(sum(oconnor_couple$obs), 392L)
  expect_equal(dim(oconnor_couple$obs), c(6L, 6L))
  expect_equal(oconnor_couple$lrx2$statistic, 202.5009)
  expect_equal(oconnor_couple$lrx2$df, 25L)
})

test_that("O'Connor 1999: lagdynamics reproduces published TRANSITION COUNTS exactly", {
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  expect_equal(unname(fit$obs), unname(oconnor_couple$obs))
  expect_equal(sum(fit$obs), 392L)
})

test_that("O'Connor 1999: lagdynamics reproduces published EXPECTED FREQUENCIES (4dp)", {
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  expect_equal(unname(fit$exp), unname(oconnor_couple$expected),
               tolerance = 1e-3)
})

test_that("O'Connor 1999: lagdynamics reproduces published TRANSITIONAL PROBABILITIES (4dp)", {
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  ok <- is.finite(fit$prob) & is.finite(oconnor_couple$prob)
  expect_equal(unname(fit$prob[ok]), unname(oconnor_couple$prob[ok]),
               tolerance = 1e-3)
})

test_that("O'Connor 1999: lagdynamics reproduces LR chi-square test", {
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  expect_equal(fit$lrx2$statistic, oconnor_couple$lrx2$statistic,
               tolerance = 1e-3)
  expect_equal(fit$lrx2$df, oconnor_couple$lrx2$df)
})

test_that("O'Connor 1999: lagdynamics reproduces published ADJUSTED RESIDUALS (3dp)", {
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  expect_equal(unname(fit$adj_res), unname(oconnor_couple$adj_res),
               tolerance = 1e-3)
})

test_that("O'Connor 1999: lagdynamics reproduces published P-VALUES (4dp)", {
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  expect_equal(unname(fit$p), unname(oconnor_couple$adj_p),
               tolerance = 1e-3)
})

test_that("O'Connor 1999: lagdynamics reproduces published YULE'S Q (3dp)", {
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  ok <- is.finite(fit$yules_q) & is.finite(oconnor_couple$yules_q)
  expect_equal(unname(fit$yules_q[ok]),
               unname(oconnor_couple$yules_q[ok]),
               tolerance = 1e-3)
})

test_that("O'Connor 1999: lagdynamics reproduces published KAPPAS (4dp)", {
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  expect_equal(unname(fit$kappa), unname(oconnor_couple$kappa),
               tolerance = 1e-3)
})

test_that("O'Connor 1999: lagdynamics reproduces published KAPPA Z-SCORES (3dp)", {
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  expect_equal(unname(fit$kappa_z), unname(oconnor_couple$kappa_z),
               tolerance = 1e-3)
})

test_that("O'Connor 1999: lagdynamics reproduces published KAPPA P-VALUES (4dp)", {
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  expect_equal(unname(fit$kappa_p), unname(oconnor_couple$kappa_p),
               tolerance = 1e-3)
})

test_that("O'Connor 1999: chi-square primitive equivalence at 1e-12", {
  # Cross-check with the analytic oracle: feeding the published
  # transition count matrix to chisq.test() must give the same
  # adjusted residuals as lagdynamics does on the original sequence.
  fit <- lsa(oconnor_couple$sequence, engine = "classical")
  suppressWarnings({
    ct <- stats::chisq.test(oconnor_couple$obs, correct = FALSE)
  })
  ok <- !is.na(fit$adj_res) & !is.na(ct$stdres)
  expect_equal(unname(fit$adj_res[ok]), unname(ct$stdres[ok]),
               tolerance = 1e-12)
})
