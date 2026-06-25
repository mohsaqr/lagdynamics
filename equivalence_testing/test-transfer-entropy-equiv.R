# Reference-equivalence battery for the experimental transfer_entropy().
# Lives here, out of tests/testthat/, so the shipped suite never validates
# against external packages (see README in this folder).
#
# Three independent oracles, all agreeing to machine precision is strong
# evidence the estimator is correct:
#   1. infotheo::condinformation  -- external package, computes I(yf;xp|yh)
#      by a different code path (transfer entropy IS a conditional MI).
#   2. direct double-sum          -- TE from first principles, a separate
#      algebraic route than lagdynamics's four-entropy decomposition.
#   3. analytic                   -- exact values for deterministic / null cases.

# --- reference helpers (deliberately NOT reusing lagdynamics internals) ----------

# Aligned (target future, target history, source past) triples for one pair.
ref_triples <- function(x, y, lag = 1, history = 1) {
  n <- length(y)
  t <- seq.int(history, n - lag)
  hist <- vapply(seq_len(history) - 1L, function(d) y[t - d], character(length(t)))
  dim(hist) <- c(length(t), history)
  data.frame(yf = y[t + lag],
             yh = apply(hist, 1L, paste, collapse = "|"),
             xp = x[t], stringsAsFactors = FALSE)
}

# Transfer entropy by the direct double sum, in bits:
#   TE = sum p(yf,yh,xp) * log2[ p(yf,yh,xp) p(yh) / (p(yh,xp) p(yf,yh)) ]
ref_te_direct <- function(d) {
  pj     <- prop.table(table(d$yf, d$yh, d$xp))
  pyh    <- prop.table(table(d$yh))
  p_yhxp <- prop.table(table(d$yh, d$xp))
  p_yfyh <- prop.table(table(d$yf, d$yh))
  idx <- which(pj > 0, arr.ind = TRUE)
  term <- function(i) {
    f <- rownames(pj)[idx[i, 1]]; h <- colnames(pj)[idx[i, 2]]
    xv <- dimnames(pj)[[3]][idx[i, 3]]; p <- pj[idx[i, 1], idx[i, 2], idx[i, 3]]
    p * log2(p * pyh[h] / (p_yhxp[h, xv] * p_yfyh[f, h]))
  }
  sum(vapply(seq_len(nrow(idx)), term, numeric(1)))
}

# Transfer entropy via infotheo's conditional mutual information (nats -> bits).
ref_te_infotheo <- function(d) {
  infotheo::condinformation(as.integer(factor(d$yf)), as.integer(factor(d$xp)),
                            as.integer(factor(d$yh)), method = "emp") / log(2)
}

te_of <- function(df, from, to) df$te[df$from == from & df$to == to]

# --- a fixed coupled dataset: X drives Y, Y has momentum, no Y -> X ---------
make_data <- function(n = 1500, k = 4, seed = 99) {
  set.seed(seed)
  states <- letters[seq_len(k)]
  x <- sample(states, n, TRUE)
  y <- Reduce(function(py, t) {
    w <- 0.5 + 3 * (states == x[t - 1]) + 1.2 * (states == py)
    sample(states, 1, prob = w / sum(w))
  }, 2:n, accumulate = TRUE, init = states[1])
  list(x = x, y = y)
}

# ---------------------------------------------------------------------------

test_that("bivariate TE matches the direct double-sum oracle (both directions)", {
  d <- make_data()
  te <- transfer_entropy(d$x, d$y, test = "none")
  expect_equal(te_of(te, "d$x", "d$y"),
               ref_te_direct(ref_triples(d$x, d$y)), tolerance = 1e-10)
  expect_equal(te_of(te, "d$y", "d$x"),
               ref_te_direct(ref_triples(d$y, d$x)), tolerance = 1e-10)
})

test_that("bivariate TE matches infotheo::condinformation", {
  skip_if_not_installed("infotheo")
  d <- make_data()
  te <- transfer_entropy(d$x, d$y, test = "none")
  expect_equal(te_of(te, "d$x", "d$y"),
               ref_te_infotheo(ref_triples(d$x, d$y)), tolerance = 1e-8)
  expect_equal(te_of(te, "d$y", "d$x"),
               ref_te_infotheo(ref_triples(d$y, d$x)), tolerance = 1e-8)
})

test_that("history = 2 and lag = 2 match both oracles", {
  skip_if_not_installed("infotheo")
  d <- make_data(n = 2500)
  for (cfg in list(c(h = 2, l = 1), c(h = 1, l = 2), c(h = 2, l = 2))) {
    te <- transfer_entropy(d$x, d$y, history = cfg["h"], lag = cfg["l"], test = "none")
    trip <- ref_triples(d$x, d$y, lag = cfg["l"], history = cfg["h"])
    expect_equal(te_of(te, "d$x", "d$y"), ref_te_direct(trip), tolerance = 1e-10,
                 label = sprintf("direct h=%d l=%d", cfg["h"], cfg["l"]))
    expect_equal(te_of(te, "d$x", "d$y"), ref_te_infotheo(trip), tolerance = 1e-8,
                 label = sprintf("infotheo h=%d l=%d", cfg["h"], cfg["l"]))
  }
})

test_that("deterministic copy Y_t = X_{t-1} gives TE = H(Y_future | Y_now) exactly", {
  set.seed(5)
  states <- letters[1:4]; n <- 1500
  x <- sample(states, n, TRUE)
  y <- c(states[1], x[-n])                       # y_t = x_{t-1}, exactly
  te <- transfer_entropy(x, y, test = "none")
  trip <- ref_triples(x, y)
  H <- function(...) { p <- prop.table(table(...)); p <- p[p > 0]; -sum(p * log2(p)) }
  H_cond <- H(trip$yf, trip$yh) - H(trip$yh)     # H(Y_future | Y_now)
  expect_equal(te_of(te, "x", "y"), H_cond, tolerance = 1e-12)
})

test_that("independent series: effective TE ~ 0 and not significant", {
  set.seed(3)
  x <- sample(letters[1:4], 4000, TRUE)
  y <- sample(letters[1:4], 4000, TRUE)
  te <- transfer_entropy(x, y, test = "surrogate", R = 299, seed = 1)
  expect_true(all(abs(te$te_effective) < 5e-3))
  expect_true(any(te$p > 0.05))
})

test_that("TE is non-negative and normalised TE equals te / H(future|history)", {
  d <- make_data()
  te <- transfer_entropy(d$x, d$y, test = "none")
  expect_true(all(te$te >= -1e-12))
  trip <- ref_triples(d$x, d$y)
  H <- function(...) { p <- prop.table(table(...)); p <- p[p > 0]; -sum(p * log2(p)) }
  leftover <- H(trip$yf, trip$yh) - H(trip$yh)
  expect_equal(te_of(te, "d$x", "d$y") / leftover,
               te$te_normalised[te$from == "d$x" & te$to == "d$y"], tolerance = 1e-10)
})

test_that("state-flow network mode matches the direct oracle on indicator channels", {
  d <- make_data(n = 2000)
  net <- transfer_entropy(d$y, test = "none")               # single sequence -> state network
  ind <- function(v, s) as.character(as.integer(v == s))
  states <- sort(unique(d$y))
  for (s in states) for (r in setdiff(states, s)) {
    ref <- ref_te_direct(ref_triples(ind(d$y, s), ind(d$y, r)))
    # absolute tolerance: some edges have TE ~ 1e-9, where a relative check
    # is meaningless but the two routes still agree to ~1e-14 in absolute terms.
    expect_true(abs(te_of(net, s, r) - ref) < 1e-9,
                info = sprintf("edge %s -> %s", s, r))
  }
})
