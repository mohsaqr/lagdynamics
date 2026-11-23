# transfer_entropy(): directed information flow between two categorical streams, in bits. Three references:
# the triple sum over the joint distribution written from the definition; infotheo (transfer entropy is the
# conditional mutual information I(target future; source past | target past)); and exact values for a
# deterministic copy. The surrogate test (te_effective and p) is compared at the end with the source permutation
# written out, a Monte Carlo summary: no established implementation of a surrogate test for discrete transfer
# entropy is installed on the ci machine (RTransferEntropy, IDTxl, pyinform, JIDT are not), so that evidence is
# a computation of this file (Indirect), not an outside program.

streams <- ref_streams()

#' Both directions of a pair of streams by the triple sum, with the share of the target's remaining uncertainty
by_definition <- function(x, y, lag = 1L, history = 1L) {
  one <- function(source, target, from, to) {
    triples <- ref_triples(source, target, lag = lag, history = history)
    te <- ref_transfer_entropy(triples)
    # the target's uncertainty left after its own past: H(future | past)
    leftover <- ref_entropy(triples$future, triples$past) - ref_entropy(triples$past)
    data.frame(from = from, to = to, te = te, te_normalised = te / leftover)
  }
  rbind(one(x, y, "x", "y"), one(y, x, "y", "x"))
}

test_that("both directions equal the triple sum over the joint distribution", {
  ours <- with(streams, transfer_entropy(x = x, y = y, lag = 1, history = 1, test = "none"))
  expect_reference(ours, by_definition(streams$x, streams$y), reference = "independent: triple sum over the joint distribution, base R",
                   covers = "transfer_entropy()", quantity = "lag 1, history 1", tolerance = 1e-10, key = c("from", "to"))
})

test_that("longer history and lag equal the triple sum", {
  longer <- ref_streams(n = 2500L)
  settings <- list(c(history = 2L, lag = 1L), c(history = 1L, lag = 2L), c(history = 2L, lag = 2L))
  invisible(lapply(settings, \(setting) {
    ours <- with(longer, transfer_entropy(x = x, y = y, lag = setting[["lag"]], history = setting[["history"]], test = "none"))
    expect_reference(ours, by_definition(longer$x, longer$y, lag = setting[["lag"]], history = setting[["history"]]),
                     reference = "independent: triple sum over the joint distribution, base R", covers = "transfer_entropy()",
                     quantity = sprintf("lag %d, history %d", setting[["lag"]], setting[["history"]]), tolerance = 1e-10,
                     key = c("from", "to"))
  }))
})

test_that("both directions equal infotheo's conditional mutual information", {
  skip_without_reference("infotheo")
  bits <- function(source, target, from, to) {
    triples <- ref_triples(source, target)
    future <- as.integer(factor(triples$future))
    past <- as.integer(factor(triples$past))
    driver <- as.integer(factor(triples$driver))
    te <- infotheo::condinformation(future, driver, past, method = "emp") / log(2)
    data.frame(from = from, to = to, te = te, te_normalised = te / (infotheo::condentropy(future, past, method = "emp") / log(2)))
  }
  ours <- with(streams, transfer_entropy(x = x, y = y, lag = 1, history = 1, test = "none"))
  # infotheo works in nats and is converted to bits: its own tolerance in the standing battery.
  expect_reference(ours, rbind(bits(streams$x, streams$y, "x", "y"), bits(streams$y, streams$x, "y", "x")), reference = "infotheo",
                   covers = "transfer_entropy()", quantity = "lag 1, history 1", tolerance = 1e-8, key = c("from", "to"))
})

test_that("a deterministic copy y[t] = x[t - 1] transfers all of the target's remaining uncertainty", {
  set.seed(5)
  x <- sample(letters[1:4], 1500L, replace = TRUE)
  y <- c("a", utils::head(x, -1L))
  triples <- ref_triples(x, y)
  # TE(x -> y) = H(future | past) exactly, so the normalised value is 1
  expected <- data.frame(from = "x", to = "y", te = ref_entropy(triples$future, triples$past) - ref_entropy(triples$past), te_normalised = 1)
  expect_reference(subset(transfer_entropy(x = x, y = y, test = "none"), from == "x"), expected,
                   reference = "independent: the analytic value for a deterministic copy", covers = "transfer_entropy()",
                   quantity = "deterministic copy", tolerance = 1e-12, key = c("from", "to"))
})

test_that("one sequence: the state network equals the triple sum on indicator streams", {
  longer <- ref_streams(n = 2000L)
  states <- sort(unique(longer$y))
  indicator <- \(state) as.character(as.integer(longer$y == state))
  pairs <- subset(expand.grid(from = states, to = states, stringsAsFactors = FALSE), from != to)
  pairs$te <- mapply(\(from, to) ref_transfer_entropy(ref_triples(indicator(from), indicator(to))), pairs$from, pairs$to)
  # Absolute: some edges are of the order 1e-9, where the two routes still agree to 1e-14.
  expect_reference(transfer_entropy(longer$y, test = "none"), pairs, reference = "independent: triple sum over the joint distribution, base R",
                   covers = "transfer_entropy()", quantity = "state network", tolerance = 1e-9, key = c("from", "to"))
})

# The surrogate test of transfer_entropy(): the source's values are permuted over the pooled triples, which keeps
# both marginals and breaks the coupling; te_effective is the observed value minus the mean of the surrogates and
# p = (1 + the number of surrogates at or above the observed value) / (R + 1). Written out below on ref_triples()
# with its own seeds (three). Tolerances per number, fixed before the run: te_effective, four standard errors of the
# difference of two means of R surrogates, 4 sd sqrt(2 / R) with sd the reference surrogates' spread; p, four
# standard errors of the difference of two shares (ref_share_se()).
test_that("surrogate test: effective transfer entropy and p match the source permutation written out (Monte Carlo)", {
  surrogates <- 999L
  surrogate_test <- function(source, target, from, to) {
    triples <- ref_triples(source, target)
    observed <- ref_transfer_entropy(triples)
    null <- replicate(surrogates, ref_transfer_entropy(transform(triples, driver = sample(driver))))
    data.frame(from = from, to = to, te_effective = observed - mean(null),
               p = (1 + sum(null >= observed)) / (surrogates + 1), spread = stats::sd(null))
  }
  rows <- lapply(1:3, \(s) {
    ours <- with(streams, transfer_entropy(x = x, y = y, R = surrogates, seed = s))
    set.seed(1000 + s)
    theirs <- rbind(surrogate_test(streams$x, streams$y, "x", "y"), surrogate_test(streams$y, streams$x, "y", "x"))
    list(ours = transform(ours, seed = s), theirs = transform(theirs, seed = s))
  })
  ours <- do.call(rbind, lapply(rows, `[[`, "ours"))
  theirs <- do.call(rbind, lapply(rows, `[[`, "theirs"))
  by_permutation <- "independent: source permutation surrogates (Schreiber 2000), base R"
  expect_each_cell(ours, theirs, tolerance = 4 * theirs$spread * sqrt(2 / surrogates), key = c("seed", "from", "to"),
                   value = "te_effective", reference = by_permutation, covers = "transfer_entropy()",
                   quantity = "Monte Carlo summary: effective transfer entropy, lag 1, history 1")
  expect_each_cell(ours, theirs, tolerance = 4 * ref_share_se(theirs$p, surrogates, surrogates),
                   key = c("seed", "from", "to"), value = "p", reference = by_permutation, covers = "transfer_entropy()",
                   quantity = "Monte Carlo summary: surrogate p-value, lag 1, history 1")
})
