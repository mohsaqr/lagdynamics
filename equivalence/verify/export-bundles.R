# Exports lagdynamics comparison cases as bundles for Prova Verify (format:
# plan/24-verify.md in the prova repository), with Prova's R writer.
#
# Every side of a case is an R script, stored in the bundle and run there in its
# own R process: its console is the side's raw output and the table it leaves in
# `result` holds the side's values. The references are the ones the package is
# already checked against (equivonly/EQUIVALENCE.md, tests/testthat/): base R
# and stats (chisq.test, loglin, glm, pnorm), tna and infotheo.
#
# Run from the project root:
#   PROVA_VERIFY_R=<prova>/reporters/verify.R \
#     Rscript equivalence/verify/export-bundles.R [case ...]
source(Sys.getenv("PROVA_VERIFY_R", file.path("..", "prova", "reporters", "verify.R")))

.verify_base <- file.path("equivalence", "verify")
.verify_generator <- "lagdynamics equivalence/verify/export-bundles.R"

# ---- data -------------------------------------------------------------------

#' Load one of the package's bundled datasets from `data/`
#'
#' @param name The dataset's name.
#' @return The dataset.
.verify_dataset <- function(name) {
  env <- new.env()
  load(file.path("data", paste0(name, ".rda")), envir = env)
  get(name, envir = env)
}

#' The data files the cases read, written to a temporary folder
#'
#' `engagement.csv`: the bundled `engagement` sequences (138 rows, 3 states).
#' `regulation-no-repeats.tsv`: lag-1 counts of the bundled `group_regulation`
#' sequences after collapsing immediate repeats, so that the diagonal is
#' impossible by construction (a structural zero).
#' `three-codes.tsv`: the 3 x 3 table the package's engine tests use.
#' `coupled-streams.csv`: two coupled streams (x drives y), the generator of
#' equivonly/test-transfer-entropy-equiv.R.
#' @return A named character vector of paths.
.verify_data <- function() {
  folder <- file.path(tempdir(), "lagdynamics-verify")
  dir.create(folder, showWarnings = FALSE)
  paths <- file.path(folder, c("engagement.csv", "regulation-no-repeats.tsv", "three-codes.tsv",
                               "coupled-streams.csv"))
  names(paths) <- c("engagement", "regulation", "three", "streams")

  utils::write.csv(.verify_dataset("engagement"), paths[["engagement"]], row.names = FALSE, na = "")

  regulation <- .verify_dataset("group_regulation")
  states <- sort(unique(stats::na.omit(unlist(regulation))))
  pairs <- do.call(rbind, lapply(seq_len(nrow(regulation)), function(row) {
    events <- rle(as.character(stats::na.omit(unlist(regulation[row, ]))))$values
    data.frame(from = utils::head(events, -1L), to = utils::tail(events, -1L))
  }))
  counts <- table(factor(pairs$from, states), factor(pairs$to, states))
  stopifnot("collapsed sequences must have no self-transitions" = all(diag(counts) == 0))
  utils::write.table(unclass(counts), paths[["regulation"]], sep = "\t", quote = FALSE, col.names = NA)

  three <- matrix(c(12, 8, 5, 4, 15, 6, 3, 7, 20), 3L, 3L, byrow = TRUE,
                  dimnames = list(c("x", "y", "z"), c("x", "y", "z")))
  utils::write.table(three, paths[["three"]], sep = "\t", quote = FALSE, col.names = NA)

  set.seed(99)
  codes <- letters[seq_len(4L)]
  x <- sample(codes, 1500L, replace = TRUE)
  y <- Reduce(function(previous, t) {
    weight <- 0.5 + 3 * (codes == x[t - 1L]) + 1.2 * (codes == previous)
    sample(codes, 1L, prob = weight / sum(weight))
  }, seq.int(2L, 1500L), accumulate = TRUE, init = codes[1L])
  utils::write.csv(data.frame(x = x, y = y), paths[["streams"]], row.names = FALSE)
  paths
}

# ---- script fragments -------------------------------------------------------

.read_engagement <- "sequences <- read.csv(\"../data/engagement.csv\", na.strings = \"\")"
.read_three <- "observed <- as.matrix(read.delim(\"../data/three-codes.tsv\", row.names = 1))"
.read_regulation <- "observed <- as.matrix(read.delim(\"../data/regulation-no-repeats.tsv\", row.names = 1))"
.read_streams <- "streams <- read.csv(\"../data/coupled-streams.csv\")"

# An independent count of lagged pairs: within each row (one sequence), the
# event at position t is paired with the event at t + lag.
.count_pairs <- c(
  "states <- sort(unique(na.omit(unlist(sequences))))",
  "count_pairs <- function(lag) {",
  "  pairs <- do.call(rbind, lapply(seq_len(nrow(sequences)), function(row) {",
  "    events <- as.character(na.omit(unlist(sequences[row, ])))",
  "    data.frame(from = head(events, -lag), to = tail(events, -lag))",
  "  }))",
  "  table(from = factor(pairs$from, states), to = factor(pairs$to, states))",
  "}")
.count_lag1 <- c(.count_pairs,
  "counts <- count_pairs(lag = 1)",
  "counts",
  "cells <- as.data.frame(counts, responseName = \"count\")",
  "transition <- paste(cells$from, \"->\", cells$to)")

# The package's tidy table of transitions, keyed by one column.
.tidy_transition <- "  transform(transition = paste(from, \"->\", to)) |>"

.lagdynamics <- function(script, key = "transition") {
  verify_side("lagdynamics", package = "lagdynamics", version = read.dcf("DESCRIPTION")[1L, "Version"],
              primary = TRUE, load = ".", key = key, script = c("library(lagdynamics)", script))
}
.base <- function(script, key = "transition") {
  verify_side("base", name = "Base R, from the definition", version = as.character(getRversion()),
              authority = "independent", key = key, script = script)
}
.stats <- function(script, key = "transition", name = "stats") {
  verify_side("stats", name = name, package = "stats", authority = "open-source", key = key, script = script)
}

# ---- cases ------------------------------------------------------------------

#' Every case: its words, data, tolerance and the script of each side
#'
#' Each tolerance was fixed before the case was first run; the comment beside
#' it says why that value.
#' @return A named list of cases.
.verify_cases <- function() {
  list(
    lag1_counts_probabilities = list(
      title = "Lag-1 transition counts and probabilities",
      description = "How often each engagement state follows each other state one step later, and the row-conditional probability of that transition, for the 138 bundled engagement sequences. Counted by lsa(), by a plain table of consecutive pairs, and by tna (frequency and relative models).",
      covers = c("lsa()", "transitions()"),
      requires = "tna", data = "engagement",
      # Counts are integers and a probability is one division of two integers:
      # the sides can differ by rounding of the last bit only.
      tolerance = 1e-12, measure = "absolute",
      sides = list(
        .lagdynamics(c(
          .read_engagement,
          "fit <- lsa(sequences, lag = 1)",
          "fit",
          "result <- transitions(fit) |>",
          .tidy_transition,
          "  subset(select = c(transition, count, prob))",
          "result")),
        .base(c(
          .read_engagement, .count_lag1,
          "result <- data.frame(transition = transition, count = cells$count,",
          "                     prob = as.vector(counts / rowSums(counts)))",
          "result")),
        verify_side("tna", package = "tna", authority = "open-source", key = "transition", script = c(
          "library(tna)",
          .read_engagement,
          "frequencies <- ftna(sequences)$weights",
          "probabilities <- tna(sequences)$weights",
          "frequencies",
          "probabilities",
          "states <- rownames(frequencies)",
          "result <- data.frame(transition = paste(states[row(frequencies)], \"->\", states[col(frequencies)]),",
          "                     count = as.vector(frequencies), prob = as.vector(probabilities))",
          "result")))),

    adjusted_residuals = list(
      title = "Expected counts, adjusted residuals and their p-values",
      description = "For the lag-1 table of the engagement sequences: the count expected under independence, the adjusted (standardised) residual of Haberman (1973) and its two-sided normal p-value. lsa() against stats::chisq.test()$expected and $stdres on the same table, and against the formulas written out.",
      covers = c("lsa()", "transitions()", "lsa_classical()"),
      data = "engagement",
      # An analytic identity (the tolerance of equivonly/EQUIVALENCE.md for
      # residuals). Relative, because the p-values run down to 1e-104.
      tolerance = 1e-10, measure = "relative",
      sides = list(
        .lagdynamics(c(
          .read_engagement,
          "fit <- lsa(sequences, lag = 1, engine = \"classical\", alternative = \"two.sided\")",
          "result <- transitions(fit) |>",
          .tidy_transition,
          "  subset(select = c(transition, expected, adj_res, p))",
          "result")),
        .stats(name = "stats::chisq.test()", c(
          .read_engagement, .count_lag1,
          "test <- chisq.test(counts, correct = FALSE)",
          "test",
          "result <- data.frame(transition = transition, expected = as.vector(test$expected),",
          "                     adj_res = as.vector(test$stdres),",
          "                     p = as.vector(2 * pnorm(-abs(test$stdres))))",
          "result")),
        .base(c(
          .read_engagement, .count_lag1,
          "observed <- unclass(counts)",
          "total <- sum(observed)",
          "row_share <- rowSums(observed) / total",
          "column_share <- colSums(observed) / total",
          "expected <- total * outer(row_share, column_share)",
          "# Haberman (1973): z = (O - E) / sqrt(E (1 - row share) (1 - column share))",
          "z <- (observed - expected) / sqrt(expected * outer(1 - row_share, 1 - column_share))",
          "result <- data.frame(transition = transition, expected = as.vector(expected),",
          "                     adj_res = as.vector(z), p = as.vector(2 * pnorm(abs(z), lower.tail = FALSE)))",
          "result")))),

    independence_tests = list(
      title = "Tablewise tests of independence",
      description = "The likelihood-ratio (G2) and Pearson (X2) chi-square statistics of the lag-1 engagement table, their degrees of freedom and p-values. tests() against stats::chisq.test() with G2 from its definition, and against stats::loglin().",
      covers = c("tests()", "lsa()"),
      data = "engagement",
      # Closed-form statistics on the same table; relative because the
      # p-values are of the order 1e-132.
      tolerance = 1e-10, measure = "relative",
      sides = list(
        .lagdynamics(key = "test", c(
          .read_engagement,
          "fit <- lsa(sequences, lag = 1)",
          "result <- tests(fit)",
          "result")),
        .stats(key = "test", name = "stats::chisq.test()", c(
          .read_engagement, .count_lag1,
          "pearson <- chisq.test(counts, correct = FALSE)",
          "pearson",
          "g2 <- 2 * sum(counts * log(counts / pearson$expected))",
          "df <- unname(pearson$parameter)",
          "result <- data.frame(test = c(\"lrx2\", \"x2\"), statistic = c(g2, unname(pearson$statistic)),",
          "                     df = c(df, df),",
          "                     p = c(pchisq(g2, df, lower.tail = FALSE), pearson$p.value))",
          "result")),
        verify_side("loglin", name = "stats::loglin()", package = "stats", authority = "open-source",
                    key = "test", script = c(
          .read_engagement, .count_lag1,
          "model <- loglin(counts, margin = list(1, 2), fit = TRUE, eps = 1e-10, print = FALSE)",
          "model$lrt",
          "model$pearson",
          "model$df",
          "result <- data.frame(test = c(\"lrx2\", \"x2\"), statistic = c(model$lrt, model$pearson),",
          "                     df = c(model$df, model$df),",
          "                     p = pchisq(c(model$lrt, model$pearson), model$df, lower.tail = FALSE))",
          "result")))),

    association_indices = list(
      title = "Yule's Q, lift and column-conditional probability",
      description = "Three per-transition association indices of the lag-1 engagement table: Yule's Q of the 2 x 2 collapse of each cell, lift (observed over expected), and the probability of the antecedent given the consequent. transitions() against the definitions written out in base R.",
      covers = c("transitions()", "lsa()"),
      data = "engagement",
      # Ratios of small integer expressions (tolerance of the package's own
      # Yule's Q identity test is 1e-12).
      tolerance = 1e-12, measure = "absolute",
      sides = list(
        .lagdynamics(c(
          .read_engagement,
          "fit <- lsa(sequences, lag = 1)",
          "result <- transitions(fit) |>",
          .tidy_transition,
          "  subset(select = c(transition, yules_q, lift, prob_col))",
          "result")),
        .base(c(
          .read_engagement, .count_lag1,
          "observed <- unclass(counts)",
          "total <- sum(observed)",
          "row_total <- rowSums(observed)[row(observed)]",
          "column_total <- colSums(observed)[col(observed)]",
          "# The 2 x 2 collapse of a cell: n_a = given then target, n_b = given then other,",
          "# n_c = other then target, n_d = other then other. Q = (ad - bc) / (ad + bc).",
          "n_a <- as.vector(observed)",
          "n_b <- row_total - n_a",
          "n_c <- column_total - n_a",
          "n_d <- total - row_total - column_total + n_a",
          "result <- data.frame(transition = transition, yules_q = (n_a * n_d - n_b * n_c) / (n_a * n_d + n_b * n_c),",
          "                     lift = n_a / (row_total * column_total / total), prob_col = n_a / column_total)",
          "result")))),

    structural_zeros = list(
      title = "Structural zeros: expected counts and adjusted residuals",
      description = "A 9 x 9 lag-1 table of the group-regulation sequences with immediate repeats collapsed, so that a state cannot follow itself. lsa(loops = FALSE) fits quasi-independence by iterative proportional fitting; the reference is the same model as a Poisson regression on the 72 possible cells (stats::glm), whose fitted values are the expected counts and whose standardised Pearson residuals are the adjusted residuals.",
      covers = c("lsa()", "lsa_ipf()", "transitions()"),
      data = "regulation",
      # Two iterative fits of the same model (IPF to a margin difference of
      # 1e-8; IRLS): the tolerance of equivonly/EQUIVALENCE.md for
      # structural-zero expected frequencies.
      tolerance = 1e-6, measure = "absolute",
      sides = list(
        .lagdynamics(c(
          .read_regulation,
          "fit <- lsa(observed, lag = 1, loops = FALSE)",
          "fit",
          "result <- transitions(fit) |>",
          "  subset(from != to) |>",
          .tidy_transition,
          "  subset(select = c(transition, expected, adj_res))",
          "result")),
        .stats(name = "stats::glm() (quasi-independence)", c(
          .read_regulation,
          "cells <- as.data.frame(as.table(observed), responseName = \"count\")",
          "names(cells) <- c(\"from\", \"to\", \"count\")",
          "possible <- subset(cells, from != to)",
          "model <- glm(count ~ from + to, family = poisson, data = possible,",
          "             control = glm.control(epsilon = 1e-14, maxit = 100))",
          "model",
          "result <- data.frame(transition = paste(possible$from, \"->\", possible$to),",
          "                     expected = fitted(model), adj_res = rstandard(model, type = \"pearson\"))",
          "result")))),

    higher_lags = list(
      title = "Lags 1 to 3: counts, probabilities and adjusted residuals",
      description = "The engagement sequences at lags 1, 2 and 3: each event paired with the event one, two and three steps later in the same sequence. lsa_lags() against a plain table of lagged pairs with stats::chisq.test()$stdres on each lag's table.",
      covers = c("lsa_lags()", "transitions()"),
      data = "engagement",
      # Counts are exact; residuals are the analytic identity of the
      # adjusted_residuals case (1e-10 there), here on an absolute scale
      # because counts, probabilities and residuals share the table.
      tolerance = 1e-10, measure = "absolute",
      sides = list(
        .lagdynamics(c(
          .read_engagement,
          "fits <- lsa_lags(sequences, lags = 1:3)",
          "fits",
          "result <- transitions(fits) |>",
          "  transform(transition = sprintf(\"lag %d: %s -> %s\", lag, from, to)) |>",
          "  subset(select = c(transition, count, prob, adj_res))",
          "result")),
        .stats(name = "Base R and stats::chisq.test()", c(
          .read_engagement, .count_pairs,
          "result <- do.call(rbind, lapply(1:3, function(lag) {",
          "  counts <- count_pairs(lag)",
          "  cells <- as.data.frame(counts, responseName = \"count\")",
          "  data.frame(transition = sprintf(\"lag %d: %s -> %s\", lag, cells$from, cells$to),",
          "             count = cells$count, prob = as.vector(prop.table(counts, margin = 1)),",
          "             adj_res = as.vector(chisq.test(counts, correct = FALSE)$stdres))",
          "}))",
          "result")))),

    alternative_engines = list(
      title = "Two-cell, bidirectional and dominance engines",
      description = "The four non-classical engines on the 3 x 3 table of the package's engine tests. Two-cell: odds ratio, log odds ratio, its Wald standard error and z of each cell's 2 x 2 collapse (Bakeman and Gottman 1997). Bidirectional: Haberman residuals of the symmetrised table O + t(O), by stats::chisq.test(). Parallel dominance: (O - t(O)) / sqrt(E + t(E)). Non-parallel dominance: (O - t(O)) / sqrt(O + t(O)).",
      covers = c("lsa()", "lsa_two_cell()", "lsa_bidirectional()", "lsa_parallel_dominance()",
                 "lsa_nonparallel_dominance()", "transitions()"),
      data = "three",
      # Closed-form expressions (the engine tests use 1e-12 for the hand
      # values and 1e-10 against chisq.test); the looser of the two.
      tolerance = 1e-10, measure = "absolute",
      sides = list(
        .lagdynamics(c(
          .read_three,
          "two_cell <- lsa(observed, engine = \"two_cell\") |> transitions() |>",
          "  transform(transition = paste(from, \"->\", to), two_cell_z = adj_res) |>",
          "  subset(select = c(transition, odds_ratio, log_or, log_or_se, two_cell_z))",
          "bidirectional <- lsa(observed, engine = \"bidirectional\") |> transitions() |>",
          "  transform(transition = paste(from, \"->\", to), bidirectional_z = adj_res) |>",
          "  subset(select = c(transition, bidirectional_z))",
          "parallel <- lsa(observed, engine = \"parallel_dominance\") |> transitions() |>",
          "  transform(transition = paste(from, \"->\", to), parallel_dominance_z = adj_res) |>",
          "  subset(select = c(transition, parallel_dominance_z))",
          "nonparallel <- lsa(observed, engine = \"nonparallel_dominance\") |> transitions() |>",
          "  transform(transition = paste(from, \"->\", to), nonparallel_dominance_z = adj_res) |>",
          "  subset(select = c(transition, nonparallel_dominance_z))",
          "result <- Reduce(merge, list(two_cell, bidirectional, parallel, nonparallel))",
          "result")),
        .stats(name = "Base R and stats::chisq.test()", c(
          .read_three,
          "observed",
          "total <- sum(observed)",
          "row_total <- rowSums(observed)[row(observed)]",
          "column_total <- colSums(observed)[col(observed)]",
          "# Two-cell: the 2 x 2 collapse of each cell (no cell is zero here).",
          "n_a <- as.vector(observed)",
          "n_b <- row_total - n_a",
          "n_c <- column_total - n_a",
          "n_d <- total - row_total - column_total + n_a",
          "log_or <- log((n_a * n_d) / (n_b * n_c))",
          "log_or_se <- sqrt(1 / n_a + 1 / n_b + 1 / n_c + 1 / n_d)",
          "# Bidirectional: Haberman residuals of the symmetrised table.",
          "symmetric <- chisq.test(observed + t(observed), correct = FALSE)",
          "# Dominance: the difference between the two directions of a pair.",
          "expected <- outer(rowSums(observed), colSums(observed)) / total",
          "difference <- observed - t(observed)",
          "states <- rownames(observed)",
          "result <- data.frame(",
          "  transition = paste(states[row(observed)], \"->\", states[col(observed)]),",
          "  odds_ratio = exp(log_or), log_or = log_or, log_or_se = log_or_se,",
          "  two_cell_z = log_or / log_or_se,",
          "  bidirectional_z = as.vector(symmetric$stdres),",
          "  parallel_dominance_z = as.vector(difference / sqrt(expected + t(expected))),",
          "  nonparallel_dominance_z = as.vector(difference / sqrt(observed + t(observed))))",
          "result")))),

    transfer_entropy = list(
      title = "Transfer entropy between two streams",
      description = "Two coupled streams of 1500 events on four codes, where x drives y and y does not drive x. Transfer entropy in bits in both directions (lag 1, history 1) and its share of the target's remaining uncertainty. transfer_entropy() against infotheo (transfer entropy is the conditional mutual information I(y future; x past | y past)) and against the direct triple sum over the joint distribution.",
      covers = "transfer_entropy()",
      requires = "infotheo", data = "streams",
      # The tolerance of equivonly/test-transfer-entropy-equiv.R against
      # infotheo (which works in nats and is converted to bits).
      tolerance = 1e-8, measure = "absolute",
      sides = list(
        .lagdynamics(key = "direction", c(
          .read_streams,
          "result <- with(streams, transfer_entropy(x = x, y = y, lag = 1, history = 1, test = \"none\")) |>",
          "  transform(direction = paste(from, \"->\", to)) |>",
          "  subset(select = c(direction, te, te_normalised))",
          "result")),
        verify_side("infotheo", package = "infotheo", authority = "open-source", key = "direction", script = c(
          "library(infotheo)",
          .read_streams,
          "bits <- function(source, target) {",
          "  n <- length(target)",
          "  future <- as.integer(factor(target[-1]))",
          "  past <- as.integer(factor(target[-n]))",
          "  driver <- as.integer(factor(source[-n]))",
          "  te <- condinformation(future, driver, past, method = \"emp\") / log(2)",
          "  c(te = te, te_normalised = te / (condentropy(future, past, method = \"emp\") / log(2)))",
          "}",
          "values <- rbind(bits(streams$x, streams$y), bits(streams$y, streams$x))",
          "result <- data.frame(direction = c(\"x -> y\", \"y -> x\"), values)",
          "result")),
        .base(key = "direction", c(
          .read_streams,
          "# TE = sum p(f, h, s) log2[ p(f, h, s) p(h) / (p(h, s) p(f, h)) ], with f the target's",
          "# next event, h its present event and s the source's present event.",
          "bits <- function(source, target) {",
          "  n <- length(target)",
          "  joint <- prop.table(table(f = target[-1], h = target[-n], s = source[-n]))",
          "  cells <- subset(as.data.frame(joint, responseName = \"p\"), p > 0)",
          "  p_h <- as.vector(apply(joint, \"h\", sum)[as.character(cells$h)])",
          "  p_hs <- apply(joint, c(\"h\", \"s\"), sum)[cbind(as.character(cells$h), as.character(cells$s))]",
          "  p_fh <- apply(joint, c(\"f\", \"h\"), sum)[cbind(as.character(cells$f), as.character(cells$h))]",
          "  te <- sum(cells$p * log2(cells$p * p_h / (p_hs * p_fh)))",
          "  # The target's uncertainty left after its own past: H(f | h).",
          "  pairs <- apply(joint, c(\"f\", \"h\"), sum)",
          "  leftover <- -sum(ifelse(pairs > 0, pairs * log2(pairs / rep(colSums(pairs), each = nrow(pairs))), 0))",
          "  c(te = te, te_normalised = te / leftover)",
          "}",
          "values <- rbind(bits(streams$x, streams$y), bits(streams$y, streams$x))",
          "result <- data.frame(direction = c(\"x -> y\", \"y -> x\"), values)",
          "result"))))
  )
}

# ---- export -----------------------------------------------------------------

#' Export one case
#'
#' A case whose reference package is not installed is skipped with a classed
#' message (`lagdynamics_verify_skipped`) and returns a row of `NA` counts.
#' @param name The case's name.
#' @param data The data files, from `.verify_data()`.
#' @return One row: the bundle, its sides and how its compared values came out.
export_case <- function(name, data = .verify_data()) {
  cases <- .verify_cases()
  stopifnot("`name` must be one exported case" = length(name) == 1L && name %in% names(cases))
  case <- cases[[name]]
  bundle <- file.path(.verify_base, "bundles", name)
  absent <- Filter(function(package) !requireNamespace(package, quietly = TRUE), case$requires %||% character())
  if (length(absent)) {
    message(messageCondition(
      sprintf("Case `%s` skipped: reference package %s is not installed.\n", name, toString(absent)),
      class = "lagdynamics_verify_skipped"))
    return(data.frame(bundle = bundle, sides = NA_integer_, compared = NA_integer_, agree = NA_integer_,
                      disagree = NA_integer_, missing = NA_integer_))
  }
  labels <- c(engagement = "Engagement sequences", regulation = "Lag-1 counts, repeats collapsed",
              three = "3 x 3 transition table", streams = "Two coupled streams")
  files <- stats::setNames(data[case$data], labels[case$data])
  # verify_run() drops a script's repeated lines when it runs it (setdiff()), so a
  # script with two identical lines would not run as it is shown.
  repeated <- vapply(case$sides, function(side) anyDuplicated(side$.script) > 0L, logical(1))
  if (any(repeated)) {
    stop(errorCondition(sprintf("Case `%s` has a script with a repeated line.", name),
                        class = "lagdynamics_verify_repeated_line", call = NULL))
  }
  verify_export(bundle = bundle, title = case$title, description = case$description, sides = case$sides,
                data = files, tolerance = case$tolerance, measure = case$measure, covers = case$covers,
                generator = .verify_generator)
}

requested <- commandArgs(trailingOnly = TRUE)
if (length(requested) == 0L) requested <- names(.verify_cases())
.data_files <- .verify_data()
print(do.call(rbind, lapply(requested, export_case, data = .data_files)), row.names = FALSE)
