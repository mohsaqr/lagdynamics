# Reference computations shared by the equivalence tests, written from the definitions in base R. Nothing
# here calls lagdynamics: each function is one side of a comparison, and the package is the other.

#' The sequences of a wide table: one per row, missing events dropped, empty rows left out
#'
#' @param wide A data frame or matrix, one sequence per row.
#' @return A list of character vectors.
ref_sequences <- function(wide) {
  rows <- lapply(seq_len(nrow(wide)), \(row) as.character(stats::na.omit(unlist(wide[row, ], use.names = FALSE))))
  Filter(length, rows)
}

#' The sequences of a long event log: one per actor, in time order
#'
#' @param log A data frame.
#' @param actor,action,time Column names.
#' @return A named list of character vectors.
ref_log_sequences <- function(log, actor, action, time) {
  ordered <- log[order(log[[actor]], log[[time]]), ]
  split(as.character(ordered[[action]]), ordered[[actor]])
}

#' Count lagged pairs: each event with the event `lag` steps later in the same sequence
#'
#' @param sequences A list of character vectors.
#' @param lag The lag.
#' @param states The states, in the order of the table.
#' @return A table `from` x `to`.
ref_pairs <- function(sequences, lag = 1L, states = sort(unique(unlist(sequences)))) {
  long <- Filter(\(events) length(events) > lag, sequences)
  from <- unlist(lapply(long, \(events) utils::head(events, -lag)))
  to <- unlist(lapply(long, \(events) utils::tail(events, -lag)))
  table(from = factor(from, states), to = factor(to, states))
}

#' A square table as one row per cell
#'
#' @param cells A matrix or table with the states as dimnames.
#' @param name The name of the value column.
#' @return A data frame `from`, `to`, `<name>`.
ref_long <- function(cells, name) {
  cells <- unclass(as.matrix(cells))
  states_from <- rownames(cells)
  states_to <- colnames(cells)
  out <- data.frame(from = states_from[as.vector(row(cells))], to = states_to[as.vector(col(cells))],
                    value = as.vector(cells), stringsAsFactors = FALSE)
  names(out)[[3L]] <- name
  out
}

#' Several square tables side by side, one row per cell
#'
#' @param ... Named matrices of one shape.
#' @return A data frame `from`, `to` and one column per matrix.
ref_cells <- function(...) {
  tables <- list(...)
  Reduce(\(left, right) merge(left, right, by = c("from", "to"), sort = FALSE),
         Map(ref_long, tables, names(tables)))
}

#' The same cells keyed as as.data.frame() of a fit keys them
#'
#' as.data.frame() of a fit is its edge list: `from` and `to` are node numbers there, and the state names
#' are in `from_label` and `to_label`.
#' @param cells A data frame with `from` and `to`.
#' @return The data frame with `from_label` and `to_label` in their place.
ref_labelled <- function(cells) {
  names(cells)[match(c("from", "to"), names(cells))] <- c("from_label", "to_label")
  cells
}

#' The 2 x 2 collapse of every cell of a table
#'
#' For the cell (given, target): a = given then target, b = given then another, c = another then target,
#' d = another then another.
#' @param observed A square matrix of counts.
#' @return A list of four matrices `a`, `b`, `c`, `d`.
ref_collapse <- function(observed) {
  observed <- unclass(as.matrix(observed))
  total <- sum(observed)
  row_total <- rowSums(observed)[as.vector(row(observed))]
  column_total <- colSums(observed)[as.vector(col(observed))]
  shaped <- \(values) matrix(values, nrow(observed), dimnames = dimnames(observed))
  list(a = observed, b = shaped(row_total - observed), c = shaped(column_total - observed),
       d = shaped(total - row_total - column_total + observed))
}

#' Lag-1 counts of sequences after collapsing immediate repeats
#'
#' A state then cannot follow itself, so the diagonal is impossible by construction.
#' @param wide A data frame, one sequence per row.
#' @return A matrix of counts with a zero diagonal.
ref_counts_without_repeats <- function(wide) {
  collapsed <- lapply(ref_sequences(wide), \(events) rle(events)$values)
  unclass(ref_pairs(collapsed, lag = 1L))
}

#' The 3 x 3 table of the engine tests (no zero cell)
#'
#' @return A matrix of counts.
ref_three <- function() {
  matrix(c(12, 8, 5, 4, 15, 6, 3, 7, 20), 3L, 3L, byrow = TRUE, dimnames = list(c("x", "y", "z"), c("x", "y", "z")))
}

#' Two coupled streams: x drives y, y has momentum, y does not drive x
#'
#' @param n Length of the streams.
#' @param k Number of codes.
#' @param seed The seed.
#' @return A data frame `x`, `y`.
ref_streams <- function(n = 1500L, k = 4L, seed = 99L) {
  set.seed(seed)
  codes <- letters[seq_len(k)]
  x <- sample(codes, n, replace = TRUE)
  y <- Reduce(\(previous, t) {
    weight <- 0.5 + 3 * (codes == x[t - 1L]) + 1.2 * (codes == previous)
    sample(codes, 1L, prob = weight / sum(weight))
  }, seq.int(2L, n), accumulate = TRUE, init = codes[1L])
  data.frame(x = x, y = y, stringsAsFactors = FALSE)
}

#' The aligned triples of transfer entropy from a source to a target
#'
#' @param source,target Character vectors of one length.
#' @param lag,history The lag and the length of the target's history.
#' @return A data frame: `future` (target at t + lag), `past` (the target's last `history` events up to t),
#'   `driver` (source at t).
ref_triples <- function(source, target, lag = 1L, history = 1L) {
  n <- length(target)
  at <- seq.int(history, n - lag)
  lags <- lapply(seq_len(history) - 1L, \(back) target[at - back])
  data.frame(future = target[at + lag], past = do.call(paste, c(lags, sep = "|")), driver = source[at],
             stringsAsFactors = FALSE)
}

#' Shannon entropy of the joint distribution of some variables, in bits
#'
#' @param ... Vectors of one length.
#' @return A number.
ref_entropy <- function(...) {
  share <- prop.table(table(...))
  share <- share[share > 0]
  -sum(share * log2(share))
}

#' Transfer entropy by the triple sum over the joint distribution, in bits
#'
#' TE = sum p(f, h, s) log2[ p(f, h, s) p(h) / (p(h, s) p(f, h)) ].
#' @param triples The table of ref_triples().
#' @return A number.
ref_transfer_entropy <- function(triples) {
  joint <- prop.table(table(f = triples$future, h = triples$past, s = triples$driver))
  cells <- subset(as.data.frame(joint, responseName = "p", stringsAsFactors = FALSE), p > 0)
  p_h <- apply(joint, "h", sum)[cells$h]
  p_hs <- apply(joint, c("h", "s"), sum)[cbind(cells$h, cells$s)]
  p_fh <- apply(joint, c("f", "h"), sum)[cbind(cells$f, cells$h)]
  sum(cells$p * log2(cells$p * p_h / (p_hs * p_fh)))
}

#' The sequences of a long event log as a wide table, with each sequence's group
#'
#' One row per actor (in the order of split()), padded with NA to the longest sequence; for tna::group_tna()
#' and tna::group_ftna(), which take a wide table and a group per row.
#' @param log A data frame.
#' @param actor,action,time,group Column names.
#' @return A list: `wide` (a data frame) and `group` (a character vector, one per row of `wide`).
ref_group_wide <- function(log, actor, action, time, group) {
  sequences <- ref_log_sequences(log, actor = actor, action = action, time = time)
  group_of <- tapply(log[[group]], log[[actor]], \(level) level[[1L]])
  longest <- max(lengths(sequences))
  wide <- as.data.frame(do.call(rbind, lapply(sequences, \(events) c(events, rep(NA, longest - length(events))))),
                        stringsAsFactors = FALSE)
  list(wide = wide, group = unname(group_of[names(sequences)]))
}

#' A weight matrix or initial vector of each group of a tna group model, one row per cell or state
#'
#' @param model A `group_tna` object (tna::group_tna() or tna::group_ftna()).
#' @param part `"weights"` or `"inits"`.
#' @param name The name of the value column.
#' @return A data frame `group`, `from`, `to`, `<name>` (weights) or `group`, `state`, `<name>` (inits).
ref_tna_groups <- function(model, part, name) {
  do.call(rbind, lapply(names(model), \(level) {
    values <- model[[level]][[part]]
    if (identical(part, "weights")) return(cbind(group = level, ref_long(values, name)))
    out <- data.frame(group = level, state = model[[level]]$labels, value = unname(values), stringsAsFactors = FALSE)
    names(out)[[3L]] <- name
    out
  }))
}

#' One comparison per cell of a Monte Carlo table, each at its own tolerance
#'
#' expect_reference() takes one tolerance for a table. A Monte Carlo summary needs one per number: four standard
#' errors of that number, from the reference's own draws. Each row of `expected` becomes one comparison, whose
#' quantity names its keys; a key on one side only is compared with NA and fails as missing.
#' @param actual,expected Data frames with the `key` columns and the column `value`.
#' @param tolerance One tolerance per row of `expected`, fixed before the comparison (NA where `expected` is NA).
#' @param key,value Column names.
#' @param quantity What the numbers are; the keys of each row are appended.
#' @param ... Passed to expect_reference(): reference, covers, note, authority.
#' @return NULL, invisibly.
expect_each_cell <- function(actual, expected, tolerance, key, value, quantity, ...) {
  undefined <- is.na(expected[[value]])
  stopifnot("one tolerance per row of `expected`" = length(tolerance) == nrow(expected),
            "tolerances must be finite and non-negative where the reference has a value" =
              all(undefined | (is.finite(tolerance) & tolerance >= 0)))
  # where the reference has no value (an untested cell) nothing is judged: NA on both sides is not comparable, a
  # value on one side only is missing
  tolerance[undefined] <- 0
  expected <- data.frame(expected[c(key, value)], tolerance_ = tolerance)
  names(expected)[match(value, names(expected))] <- "reference_"
  actual <- actual[c(key, value)]
  names(actual)[match(value, names(actual))] <- "package_"
  matched <- merge(expected, actual, by = key, all = TRUE, sort = TRUE)
  labels <- do.call(paste, c(lapply(matched[key], as.character), sep = ", "))
  invisible(Map(\(ours, theirs, allowed, label) {
    expect_reference(ours, theirs, quantity = sprintf("%s: %s", quantity, label),
                     tolerance = if (is.na(allowed)) 0 else allowed, ...)
  }, matched$package_, matched$reference_, matched$tolerance_, labels))
  invisible(NULL)
}

#' Standard error of the difference of two Monte Carlo shares, from the reference's share
#'
#' Each share is a mean of Bernoulli draws: the package's over `n_package` draws, the reference's over
#' `n_reference`. The variance p (1 - p) is taken at the reference's share, held at least 1 / n_reference away
#' from 0 and 1: a share of exactly 0 or 1 in the reference's draws does not prove a variance of 0.
#' @param share The reference's shares.
#' @param n_package,n_reference The numbers of draws.
#' @return One standard error per share.
ref_share_se <- function(share, n_package, n_reference) {
  bounded <- pmin(pmax(share, 1 / n_reference), 1 - 1 / n_reference)
  sqrt(bounded * (1 - bounded) * (1 / n_package + 1 / n_reference))
}
