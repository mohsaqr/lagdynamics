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
