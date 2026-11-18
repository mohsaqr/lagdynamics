sequences <- read.csv("../data/engagement.csv", na.strings = "")
states <- sort(unique(na.omit(unlist(sequences))))
count_pairs <- function(lag) {
  pairs <- do.call(rbind, lapply(seq_len(nrow(sequences)), function(row) {
    events <- as.character(na.omit(unlist(sequences[row, ])))
    data.frame(from = head(events, -lag), to = tail(events, -lag))
  }))
  table(from = factor(pairs$from, states), to = factor(pairs$to, states))
}
counts <- count_pairs(lag = 1)
counts
cells <- as.data.frame(counts, responseName = "count")
transition <- paste(cells$from, "->", cells$to)
observed <- unclass(counts)
total <- sum(observed)
row_total <- rowSums(observed)[row(observed)]
column_total <- colSums(observed)[col(observed)]
# The 2 x 2 collapse of a cell: n_a = given then target, n_b = given then other,
# n_c = other then target, n_d = other then other. Q = (ad - bc) / (ad + bc).
n_a <- as.vector(observed)
n_b <- row_total - n_a
n_c <- column_total - n_a
n_d <- total - row_total - column_total + n_a
result <- data.frame(transition = transition, yules_q = (n_a * n_d - n_b * n_c) / (n_a * n_d + n_b * n_c),
                     lift = n_a / (row_total * column_total / total), prob_col = n_a / column_total)
result
