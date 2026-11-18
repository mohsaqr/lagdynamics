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
row_share <- rowSums(observed) / total
column_share <- colSums(observed) / total
expected <- total * outer(row_share, column_share)
# Haberman (1973): z = (O - E) / sqrt(E (1 - row share) (1 - column share))
z <- (observed - expected) / sqrt(expected * outer(1 - row_share, 1 - column_share))
result <- data.frame(transition = transition, expected = as.vector(expected),
                     adj_res = as.vector(z), p = as.vector(2 * pnorm(abs(z), lower.tail = FALSE)))
result
