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
test <- chisq.test(counts, correct = FALSE)
test
result <- data.frame(transition = transition, expected = as.vector(test$expected),
                     adj_res = as.vector(test$stdres),
                     p = as.vector(2 * pnorm(-abs(test$stdres))))
result
