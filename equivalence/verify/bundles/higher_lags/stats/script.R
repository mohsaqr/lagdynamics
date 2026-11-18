sequences <- read.csv("../data/engagement.csv", na.strings = "")
states <- sort(unique(na.omit(unlist(sequences))))
count_pairs <- function(lag) {
  pairs <- do.call(rbind, lapply(seq_len(nrow(sequences)), function(row) {
    events <- as.character(na.omit(unlist(sequences[row, ])))
    data.frame(from = head(events, -lag), to = tail(events, -lag))
  }))
  table(from = factor(pairs$from, states), to = factor(pairs$to, states))
}
result <- do.call(rbind, lapply(1:3, function(lag) {
  counts <- count_pairs(lag)
  cells <- as.data.frame(counts, responseName = "count")
  data.frame(transition = sprintf("lag %d: %s -> %s", lag, cells$from, cells$to),
             count = cells$count, prob = as.vector(prop.table(counts, margin = 1)),
             adj_res = as.vector(chisq.test(counts, correct = FALSE)$stdres))
}))
result
