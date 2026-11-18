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
pearson <- chisq.test(counts, correct = FALSE)
pearson
g2 <- 2 * sum(counts * log(counts / pearson$expected))
df <- unname(pearson$parameter)
result <- data.frame(test = c("lrx2", "x2"), statistic = c(g2, unname(pearson$statistic)),
                     df = c(df, df),
                     p = c(pchisq(g2, df, lower.tail = FALSE), pearson$p.value))
result
