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
model <- loglin(counts, margin = list(1, 2), fit = TRUE, eps = 1e-10, print = FALSE)
model$lrt
model$pearson
model$df
result <- data.frame(test = c("lrx2", "x2"), statistic = c(model$lrt, model$pearson),
                     df = c(model$df, model$df),
                     p = pchisq(c(model$lrt, model$pearson), model$df, lower.tail = FALSE))
result
