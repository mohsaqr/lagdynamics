library(lagdynamics)
sequences <- read.csv("../data/engagement.csv", na.strings = "")
fits <- lsa_lags(sequences, lags = 1:3)
fits
result <- transitions(fits) |>
  transform(transition = sprintf("lag %d: %s -> %s", lag, from, to)) |>
  subset(select = c(transition, count, prob, adj_res))
result
