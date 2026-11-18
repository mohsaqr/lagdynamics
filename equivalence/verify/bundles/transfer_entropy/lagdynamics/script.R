library(lagdynamics)
streams <- read.csv("../data/coupled-streams.csv")
result <- with(streams, transfer_entropy(x = x, y = y, lag = 1, history = 1, test = "none")) |>
  transform(direction = paste(from, "->", to)) |>
  subset(select = c(direction, te, te_normalised))
result
