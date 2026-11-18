library(lagdynamics)
sequences <- read.csv("../data/engagement.csv", na.strings = "")
fit <- lsa(sequences, lag = 1, engine = "classical", alternative = "two.sided")
result <- transitions(fit) |>
  transform(transition = paste(from, "->", to)) |>
  subset(select = c(transition, expected, adj_res, p))
result
