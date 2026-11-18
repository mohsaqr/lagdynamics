library(lagdynamics)
sequences <- read.csv("../data/engagement.csv", na.strings = "")
fit <- lsa(sequences, lag = 1)
fit
result <- transitions(fit) |>
  transform(transition = paste(from, "->", to)) |>
  subset(select = c(transition, count, prob))
result
