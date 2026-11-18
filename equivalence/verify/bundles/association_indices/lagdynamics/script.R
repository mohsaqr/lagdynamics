library(lagdynamics)
sequences <- read.csv("../data/engagement.csv", na.strings = "")
fit <- lsa(sequences, lag = 1)
result <- transitions(fit) |>
  transform(transition = paste(from, "->", to)) |>
  subset(select = c(transition, yules_q, lift, prob_col))
result
