library(lagdynamics)
sequences <- read.csv("../data/engagement.csv", na.strings = "")
fit <- lsa(sequences, lag = 1)
result <- tests(fit)
result
