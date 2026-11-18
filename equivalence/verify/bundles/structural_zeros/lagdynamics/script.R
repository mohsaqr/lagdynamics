library(lagdynamics)
observed <- as.matrix(read.delim("../data/regulation-no-repeats.tsv", row.names = 1))
fit <- lsa(observed, lag = 1, loops = FALSE)
fit
result <- transitions(fit) |>
  subset(from != to) |>
  transform(transition = paste(from, "->", to)) |>
  subset(select = c(transition, expected, adj_res))
result
