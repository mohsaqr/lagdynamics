library(lagdynamics)
observed <- as.matrix(read.delim("../data/three-codes.tsv", row.names = 1))
two_cell <- lsa(observed, engine = "two_cell") |> transitions() |>
  transform(transition = paste(from, "->", to), two_cell_z = adj_res) |>
  subset(select = c(transition, odds_ratio, log_or, log_or_se, two_cell_z))
bidirectional <- lsa(observed, engine = "bidirectional") |> transitions() |>
  transform(transition = paste(from, "->", to), bidirectional_z = adj_res) |>
  subset(select = c(transition, bidirectional_z))
parallel <- lsa(observed, engine = "parallel_dominance") |> transitions() |>
  transform(transition = paste(from, "->", to), parallel_dominance_z = adj_res) |>
  subset(select = c(transition, parallel_dominance_z))
nonparallel <- lsa(observed, engine = "nonparallel_dominance") |> transitions() |>
  transform(transition = paste(from, "->", to), nonparallel_dominance_z = adj_res) |>
  subset(select = c(transition, nonparallel_dominance_z))
result <- Reduce(merge, list(two_cell, bidirectional, parallel, nonparallel))
result
