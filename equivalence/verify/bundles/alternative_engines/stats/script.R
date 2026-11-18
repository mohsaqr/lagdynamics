observed <- as.matrix(read.delim("../data/three-codes.tsv", row.names = 1))
observed
total <- sum(observed)
row_total <- rowSums(observed)[row(observed)]
column_total <- colSums(observed)[col(observed)]
# Two-cell: the 2 x 2 collapse of each cell (no cell is zero here).
n_a <- as.vector(observed)
n_b <- row_total - n_a
n_c <- column_total - n_a
n_d <- total - row_total - column_total + n_a
log_or <- log((n_a * n_d) / (n_b * n_c))
log_or_se <- sqrt(1 / n_a + 1 / n_b + 1 / n_c + 1 / n_d)
# Bidirectional: Haberman residuals of the symmetrised table.
symmetric <- chisq.test(observed + t(observed), correct = FALSE)
# Dominance: the difference between the two directions of a pair.
expected <- outer(rowSums(observed), colSums(observed)) / total
difference <- observed - t(observed)
states <- rownames(observed)
result <- data.frame(
  transition = paste(states[row(observed)], "->", states[col(observed)]),
  odds_ratio = exp(log_or), log_or = log_or, log_or_se = log_or_se,
  two_cell_z = log_or / log_or_se,
  bidirectional_z = as.vector(symmetric$stdres),
  parallel_dominance_z = as.vector(difference / sqrt(expected + t(expected))),
  nonparallel_dominance_z = as.vector(difference / sqrt(observed + t(observed))))
result
