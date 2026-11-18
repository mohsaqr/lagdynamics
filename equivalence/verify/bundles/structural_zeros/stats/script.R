observed <- as.matrix(read.delim("../data/regulation-no-repeats.tsv", row.names = 1))
cells <- as.data.frame(as.table(observed), responseName = "count")
names(cells) <- c("from", "to", "count")
possible <- subset(cells, from != to)
model <- glm(count ~ from + to, family = poisson, data = possible,
             control = glm.control(epsilon = 1e-14, maxit = 100))
model
result <- data.frame(transition = paste(possible$from, "->", possible$to),
                     expected = fitted(model), adj_res = rstandard(model, type = "pearson"))
result
