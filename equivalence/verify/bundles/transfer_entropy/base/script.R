streams <- read.csv("../data/coupled-streams.csv")
# TE = sum p(f, h, s) log2[ p(f, h, s) p(h) / (p(h, s) p(f, h)) ], with f the target's
# next event, h its present event and s the source's present event.
bits <- function(source, target) {
  n <- length(target)
  joint <- prop.table(table(f = target[-1], h = target[-n], s = source[-n]))
  cells <- subset(as.data.frame(joint, responseName = "p"), p > 0)
  p_h <- as.vector(apply(joint, "h", sum)[as.character(cells$h)])
  p_hs <- apply(joint, c("h", "s"), sum)[cbind(as.character(cells$h), as.character(cells$s))]
  p_fh <- apply(joint, c("f", "h"), sum)[cbind(as.character(cells$f), as.character(cells$h))]
  te <- sum(cells$p * log2(cells$p * p_h / (p_hs * p_fh)))
  # The target's uncertainty left after its own past: H(f | h).
  pairs <- apply(joint, c("f", "h"), sum)
  leftover <- -sum(ifelse(pairs > 0, pairs * log2(pairs / rep(colSums(pairs), each = nrow(pairs))), 0))
  c(te = te, te_normalised = te / leftover)
}
values <- rbind(bits(streams$x, streams$y), bits(streams$y, streams$x))
result <- data.frame(direction = c("x -> y", "y -> x"), values)
result
