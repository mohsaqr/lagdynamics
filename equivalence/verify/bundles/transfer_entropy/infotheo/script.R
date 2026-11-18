library(infotheo)
streams <- read.csv("../data/coupled-streams.csv")
bits <- function(source, target) {
  n <- length(target)
  future <- as.integer(factor(target[-1]))
  past <- as.integer(factor(target[-n]))
  driver <- as.integer(factor(source[-n]))
  te <- condinformation(future, driver, past, method = "emp") / log(2)
  c(te = te, te_normalised = te / (condentropy(future, past, method = "emp") / log(2)))
}
values <- rbind(bits(streams$x, streams$y), bits(streams$y, streams$x))
result <- data.frame(direction = c("x -> y", "y -> x"), values)
result
