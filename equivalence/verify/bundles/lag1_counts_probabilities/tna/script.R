library(tna)
sequences <- read.csv("../data/engagement.csv", na.strings = "")
frequencies <- ftna(sequences)$weights
probabilities <- tna(sequences)$weights
frequencies
probabilities
states <- rownames(frequencies)
result <- data.frame(transition = paste(states[row(frequencies)], "->", states[col(frequencies)]),
                     count = as.vector(frequencies), prob = as.vector(probabilities))
result
