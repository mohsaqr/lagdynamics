# Build lagdynamics's `ai_long` data set from the Nestimate package.
# Run with: Rscript data-raw/ai_long.R

stopifnot(requireNamespace("Nestimate", quietly = TRUE))

data("ai_long", package = "Nestimate")

stopifnot(
  is.data.frame(ai_long),
  nrow(ai_long) == 8551L,
  all(c("message_id", "project", "session_id", "timestamp", "session_date",
        "code", "cluster", "code_order", "order_in_session") %in%
        names(ai_long)),
  length(unique(ai_long$session_id)) == 428L,
  length(unique(ai_long$code)) == 8L
)

save(ai_long, file = "data/ai_long.rda", compress = "bzip2")

