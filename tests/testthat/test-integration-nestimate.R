test_that("bundled Nestimate-derived data stay compatible with Nestimate", {
  skip_if_not_installed("Nestimate")

  expect_identical(as.data.frame(Nestimate::trajectories,
                                 stringsAsFactors = FALSE), engagement)

  ns_data <- new.env(parent = emptyenv())
  utils::data("ai_long", package = "Nestimate", envir = ns_data)
  if (!exists("ai_long", envir = ns_data, inherits = FALSE)) {
    skip("Installed Nestimate does not ship ai_long.")
  }
  expect_named(ns_data$ai_long, names(ai_long))

  fit <- lsa(ns_data$ai_long, actor = "project", session = "session_id",
             action = "code", order = "order_in_session")
  expect_s3_class(fit, "lsa")
  expect_equal(fit$data$n_sequences, length(unique(ns_data$ai_long$session_id)))
  expect_equal(fit$data$n_events, nrow(ns_data$ai_long))
})

test_that("lsa ingests Nestimate prepared sequence containers", {
  skip_if_not_installed("Nestimate")

  prepared <- structure(
    list(sequence_data = as.data.frame(Nestimate::trajectories,
                                       stringsAsFactors = FALSE)),
    class = "nestimate_data"
  )

  from_prepared <- lsa(prepared)
  from_matrix <- lsa(Nestimate::trajectories)

  expect_s3_class(from_prepared, "lsa")
  expect_equal(from_prepared$obs, from_matrix$obs)
  expect_equal(from_prepared$data$n_sequences, from_matrix$data$n_sequences)
  expect_identical(from_prepared$data$labels, from_matrix$data$labels)
})

test_that("lsa ingests Nestimate network objects", {
  skip_if_not_installed("Nestimate")

  net <- Nestimate::build_network(ai_long, method = "tna",
                                  actor = "project",
                                  session = "session_id",
                                  action = "code",
                                  order = "order_in_session")

  from_net <- lsa(net)
  from_log <- lsa(ai_long, actor = "project", session = "session_id",
                  action = "code", order = "order_in_session")

  expect_s3_class(from_net, "lsa")
  expect_equal(from_net$obs, from_log$obs)
  expect_equal(from_net$data$n_sequences, from_log$data$n_sequences)
  expect_identical(from_net$data$labels, from_log$data$labels)
})
