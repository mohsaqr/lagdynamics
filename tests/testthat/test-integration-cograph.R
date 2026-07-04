test_that("lsa fits satisfy the cograph network protocol", {
  skip_if_not_installed("cograph")

  fit <- lsa(engagement)

  expect_s3_class(fit, "cograph_network")
  expect_named(fit$nodes, c("id", "label", "name", "outgoing", "incoming"))
  expect_true(all(c("from", "to", "weight") %in% names(fit$edges)))
  expect_type(fit$edges$from, "integer")
  expect_type(fit$edges$to, "integer")
  expect_equal(sort(unique(c(fit$edges$from, fit$edges$to))), fit$nodes$id)
  expect_identical(fit$weights, fit$obs)

  expect_identical(cograph::get_nodes(fit), fit$nodes)
  expect_identical(cograph::get_edges(fit), fit$edges)

  flat <- cograph::to_df(fit)
  expect_named(flat, c("from", "to", "weight"))
  expect_equal(nrow(flat), nrow(fit$edges))
  expect_equal(sum(flat$weight), sum(fit$obs))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(cograph::splot(fit), "cograph_network")
})
