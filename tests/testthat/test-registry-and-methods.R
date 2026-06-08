# Coverage for the public engine registry and the inference-object
# print / as.data.frame S3 methods. These are exported surfaces that
# were previously exercised only indirectly.

# --- engine registry ---------------------------------------------------

test_that("register / get / list / unregister round-trips", {
  nm <- "test_custom_engine_xyz"
  on.exit(if (nm %in% list_lsa_engines()$name) unregister_lsa_engine(nm),
          add = TRUE)

  fn <- function(transitions, ...) .engine_classical(transitions, ...)
  expect_invisible(register_lsa_engine(nm, fn, "A test engine"))

  reg <- list_lsa_engines()
  expect_s3_class(reg, "data.frame")
  expect_true(nm %in% reg$name)

  entry <- get_lsa_engine(nm)
  expect_identical(entry$name, nm)
  expect_true(is.function(entry$fn))

  # The custom engine actually drives lsa().
  fit <- lsa(c("a", "b", "a", "b", "a"), engine = nm)
  expect_s3_class(fit, "lsa")

  expect_invisible(unregister_lsa_engine(nm))
  expect_false(nm %in% list_lsa_engines()$name)
})

test_that("register_lsa_engine() rejects a fn without a transitions arg", {
  expect_error(
    register_lsa_engine("bad_engine", function(x) x, "no transitions arg"),
    "transitions"
  )
})

test_that("get_lsa_engine() errors on an unknown engine", {
  expect_error(get_lsa_engine("definitely_not_registered"),
               "not registered")
})

test_that("get_lsa_engine() errors when a required package is missing", {
  nm <- "needs_missing_pkg"
  on.exit(if (nm %in% list_lsa_engines()$name) unregister_lsa_engine(nm),
          add = TRUE)
  register_lsa_engine(nm, function(transitions, ...) transitions,
                      "needs a missing package",
                      requires = "a.package.that.does.not.exist")
  expect_error(get_lsa_engine(nm), "not installed")
})

test_that("unregister_lsa_engine() errors on an unknown engine", {
  expect_error(unregister_lsa_engine("never_registered"), "not registered")
})

test_that("built-in engines are all registered", {
  reg <- list_lsa_engines()
  expect_true(all(c("classical", "two_cell", "bidirectional",
                    "parallel_dominance", "nonparallel_dominance")
                  %in% reg$name))
})

# --- inference-object print + as.data.frame ----------------------------

test_that("inference objects print and coerce to their edge frames", {
  set.seed(101)
  fit <- lsa(list(c("a", "b", "a", "b", "a"), c("b", "a", "b", "a", "b")),
             engine = "classical")

  bs <- bootstrap_lsa(fit, R = 20)
  pm <- permute_lsa(fit, R = 20)
  st <- stability_lsa(fit, R = 20, proportion = 0.5)

  expect_output(print(bs), "<lsa_bootstrap>")
  expect_output(print(pm), "<lsa_permutation>")
  expect_output(print(st), "<lsa_stability>")

  expect_identical(as.data.frame(bs), bs$edges)
  expect_identical(as.data.frame(pm), pm$edges)
  expect_identical(as.data.frame(st), st$edges)
  expect_s3_class(as.data.frame(bs), "data.frame")
})
