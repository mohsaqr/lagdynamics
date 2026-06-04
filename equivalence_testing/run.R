#!/usr/bin/env Rscript
# Run the equivalence battery against the current source tree.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
testthat::test_dir("equivalence_testing", stop_on_failure = TRUE)
