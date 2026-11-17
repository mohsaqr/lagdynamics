Sys.setenv(OMP_THREAD_LIMIT = 2)
library(testthat)
library(lagdynamics)

test_check("lagdynamics")
