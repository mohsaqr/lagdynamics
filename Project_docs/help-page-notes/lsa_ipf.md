# Implementation notes removed from `?lsa_ipf` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/ipf.R); it is kept here for maintainers.

- **Algorithm** (Wickens, 1989, pp. 107-112): the expected table starts
  as the 0/1 `structure` matrix and alternates row and column scaling
  (scale factor 0 for an all-zero row or column) until the largest
  absolute difference between fitted and observed margins is below
  `tol`. Validated against `stats::loglin()` (tests/testthat/test-ipf.R).
- **Non-convergence** is reported only in the `converged` field; no
  warning is raised by `lsa_ipf()` itself (the classical engine warns).
- **Validation**: a non-square or non-numeric `obs` and an invalid
  `structure` fail a plain `stopifnot()` (unclassed).
