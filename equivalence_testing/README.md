# Equivalence testing (does NOT ship with the package)

This folder holds the cross-validation / published-oracle test battery:
it checks lagseq's numbers against independently published worked
examples and base-R primitives. It is intentionally kept **out of the
package's own test suite** (`tests/testthat/`) and excluded from the
build (see `.Rbuildignore`), so the shipped package does not depend on
or reference external tools.

These tests use the package's exported oracle datasets
(`oconnor_couple`, `kg_logs`, `kg_lsa_oracle`, `qi2026_grandmother`,
`imdb_genres`) and base-R checks (`stats::chisq.test`, `stats::loglin`).

## Run

```r
# from the package root
pkgload::load_all(".", quiet = TRUE)
testthat::test_dir("equivalence_testing")
```
