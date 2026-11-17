# Implementation notes removed from `?lagdynamics-package` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/lagdynamics-package.R); it is kept here for maintainers.

- **Design principles section (removed)**: numerical methods are
  implemented from primary literature and cross-validated against base-R
  primitives and hand-formula identities (`stats::loglin()`,
  `stats::chisq.test()`, `stats::pchisq()`); see
  `tests/testthat/test-references.R` and `inst/REFERENCES.md`.
- **Roadmap section (removed)**: stationarity tests (`stationarity_lsa()`)
  are planned but not exported. The page no longer advertises unexported
  functions.
- **Reference dropped**: "Wampold, B. E. (1982). Sequential analysis:
  Stationarity, statistical procedures, and applications. Psychological
  Bulletin, 92(2), 380-387." It was cited only for the unexported
  stationarity roadmap, is not cited in `inst/REFERENCES.md`, and its
  bibliographic details could not be confirmed. Verify before restoring.
- **Reference title completed**: Bakeman & Quera (1995) is cited with its
  full title, *Analyzing interaction: Sequential analysis with SDIS and
  GSEQ* (the package page had a truncated title).
