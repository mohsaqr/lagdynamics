# Implementation notes removed from `?lsa_transitions` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/lsa-transitions.R); it is kept here for maintainers.

- **Counting** (`.count_transitions()`): vectorised pair index
  `(from - 1) * K + to` tabulated with `tabulate()` within sequences;
  `inst/REFERENCES.md` section 1.
- **lag = 0** pairs each event with itself (a degenerate diagonal).
- **Validation of `lag`** is a plain `stopifnot()` (unclassed error).
- **Elements of the object** (formerly in `@return`): `obs` (`K x K`
  with dimnames), `row_totals`, `col_totals`, `n_transitions`, `lag`,
  `labels`, `edges` (the table returned by `as.data.frame()`).
