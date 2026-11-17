# Implementation notes removed from `?lsa_to_tna` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/tna-bridge.R); it is kept here for maintainers.

- **Name**: deliberately `lsa_to_tna()`, not `as_tna()`, to avoid an
  export-name collision with sibling packages.
- **Object built**: a list `weights`, `inits`, `labels`, `data` with
  attributes `type` (`"frequency"` for counts, else `"relative"`) and
  `scaling = character(0)`, class `"tna"`. `data` is a `tna_seq_data`
  matrix rebuilt from the events (`.lsa_seqdata_matrix()`), with
  `alphabet`, `labels` and `grDevices::rainbow()` colours; `NULL` for
  transition-matrix fits.
- **Weights**: non-finite cells become 0 and structural-zero cells are
  forced to 0 (`.lsa_weight_matrix()`); for `"adj_res"` negative
  residuals are set to 0.
- **adj_res / lift**: `inits` and `data` are `NULL`. Printing such an
  object with `tna`'s print method currently fails (see the report of
  2026-11-17: `print(lsa_to_tna(lsa(engagement), weights = "adj_res"))`
  errors with "attempt to set an attribute on NULL").
