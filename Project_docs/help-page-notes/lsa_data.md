# Implementation notes removed from `?lsa_data` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/lsa-data.R); it is kept here for maintainers.

- **Matrix detection** (`.is_transition_matrix()`): a numeric matrix with
  `nrow == ncol >= 2` is read as transition counts. A data frame is never
  read as counts, even if square.
- **Sequence objects** are detected before the shape heuristics
  (`.is_seq_object()`, `.sequences_from_object()` in R/ingest-objects.R)
  so that a fitted object, which is a list, is not misread as a list of
  sequences. A `tna` built from a bare matrix (no retained sequences)
  errors, because counts cannot be recovered from probability weights.
- **Empty sequences** (all `NA`) are dropped from lists and wide tables.
- **Integer codes** must be finite whole numbers >= 1; with `labels`
  they index into it and the maximum code must not exceed its length.
- **Elements of the object** (formerly in `@return`): `events` (integer
  codes, `NULL` for a matrix), `seq_id`, `labels`, `n_states`,
  `n_sequences`, `n_events`, `transitions_per_seq` (lag-1 transitions
  per sequence), `source` (`"events"` / `"transitions"`), `obs_input`
  (the matrix when `source = "transitions"`).
