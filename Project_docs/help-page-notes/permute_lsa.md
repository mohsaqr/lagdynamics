# Implementation notes removed from `?permute_lsa` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/permute_lsa.R); it is kept here for
maintainers.

- **p-value**: `p_perm = (1 + #{|z*| >= |z_obs|}) / (1 + n_finite)`,
  with `n_finite` the number of finite permuted residuals for that cell
  (not a flat `R`; the old help text gave `(1 + R)`, true only when
  every replicate is finite). Exceedances are counted with
  `na.rm = TRUE`, so `NA` replicates (structural zeros, zero-margin
  cells in a permuted table) are neither exceedances nor
  non-exceedances.
- **Untestable cells**: a cell with a non-finite observed residual, or
  with no finite replicate, gets `p_perm = NA`. Without this guard an
  `NA` observed residual would give 0 exceedances and a spurious
  `p = 1 / (R + 1)`.
- **Within-sequence shuffle** (`.shuffle_within_sequences()`): for each
  sequence's positions `sp`, `events[sp] <- events[sp[sample.int(length(sp))]]`.
  Sources are drawn from `sp`, not `1..length(sp)` (an earlier bug copied
  one sequence's content into another); `sample.int()` avoids the
  `sample(x)` length-1 quirk. Sequences of length 1 are left alone.
- **Refit**: `.refit_from_events()` with the fit's recipe and the
  original `seq_id`, so sequence boundaries are unchanged even with
  `within_sequence = FALSE` (only the events move).
- **shuffles replay**: must be a list of length `>= R` (extras dropped);
  every element must be a permutation of `seq_len(T)` (checked by
  `.is_permutation()`: length `T`, integer-coercible, no `NA`, a
  bijection). The first offending index is named in the error.
- **Single sequence**: with `within_sequence = TRUE` and one sequence, a
  warning says within-sequence and global shuffling coincide.
- **significant**: `is.finite(p_perm) & p_perm < fit$params$alpha`.
- **perm_adj_res**: `R x K^2`, cells in `as.vector(adj_res)` order.
