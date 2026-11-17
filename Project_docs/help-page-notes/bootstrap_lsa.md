# Implementation notes removed from `?bootstrap_lsa` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/bootstrap_lsa.R); it is kept here for
maintainers.

- **Recipe refit**: every resample is refitted by `.refit_from_events()`
  with the snapshot in `fit$params` (engine, lag, alternative,
  structural zeros, engine params). It builds a bare `lsa_data` directly
  rather than going through `lsa_data()`, for speed, and passes
  `n_events = length(events)` and the event column totals (needed by the
  kappa statistics).
- **Sequence level**: each resample draws `S` sequence indices with
  replacement from `seq_len(S)` (`sample.int(S, R * S, replace = TRUE)`
  filled into an `R x S` matrix) and concatenates the drawn sequences,
  each kept as its own sequence.
- **Event level** (`.stationary_indices()`): block lengths are
  `1 + rgeom(1, 1 / mean_block)`, start positions uniform on `1..T`,
  blocks wrap around the end of the stream, and blocks are concatenated
  until the resample has exactly `T` events. All events are one sequence
  in the refit. `block_length` is coerced with `as.integer()`.
- **Single-sequence fallback**: with one sequence and `level = "auto"`
  the warning is "Only one sequence available; using event-level
  stationary block bootstrap ..."; with `level = "sequence"` it is
  "... falling back to event-level stationary block bootstrap." The
  warning is raised after the indices are drawn, before refitting.
- **Replay hook** (`.validate_boot_indices()`): `indices` must be a
  numeric matrix with at least `R` rows (extra rows are dropped), finite
  whole numbers only, exactly `S` columns (sequence level) or `T` columns
  (event level), values in `1..S` or `1..T`. The event-level matrix holds
  the fully expanded positions (the form `.stationary_indices()`
  produces), not block starts, so `indices_used` can be fed back in for
  bit-identical replay across sessions or languages.
- **Summaries** (`.summarize_bootstrap()`): means, SDs and percentile
  quantiles (`stats::quantile()`, default type 7) per cell, all with
  `na.rm = TRUE`. Replicates that are `NA` for a cell (structural zeros,
  rows whose margin is zero in the resample) are excluded from that
  cell's summary; a cell that is `NA` in every replicate stays `NA`.
- **Bootstrap p-value**: `adj_res_p_boot = 2 * min(P(z* <= 0), P(z* >= 0))`
  over the finite replicates. Only the residual gets a p-value; a
  p-value "around zero" is not meaningful for counts or probabilities.
- **Stable flag**: `adj_res_stable` is `TRUE` when both interval ends are
  finite, non-zero and of the same sign.
- **Per-resample matrices**: `boot_obs`, `boot_adj_res`, `boot_prob`,
  `boot_yulesq` are `R x K^2` with cells in `as.vector(obs)`
  (column-major) order.
- **Parallel** (`.run_parallel()`, shared with permute/stability/
  reliability/compare): `parallel::mclapply()` on Unix; on Windows a
  PSOCK cluster seeded with `parallel::clusterSetRNGStream()` from the
  caller's RNG state, so `set.seed(); bootstrap_lsa(parallel = TRUE)` is
  reproducible there as well. Per-replicate progress messages (`verbose`) are only
  emitted in the sequential path.
