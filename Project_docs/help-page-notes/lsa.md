# Implementation notes removed from `?lsa` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/lsa.R); it is kept here for maintainers.

- **lag = 0** pairs each event with itself: a degenerate diagonal for
  single-stream event data. Genuine co-occurrence needs concurrent codes,
  which are not supported.
- **Structural zeros**: when any cell is forbidden the classical engine
  fits the expected table with `lsa_ipf()` and computes Christensen's
  design-matrix residuals (`inst/REFERENCES.md` sections 2.2 and 4.2). A
  non-converged IPF and a singular design matrix raise unclassed
  warnings (R/engine-classical.R).
- **`loops = FALSE`** zeros the diagonal of a supplied `structural_zeros`
  matrix as well (`.resolve_structural_zeros()`). The matrix gets the
  state labels as dimnames.
- **Grouped fits** (`.lsa_grouped()`): sequences are partitioned by group
  and each group is refitted with `labels = global_labels`, so every
  group has a full `K x K` table even if it never visits some state.
  Empty factor levels of `group` are dropped. The result carries
  `levels`, `group_sizes`, `labels`, `engine` and `call` attributes.
- **Long-format group column**: `.prepare_long()` derives one group label
  per recovered sequence and returns it as an attribute that drives the
  grouped fit.
- **actor x session overflow**: with both `actor` and `session`, each
  sequence is keyed by `interaction()`, an integer index over the full
  cross of distinct actors and sessions. It overflows once that product
  exceeds `.Machine$integer.max` (about 2.1e9); overflowing keys become
  `NA` and those sequences are dropped with only R's warning "NAs produced
  by integer overflow". The help page now only advises passing a globally
  unique session id as `actor`.
- **Slots of the `lsa` object** (formerly listed in `@return`):
  `edges` (backs `transitions()`, plus cograph protocol columns `from`,
  `to` integer ids, `from_label`, `to_label`, `edge`, `weight`);
  `nodes` (`id, label, name, outgoing, incoming`); `K x K` matrices
  `obs, exp, prob, prob_col, adj_res, p, yules_q, kappa, kappa_z,
  kappa_p`; `lrx2`, `x2` (lists `statistic, df, p`, `NULL` for engines
  without an expected table); `inits` (named vector, `NULL` for
  transition-matrix input); `weights` (`K x K`, equal to `obs`);
  `directed` (`FALSE` for the bidirectional engine); `method` (engine
  name, read by the cograph protocol, also in `params$engine`); `data`
  (the `lsa_data`); `params` (recipe snapshot used by the resampling
  verbs); `meta` (source, IPF info, version, call).
- **`sign` column**: `"over"`, `"under"` or `"expected"` from `count`
  versus `expected` (lift > 1, < 1, = 1). For the bidirectional engine the
  residual is computed on the symmetrised table but `count`/`expected`
  are directional, so `sign` can disagree with the sign of `adj_res`
  (engagement: Active -> Active has `adj_res = 30.6` and `sign = "under"`).
- **Engine wrappers** were previously documented under `@details`; they
  are `lsa()` with `engine` fixed.
