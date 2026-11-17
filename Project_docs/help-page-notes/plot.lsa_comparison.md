# Implementation notes removed from `?plot.lsa_comparison` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/plot-comparison.R); it is kept here for maintainers.

- **Barrel details**: bar ends carry the observed count; the bar of the
  group with the higher value gets a border that darkens with the size of
  the difference; the centre chip is starred when significant.
- **Palette**: vcd / Nestimate mosaic convention, `.cmp_low =
  "#D33F6A"` (avoided), `.cmp_mid = "#F7F7F7"`, `.cmp_high = "#4A6FE3"`
  (over-represented). The log odds ratio fill is `.lsa_log_or()`.
- **rank = "effect"** considers only tested cells (finite `p_perm`).
- **Pairwise objects**: one facet per pair; the ranking is shared across
  facets so rows line up.
- The heatmap subtitle states the measure, `R`, the number of
  significant edges and `adjust`.
