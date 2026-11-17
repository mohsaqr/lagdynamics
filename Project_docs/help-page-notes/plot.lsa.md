# Implementation notes removed from `?plot.lsa` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/plot-lsa.R); it is kept here for maintainers.

- **Heatmap worker** `.plot_lsa_heatmap()` (ggplot2, styled after
  pathtree::plot_pathways()): rows (sources) ordered by outgoing volume,
  most active at the top; columns in state order; cells with zero count
  are grey (`na.value = "grey85"`) and unlabelled.
- **Bold cell**: per row, the largest value (largest absolute residual
  for `which = "residuals"`, where it must also be significant); other
  cells labelled in grey. Labels have 0 decimals for counts, 2 otherwise.
- **Colours**: residuals use `scale_fill_gradient2()` with
  `.lsa_div_low = "#2166AC"` (negative, blue), `.lsa_div_mid =
  "#F7F7F7"`, `.lsa_div_high = "#B2182B"` (positive, red), symmetric
  limits; probabilities use `#DEEBF7 -> #08519C` on `[0, 1]`. These are
  not the Okabe-Ito / `#D33F6A`-`#4A6FE3` palette of the house rules,
  and the residual direction is opposite to the network and comparison
  plots (blue = over-represented there).
- **Grouped fits**: `combined = FALSE` prints one ggplot per group (or
  draws one base plot per group with `title()`); `combined = TRUE` tiles
  ggplots with grid viewports or base plots with `par(mfrow)`, restored
  on exit. A missing ggplot2 raises `lagdynamics_missing_package`
  (unreachable while ggplot2 is in Imports).
