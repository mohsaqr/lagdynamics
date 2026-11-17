# Implementation notes removed from `?plot_transitions` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/plot-transitions.R); it is kept here for maintainers.

- **Probability / count networks** use `cograph::splot(tna_styling =
  TRUE)`; the donut ring holds `inits`; `edge_cutoff = 0.05` is the
  default for probabilities only. `node_fill` is passed only when it
  differs from `"white"`.
- **Signed networks** (residuals, Yule's Q): positive edges `.cmp_high =
  "#4A6FE3"`, negative `.cmp_low = "#D33F6A"`; avoided edges dashed
  (`edge_style = 2`) with a halo drawn through cograph's CI band
  (`edge_ci` = |weight| / max |negative weight|, colour `.cmp_low`,
  alpha 0.18, scale 2). Lift edges are a single `#4A6FA5`. Node border
  `"steelblue"`, width 1.1. All overridable through `...`.
- **top**: fraction keeps `ceiling(top * n_edges)`, integer keeps
  `min(top, n_edges)`; ties broken by `order()`.
- **significant** uses `fit$p < fit$params$alpha`.
- A missing cograph raises `lagdynamics_missing_package` (unreachable
  while cograph is in Imports).
