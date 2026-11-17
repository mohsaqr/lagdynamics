# Plot a lag sequential analysis fit

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws a lag
sequential analysis fit as a residual heatmap, a lag transition network,
a chord diagram or a polar sunburst, chosen with `type`; each display
suits a different question. A residual heatmap lays out the full table
of antecedent by consequent states, coloured by the adjusted residual,
and shows at a glance which transitions are over-represented and which
are avoided. A lag transition network draws the states as nodes and the
transitions as directed edges, which shows the paths the process takes.
A chord diagram shows how the transitions out of each state are shared
among the targets, and a polar sunburst shows the same per-source
breakdown on a ring.

The heatmap (default) shows the adjusted residuals, with significant
cells in bold; `which` switches it to transition probabilities, observed
counts or expected counts. The other types are drawn by
[`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md),
[`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md)
and
[`plot_polar()`](https://pak.dynasite.org/lagdynamics/reference/plot_polar.md),
whose arguments can be passed through `...`. A grouped fit is drawn once
per group, or in one figure with `combined = TRUE`.

For the uncertainty of the residuals use
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)
on a bootstrap; for differences between groups,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the result
of
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
(see
[`plot.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md)).

## Usage

``` r
# S3 method for class 'lsa'
plot(x, type = c("heatmap", "network", "chord", "sunburst"), ...)

# S3 method for class 'lsa_group'
plot(
  x,
  type = c("heatmap", "network", "chord", "sunburst"),
  combined = FALSE,
  ...
)
```

## Arguments

- x:

  An `lsa` fit from
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or a
  grouped fit from `lsa(..., group = )`.

- type:

  `"heatmap"` (default), `"network"`, `"chord"` or `"sunburst"`.

- ...:

  Passed to the chosen display. For `"heatmap"`: `which`, one of
  `"residuals"` (default), `"prob"`, `"count"` or `"expected"`. For the
  other types, the arguments of
  [`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md),
  [`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md)
  or
  [`plot_polar()`](https://pak.dynasite.org/lagdynamics/reference/plot_polar.md).

- combined:

  For a grouped fit: `FALSE` (default) draws one figure per group;
  `TRUE` tiles the groups in one figure.

## Value

For `"heatmap"` and `"sunburst"`, a `ggplot` object, drawn when printed.
For `"network"` and `"chord"`, the value of
[`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md)
or
[`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md),
invisibly, with the plot drawn on the current device. For a grouped fit,
a list with one such value per group, invisibly, after drawing.

## Conditions

Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`
when the fit holds no finite values of the quantity requested for the
heatmap.

## See also

[`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md),
[`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md),
[`plot_polar()`](https://pak.dynasite.org/lagdynamics/reference/plot_polar.md),
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md),
[`plot.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md),
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)

## Examples

``` r
fit <- lsa(engagement)
plot(fit)

plot(fit, which = "prob")

plot(fit, type = "network")
#> Registered S3 method overwritten by 'cograph':
#>   method             from
#>   plot.tna_bootstrap tna 

plot(fit, type = "chord")


gfit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
            time = "Time", group = "Achiever")
plot(gfit, combined = TRUE)

```
