# Plot a lag transition network

`plot_transitions()` draws a lag transition network with cograph: the
states are nodes and the transitions between them directed edges. Which
quantity weights the edges decides what the network shows. Weighted by
adjusted residuals it is a residual network: edges show where a
transition is more frequent (over-represented) or less frequent
(avoided) than independence predicts, the reading of lag sequential
analysis. Weighted by transition probabilities it is the network of
Transition Network Analysis, which shows where the process goes next,
whether or not that is more than chance.

With `weights = "residuals"` (default), over-represented transitions are
solid blue edges and avoided transitions dashed red edges; `"yules_q"`
draws the same signed network on the bounded Yule's Q scale, and
`"lift"` the ratio of observed to expected counts. With `"prob"` (or its
aliases `"tna"` and `"relative"`) or `"count"`, the network is drawn in
the style of Transition Network Analysis, each node ringed by its
initial-state probability; probability edges below 0.05 are hidden
unless `edge_cutoff` is passed. Dense networks are thinned with
`significant` or, better at large sample sizes where nearly every cell
is significant, with `top`.

The same display is `plot(fit, type = "network")`. For the full table
use the residual heatmap of
[`plot.lsa()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md);
for the per-source breakdown,
[`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md)
or
[`plot_polar()`](https://pak.dynasite.org/lagdynamics/reference/plot_polar.md).

## Usage

``` r
plot_transitions(
  fit,
  weights = c("residuals", "tna", "relative", "count", "prob", "lift", "yules_q"),
  significant = FALSE,
  top = NULL,
  decimals = 1,
  node_fill = "white",
  edge_labels = TRUE,
  ...
)
```

## Arguments

- fit:

  An `lsa` fit from
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

- weights:

  Edge weight: `"residuals"` (default), `"tna"`, `"relative"`, `"prob"`
  (the three are the transition probability), `"count"`, `"lift"` or
  `"yules_q"`.

- significant:

  Keep only transitions whose residual p-value is below the fit's
  `alpha` (default `FALSE`).

- top:

  Keep only the strongest edges by absolute weight: a fraction below 1
  keeps that share of the edges, a number of 1 or more keeps that many.
  Default `NULL` (all edges). Applied after `significant`.

- decimals:

  Decimal places of the edge labels (default `1`).

- node_fill:

  Node colour (default `"white"`). The probability and count networks
  colour nodes by state unless `node_fill` is set.

- edge_labels:

  Show the edge weights as labels (default `TRUE`).

- ...:

  Passed to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html),
  for example `layout`, `node_shape`, `edge_cutoff` or `curvature`.

## Value

The network drawn by
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html),
invisibly; the plot is drawn on the current device.

## See also

[`plot.lsa()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md),
[`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md),
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md)

## Examples

``` r
fit <- lsa(engagement)
plot_transitions(fit)

plot_transitions(fit, weights = "prob")


fit9 <- lsa(group_regulation)
plot_transitions(fit9, top = 12)

plot_transitions(fit9, weights = "yules_q", significant = TRUE)

```
