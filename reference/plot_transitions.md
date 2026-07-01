# Plot the Transition Network

Draws the directed transition network of an `lsa` fit with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).
Pick the edge weight with `weights`; optionally keep only significant
edges. Nodes are white and edges are labelled by default. Returns the
`cograph` network invisibly.

## Usage

``` r
plot_transitions(
  fit,
  weights = c("residuals", "count", "prob", "lift", "yules_q"),
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
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md).

- weights:

  Which matrix becomes the edge weight, and how it is drawn:

  - `"residuals"` (default) – a **residual network** (not a transition
    one): adjusted residuals coloured by sign on the TNA / Nestimate
    convention, **blue = more** (over-represented) solid and **red =
    less** (avoided) dashed with a soft halo.

  - `"prob"` / `"count"` – the familiar **transition network** of
    Transition Network Analysis (TNA), drawn with
    `cograph::splot(tna_styling = TRUE)`: cograph's own TNA styling
    (coloured nodes, weighted directed edges) plus a donut ring per node
    carrying its initial-state probability, and edges labelled with the
    transition probability (`"prob"`) or observed count (`"count"`). For
    `"prob"`, edges below `0.05` are dropped by default so weak
    transitions do not clutter the plot (override with `edge_cutoff`).

  - `"lift"` – observed / expected, drawn in a single neutral colour
    with magnitude carried by edge width.

  - `"yules_q"` – a **signed association network**: Yule's Q on a fixed
    `[-1, 1]` scale, coloured by sign like the residual network (blue
    over-represented, red avoided) but bounded and not growing with
    sample size.

- significant:

  Logical. Keep only edges whose adjusted-residual p-value is below the
  fit's alpha; weaker cells are set to 0 (no edge). Default `FALSE`.
  Note that at large sample sizes almost every cell is significant, so
  this is a weak visual filter – prefer `top` (effect-size pruning) to
  declutter a dense residual network.

- top:

  Numeric or `NULL`. Keep only the strongest edges by absolute weight
  (applied after `significant`); the rest are set to 0. A fraction
  `0 < top < 1` keeps that **proportion** of the present edges
  (`top = 0.5` -\> the strongest half); a value `top >= 1` keeps that
  many edges (`top = 12` -\> the 12 strongest). The legible way to thin
  a dense residual network: it prunes by effect size
  (`|adjusted residual|`) rather than by p-value. `NULL` (default) keeps
  every edge. Applies to every view; for the probability network it
  composes with the default `edge_cutoff = 0.05`.

- decimals:

  Number of decimal places for edge labels. Default `1`.

- node_fill:

  Node fill colour. Default `"white"`; the probability / count networks
  use a per-state palette instead unless `node_fill` is set explicitly.

- edge_labels:

  Logical (or a label vector). Show edge weights as labels. Default
  `TRUE`.

- ...:

  Passed to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
  (e.g. `node_shape`, `layout`, `edge_cutoff`, `curvature`).

## Value

The `cograph_network` object, invisibly (drawn as a side effect).

## See also

[`plot.lsa()`](https://mohsaqr.github.io/lagdynamics/reference/plot.lsa.md)
(heatmap),
[`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md),
[`transition_probabilities()`](https://mohsaqr.github.io/lagdynamics/reference/transition_probabilities.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lsa(group_regulation)
plot_transitions(fit)                                   # residual network
plot_transitions(fit, weights = "prob")                 # probabilities
plot_transitions(fit, weights = "residuals",            # residual network,
                 significant = TRUE)                     #   significant only
plot_transitions(fit, top = 12)                         # 12 strongest edges
plot_transitions(fit, top = 0.5)                        # strongest 50%
plot_transitions(fit, decimals = 2)                     # 2-dp edge labels
plot_transitions(fit, node_shape = "square")            # splot passthrough
} # }
```
