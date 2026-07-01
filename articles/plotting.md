# Plotting lag-sequential models

`lagdynamics` has one plotting verb, `plot(fit, type = )`, plus
dedicated plots for the resampling and comparison results. Every view of
a fit displays the same quantity – the adjusted residual, the departure
of a transition from chance – and differs only in geometry. The single
fit used throughout is the bundled `engagement` data (138 students,
three engagement states).

| view | call | backend |
|----|----|----|
| residual heatmap | `plot(fit)` | ggplot2 |
| residual network | `plot(fit, type = "network")` | cograph |
| transition (TNA) network | `plot(fit, type = "network", weights = "prob")` | cograph |
| chord diagram | `plot(fit, type = "chord")` | cograph |
| polar sunburst | `plot(fit, type = "sunburst")` | ggplot2 |
| uncertainty forest | `plot(bootstrap_lsa(fit))`, `plot(certainty_lsa(fit))` | ggplot2 |
| group barrel | `plot(compare_lsa(g))`, `plot(bayes_compare_lsa(g))` | ggplot2 |

**Colour.** Two conventions are used, by purpose. The residual
*heatmap*, *chord*, and *sunburst* use the residual diverging scale
(warm = over- represented, cool = avoided). The *network* and *group
comparison* plots follow the wider Transition Network Analysis (TNA)
convention – **blue = more than chance**, **red = less** (avoided edges
dashed) – so an `lsa` network reads like any other transition network.

## Residual heatmap

The default. Rows are the current state, columns the next; colour is the
adjusted residual. `which` selects the matrix.

``` r

plot(fit)                       # adjusted residuals (default)
```

![](plotting_files/figure-html/heatmap-1.png)

``` r

plot(fit, which = "prob")       # transition probabilities
```

![](plotting_files/figure-html/heatmap-2.png)

## Residual network

The same residuals as a directed graph: **blue** edges (solid) are
over-represented, **red** edges (dashed, with a soft halo) are avoided.

``` r

plot(fit, type = "network")
```

![](plotting_files/figure-html/net-residual-1.png)

## Transition network (a TNA model)

Weight the edges by probability and the view becomes the familiar
transition network, drawn in the Transition Network Analysis (TNA) style
by
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html):
coloured nodes, a per-node initial-probability ring, and weighted
directed edges. This shows what happens, where the residual network
shows what is surprising.

``` r

plot(fit, type = "network", weights = "prob")
```

![](plotting_files/figure-html/net-transition-1.png)

## Chord and sunburst

A chord diagram of the transition flow, and a polar sunburst of each
state’s outgoing distribution.

``` r

plot(fit, type = "chord")
```

![](plotting_files/figure-html/chord-1.png)

``` r

plot(fit, type = "sunburst")                 # rose (default)
```

![](plotting_files/figure-html/sunburst-1.png)

``` r

plot(fit, type = "sunburst", style = "wedge")  # frequency wedges
```

![](plotting_files/figure-html/sunburst-2.png)

## Uncertainty forests

A fitted edge is one estimate; the forest shows its interval. Both the
resampling bootstrap and the analytic certainty plot as a circular
forest of per-edge intervals.

``` r

plot(bootstrap_lsa(fit, R = 200))            # resampling CIs
```

![](plotting_files/figure-html/forest-1.png)

``` r

plot(certainty_lsa(fit))                     # analytic credible intervals
```

![](plotting_files/figure-html/forest-2.png)

## Group comparison

For a real grouping we use the bundled long event log
`group_regulation_long`, and fit one model per achievement group.

``` r

gfit <- lsa(group_regulation_long, actor = "Actor",
            action = "Action", time = "Time", group = "Achiever")
```

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a comparison
draws a back-to-back **barrel**: each group’s bar runs to one side, the
higher group’s bar is bordered (darker for a larger difference), and the
centre chip carries the difference p-value.

``` r

cmp <- compare_lsa(gfit, R = 500)
plot(cmp)                                    # frequency-ranked rows
```

![](plotting_files/figure-html/barrel-1.png)

``` r

plot(cmp, rank = "effect")                   # surface avoided (red) edges
```

![](plotting_files/figure-html/barrel-2.png)

``` r

plot(cmp, style = "heatmap")                 # full-grid difference
```

![](plotting_files/figure-html/barrel-3.png)

The Bayesian comparison reuses the same barrel, with the credible
posterior in place of the permutation p-value.

``` r

plot(bayes_compare_lsa(gfit, seed = 1))
```

![](plotting_files/figure-html/bayes-barrel-1.png)

## Grouped fits

Every `type` works on a grouped fit, drawing one panel per group
(`combined = FALSE`, the default) or a single tiled figure
(`combined = TRUE`).

``` r

plot(gfit)                                   # one heatmap per group
```

![](plotting_files/figure-html/grouped-1.png)![](plotting_files/figure-html/grouped-2.png)

## The worker functions

`plot(fit, type = )` is a front door over exported workers, each with
its own full argument list:

``` r

plot_transitions(fit, weights = "residuals", significant = TRUE)  # network
plot_chords(fit, width = "prob")                                  # chord
plot_polar(fit, style = "wedge")                                  # sunburst
plot_forest(bootstrap_lsa(fit), n_top = 12)                       # forest
```

## In short

``` r

plot(fit)                                   # residual heatmap (default)
plot(fit, type = "network")                 # residual network (blue = more)
plot(fit, type = "network", weights = "prob")  # transition network (a TNA model)
plot(fit, type = "chord"); plot(fit, type = "sunburst")
plot(bootstrap_lsa(fit)); plot(certainty_lsa(fit))   # uncertainty forests
plot(compare_lsa(gfit)); plot(bayes_compare_lsa(gfit))  # group barrels
plot(gfit)                                  # grouped: one panel per group
```
