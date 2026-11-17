# Plot a comparison of transitions between groups

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
result of
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
or
[`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
as a barrel chart (default) or a heatmap. When the same states are
observed in two or more groups, a comparison shows which transitions
differ between them, by how much, and whether the difference is larger
than chance, so a display of it has to show both the transitions
themselves in each group and the test of their difference.

The barrel is a back-to-back chart with one row per transition: the
first group's bar runs to the left and the second's to the right, bar
length is the transition probability (or count), bar colour is the log
odds ratio of the transition within the group (blue over-represented,
red avoided), and the centre shows the p-value of the difference, in
bold when significant. The rows are the `top_n` most frequent
transitions, or with `rank = "effect"` those with the strongest
association in either group. The heatmap draws the difference between
the groups as a grid of `from` by `to` states. With more than two
groups, each pair is drawn as a facet.

The tested differences are read as a table with `transitions(cmp)`. To
colour the ribbons of a chord diagram by a group difference, see
[`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md).

## Usage

``` r
# S3 method for class 'lsa_comparison'
plot(
  x,
  style = c("barrel", "heatmap"),
  value = c("prob", "count"),
  rank = c("frequency", "effect"),
  top_n = 12L,
  ...
)

# S3 method for class 'lsa_comparison_pairwise'
plot(
  x,
  style = c("barrel", "heatmap"),
  value = c("prob", "count"),
  rank = c("frequency", "effect"),
  top_n = 12L,
  ...
)
```

## Arguments

- x:

  The result of
  [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
  or
  [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md).

- style:

  `"barrel"` (default) or `"heatmap"`.

- value:

  For the barrel, bar length: `"prob"` (default) or `"count"`.

- rank:

  For the barrel, which transitions to show: `"frequency"` (default, the
  most frequent overall) or `"effect"` (the strongest association,
  largest absolute log odds ratio, in either group).

- top_n:

  For the barrel, the number of transitions shown (default `12`); with
  several pairs the same transitions are shown in every facet.

- ...:

  Ignored.

## Value

A `ggplot` object, drawn when printed.

## See also

[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md),
[`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md),
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md)

## Examples

``` r
gfit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
            time = "Time", group = "Achiever")
set.seed(1)
cmp <- compare_lsa(gfit, R = 100, adjust = "BH")
plot(cmp)

plot(cmp, rank = "effect")

plot(cmp, style = "heatmap")

```
