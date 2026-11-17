# Plot a chord diagram of transitions

`plot_chords()` draws the transitions of a fit as a chord diagram with
cograph: the states are arcs around a circle and each transition is a
ribbon from its source arc to its target arc. The width of a ribbon
shows how much of the traffic the transition carries, and its colour can
show whether the transition is more or less frequent than independence
predicts, so the diagram shows where the transitions out of each state
go, and which of the frequent routes are also over-represented.

Ribbon width is the count (default) or the transition probability;
ribbon colour is the adjusted residual on a diverging scale (red
over-represented, blue avoided), or the lift, probability or count on a
sequential scale. With `compare`, a second fit on the same states, each
ribbon is coloured by the difference between the two fits in the colour
quantity, which shows where two groups differ.

The same display is `plot(fit, type = "chord")`. For the full table use
the residual heatmap of
[`plot.lsa()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md);
for a network layout,
[`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md);
for a test of group differences,
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md).

## Usage

``` r
plot_chords(
  fit,
  compare = NULL,
  width = c("count", "prob"),
  color = c("residuals", "lift", "prob", "count"),
  significant = FALSE,
  self_loops = TRUE,
  alpha = 0.6,
  ...
)
```

## Arguments

- fit:

  An `lsa` fit from
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

- compare:

  Optional second `lsa` fit with the same states; ribbons are then
  coloured by `fit` minus `compare`. Default `NULL`.

- width:

  Ribbon width: `"count"` (default) or `"prob"`.

- color:

  Ribbon colour: `"residuals"` (default), `"lift"`, `"prob"` or
  `"count"`.

- significant:

  Draw only transitions whose residual p-value is below the fit's
  `alpha` (default `FALSE`); ignored with `compare`.

- self_loops:

  Draw self-transitions (default `TRUE`).

- alpha:

  Ribbon opacity (default `0.6`).

- ...:

  Passed to
  [`cograph::plot_chord()`](https://sonsoles.me/cograph/reference/plot_chord.html),
  for example `ticks`, `segment_width`, `label_size` or `title`.

## Value

The value of
[`cograph::plot_chord()`](https://sonsoles.me/cograph/reference/plot_chord.html)
(the drawn segments and chords), invisibly; the plot is drawn on the
current device.

## Conditions

Errors carry the class `lagdynamics_error` and
`lagdynamics_incompatible_fits` when `fit` and `compare` have different
states.

## See also

[`plot.lsa()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md),
[`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md),
[`plot_polar()`](https://pak.dynasite.org/lagdynamics/reference/plot_polar.md),
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)

## Examples

``` r
fit <- lsa(engagement)
plot_chords(fit)

plot_chords(fit, width = "prob", self_loops = FALSE)


# Ribbons coloured by the difference between two groups
high <- lsa(subset(group_regulation_long, Achiever == "High"),
            actor = "Actor", action = "Action", time = "Time")
low <- lsa(subset(group_regulation_long, Achiever == "Low"),
           actor = "Actor", action = "Action", time = "Time")
plot_chords(high, compare = low)

```
