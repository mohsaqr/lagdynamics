# Plot credible intervals of transition probabilities

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
result of
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
as a circular forest plot: one spoke per transition, showing the
observed transition probability and its credible interval, so the
precisely and imprecisely estimated transitions can be seen side by
side. The plot is drawn by
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)
with `metric = "prob"`; see that page for the layout and the arguments
that can be passed through `...`.

## Usage

``` r
# S3 method for class 'lsa_certainty'
plot(x, metric = "prob", ...)
```

## Arguments

- x:

  An `lsa_certainty` object from
  [`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md).

- metric:

  Interval to draw. Default `"prob"`, the only interval a certainty
  result has.

- ...:

  Passed to
  [`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md).

## Value

A `ggplot` object, drawn when printed.

## See also

[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md),
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md).

## Examples

``` r
cert <- certainty_lsa(lsa(group_regulation))
plot(cert)

```
