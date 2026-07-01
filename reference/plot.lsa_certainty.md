# Plot an Analytic-Certainty Result

Circular forest of the per-edge transition-probability credible
intervals from
[`certainty_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/certainty_lsa.md)
(delegates to
[`plot_forest()`](https://mohsaqr.github.io/lagdynamics/reference/plot_forest.md)
with `metric = "prob"`).

## Usage

``` r
# S3 method for class 'lsa_certainty'
plot(x, metric = "prob", ...)
```

## Arguments

- x:

  An `lsa_certainty` object.

- metric:

  Which credible interval to draw. Default `"prob"`.

- ...:

  Passed to
  [`plot_forest()`](https://mohsaqr.github.io/lagdynamics/reference/plot_forest.md).

## Value

A `ggplot` object (drawn when printed). Needs `ggplot2`.
