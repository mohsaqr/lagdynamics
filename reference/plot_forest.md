# Plot a circular forest of bootstrap confidence intervals

`plot_forest()` draws the bootstrap confidence intervals of a fit's
transitions as a circular forest plot, one spoke per transition: the
spoke spans the bootstrap confidence interval, a square marks the
observed estimate, and a dashed ring marks the null value. A forest plot
shows an estimate and its confidence interval for each of many
quantities, so that their size and precision can be compared at once;
here the quantities are the transitions, and the plot shows which
transitions keep their direction across resamples of the sequences and
how precisely each is estimated.

Transitions whose adjusted residual is stable across resamples are
coloured by direction (red over-represented, blue avoided) and the
others are grey. The interval can be that of the adjusted residual
(default), the count, the transition probability or Yule's Q.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a bootstrap
result draws the same display.

The intervals come from
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md);
read them as a table with `transitions(boot)`. For analytic intervals of
the transition probabilities, see
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
and its plot method.

## Usage

``` r
plot_forest(
  boot,
  metric = c("residuals", "count", "prob", "yules_q"),
  n_top = NULL,
  show_nonsig = TRUE,
  label_size = 2.6
)

# S3 method for class 'lsa_bootstrap'
plot(x, ...)
```

## Arguments

- boot:

  A bootstrap result from
  [`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md).

- metric:

  Interval to draw: `"residuals"` (default), `"count"`, `"prob"` or
  `"yules_q"`.

- n_top:

  Keep only the `n_top` transitions with the largest absolute estimate.
  Default `NULL` (all).

- show_nonsig:

  Draw transitions whose residual is not stable, in grey (default
  `TRUE`).

- label_size:

  Size of the transition labels (default `2.6`).

- x:

  A bootstrap result from
  [`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
  for [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- ...:

  Passed to `plot_forest()`, for example `metric` or `n_top`.

## Value

A `ggplot` object, drawn when printed.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a bootstrap
result returns the same object.

## Conditions

Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`
when the bootstrap holds no interval for `metric`, or when no transition
is left to draw.

## See also

[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`plot.lsa()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md),
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)

## Examples

``` r
set.seed(1)
boot <- bootstrap_lsa(lsa(engagement), R = 100)
plot_forest(boot)

plot(boot, metric = "prob")

# \donttest{
boot9 <- bootstrap_lsa(lsa(group_regulation), R = 500)
plot_forest(boot9, n_top = 20, show_nonsig = FALSE)

# }
```
