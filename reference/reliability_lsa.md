# Estimate the split-half reliability of a lag transition network

`reliability_lsa()` estimates the split-half reliability of a lag
transition network from one data set. A lag transition network is
reliable when two independent samples of the same population give
similar networks. Split-half reliability follows the resampling approach
that Epskamp, Borsboom and Fried (2018) describe for the accuracy of
psychological networks: the sequences are divided at random into two
halves, a network is estimated on each half, and the two networks are
compared transition by transition. A high correlation between the halves
means that the structure of the network does not depend on which
sequences happened to be sampled.

The random split is repeated `R` times. Each split assigns half of the
sequences to one half and the rest to the other, refits the model on
each with the same engine and settings, and correlates the two halves'
transition probabilities, counts or adjusted residuals over all
transitions. `reliability_lsa()` reports the mean, standard deviation
and 95% interval (2.5% and 97.5% quantiles) of the split-half
correlations, and the number of failed replicates, in which a half had
no transitions or a half's values did not vary. A grouped fit gives one
result per group.

Reliability describes the network as a whole. For the precision of each
transition use
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md)
or
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md);
for whether each transition stays significant when sequences are
dropped, use
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md).

## Usage

``` r
reliability_lsa(fit, ...)

# S3 method for class 'lsa'
reliability_lsa(
  fit,
  R = 100L,
  weights = c("prob", "count", "adj_res"),
  method = c("pearson", "spearman"),
  parallel = FALSE,
  n_cores = NULL,
  verbose = FALSE,
  ...
)

# S3 method for class 'lsa_group'
reliability_lsa(fit, ...)

# S3 method for class 'lsa_reliability_group'
print(x, ...)

# S3 method for class 'lsa_reliability'
print(x, ...)
```

## Arguments

- fit:

  An `lsa` fit from
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or a
  grouped fit from `lsa(..., group = )`, built from at least two
  sequences (not from a transition matrix).

- ...:

  For a grouped fit, the arguments above, passed on to each group;
  otherwise not used.

- R:

  Number of random splits. Default `100`.

- weights:

  Transition values correlated between the halves: `"prob"` (default,
  transition probabilities), `"count"` or `"adj_res"` (adjusted
  residuals).

- method:

  Correlation: `"pearson"` (default) or `"spearman"`.

- parallel:

  Run the splits on several cores. Default `FALSE`.

- n_cores:

  Number of cores when `parallel = TRUE`. Default `NULL`, all cores but
  one.

- verbose:

  Report progress every 100 splits. Default `FALSE`.

- x:

  An `lsa_reliability` or `lsa_reliability_group` object, for the
  [`print()`](https://rdrr.io/r/base/print.html) methods.

## Value

An object of class `c("lsa_reliability", "list")`. Read it with
[as.data.frame()](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_reliability.md):
`as.data.frame(x)` gives one row per split with its correlation (`NA`
for a failed replicate), and `as.data.frame(x, what = "summary")` one
row with the settings, the number of failed replicates `n_failed`, and
the `mean`, `sd`, `ci_low` and `ci_high` of the correlations.
[`print()`](https://rdrr.io/r/base/print.html) shows the same summary
and returns `x` invisibly. For a grouped fit, an object of class
`c("lsa_reliability_group", "list")` with one result per group, whose
tables gain a leading `group` column.

## Conditions

Errors carry the class `lagdynamics_error` and one of:

- `lagdynamics_needs_events`: the fit was built from a transition
  matrix, so there are no sequences to split.

- `lagdynamics_bad_input`: the fit has fewer than two sequences.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods, 50*(1), 195-212.

## See also

[`as.data.frame.lsa_reliability()`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_reliability.md)
to read the result;
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md)
and
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)
for transition-level checks.

## Examples

``` r
fit <- lsa(engagement)
rel <- reliability_lsa(fit, R = 50)
rel
#> <lsa_reliability>
#>   engine:        classical
#>   replicates:    50
#>   weights:       prob
#>   method:        pearson
#>   n sequences:   136
#>   split-half r:  0.965  (sd = 0.029)
#>   95% CI:        [0.896, 0.997]
#>   failed:        0 of 50 replicates
as.data.frame(rel, what = "summary")
#>   weights  method  R n_failed      mean         sd    ci_low   ci_high
#> 1    prob pearson 50        0 0.9653006 0.02936169 0.8959057 0.9971183

# One result per group
g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
         time = "Time", group = "Achiever")
reliability_lsa(g, R = 20) |>
  as.data.frame(what = "summary")
#>   group weights  method  R n_failed      mean          sd    ci_low   ci_high
#> 1  High    prob pearson 20        0 0.9866562 0.007352071 0.9691266 0.9935066
#> 2   Low    prob pearson 20        0 0.9864346 0.004408260 0.9776246 0.9926271
# \donttest{
# A full analysis uses more splits
reliability_lsa(fit, R = 500)
#> <lsa_reliability>
#>   engine:        classical
#>   replicates:    500
#>   weights:       prob
#>   method:        pearson
#>   n sequences:   136
#>   split-half r:  0.968  (sd = 0.027)
#>   95% CI:        [0.895, 0.996]
#>   failed:        0 of 500 replicates
# }
```
