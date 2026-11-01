# Split-Half Reliability for an LSA Fit

Estimates the reliability of an LSA network by repeated random
split-half resampling of sequences. Each replicate draws two disjoint
halves of the sequences without replacement, refits the engine on each
half, and computes the correlation between the two half-network
edge-weight vectors. Returns the distribution of replicate correlations
plus a point summary.

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

  An `lsa` object returned by
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or
  an `lsa_group` from `lsa(..., group = )`. Must be built from
  event-level data (sequences), not from a pre-computed transition
  matrix.

- ...:

  For the `lsa_group` method, the arguments above, passed to each
  group's fit; otherwise reserved (and ignored by
  [`print()`](https://rdrr.io/r/base/print.html)).

- R:

  Integer. Number of split-half replicates. Default `100`.

- weights:

  Character. Which edge matrix to correlate across halves: `"prob"`
  (default), `"count"`, or `"adj_res"`.

- method:

  Character. Correlation method: `"pearson"` (default) or `"spearman"`.

- parallel:

  Logical. Use multi-core resampling. Default `FALSE`.

- n_cores:

  Integer. Worker count when `parallel = TRUE`.

- verbose:

  Logical. Print progress every 100 replicates.

- x:

  An `lsa_reliability` or `lsa_reliability_group` object
  ([`print()`](https://rdrr.io/r/base/print.html) methods).

## Value

An object of class `c("lsa_reliability", "list")`.
[`print()`](https://rdrr.io/r/base/print.html) shows the mean split-half
correlation with its spread and interval and returns `x` invisibly;
[as.data.frame()](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_reliability.md)
gives the per-replicate correlations. Elements:

- correlations:

  Numeric vector of length `R`: the split-half correlation of each
  replicate.

- mean, sd:

  Mean and standard deviation of the finite replicate correlations.

- ci_low, ci_high:

  Empirical 2.5% and 97.5% quantiles.

- R, weights, method, n_sequences:

  Recipe metadata.

- fit:

  Reference to the original fit.

For an `lsa_group`, a named list of these, one per group, of class
`c("lsa_reliability_group", "list")`.

## Details

For a grouped fit (`lsa_group`), reliability is estimated separately
within each group and the per-group `lsa_reliability` objects are
returned in an `lsa_reliability_group` container.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods, 50*(1), 195-212.

## See also

[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)

## Examples

``` r
fit <- lsa(engagement, engine = "classical")
rel <- reliability_lsa(fit, R = 50)
rel
#> <lsa_reliability>
#>   engine:        classical
#>   replicates:    50
#>   weights:       prob
#>   method:        pearson
#>   n sequences:   136
#>   split-half r:  0.971  (sd = 0.029)
#>   95% CI:        [0.900, 0.995]
# \donttest{
# A full analysis uses more split-half replicates
reliability_lsa(fit, R = 500)
#> <lsa_reliability>
#>   engine:        classical
#>   replicates:    500
#>   weights:       prob
#>   method:        pearson
#>   n sequences:   136
#>   split-half r:  0.969  (sd = 0.026)
#>   95% CI:        [0.902, 0.996]
# }
```
