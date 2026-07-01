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
```

## Arguments

- fit:

  An `lsa` object returned by
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md).
  Must be built from event-level data (sequences), not from a
  pre-computed transition matrix.

- ...:

  Reserved.

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

## Value

An object of class `c("lsa_reliability", "list")` with:

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

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods, 50*(1), 195-212.

For a grouped fit (`lsa_group`), reliability is estimated separately
within each group and the per-group `lsa_reliability` objects are
returned in an `lsa_reliability_group` container with its own print
method.

## See also

[`bootstrap_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/bootstrap_lsa.md),
[`stability_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/stability_lsa.md),
[`permute_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/permute_lsa.md)

## Examples

``` r
# \donttest{
fit <- lsa(engagement, engine = "classical")
rel <- reliability_lsa(fit, R = 50)
rel
#> <lsa_reliability>
#>   engine:        classical
#>   replicates:    50
#>   weights:       prob
#>   method:        pearson
#>   n sequences:   136
#>   split-half r:  0.972  (sd = 0.020)
#>   95% CI:        [0.918, 0.994]
# }
```
