# Convert an lsa Fit to a tna Network

Convert an `lsa` fit to a `tna`-class network object usable by the `tna`
package's centrality, pruning, community, and bootstrap routines.

## Usage

``` r
lsa_to_tna(x, ...)

# S3 method for class 'lsa'
lsa_to_tna(x, weights = c("prob", "count", "adj_res", "lift"), ...)

# S3 method for class 'lsa_group'
lsa_to_tna(x, weights = c("prob", "count", "adj_res", "lift"), ...)
```

## Arguments

- x:

  An `lsa` fit from
  [`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md), or an
  `lsa_group` from `lsa(..., group = )`.

- ...:

  Method-specific arguments.

- weights:

  Character. Which matrix to expose as the tna edge weights. One of
  `"prob"` (row-normalised probabilities, default), `"count"` (raw
  observed counts), `"adj_res"` (adjusted residuals, over-representation
  network), or `"lift"` (observed / expected association strength).

## Value

For an `lsa` fit, a `tna` object with `weights`, `inits`, `labels`, and,
when available, `data` slots. For an `lsa_group`, a `group_tna` object.

## Details

This function is deliberately named `lsa_to_tna()`, not `as_tna()`, to
avoid overlapping export names with sibling packages.

## Examples

``` r
fit <- lsa(engagement)
net <- lsa_to_tna(fit, weights = "prob")
tna::centralities(net)
#> # A tibble: 3 × 10
#>   state    OutStrength InStrength ClosenessIn ClosenessOut Closeness Betweenness
#> * <fct>          <dbl>      <dbl>       <dbl>        <dbl>     <dbl>       <dbl>
#> 1 Active         0.302      0.324      0.0811       0.0779     0.100           0
#> 2 Average        0.390      0.664      0.160        0.0973     0.160           2
#> 3 Disenga…       0.517      0.221      0.0691       0.101      0.114           0
#> # ℹ 3 more variables: BetweennessRSP <dbl>, Diffusion <dbl>, Clustering <dbl>
```
