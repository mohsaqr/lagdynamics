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
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or
  an `lsa_group` from `lsa(..., group = )`.

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
if (FALSE) { # \dontrun{
fit <- lsa(engagement)
net <- lsa_to_tna(fit, weights = "prob")
tna::centralities(net)
} # }
```
