# Case-Drop Stability for an LSA Fit

Resamples the underlying sequence data **without** replacement at a
specified retention proportion, refits the engine on each subsample, and
records which edges remain significant at the recipe's alpha threshold.
Returns a per-edge "stability" proportion: the fraction of subsamples in
which the edge was significant. Edges with stability \>= `min_stable`
(default `0.95`) are flagged as robust.

## Usage

``` r
stability_lsa(
  fit,
  R = 500L,
  proportion = 0.8,
  min_stable = 0.95,
  parallel = FALSE,
  n_cores = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- fit:

  An `lsa` object returned by
  [`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md).

- R:

  Integer. Number of subsamples. Default `500`.

- proportion:

  Numeric in (0, 1). Fraction of cases retained per subsample. Default
  `0.8`.

- min_stable:

  Numeric in (0, 1). Stability threshold for the `stable` flag in the
  output edge frame. Default `0.95`.

- parallel:

  Logical. Use multi-core resampling. Default `FALSE`.

- n_cores:

  Integer. Worker count when `parallel = TRUE`.

- verbose:

  Logical. Print progress every 100 replicates.

- ...:

  Reserved.

## Value

An object of class `c("lsa_stability", "list")` with:

- edges:

  Tidy per-edge data frame with `observed_sig` (whether the cell was
  significant in the original fit), `stability` (fraction of subsamples
  in which the cell was significant), and `stable`
  (`stability >= min_stable`).

- stability_matrix:

  `R x K^2` 0/1 matrix recording per-cell significance across
  replicates.

- R, proportion, min_stable:

  Recipe metadata.

- fit:

  Reference to the original fit.

## See also

[`bootstrap_lsa()`](https://saqr.me/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://saqr.me/lagdynamics/reference/permute_lsa.md)

## Examples

``` r
# \donttest{
fit <- lsa(engagement, engine = "classical")
st <- stability_lsa(fit, R = 100)
head(as.data.frame(st))
#>         from      to observed_sig stability stable
#> 1     Active  Active         TRUE      1.00   TRUE
#> 2    Average  Active         TRUE      1.00   TRUE
#> 3 Disengaged  Active         TRUE      1.00   TRUE
#> 4     Active Average         TRUE      1.00   TRUE
#> 5    Average Average         TRUE      1.00   TRUE
#> 6 Disengaged Average        FALSE      0.34  FALSE
# }
```
