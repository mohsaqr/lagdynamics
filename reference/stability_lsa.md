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
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

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

An object of class `c("lsa_stability", "list")`. Read the per-edge table
with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md);
see
[lsa-inference-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)
for [`print()`](https://rdrr.io/r/base/print.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).
Elements:

- edges:

  Tidy data frame, one row per `from -> to` cell: `from`, `to`,
  `observed_sig` (whether the cell was significant in the original fit),
  `stability` (fraction of subsamples in which the cell was
  significant), and `stable` (`stability >= min_stable`).

- stability_matrix:

  `R x K^2` 0/1 matrix recording per-cell significance across
  replicates.

- R, proportion, min_stable:

  Recipe metadata.

- fit:

  Reference to the original fit.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods, 50*(1), 195-212.

## See also

[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)

## Examples

``` r
fit <- lsa(engagement, engine = "classical")
st <- stability_lsa(fit, R = 100)
st
#> <lsa_stability>
#>   engine:        classical
#>   replicates:    100
#>   proportion:    80%
#>   min stable:    95%
#>   stable edges:  7 of 9 (>= 95% across replicates)
transitions(st, significant = TRUE)
#>         from         to observed_sig stability stable
#> 1     Active     Active         TRUE         1   TRUE
#> 2    Average     Active         TRUE         1   TRUE
#> 3 Disengaged     Active         TRUE         1   TRUE
#> 4     Active    Average         TRUE         1   TRUE
#> 5    Average    Average         TRUE         1   TRUE
#> 6     Active Disengaged         TRUE         1   TRUE
#> 7 Disengaged Disengaged         TRUE         1   TRUE
# \donttest{
# A full analysis uses more subsamples
transitions(stability_lsa(fit, R = 500), sort = "strength")
#>         from         to observed_sig stability stable
#> 1     Active     Active         TRUE     1.000   TRUE
#> 2    Average     Active         TRUE     1.000   TRUE
#> 3 Disengaged     Active         TRUE     1.000   TRUE
#> 4     Active    Average         TRUE     1.000   TRUE
#> 5    Average    Average         TRUE     1.000   TRUE
#> 6     Active Disengaged         TRUE     1.000   TRUE
#> 7 Disengaged Disengaged         TRUE     1.000   TRUE
#> 8 Disengaged    Average        FALSE     0.270  FALSE
#> 9    Average Disengaged        FALSE     0.012  FALSE
# }
```
