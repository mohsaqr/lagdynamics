# Convert a lag sequential analysis fit to a transition network for tna

`lsa_to_tna()` converts a lag sequential analysis fit to a transition
network for the tna package: a `tna` network from an `lsa` fit, or a
`group_tna` from a grouped fit, so that a network found by LSA can be
analysed with the network tools of tna. Lag sequential analysis and
Transition Network Analysis (TNA) read the same transitions from two
angles. LSA tests which transitions are more or less frequent than
independence predicts; TNA treats the transition probabilities as a
weighted network and studies its structure with centralities,
communities, cliques and resampling-based pruning.

The edge weights are chosen by `weights`: the transition probabilities
(the TNA model), the counts, the positive adjusted residuals (the
over-representation network) or the lift. With probabilities or counts
the network also carries the initial-state probabilities and the
sequences, so the resampling functions of tna can be used on it. With
residuals or lift it carries neither, because tna would re-estimate
probabilities from the sequences; resample with
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md)
instead.

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
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or a
  grouped fit from `lsa(..., group = )`.

- ...:

  Passed between methods.

- weights:

  Edge weights: `"prob"` (default, transition probabilities), `"count"`
  (observed counts), `"adj_res"` (adjusted residuals, negative values
  set to 0) or `"lift"` (observed over expected counts).

## Value

A `tna` object (for an `lsa` fit) or a `group_tna` object (one `tna` per
group, for a grouped fit), read with the functions of the tna package.

## Conditions

Errors carry the class `lagdynamics_error` and
`lagdynamics_missing_package` when tna is not installed. With
`weights = "adj_res"` or `"lift"` and a fit that has sequences, a
warning of class `lagdynamics_tna_no_data` notes that the network
carries no sequence data.

## See also

[`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md)
and
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
for the same parameters as tidy tables;
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

## Examples

``` r
fit <- lsa(engagement)
net <- lsa_to_tna(fit, weights = "prob")
tna::centralities(net)
#> # A tibble: 3 × 11
#>   state    OutStrength InStrength ClosenessIn ClosenessOut Closeness Betweenness
#> * <fct>          <dbl>      <dbl>       <dbl>        <dbl>     <dbl>       <dbl>
#> 1 Active         0.302      0.324      0.0811       0.0779     0.100           0
#> 2 Average        0.390      0.664      0.160        0.0973     0.160           2
#> 3 Disenga…       0.517      0.221      0.0691       0.101      0.114           0
#> # ℹ 4 more variables: BetweennessRSP <dbl>, Diffusion <dbl>, Clustering <dbl>,
#> #   PageRank <dbl>
```
