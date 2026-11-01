# Print, Summarise and Tidy an LSA Fit

Methods for the objects returned by
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).
[`print()`](https://rdrr.io/r/base/print.html) gives a compact model
overview: engine, sizes, the tablewise independence test, the strongest
over-represented transitions, and the initial-state distribution.
[`summary()`](https://rdrr.io/r/base/summary.html) prints that overview
followed by the node activity and the observed, expected, probability
and residual matrices, and returns a one-row tidy overview.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the full edge table, including the integer ids and `cograph_network`
protocol columns that
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
leaves out. [`nobs()`](https://rdrr.io/r/stats/nobs.html) is the number
of transitions the model was fitted to.

## Usage

``` r
# S3 method for class 'lsa_group'
print(x, ...)

# S3 method for class 'lsa'
print(x, ...)

# S3 method for class 'lsa'
summary(object, ...)

# S3 method for class 'lsa'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_group'
summary(object, ...)

# S3 method for class 'lsa'
nobs(object, ...)
```

## Arguments

- x, object:

  An `lsa` fit from
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or
  (for the `lsa_group` methods) a grouped fit from `lsa(..., group = )`.

- ...:

  Ignored.

- row.names:

  Optional row names for the returned data frame (`as.data.frame.lsa()`
  only); ignored by the grouped method.

- optional:

  Ignored; present for consistency with the generic.

## Value

- [`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

- `summary.lsa()` prints the extended view and returns, invisibly, a
  one-row `data.frame` with columns `engine`, `lag`, `n_states`,
  `n_sequences`, `n_events`, `n_transitions`, `n_significant`, `alpha`,
  `lrx2`, `lrx2_df`, `lrx2_p`, `x2`, `x2_df`, `x2_p` (`NA` where the
  engine computed no tablewise test or the input was a transition
  matrix).

- `summary.lsa_group()` returns that frame visibly with one row per
  group and a leading `group` column; it prints nothing else.

- `as.data.frame.lsa()` returns a `data.frame` with one row per
  `from -> to` cell: integer `from`/`to` ids, `from_label`, `to_label`,
  the statistic columns of
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
  and the protocol columns `edge` and `weight`.
  `as.data.frame.lsa_group()` stacks these with a leading `group`
  column.

- [`nobs()`](https://rdrr.io/r/stats/nobs.html) returns a single number,
  the total transition count.

## Details

For reading results, prefer the verbs:
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
for the per-transition table,
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md) and
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
for the other parts of the fit.

## See also

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)

## Examples

``` r
fit <- lsa(engagement)
print(fit)
#> Lag Sequential Analysis  -  classical  (lag 1, directed)
#>   3 states | 1734 transitions | 1870 events | 136 sequences
#>   states: Active, Average, Disengaged
#>   independence: G² = 618.3, df = 4, p <2e-16
#> 
#>   Significant transitions (p < 0.05): 7 of 9
#>   strongest over-represented (of 3):
#>     Active -> Active          z =  +21.7  ***
#>     Disengaged -> Disengaged  z =  +15.4  ***
#>     Average -> Average        z =  +12.5  ***
#> 
#>   Initial states:
#>     Active     0.382  ████████████████████████
#>     Average    0.368  ███████████████████████
#>     Disengaged 0.250  ████████████████
summary(fit)
#> Lag Sequential Analysis
#> =======================
#> Lag Sequential Analysis  -  classical  (lag 1, directed)
#>   3 states | 1734 transitions | 1870 events | 136 sequences
#>   states: Active, Average, Disengaged
#>   independence: G² = 618.3, df = 4, p <2e-16
#> 
#>   Significant transitions (p < 0.05): 7 of 9
#>   strongest over-represented (of 3):
#>     Active -> Active          z =  +21.7  ***
#>     Disengaged -> Disengaged  z =  +15.4  ***
#>     Average -> Average        z =  +12.5  ***
#> 
#>   Initial states:
#>     Active     0.382  ████████████████████████
#>     Average    0.368  ███████████████████████
#>     Disengaged 0.250  ████████████████
#> 
#> Node activity (share of transitions):
#>     Active      out 0.379 ██████████  in 0.375 ██████████
#>     Average     out 0.433 ████████████  in 0.440 ████████████
#>     Disengaged  out 0.187 █████  in 0.185 █████
#> 
#> Observed counts (obs):
#>            Active Average Disengaged
#> Active        459     176         23
#> Average       153     458        140
#> Disengaged     39     129        157
#> 
#> Expected counts (exp):
#>             Active Average Disengaged
#> Active     247.035 289.535    121.430
#> Average    281.950 330.457    138.593
#> Disengaged 122.016 143.007     59.977
#> 
#> Transitional probabilities (prob):
#>            Active Average Disengaged
#> Active      0.698   0.267      0.035
#> Average     0.204   0.610      0.186
#> Disengaged  0.120   0.397      0.483
#> 
#> Adjusted residuals (adj_res):
#>             Active Average Disengaged
#> Active      21.663 -11.319    -12.557
#> Average    -12.906  12.453      0.176
#> Disengaged -10.549  -1.736     15.390
nobs(fit)
#> [1] 1734
head(as.data.frame(fit))
#>   from to from_label to_label lag count expected      prob   prob_col
#> 1    1  1     Active   Active   1   459 247.0346 0.6975684 0.70506912
#> 2    2  1    Average   Active   1   153 281.9498 0.2037284 0.23502304
#> 3    3  1 Disengaged   Active   1    39 122.0156 0.1200000 0.05990783
#> 4    1  2     Active  Average   1   176 289.5352 0.2674772 0.23066841
#> 5    2  2    Average  Average   1   458 330.4573 0.6098535 0.60026212
#> 6    3  2 Disengaged  Average   1   129 143.0075 0.3969231 0.16906946
#>      adj_res             p    yules_q      kappa    kappa_z      kappa_p
#> 1  21.662811 4.602765e-104  0.8278779  0.4438330  19.187756 4.684307e-82
#> 2 -12.906038  4.161672e-38 -0.6010580 -0.4994042 -14.696262 6.811736e-49
#> 3 -10.549486  5.107575e-26 -0.6983917 -0.7069460 -11.462654 2.031922e-30
#> 4 -11.319132  1.055115e-29 -0.5335259 -0.4241513 -12.481749 9.390092e-36
#> 5  12.452615  1.353374e-35  0.5530510  0.2274894   9.834808 7.971952e-23
#> 6  -1.736461  8.248229e-02 -0.1083175 -0.1618196  -2.964801 3.028788e-03
#>        lift  sign significant                  edge weight
#> 1 1.8580393  over        TRUE      Active -> Active    459
#> 2 0.5426497 under        TRUE     Average -> Active    153
#> 3 0.3196313 under        TRUE  Disengaged -> Active     39
#> 4 0.6078709 under        TRUE     Active -> Average    176
#> 5 1.3859581  over        TRUE    Average -> Average    458
#> 6 0.9020506 under       FALSE Disengaged -> Average    129

gfit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
            time = "Time", group = "Achiever")
gfit
#> <lsa_group>
#>   engine:    classical
#>   states:    9 (adapt, cohesion, consensus, coregulate, discuss, emotion, monitor, plan, synthesis)
#>   groups:    2
#>     - High:        1000 sequences
#>     - Low:         1000 sequences
summary(gfit)
#>   group    engine lag n_states n_sequences n_events n_transitions n_significant
#> 1  High classical   1        9        1000    13721         12721            66
#> 2   Low classical   1        9        1000    13812         12812            67
#>   alpha     lrx2 lrx2_df lrx2_p       x2 x2_df x2_p
#> 1  0.05 6146.793      64      0 6806.682    64    0
#> 2  0.05 7714.714      64      0 9197.143    64    0
```
