# Print, summarise and tidy a lag sequential analysis fit

[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) describe a fitted lag
sequential analysis, an `lsa` object from
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md). A fit
has three levels of description: the model as a whole (its size and the
tablewise test of independence), the states (how often each is left and
entered, and how often a sequence starts in it), and the individual
transitions. These methods give the first two levels;
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
gives the third.

[`print()`](https://rdrr.io/r/base/print.html) shows a compact overview:
the engine and lag, the numbers of states, transitions, events and
sequences, the likelihood-ratio test of independence, the number of
significant transitions out of the cells tested (a structural zero, such
as a self-transition with `loops = FALSE`, has no p-value and is not
counted), the strongest over-represented transitions and the
initial-state distribution.
[`summary()`](https://rdrr.io/r/base/summary.html) returns a
`summary.lsa` object holding two tables, a one-row overview of the fit
and one row per state, which are read with
`as.data.frame(summary(fit), what = )`. For a grouped fit,
[`summary()`](https://rdrr.io/r/base/summary.html) returns one overview
row per group.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on a fit
returns the full edge table, including the integer state ids and the
network columns that
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
leaves out. [`nobs()`](https://rdrr.io/r/stats/nobs.html) is the number
of transitions the model was fitted to.

For the per-transition, per-state and per-test tables, prefer the
reading verbs
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md) and
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md).

## Usage

``` r
# S3 method for class 'lsa_group'
print(x, ...)

# S3 method for class 'lsa'
print(x, ...)

# S3 method for class 'lsa'
summary(object, ...)

# S3 method for class 'summary.lsa'
print(x, ...)

# S3 method for class 'summary.lsa'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("overview", "nodes"),
  ...
)

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
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), or a
  grouped fit from `lsa(..., group = )`. For
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on a
  summary, a `summary.lsa` object.

- ...:

  Ignored.

- row.names:

  Optional row names for the returned data frame; ignored by the grouped
  method.

- optional:

  Ignored; present for consistency with the generic.

- what:

  For [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on
  a summary: `"overview"` (default, the one-row overview of the fit) or
  `"nodes"` (one row per state).

## Value

- [`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

- [`summary()`](https://rdrr.io/r/base/summary.html) on an `lsa` fit
  returns an object of class `summary.lsa` holding two tables. The
  overview has one row and the columns `engine`, `lag`, `n_states`,
  `n_sequences`, `n_events`, `n_transitions`, `n_significant`, `alpha`,
  `lrx2`, `lrx2_df`, `lrx2_p`, `x2`, `x2_df`, `x2_p` (`NA` where the
  engine computed no tablewise test or the input was a transition
  matrix). The node table has one row per state and the columns `state`,
  `outgoing`, `incoming`, `out_share`, `in_share`, `init_prob`. Read
  them with `as.data.frame(s)` and `as.data.frame(s, what = "nodes")`.

- [`summary()`](https://rdrr.io/r/base/summary.html) on a grouped fit
  returns a `data.frame` with one row per group: a leading `group`
  column followed by the overview columns.

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on an
  `lsa` fit returns a `data.frame` with one row per `from -> to` cell:
  integer `from` and `to` ids, `from_label`, `to_label`, the statistic
  columns of
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
  and the network columns `edge` and `weight`. On a grouped fit the
  tables are stacked under a leading `group` column.

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on a
  summary returns the overview or the node table.

- [`nobs()`](https://rdrr.io/r/stats/nobs.html) returns a single number,
  the total transition count.

## See also

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md);
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md) and
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
for the reading verbs.

## Examples

``` r
fit <- lsa(engagement)
fit
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
s <- summary(fit)
s
#> Lag sequential analysis: classical engine, lag 1
#>   3 states, 1734 transitions, 1870 events, 136 sequences
#>   significant transitions (p < 0.05): 7 of 9
#>   likelihood-ratio test: G² = 618.31, df = 4, p <2e-16
#>   Pearson test: χ² = 629.20, df = 4, p <2e-16
#> 
#>   State activity (share of transitions out and in):
#>     Active      out 0.379 ██████████    in 0.375 ██████████
#>     Average     out 0.433 ████████████  in 0.440 ████████████
#>     Disengaged  out 0.187 █████         in 0.185 █████
as.data.frame(s)
#>      engine lag n_states n_sequences n_events n_transitions n_significant alpha
#> 1 classical   1        3         136     1870          1734             7  0.05
#>       lrx2 lrx2_df       lrx2_p       x2 x2_df          x2_p
#> 1 618.3053       4 1.69152e-132 629.2022     4 7.406567e-135
as.data.frame(s, what = "nodes")
#>        state outgoing incoming out_share  in_share init_prob
#> 1     Active      658      651 0.3794694 0.3754325 0.3823529
#> 2    Average      751      763 0.4331027 0.4400231 0.3676471
#> 3 Disengaged      325      320 0.1874279 0.1845444 0.2500000
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
