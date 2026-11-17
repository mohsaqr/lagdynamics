# Fit a lag sequential analysis model

`lsa()` fits a lag sequential analysis model to event sequences, wide
sequence tables, transition-count matrices, long event logs, or sequence
objects from other packages. Lag sequential analysis (Sackett, 1979;
Bakeman & Gottman, 1997) tests whether, in a stream of coded events, one
state follows another more or less often than chance predicts.

The transitions observed at a given lag are counted into a table of
antecedent (`from`) by consequent (`to`) states; an expected table is
derived under the hypothesis that the two are independent; and each cell
is tested with an adjusted residual (Haberman, 1979), a standard normal
deviate that is positive when the transition is over-represented and
negative when it is avoided. The significant cells define a lag
transition network, and the tablewise likelihood-ratio and Pearson tests
show whether there is any sequential dependence at all. For every cell,
`lsa()` estimates the observed and expected counts, the transition
probability, the adjusted residual and its p-value, Yule's Q, kappa and
lift. The `engine` selects the test: the classical analysis of Bakeman
and Quera (1995) by default; a two-cell (2 x 2) test with odds ratios;
Sackett's bidirectional test on the symmetrised table; or the parallel
and non-parallel dominance tests of Sackett (1979) and Wampold (1984).
Self-transitions or any other cells that cannot occur by construction
are declared as structural zeros, and the expected table is then fitted
by iterative proportional fitting (Christensen, 1997). With `group`, one
model is fitted per group on a shared state set so the groups can be
compared.

It returns an `lsa` object, or with `group` an `lsa_group` object
holding one fit per group.

The fit is read with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md),
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
and
[`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md),
and drawn with
[`plot.lsa()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md).
To fit several lags at once use
[`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md);
to assess the transitions use
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
or
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md);
to compare groups use
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md).
The engine wrappers `lsa_classical()`, `lsa_two_cell()`,
`lsa_bidirectional()`, `lsa_parallel_dominance()` and
`lsa_nonparallel_dominance()` are `lsa()` with `engine` fixed and return
the same object.

## Usage

``` r
lsa(
  data,
  lag = 1,
  engine = "classical",
  alternative = c("two.sided", "greater", "less"),
  alpha = 0.05,
  loops = TRUE,
  structural_zeros = NULL,
  labels = NULL,
  group = NULL,
  actor = NULL,
  action = NULL,
  time = NULL,
  order = NULL,
  session = NULL,
  time_threshold = 900,
  custom_format = NULL,
  is_unix_time = FALSE,
  unix_time_unit = "seconds",
  params = list(),
  ...
)

lsa_classical(data, ...)

lsa_two_cell(data, ...)

lsa_bidirectional(data, ...)

lsa_parallel_dominance(data, ...)

lsa_nonparallel_dominance(data, ...)
```

## Arguments

- data:

  Sequences in any form accepted by
  [`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md):
  a vector, a list of sequences, a wide table (one row per sequence), a
  square transition-count matrix, or a `tna`, `group_tna`,
  `tna_seq_data`, `tna_data`, `nestimate_data`, `netobject` or `stslist`
  object. With `action`, a long event log (one row per event). `NA` and
  empty strings are missing values, not a state; to model missingness as
  a state, recode it first.

- lag:

  Integer lag (default `1`). Positive lags pair each event with the one
  `lag` steps later, negative lags with the one `|lag|` steps earlier. A
  transition-count matrix supports `lag = 1` only.

- engine:

  Name of a registered engine (default `"classical"`); the built-in
  engines are listed by
  [`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md).

- alternative:

  Alternative hypothesis for the residual and kappa p-values:
  `"two.sided"` (default), `"greater"` or `"less"`.

- alpha:

  Significance level that flags a transition as significant (default
  `0.05`).

- loops:

  Keep self-transitions in the model (default `TRUE`). `FALSE` declares
  every self-transition a structural zero.

- structural_zeros:

  Optional `K x K` 0/1 matrix in which `0` marks a transition that
  cannot occur; default `NULL` (none). Combined with `loops = FALSE`,
  the diagonal is also forbidden. Supported by the classical engine.

- labels:

  Optional character vector of state names. Default `NULL` takes them
  from the data (sorted values of character input).

- group:

  Optional grouping: one value per sequence, or, with a long event log,
  the name of a column that is constant within each sequence. Default
  `NULL` (one model).

- actor, action, time, order, session:

  Column names for a long event log (all default `NULL`). Supplying
  `action` (the state column) switches to long-format input. Rows are
  grouped into sequences by `actor`, `session`, or both, ordered by
  `order` or else by `time`.

- time_threshold:

  Gap in seconds that starts a new sequence when `time` is given and
  `session` is not (default `900`).

- custom_format:

  Optional [`strptime()`](https://rdrr.io/r/base/strptime.html) format
  for the `time` column (default `NULL`).

- is_unix_time:

  Read the `time` column as a Unix epoch (default `FALSE`).

- unix_time_unit:

  Unit of the epoch: `"seconds"` (default), `"milliseconds"` or
  `"microseconds"`.

- params:

  Named list of engine-specific parameters (default
  [`list()`](https://rdrr.io/r/base/list.html)).

- ...:

  Engine-specific parameters, merged into `params`. For the engine
  wrappers, any argument of `lsa()` other than `engine`.

## Value

An object of class `c("lsa", "cograph_network")`. Read it with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
(one row per `from -> to` pair: `from`, `to`, `lag`, `count`,
`expected`, `prob`, `prob_col`, `adj_res`, `p`, `yules_q`, `kappa`,
`kappa_z`, `kappa_p`, `lift`, `sign`, `significant`, plus any
engine-specific columns),
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md)
(one row per state),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md)
(one row per tablewise test),
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
(one row per state),
[`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md)
and [`summary()`](https://rdrr.io/r/base/summary.html) (see
[lsa-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)).
The object also stores the per-cell quantities as `K x K` matrices, the
canonical
[`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md)
input and the parameters of the call, for the inference and plotting
functions.

With `group`, an object of class `c("lsa_group", "list")`: one `lsa` fit
per group level. The reading verbs return the same tables with a leading
`group` column, and [`summary()`](https://rdrr.io/r/base/summary.html)
returns one row per group.

## Long-format input

A long event log is turned into sequences before the model is fitted:
one sequence per `actor`, per `session`, or per actor-session pair, or a
single sequence when neither is given. Events are ordered within a
sequence by `order`, or else by `time`; when `time` is given without
`session`, a gap longer than `time_threshold` seconds starts a new
sequence. If the session identifier is already unique across actors,
pass it as `actor` and leave `session` unset. Check the number of
sequences reported by `summary(fit)` against the number expected.

## Conditions

Errors carry the class `lagdynamics_error` and one of:

- `lagdynamics_bad_input`: invalid data, labels, group, structural-zero
  matrix or long-format column.

- `lagdynamics_unknown_engine`: `engine` is not registered.

- `lagdynamics_missing_package`: the engine needs a package that is not
  installed.

- `lagdynamics_no_transitions`: the data contain no transition at the
  requested lag.

- `lagdynamics_unsupported_option`: the engine does not support
  structural zeros.

- `lagdynamics_ipf_infeasible`: the structural zeros leave a state with
  observed transitions but no permitted cell.

- `lagdynamics_needs_events`: `group` was given with a transition-count
  matrix, which cannot be split.

A warning of class `lagdynamics_structural_zero_counts` is raised when
the data contain transitions in cells declared structural zeros: the row
and column totals then include transitions the model declares
impossible, and the tablewise test and residuals are not interpretable.

## References

Bakeman, R., & Gottman, J. M. (1997). *Observing interaction: An
introduction to sequential analysis* (2nd ed.). Cambridge University
Press.

Bakeman, R., & Quera, V. (1995). *Analyzing interaction: Sequential
analysis with SDIS and GSEQ*. Cambridge University Press.

Christensen, R. (1997). *Log-linear models and logistic regression* (2nd
ed.). Springer.

Haberman, S. J. (1979). *Analysis of qualitative data: Volume 2, New
developments*. Academic Press.

Sackett, G. P. (1979). The lag sequential analysis of contingency and
cyclicity in behavioral interaction research. In J. D. Osofsky (Ed.),
*Handbook of infant development* (pp. 623-649). Wiley.

Wampold, B. E. (1984). Tests of dominance in sequential categorical
data. *Psychological Bulletin, 96*(2), 424-429.

The formula-by-formula source map is installed with the package at
`system.file("REFERENCES.md", package = "lagdynamics")`.

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md),
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
to read the fit;
[`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)
for several lags;
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)
and
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
for inference;
[`plot.lsa()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md);
[`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md)
and
[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md)
for engines;
[`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md)
for the accepted input forms.

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
transitions(fit, significant = TRUE, sort = "strength")
#>         from         to lag count  expected       prob   prob_col   adj_res
#> 1     Active     Active   1   459 247.03460 0.69756839 0.70506912  21.66281
#> 2 Disengaged Disengaged   1   157  59.97693 0.48307692 0.49062500  15.39047
#> 3    Average     Active   1   153 281.94983 0.20372836 0.23502304 -12.90604
#> 4     Active Disengaged   1    23 121.43022 0.03495441 0.07187500 -12.55693
#> 5    Average    Average   1   458 330.45732 0.60985353 0.60026212  12.45261
#> 6     Active    Average   1   176 289.53518 0.26747720 0.23066841 -11.31913
#> 7 Disengaged     Active   1    39 122.01557 0.12000000 0.05990783 -10.54949
#>               p    yules_q      kappa    kappa_z      kappa_p      lift  sign
#> 1 4.602765e-104  0.8278779  0.4438330  19.187756 4.684307e-82 1.8580393  over
#> 2  1.896612e-53  0.7544102  0.3135556  13.555614 7.339504e-42 2.6176731  over
#> 3  4.161672e-38 -0.6010580 -0.4994042 -14.696262 6.811736e-49 0.5426497 under
#> 4  3.642124e-36 -0.8264798 -0.8271733 -13.412058 5.139269e-41 0.1894092 under
#> 5  1.353374e-35  0.5530510  0.2274894   9.834808 7.971952e-23 1.3859581  over
#> 6  1.055115e-29 -0.5335259 -0.4241513 -12.481749 9.390092e-36 0.6078709 under
#> 7  5.107575e-26 -0.6983917 -0.7069460 -11.462654 2.031922e-30 0.3196313 under
#>   significant
#> 1        TRUE
#> 2        TRUE
#> 3        TRUE
#> 4        TRUE
#> 5        TRUE
#> 6        TRUE
#> 7        TRUE
tests(fit)
#>   test statistic df             p
#> 1 lrx2  618.3053  4 1.691520e-132
#> 2   x2  629.2022  4 7.406567e-135

# The same model through an engine wrapper, at lag 2
lsa_classical(engagement, lag = 2)
#> Lag Sequential Analysis  -  classical  (lag 2, directed)
#>   3 states | 1603 transitions | 1870 events | 136 sequences
#>   states: Active, Average, Disengaged
#>   independence: G² = 486.6, df = 4, p <2e-16
#> 
#>   Significant transitions (p < 0.05): 8 of 9
#>   strongest over-represented (of 3):
#>     Active -> Active          z =  +18.7  ***
#>     Disengaged -> Disengaged  z =  +15.2  ***
#>     Average -> Average        z =  +11.6  ***
#> 
#>   Initial states:
#>     Active     0.382  ████████████████████████
#>     Average    0.368  ███████████████████████
#>     Disengaged 0.250  ████████████████

# A long event log, one model per achievement group
gfit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
            time = "Time", group = "Achiever")
summary(gfit)
#>   group    engine lag n_states n_sequences n_events n_transitions n_significant
#> 1  High classical   1        9        1000    13721         12721            66
#> 2   Low classical   1        9        1000    13812         12812            67
#>   alpha     lrx2 lrx2_df lrx2_p       x2 x2_df x2_p
#> 1  0.05 6146.793      64      0 6806.682    64    0
#> 2  0.05 7714.714      64      0 9197.143    64    0
```
