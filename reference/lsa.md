# Lag Sequential Analysis

Fits a lag sequential analysis (LSA) on categorical event sequence data
using a registered engine. Returns a tidy S3 object with named slots for
observed/expected/probability/residual matrices and a long-format edge
table suitable for transition-network visualization.

## Usage

``` r
lsa_classical(data, ...)

lsa_two_cell(data, ...)

lsa_bidirectional(data, ...)

lsa_parallel_dominance(data, ...)

lsa_nonparallel_dominance(data, ...)

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
```

## Arguments

- data:

  Sequence input (any form accepted by
  [`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md)),
  *or* a raw long-format event-log `data.frame` when the `actor` /
  `action` arguments are supplied (see below). Accepted already-
  sequenced forms include vectors, lists of sequences, wide
  matrices/data.frames, transition-count matrices, and sequence- bearing
  objects (`tna`, `group_tna`, `tna_seq_data`, `tna_data`,
  `nestimate_data`, `netobject`, `stslist`); see
  [`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md).
  `NA` and empty-string cells are treated as missingness, not as a
  state: they are dropped wherever they occur and no transition is
  counted into or out of them. To model missingness as its own state,
  recode it (e.g. `NA -> "missing"`) before calling `lsa()`.

- ...:

  Additional engine-specific parameters (merged into `params`). For the
  engine wrappers (`lsa_classical()`, `lsa_two_cell()`,
  `lsa_bidirectional()`, `lsa_parallel_dominance()`,
  `lsa_nonparallel_dominance()`), any argument of `lsa()` other than
  `engine`.

- lag:

  Integer. The transition lag. Default `1`. Positive lags count
  successors (state at `t -> t + lag`); negative lags count predecessors
  (what occurred `|lag|` steps before); `0` pairs each event with itself
  (degenerate for single-stream event data – genuine co-occurrence needs
  concurrent codes, not yet supported). Pre-computed transition-matrix
  input supports `lag = 1` only. To analyse several lags at once, see
  [`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md).

- engine:

  Character scalar. The engine name, registered via
  [`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md).
  Built-in engines: `"classical"`, `"two_cell"`, `"bidirectional"`,
  `"parallel_dominance"`, `"nonparallel_dominance"`. Default
  `"classical"`.

- alternative:

  Character scalar. The alternative hypothesis for adjusted-residual and
  kappa p-values: one of `"two.sided"` (default), `"greater"`, or
  `"less"`.

- alpha:

  Numeric. Significance threshold used to mark edges as significant (the
  `significant` column of
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)).
  Default `0.05`.

- loops:

  Logical. Keep self-transitions (the diagonal)? **Default `TRUE`.** Set
  `loops = FALSE` to forbid every self-transition – the common reason to
  exclude cells – without building a matrix by hand.

- structural_zeros:

  Optional `K x K` 0/1 matrix for an *arbitrary* forbidden-cell pattern,
  where `0` marks a forbidden (structural-zero) cell and `1` an
  estimable one. **Default `NULL`: every cell is part of the model.**
  Combines with `loops`: `loops = FALSE` also zeros the diagonal of a
  supplied matrix. When any cell is forbidden the engine switches to
  iterative proportional fitting and Christensen's design-matrix
  residuals (see `inst/REFERENCES.md`, sections 2.2 and 4.2).

- labels:

  Optional character vector of state labels.

- group:

  Optional grouping for a multi-group fit. Either a vector with one
  entry per input sequence (length `n_sequences`), or — for
  **long-format** input (see `actor`/`action`) — the **name of a
  grouping column** in the log, which must be constant within each
  actor/session so each recovered sequence maps to one group. Sequences
  are partitioned by group and a separate `lsa` fit is built for each.
  All group fits share one global label set (derived from the full data)
  so their `K x K` matrices are directly comparable, even when a group
  never visits some state. Returns an `lsa_group` object (a named list
  of `lsa` fits). Requires event-level input; a pre-computed transition
  matrix cannot be split by group. Default `NULL` (single-group fit).

- actor, action, time, order, session:

  Column names (each a single string) for **long-format** event-log
  input. Supplying `action` switches `lsa()` into long-format mode: the
  raw log in `data` is sequenced by grouping rows per `actor`, per
  `session`, per `actor x session`, or as one global sequence when
  neither boundary is supplied. Events are ordered within each sequence
  by `order` if given else by `time`, and – when `time` is given and no
  `session` column is – a new sequence starts whenever the gap between
  consecutive events exceeds `time_threshold` seconds. All `NULL` by
  default; `action` is required in long-format mode (without it the
  input is taken as already sequenced). See the *Long-format input*
  section for a limit on combining `actor` and `session`.

- time_threshold:

  Numeric. Maximum gap in seconds between consecutive events before a
  new session is started in long-format mode. Default `900` (15
  minutes). Ignored unless `time` is given and `session` is not.

- custom_format:

  Optional `strptime` format string for parsing the `time` column (e.g.
  `"%Y-%m-%d %H:%M:%S"`). Default `NULL` (native date/time classes and
  ISO strings are parsed directly).

- is_unix_time:

  Logical. Treat the `time` column as a Unix epoch. Default `FALSE`.

- unix_time_unit:

  Character. Unit of the Unix epoch when `is_unix_time = TRUE`:
  `"seconds"` (default), `"milliseconds"`, or `"microseconds"`.

- params:

  Optional named list of engine-specific parameters forwarded to the
  engine function.

## Value

An object of class `c("lsa", "cograph_network")`. Read it with the verbs
rather than by reaching into slots:
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
for the tidy edge table,
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md),
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md),
and [`summary()`](https://rdrr.io/r/base/summary.html) for the other
results, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)/[`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md)
to draw it. Every number a verb returns is backed by these slots:

- edges:

  The tidy one-row-per-transition frame that backs
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
  (with extra `cograph_network` protocol columns).

- nodes:

  Data frame backing
  [`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md):
  `id, label, name, outgoing, incoming`.

- obs, exp, prob, prob_col, adj_res, p, yules_q, kappa, kappa_z,
  kappa_p:

  The same per-cell quantities as `edges`, in `K x K` matrix form
  (`prob` is row-conditional P(to \| from), `prob_col`
  column-conditional P(from \| to)). Convenient for matrix algebra; not
  the primary interface.

- lrx2, x2:

  Lists `(statistic, df, p)` backing
  [`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md):
  the tablewise likelihood-ratio (G^2) and Pearson chi-square tests of
  independence; `NULL` for engines without an expected table.

- inits:

  Named numeric vector backing
  [`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
  (proportion of sequences starting in each state, sums to 1); `NULL`
  for transition-matrix input.

- weights:

  `K x K` matrix used as the default edge weight for plotting. Equal to
  `obs` (counts) by default.

- directed:

  Logical scalar; `TRUE` for directed engines, `FALSE` for
  `bidirectional`.

- method:

  Engine name (the slot the `cograph_network` protocol reads). Also
  recorded in `params$engine`.

- data:

  The canonical `lsa_data` object (events + seq_id).

- params:

  Immutable snapshot of all parameters used (recipe), including
  `params$engine`.

- meta:

  List with source, IPF info, version, and call.

When `group` is supplied, returns an object of class
`c("lsa_group", "list")`: a named list of `lsa` fits (one per group
level) carrying `levels`, `group_sizes`, `labels`, and `engine`
attributes. Downstream verbs
([`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md),
[`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md),
etc.) dispatch on it and return grouped results.

## Details

The engine wrappers `lsa_classical()`, `lsa_two_cell()`,
`lsa_bidirectional()`, `lsa_parallel_dominance()` and
`lsa_nonparallel_dominance()` are `lsa()` with `engine` fixed to the
matching built-in engine; they return the same object.

## Long-format input

When both `actor` and `session` are given, each sequence is keyed by the
actor x session combination, built with
[`interaction()`](https://rdrr.io/r/base/interaction.html). That key is
an integer index over the full cross of distinct actors and distinct
sessions, so it overflows once the number of distinct actors times the
number of distinct sessions exceeds `.Machine$integer.max` (about
2.1e9). The overflowing keys become `NA` and those sequences are
dropped, with only R's warning "NAs produced by integer overflow". Check
that the sequence count reported by `summary(fit)` matches the number of
distinct sessions in the log. When the session identifier is already
unique across actors, pass it as `actor` and leave `session` `NULL`: no
cross is formed and every sequence is kept.

## References

Bakeman, R., & Quera, V. (1995). *Analyzing interaction: Sequential
analysis with SDIS and GSEQ*. Cambridge University Press.

Bakeman, R., & Gottman, J. M. (1997). *Observing interaction: An
introduction to sequential analysis* (2nd ed.). Cambridge University
Press.

Christensen, R. (1997). *Log-linear models and logistic regression* (2nd
ed.). Springer.

Haberman, S. J. (1979). *Analysis of qualitative data: Volume 2, New
developments*. Academic Press.

Sackett, G. P. (1979). The lag sequential analysis of contingency and
cyclicity in behavioral interaction research. In J. D. Osofsky (Ed.),
*Handbook of infant development* (pp. 623-649). Wiley.

Wampold, B. E. (1984). Tests of dominance in sequential categorical
data. *Psychological Bulletin, 96*(2), 424-429.

## See also

[`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md),
[`lsa_transitions()`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md),
[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md),
[`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md)

## Examples

``` r
seq <- c("Question", "Explain", "Agree",
         "Question", "Explain", "Elaborate",
         "Agree", "Question", "Explain")
fit <- lsa(seq, engine = "classical")
fit
#> Lag Sequential Analysis  -  classical  (lag 1, directed)
#>   4 states | 8 transitions | 9 events | 1 sequences
#>   states: Agree, Elaborate, Explain, Question
#>   independence: G² = 18.4, df = 9, p 0.0312
#> 
#>   Significant transitions (p < 0.05): 2 of 16
#>   strongest over-represented (of 2):
#>     Agree -> Question    z =   +2.8  ** 
#>     Question -> Explain  z =   +2.8  ** 
#> 
#>   Initial states:
#>     Question  1.000  ████████████████████████
#>     Agree     0.000  
#>     Elaborate 0.000  
#>     Explain   0.000  
transitions(fit)
#>         from        to lag count expected prob prob_col    adj_res           p
#> 1      Agree     Agree   1     0    0.500  0.0      0.0 -0.9428090 0.345778586
#> 2  Elaborate     Agree   1     1    0.250  1.0      0.5  1.8516402 0.064077506
#> 3    Explain     Agree   1     1    0.500  0.5      0.5  0.9428090 0.345778586
#> 4   Question     Agree   1     0    0.750  0.0      0.0 -1.2649111 0.205903211
#> 5      Agree Elaborate   1     0    0.250  0.0      0.0 -0.6172134 0.537093978
#> 6  Elaborate Elaborate   1     0    0.125  0.0      0.0 -0.4040610 0.686167850
#> 7    Explain Elaborate   1     1    0.250  0.5      1.0  1.8516402 0.064077506
#> 8   Question Elaborate   1     0    0.375  0.0      0.0 -0.8280787 0.407625948
#> 9      Agree   Explain   1     0    0.750  0.0      0.0 -1.2649111 0.205903211
#> 10 Elaborate   Explain   1     0    0.375  0.0      0.0 -0.8280787 0.407625948
#> 11   Explain   Explain   1     0    0.750  0.0      0.0 -1.2649111 0.205903211
#> 12  Question   Explain   1     3    1.125  1.0      1.0  2.8284271 0.004677735
#> 13     Agree  Question   1     2    0.500  1.0      1.0  2.8284271 0.004677735
#> 14 Elaborate  Question   1     0    0.250  0.0      0.0 -0.6172134 0.537093978
#> 15   Explain  Question   1     0    0.500  0.0      0.0 -0.9428090 0.345778586
#> 16  Question  Question   1     0    0.750  0.0      0.0 -1.2649111 0.205903211
#>       yules_q kappa    kappa_z     kappa_p     lift  sign significant
#> 1  -1.0000000 -1.00 -0.8081220 0.419020334 0.000000 under       FALSE
#> 2   1.0000000  1.00  1.8708287 0.061368829 4.000000  over       FALSE
#> 3   0.6666667  0.25  0.5345225 0.592980098 2.000000  over       FALSE
#> 4  -1.0000000 -1.00 -1.0690450 0.285049407 0.000000 under       FALSE
#> 5  -1.0000000 -1.00 -0.5345225 0.592980098 0.000000 under       FALSE
#> 6  -1.0000000 -1.00 -0.3535534 0.723673610 0.000000 under       FALSE
#> 7   1.0000000  1.00  1.4142136 0.157299207 4.000000  over       FALSE
#> 8  -1.0000000 -1.00 -0.7071068 0.479500122 0.000000 under       FALSE
#> 9  -1.0000000 -1.00 -1.0690450 0.285049407 0.000000 under       FALSE
#> 10 -1.0000000 -1.00 -0.7071068 0.479500122 0.000000 under       FALSE
#> 11 -1.0000000 -1.00 -1.4142136 0.157299207 0.000000 under       FALSE
#> 12  1.0000000  1.00  2.8284271 0.004677735 2.666667  over        TRUE
#> 13  1.0000000  1.00  2.1380899 0.032509445 4.000000  over        TRUE
#> 14 -1.0000000 -1.00 -0.7071068 0.479500122 0.000000 under       FALSE
#> 15 -1.0000000 -1.00 -1.4142136 0.157299207 0.000000 under       FALSE
#> 16 -1.0000000 -1.00 -1.4142136 0.157299207 0.000000 under       FALSE

# The same model through an engine wrapper
lsa_classical(seq)
#> Lag Sequential Analysis  -  classical  (lag 1, directed)
#>   4 states | 8 transitions | 9 events | 1 sequences
#>   states: Agree, Elaborate, Explain, Question
#>   independence: G² = 18.4, df = 9, p 0.0312
#> 
#>   Significant transitions (p < 0.05): 2 of 16
#>   strongest over-represented (of 2):
#>     Agree -> Question    z =   +2.8  ** 
#>     Question -> Explain  z =   +2.8  ** 
#> 
#>   Initial states:
#>     Question  1.000  ████████████████████████
#>     Agree     0.000  
#>     Elaborate 0.000  
#>     Explain   0.000  

# Long-format event log with a recorded grouping column
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
