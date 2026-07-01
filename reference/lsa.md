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
  [`lsa_data()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_data.md)),
  *or* a raw long-format event-log `data.frame` when the `actor` /
  `action` arguments are supplied (see below). Accepted already-
  sequenced forms include vectors, lists of sequences, wide
  matrices/data.frames, transition-count matrices, and sequence- bearing
  objects (`tna`, `group_tna`, `nestimate_data`, `stslist`). `NA` and
  empty-string cells are treated as missingness, not as a state: they
  are dropped wherever they occur and no transition is counted into or
  out of them. To model missingness as its own state, recode it (e.g.
  `NA -> "missing"`) before calling `lsa()`.

- ...:

  Additional engine-specific parameters (merged into `params`).

- lag:

  Integer. The transition lag. Default `1`. Positive lags count
  successors (state at `t -> t + lag`); negative lags count predecessors
  (what occurred `|lag|` steps before); `0` pairs each event with itself
  (degenerate for single-stream event data – genuine co-occurrence needs
  concurrent codes, not yet supported). Pre-computed transition-matrix
  input supports `lag = 1` only. To analyse several lags at once, see
  [`lsa_lags()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_lags.md).

- engine:

  Character scalar. The engine name, registered via
  [`register_lsa_engine()`](https://mohsaqr.github.io/lagdynamics/reference/register_lsa_engine.md).
  Built-in engines: `"classical"`, `"two_cell"`, `"bidirectional"`,
  `"parallel_dominance"`, `"nonparallel_dominance"`. Default
  `"classical"`.

- alternative:

  Character scalar. The alternative hypothesis for adjusted-residual and
  kappa p-values: one of `"two.sided"` (default), `"greater"`, or
  `"less"`.

- alpha:

  Numeric. Significance threshold used to mark edges as significant in
  `fit$edges$significant`. Default `0.05`.

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
  residuals (see `inst/REFERENCES.md` §2.2, §4.2).

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
  input. Supplying `action` (and `actor`) switches `lsa()` into
  long-format mode: the raw log in `data` is sequenced into event
  sequences by grouping rows per `actor` (optionally crossed with an
  explicit `session` id), ordering within each group by `order` if given
  else by `time`, and – when `time` is given and no `session` column is
  – starting a new session whenever the gap between consecutive events
  exceeds `time_threshold` seconds. All `NULL` by default (input is
  taken as already-sequenced). Cannot be combined with `group`.

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
[`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md)
for the tidy edge table,
[`nodes()`](https://mohsaqr.github.io/lagdynamics/reference/nodes.md),
[`tests()`](https://mohsaqr.github.io/lagdynamics/reference/tests.md),
[`initial()`](https://mohsaqr.github.io/lagdynamics/reference/initial.md),
and [`summary()`](https://rdrr.io/r/base/summary.html) for the other
results, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)/[`plot_transitions()`](https://mohsaqr.github.io/lagdynamics/reference/plot_transitions.md)
to draw it. Every number a verb returns is backed by these slots:

- edges:

  The tidy one-row-per-transition frame that backs
  [`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md)
  (with extra `cograph_network` protocol columns).

- nodes:

  Data frame backing
  [`nodes()`](https://mohsaqr.github.io/lagdynamics/reference/nodes.md):
  `id, label, name, outgoing, incoming`.

- obs, exp, prob, prob_col, adj_res, p, yules_q, kappa, kappa_z,
  kappa_p:

  The same per-cell quantities as `edges`, in `K x K` matrix form
  (`prob` is row-conditional P(to \| from), `prob_col`
  column-conditional P(from \| to)). Convenient for matrix algebra; not
  the primary interface.

- lrx2, x2:

  Lists `(statistic, df, p)` backing
  [`tests()`](https://mohsaqr.github.io/lagdynamics/reference/tests.md):
  the tablewise likelihood-ratio (G^2) and Pearson chi-square tests of
  independence; `NULL` for engines without an expected table.

- inits:

  Named numeric vector backing
  [`initial()`](https://mohsaqr.github.io/lagdynamics/reference/initial.md)
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
([`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md),
[`transition_probabilities()`](https://mohsaqr.github.io/lagdynamics/reference/transition_probabilities.md),
[`reliability_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/reliability_lsa.md),
etc.) dispatch on it and return grouped results.

## See also

[`lsa_data()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_data.md),
[`lsa_transitions()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_transitions.md),
[`register_lsa_engine()`](https://mohsaqr.github.io/lagdynamics/reference/register_lsa_engine.md),
[`list_lsa_engines()`](https://mohsaqr.github.io/lagdynamics/reference/list_lsa_engines.md)

## Examples

``` r
seq <- c("Question", "Explain", "Agree",
         "Question", "Explain", "Elaborate",
         "Agree", "Question", "Explain")
fit <- lsa(seq, engine = "classical")
fit
#> Lag Sequential Analysis  —  classical  (lag 1, directed)
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
head(fit$edges)
#>   from to from_label  to_label lag count expected prob prob_col    adj_res
#> 1    1  1      Agree     Agree   1     0    0.500  0.0      0.0 -0.9428090
#> 2    2  1  Elaborate     Agree   1     1    0.250  1.0      0.5  1.8516402
#> 3    3  1    Explain     Agree   1     1    0.500  0.5      0.5  0.9428090
#> 4    4  1   Question     Agree   1     0    0.750  0.0      0.0 -1.2649111
#> 5    1  2      Agree Elaborate   1     0    0.250  0.0      0.0 -0.6172134
#> 6    2  2  Elaborate Elaborate   1     0    0.125  0.0      0.0 -0.4040610
#>            p    yules_q kappa    kappa_z    kappa_p lift  sign significant
#> 1 0.34577859 -1.0000000 -1.00 -0.8081220 0.41902033    0 under       FALSE
#> 2 0.06407751  1.0000000  1.00  1.8708287 0.06136883    4  over       FALSE
#> 3 0.34577859  0.6666667  0.25  0.5345225 0.59298010    2  over       FALSE
#> 4 0.20590321 -1.0000000 -1.00 -1.0690450 0.28504941    0 under       FALSE
#> 5 0.53709398 -1.0000000 -1.00 -0.5345225 0.59298010    0 under       FALSE
#> 6 0.68616785 -1.0000000 -1.00 -0.3535534 0.72367361    0 under       FALSE
#>                     edge weight
#> 1         Agree -> Agree      0
#> 2     Elaborate -> Agree      1
#> 3       Explain -> Agree      1
#> 4      Question -> Agree      0
#> 5     Agree -> Elaborate      0
#> 6 Elaborate -> Elaborate      0
```
