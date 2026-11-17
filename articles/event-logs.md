# Event logs and interoperability

Lag sequential analysis counts transitions between consecutive events of
the same sequence. What counts as a sequence therefore decides what
counts as a transition: two actions by different students are not a
transition, and neither are the last action of one working session and
the first action of the next one a day later. With data that arrive as
an event log, one row per recorded action, the analysis starts by
recovering the sequences: which rows belong together, in what order, and
where a sequence ends.

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
performs that step itself. It accepts wide sequence data, a long event
log described by a small set of column arguments, and the sequence
objects of the `tna` and `Nestimate` packages, and it returns a fit
whose transition probabilities, initial probabilities, and network can
be handed back to those packages. This guide covers each input form, the
rules that split a log into sequences, grouped logs, and the export of a
fit. The column grammar is the one used by
[`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html)
and
[`Nestimate::build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.html),
so a log prepared for one of them needs no change.

``` r

library(lagdynamics)
```

## Data

Three bundled data sets cover the input forms.

- `group_regulation_long`: an event log of collaborative learning,
  27,533 rows. **Actor:** `Actor` (2,000 students). **Action:**
  `Action`, nine regulation codes. **Time:** `Time`, a POSIXct time
  stamp. **Group:** `Achiever` (`High` or `Low`), fixed for each actor.
- `ai_long`: an event log of coded AI behaviours in human-AI programming
  sessions, 8,551 rows. **Actor:** `project` (34 projects). **Session:**
  `session_id` (428 sessions). **Action:** `code`, eight behaviour
  codes. **Order:** `order_in_session`, the position of the coded
  message within its session.
- `group_regulation` and `engagement`: wide data, one row per sequence
  and one column per time point.

``` r

head(group_regulation_long)
#>   Actor Achiever Group Course                Time    Action
#> 1     1     High     1      A 2025-01-01 08:27:07  cohesion
#> 2     1     High     1      A 2025-01-01 08:35:20 consensus
#> 3     1     High     1      A 2025-01-01 08:42:18   discuss
#> 4     1     High     1      A 2025-01-01 08:50:00 synthesis
#> 5     1     High     1      A 2025-01-01 08:52:25     adapt
#> 6     1     High     1      A 2025-01-01 08:57:31 consensus
head(ai_long)
#>   message_id   project   session_id  timestamp session_date        code cluster
#> 1       3441 Project_7 0086cabebd15 1772661600   2026-03-05    Delegate  Action
#> 2       3441 Project_7 0086cabebd15 1772661600   2026-03-05        Plan  Repair
#> 3       3443 Project_7 0086cabebd15 1772661600   2026-03-05     Execute  Action
#> 4       3443 Project_7 0086cabebd15 1772661600   2026-03-05        Plan  Repair
#> 5       3445 Project_7 0086cabebd15 1772661600   2026-03-05     Execute  Action
#> 6       3447 Project_7 0086cabebd15 1772661600   2026-03-05 Investigate  Action
#>   code_order order_in_session
#> 1          1                5
#> 2          2                6
#> 3          1                8
#> 4          2                9
#> 5          1               12
#> 6          1               14
```

## Wide sequence data

Wide data need no column arguments: each row is a sequence and the
columns are read in order. `group_regulation` holds the same 2,000
sequences as the long log above, already assembled.

``` r

fit_wide <- lsa(group_regulation)
as.data.frame(summary(fit_wide))
#>      engine lag n_states n_sequences n_events n_transitions n_significant alpha  lrx2
#> 1 classical   1        9        2000    27533         25533            72  0.05 13204
#>   lrx2_df lrx2_p    x2 x2_df x2_p
#> 1      64      0 15435    64    0
```

The wide data give 2,000 sequences with 25,533 transitions among nine
states. Sequences differ in length, and the missing cells that pad the
shorter rows are dropped, so no transition is counted into or out of
them.

## Long event logs

A long log is described by naming its columns, as strings:

| Argument | Role |
|----|----|
| `action` | The state of each event. Required; naming it switches to long-format input. |
| `actor` | Who produced the event. Each actor’s events form a sequence. |
| `session` | A working session. Each session, or each actor-session pair when both are given, forms a sequence. |
| `time` | A time stamp that orders the events (POSIXct, a parseable string, or Unix time with `is_unix_time = TRUE`). |
| `order` | An explicit position that orders the events; used in preference to `time`. |
| `group` | A grouping column, fixed within each sequence, for one model per group. |

The regulation log has an actor, an action, and a time stamp.

``` r

fit_long <- lsa(group_regulation_long, actor = "Actor", action = "Action",
                time = "Time")
as.data.frame(summary(fit_long))
#>      engine lag n_states n_sequences n_events n_transitions n_significant alpha  lrx2
#> 1 classical   1        9        2000    27533         25533            72  0.05 13204
#>   lrx2_df lrx2_p    x2 x2_df x2_p
#> 1      64      0 15435    64    0
```

The log is assembled into 2,000 sequences, one per actor, with the same
25,533 transitions as the wide data. The two fits are the same model:

``` r

isTRUE(all.equal(transitions(fit_long), transitions(fit_wide),
                 check.attributes = FALSE))
#> [1] TRUE
```

## Where a sequence ends

When `time` is given without `session`, a sequence also ends wherever
two consecutive events of an actor are more than `time_threshold`
seconds apart, 900 seconds (15 minutes) by default. The threshold is an
analytic choice: it states how long a pause has to be before the next
action no longer counts as a response to the previous one. A shorter
threshold splits the same log into more, shorter sequences.

``` r

fit_5min <- lsa(group_regulation_long, actor = "Actor", action = "Action",
                time = "Time", time_threshold = 300)
as.data.frame(summary(fit_5min))
#>      engine lag n_states n_sequences n_events n_transitions n_significant alpha lrx2
#> 1 classical   1        9       12305    27533         15228            69  0.05 7927
#>   lrx2_df lrx2_p   x2 x2_df x2_p
#> 1      64      0 9222    64    0
```

With a five-minute threshold the 2,000 actor logs break into 12,305
sequences, and the number of transitions falls from 25,533 to 15,228:
every pause longer than five minutes removes the transition that spanned
it. The independence test remains decisive (G² = 7,927 on 64 degrees of
freedom) and 69 transitions stay significant, against 72 with the
default, so the structure of the process does not hinge on the
threshold. Comparing the fit at two thresholds in this way is the
sensitivity check the choice needs.

Without `actor` or `session`, the whole log is one stream ordered by
`time`. The regulation log interleaves 2,000 students working at the
same time, so a single stream mixes their actions into transitions that
never happened:

``` r

fit_stream <- lsa(group_regulation_long, action = "Action", time = "Time")
as.data.frame(summary(fit_stream))
#>      engine lag n_states n_sequences n_events n_transitions n_significant alpha lrx2
#> 1 classical   1        9           1    27533         27532             1  0.05 46.2
#>   lrx2_df lrx2_p x2 x2_df  x2_p
#> 1      64  0.954 46    64 0.957
```

The stream gives one sequence of 27,532 transitions, and its
independence test is not significant (G² = 46.2, df = 64, p = 0.954):
interleaving the students destroys the sequential structure. The `actor`
argument is what makes the regulation log analysable.

## Sessions and explicit order

The AI log records sessions within projects and an explicit position of
each coded message within its session. Naming `actor` and `session`
makes each project-session pair a sequence, and `order` sets the order
of events within it. Several codes can share a message, and the position
column keeps them in the order they were coded.

``` r

fit_ai <- lsa(ai_long, actor = "project", session = "session_id",
              action = "code", order = "order_in_session")
fit_ai
#> Lag Sequential Analysis  -  classical  (lag 1, directed)
#>   8 states | 8123 transitions | 8551 events | 428 sequences
#>   states: Ask, Delegate, Execute, Explain, Investigate, Plan, Repair, Report
#>   independence: G² = 2168.1, df = 49, p <2e-16
#> 
#>   Significant transitions (p < 0.05): 40 of 64
#>   strongest over-represented (of 19):
#>     Investigate -> Plan  z =  +33.4  ***
#>     Delegate -> Plan     z =  +18.1  ***
#>     Ask -> Explain       z =  +16.4  ***
#>     Execute -> Execute   z =  +16.1  ***
#>     Explain -> Report    z =  +14.0  ***
#>     ... and 14 more
#> 
#>   Initial states:
#>     Investigate 0.715  ████████████████████████
#>     Delegate    0.185  ██████
#>     Execute     0.058  ██
#>     Plan        0.035  █
#>     Ask         0.005  
#>     Explain     0.002  
#>     Repair      0.000  
#>     Report      0.000
```

The log yields 428 sequences, one per session, with 8,123 transitions
among eight codes. Forty of the 64 transitions depart from chance, and
the strongest favoured transition is `Investigate -> Plan` (z = 33.4).
Most sessions open with `Investigate` (initial probability 0.715).

## Large logs: the actor-session key

When both `actor` and `session` are given, the sequence key is built
over the full cross of distinct actors and distinct sessions. The size
of that cross is held in an R integer, so once the number of distinct
actors times the number of distinct sessions passes 2,147,483,647, the
largest integer R can hold, the key overflows. The affected events lose
their sequence and are dropped, and the only sign is an “NAs produced by
integer overflow” warning. A log of 50,000 actors and 50,000 sessions is
enough to trigger it.

On a large log, compare the number of sequences in the summary with the
number of distinct sessions. When the session identifiers are already
unique across actors, as in `ai_long`, pass the session column as
`actor` and omit `session`; the key is then built over the sessions
alone.

``` r

fit_ai_sessions <- lsa(ai_long, actor = "session_id", action = "code",
                       order = "order_in_session")
as.data.frame(summary(fit_ai_sessions))
#>      engine lag n_states n_sequences n_events n_transitions n_significant alpha lrx2
#> 1 classical   1        8         428     8551          8123            40  0.05 2168
#>   lrx2_df lrx2_p   x2 x2_df x2_p
#> 1      49      0 2485    49    0
isTRUE(all.equal(transitions(fit_ai_sessions), transitions(fit_ai),
                 check.attributes = FALSE))
#> [1] TRUE
```

Keying on the session alone recovers the same 428 sequences and an
identical transition table, so the two calls are interchangeable here,
and only the second one is safe on a log too large for the cross.

## Grouped long logs

Adding `group` fits one model per group from the same log, under one
shared set of states. The grouping column must take a single value
within each recovered sequence; a sequence that spans two groups is an
error.

``` r

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
#>   group    engine lag n_states n_sequences n_events n_transitions n_significant alpha
#> 1  High classical   1        9        1000    13721         12721            66  0.05
#> 2   Low classical   1        9        1000    13812         12812            67  0.05
#>   lrx2 lrx2_df lrx2_p   x2 x2_df x2_p
#> 1 6147      64      0 6807    64    0
#> 2 7715      64      0 9197    64    0
```

The high and the low achievers contribute 1,000 sequences each, with
12,721 and 12,812 transitions. Every reading verb on a grouped fit
returns a leading `group` column, and the grouped fit is the input of
the group comparisons covered in [the group-comparison
vignette](https://pak.dynasite.org/lagdynamics/articles/group-comparison.md).

``` r

transitions(gfit, significant = TRUE, sort = "strength") |> head(4)
#>   group      from        to lag count expected   prob prob_col adj_res         p yules_q
#> 1  High   emotion  cohesion   1   505    114.0 0.3258  0.53953    40.6  0.00e+00   0.847
#> 2  High consensus consensus   1   283    917.0 0.0834  0.08229   -28.6 3.68e-180  -0.698
#> 3  High   discuss synthesis   1   213     46.1 0.1063  0.72696    27.1 1.79e-161   0.881
#> 4  High   discuss      plan   1    25    458.8 0.0125  0.00858   -25.1 2.35e-139  -0.934
#>    kappa kappa_z   kappa_p   lift  sign significant
#> 1  0.425    37.7 6.36e-311 4.4280  over        TRUE
#> 2 -0.709   -30.1 5.36e-199 0.3086 under        TRUE
#> 3  0.613    25.5 5.39e-143 4.6169  over        TRUE
#> 4 -0.949   -26.0 3.81e-149 0.0545 under        TRUE
```

## tna objects

A `tna` model keeps the sequences it was fitted from, and
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) reads
them directly. The model below is fitted by `tna` to the wide
`engagement` data and passed to
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

``` r

tna_fit <- tna::tna(engagement)
tna_fit
#> State Labels : 
#> 
#>    Active, Average, Disengaged 
#> 
#> Transition Probability Matrix :
#> 
#>            Active Average Disengaged
#> Active      0.698   0.267      0.035
#> Average     0.204   0.610      0.186
#> Disengaged  0.120   0.397      0.483
#> 
#> Initial Probabilities : 
#> 
#>     Active    Average Disengaged 
#>      0.382      0.368      0.250
fit_from_tna <- lsa(tna_fit)
isTRUE(all.equal(transitions(fit_from_tna), transitions(lsa(engagement)),
                 check.attributes = FALSE))
#> [1] TRUE
```

Loading `tna` prints a notice that it registers its own plot method for
bootstrap results in place of the one from `cograph`; the notice does
not concern the models. The fit from the `tna` model and the fit from
the original data are identical. The same holds for the output of
[`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html),
which [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
accepts as it stands. A grouped `tna` model (`group_tna`) is read as
well, but its groups are pooled into a single fit; for one model per
group, fit the log with `group =` as in the previous section.

## Nestimate objects

[`Nestimate::build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.html)
uses the same long-log grammar and keeps the prepared sequences in the
network object, which
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) reads.

``` r

net <- Nestimate::build_network(ai_long, method = "tna", actor = "project",
                                session = "session_id", action = "code",
                                order = "order_in_session")
#> Metadata aggregated per session: ties resolved by first occurrence in 'session_date' (4 sessions), 'cluster' (22 sessions)
fit_from_nestimate <- lsa(net)
isTRUE(all.equal(transitions(fit_from_nestimate), transitions(fit_ai),
                 check.attributes = FALSE))
#> [1] TRUE
```

`Nestimate` reports that two metadata columns (`session_date` and
`cluster`) vary within some sessions and keeps their first value; the
message concerns those columns only, not the sequences. The `Nestimate`
network and the direct
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) call on
the log recover the same transition table.

## From a fit to tna and other tools

A fit exposes the quantities that transition network tools consume.
[`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md)
returns the row-stochastic transition matrix as a tidy data frame, one
row per `from -> to` pair, and
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
the distribution of first states.

``` r

transition_probabilities(fit_long) |> head(9)
#>    from         to   prob
#> 1 adapt      adapt 0.0000
#> 2 adapt   cohesion 0.2731
#> 3 adapt  consensus 0.4774
#> 4 adapt coregulate 0.0216
#> 5 adapt    discuss 0.0589
#> 6 adapt    emotion 0.1198
#> 7 adapt    monitor 0.0334
#> 8 adapt       plan 0.0157
#> 9 adapt  synthesis 0.0000
initial(fit_long)
#>        state init_prob
#> 1      adapt    0.0115
#> 2   cohesion    0.0605
#> 3  consensus    0.2140
#> 4 coregulate    0.0190
#> 5    discuss    0.1755
#> 6    emotion    0.1515
#> 7    monitor    0.1440
#> 8       plan    0.2045
#> 9  synthesis    0.0195
```

The nine rows shown are the transitions out of `adapt`; their
probabilities sum to one. The most common first states are `consensus`
(0.2140) and `plan` (0.2045). For a grouped fit both tables carry a
leading `group` column.

[`lsa_to_tna()`](https://pak.dynasite.org/lagdynamics/reference/lsa_to_tna.md)
converts a fit to a `tna` model, so the `tna` toolkit (centralities,
communities, cliques, its plots) applies to a model estimated here.
`weights` chooses the edge weights; `"prob"` gives the transition
probabilities of a standard transition network.

``` r

tn <- lsa_to_tna(fit_from_tna, weights = "prob")
tna::centralities(tn, measures = c("InStrength", "OutStrength", "Betweenness"))
#> # A tibble: 3 × 4
#>   state      InStrength OutStrength Betweenness
#> * <fct>           <dbl>       <dbl>       <dbl>
#> 1 Active          0.324       0.302           0
#> 2 Average         0.664       0.390           2
#> 3 Disengaged      0.221       0.517           0
```

`Average` receives the most probability from the other states
(in-strength 0.664) and is the only state that lies on shortest paths
between the other two (betweenness 2, against 0 for `Active` and
`Disengaged`). A grouped fit converts to a `group_tna` model, one `tna`
model per group:

``` r

tna::centralities(lsa_to_tna(gfit), measures = "InStrength")
#> # A tibble: 18 × 3
#>    group state      InStrength
#>  * <chr> <fct>           <dbl>
#>  1 High  adapt           0.215
#>  2 High  cohesion        0.815
#>  3 High  consensus       2.95 
#>  4 High  coregulate      0.436
#>  5 High  discuss         1.12 
#>  6 High  emotion         1.00 
#>  7 High  monitor         0.301
#>  8 High  plan            1.27 
#>  9 High  synthesis       0.167
#> 10 Low   adapt           0.452
#> 11 Low   cohesion        0.799
#> 12 Low   consensus       2.42 
#> 13 Low   coregulate      0.689
#> 14 Low   discuss         1.21 
#> 15 Low   emotion         0.808
#> 16 Low   monitor         0.393
#> 17 Low   plan            1.15 
#> 18 Low   synthesis       0.216
```

`consensus` is the state with the most incoming probability in both
groups, with in-strength 2.95 among the high achievers and 2.42 among
the low achievers.

A fit is also a `cograph_network`, so it can be drawn directly as a
transition network with probability weights; the displays are covered in
[the plotting
guide](https://pak.dynasite.org/lagdynamics/articles/plotting.md).

## Summary

| Input | Call |
|----|----|
| Wide data, one row per sequence | `lsa(data)` |
| Long log, one sequence per actor | `lsa(log, actor = , action = , time = )` |
| Long log with sessions and explicit order | `lsa(log, actor = , session = , action = , order = )` |
| Very large log with globally unique sessions | `lsa(log, actor = "<session column>", action = , order = )` |
| One model per group | add `group =` |
| `tna` model, [`tna::prepare_data()`](https://sonsoles.me/tna/reference/prepare_data.html) output | `lsa(object)` |
| `Nestimate` network or prepared data | `lsa(object)` |
| Hand a fit to `tna` | `lsa_to_tna(fit)` |
| Transition and initial probabilities as tables | `transition_probabilities(fit)`, `initial(fit)` |
