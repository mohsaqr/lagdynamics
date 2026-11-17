# Transfer entropy between states

Transfer entropy (Schreiber, 2000) measures directed predictive
coupling: how much knowing the past of a source series $`X`$ reduces the
uncertainty about the future of a target series $`Y`$, over and above
what the target’s own past already tells. With $`Y_{t+\ell}`$ the target
$`\ell`$ steps ahead, $`Y_t^{(k)}`$ the target’s last $`k`$ values and
$`X_t`$ the source’s current value,

``` math
T_{X \to Y} = H\big(Y_{t+\ell} \mid Y_t^{(k)}\big) -
  H\big(Y_{t+\ell} \mid Y_t^{(k)}, X_t\big),
```

measured in bits. It is zero when the source adds nothing to the
target’s own history, and it is not symmetric: $`T_{X \to Y}`$ and
$`T_{Y \to X}`$ are separate quantities. Because it conditions on the
target’s past, a process with strong momentum does not inflate it the
way it inflates a plain lagged association.

Transfer entropy answers a different question from lag sequential
analysis. An adjusted residual says whether the transition $`i \to j`$
occurs more or less often than base rates predict, and its sign says
which. Transfer entropy says how much information one series carries
about the next step of another, and it has no sign: a source that makes
the target’s next state less likely is as informative as one that makes
it more likely. The two are read together, as this guide shows.
[`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md)
is marked experimental in the package.

``` r

library(lagdynamics)
```

## Data

- **`engagement`**: weekly engagement states (`Active`, `Average`,
  `Disengaged`) of 138 students over 15 weeks, one row per student; two
  rows are entirely missing.
- **`group_regulation`**: coded regulation actions in 2,000
  collaborative learning sessions, one row per session and one column
  per ordered time point (`T1` to `T26`), nine states.
- **`ai_long`**: a long log of 8,551 coded behaviours of an AI coding
  assistant in 428 “vibe-coding” sessions from 34 projects, in which a
  person builds software by instructing the assistant in natural
  language. Each row is one coded behaviour, ordered within its session
  by `order_in_session`. The column `code` holds one of eight behaviour
  codes (`Ask`, `Delegate`, `Execute`, `Explain`, `Investigate`, `Plan`,
  `Repair`, `Report`) and the column `cluster` the broader category the
  code belongs to (`Action`, `Communication`, `Repair`).

## The estimator

[`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md)
computes plug-in Shannon entropies from pooled (target future, target
history, source value) triples. Each sequence contributes the triples
that fit inside it, so the lag never crosses from one sequence into the
next, and triples with a missing value are dropped. The returned table
has one row per ordered pair with the columns:

- `te`: the transfer entropy in bits.
- `te_normalised`: `te` divided by the target’s leftover uncertainty
  $`H(Y_{t+\ell} \mid Y_t^{(k)})`$, so it is the share of what the
  target’s own past leaves unexplained that the source explains (0 to
  1).
- `te_effective`: `te` minus the mean transfer entropy of the
  surrogates, a bias-corrected value (surrogate test only).
- `p`: the surrogate p-value (surrogate test only).
- `n`: the number of pooled triples.

Rows are ordered by decreasing `te`. The arguments `lag` ($`\ell`$,
default 1) and `history` ($`k`$, default 1) set the prediction horizon
and the length of the target’s own past.

## A state network from one sequence dataset

Given one set of categorical sequences,
[`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md)
returns a directed network among the states. Each ordered pair of states
$`s \to r`$ is turned into two binary series: the source is “the unit is
in $`s`$” and the target is “the unit is in $`r`$”. The transfer entropy
of the pair is then the information that being in $`s`$ now gives about
being in $`r`$ next, beyond whether the unit is in $`r`$ now. Self-pairs
are not computed, so $`K`$ states give $`K(K-1)`$ rows. The call below
skips the significance test.

``` r

transfer_entropy(engagement, test = "none")
#>         from         to      te te_normalised    n
#> 1     Active Disengaged 0.03645       0.06002 1734
#> 2    Average Disengaged 0.03645       0.06002 1734
#> 3     Active    Average 0.00695       0.00752 1734
#> 4 Disengaged    Average 0.00695       0.00752 1734
#> 5    Average     Active 0.00479       0.00635 1734
#> 6 Disengaged     Active 0.00479       0.00635 1734
```

The largest flows end in `Disengaged`: knowing whether a student is
`Active` this week gives 0.036 bits about whether the student is
`Disengaged` next week, 6.0% of the uncertainty that the student’s own
disengagement status leaves. Flows into `Average` and `Active` are five
to eight times smaller (0.007 and 0.005 bits). The rows come in
identical pairs. This is a property of the binary decomposition with
three states: once the target’s own indicator says the student is not
`Disengaged`, the student is either `Active` or `Average`, so the two
source indicators carry the same information. The documentation of
[`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md)
notes this redundancy; with three states the network has three distinct
values, one per target, and a larger alphabet or the bivariate mode
below is needed for source-specific flows.

## Testing against surrogates

A plug-in estimate of transfer entropy is positive even when the source
and the target are unrelated, because sampling noise always produces
some apparent information. The surrogate test (`test = "surrogate"`, the
default) shuffles the source values across the pooled triples `R` times,
which breaks the coupling between source and target while keeping both
marginal distributions, and recomputes the transfer entropy each time.
The p-value is
$`(1 + \#\{T_\text{surrogate} \ge T_\text{observed}\}) / (R + 1)`$, and
the effective transfer entropy subtracts the surrogate mean. `seed`
makes the shuffles reproducible.

``` r

transfer_entropy(engagement, R = 199, seed = 2026)
#>         from         to      te te_effective te_normalised     p    n
#> 1     Active Disengaged 0.03645      0.03569       0.06002 0.005 1734
#> 2    Average Disengaged 0.03645      0.03559       0.06002 0.005 1734
#> 3     Active    Average 0.00695      0.00607       0.00752 0.005 1734
#> 4 Disengaged    Average 0.00695      0.00613       0.00752 0.005 1734
#> 5    Average     Active 0.00479      0.00396       0.00635 0.010 1734
#> 6 Disengaged     Active 0.00479      0.00397       0.00635 0.010 1734
```

The four flows into `Disengaged` and `Average` exceed all 199 surrogates
($`p = 0.005`$, the smallest attainable value), and the two flows into
`Active` exceed all but one ($`p = 0.010`$). The bias correction is
small: the effective transfer entropy of `Active -> Disengaged` is
0.0357 bits against an estimate of 0.0365, because 1,734 triples leave
little room for noise. The directed flows are real, and the flow towards
disengagement is the one of practical size.

## Longer history and longer horizon

The `history` argument conditions on more of the target’s own past.

``` r

transfer_entropy(engagement, history = 2, test = "none")
#>         from         to       te te_normalised    n
#> 1     Active Disengaged 0.033837       0.06105 1603
#> 2    Average Disengaged 0.033837       0.06105 1603
#> 3    Average     Active 0.001273       0.00180 1603
#> 4 Disengaged     Active 0.001273       0.00180 1603
#> 5     Active    Average 0.000938       0.00106 1603
#> 6 Disengaged    Average 0.000938       0.00106 1603
```

With two weeks of the target’s own status known, the flows into `Active`
and `Average` shrink to between a quarter and an eighth of their size
(from 0.0048 to 0.0013 bits into `Active`, from 0.0070 to 0.0009 into
`Average`), while the flow into `Disengaged` barely moves (0.0338 bits).
Most of what the current state said about the next active or average
week was already contained in the student’s own recent history; what it
says about the next disengaged week is not. The `lag` argument moves the
prediction horizon in the same way as in
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md);
`lag = 2` asks about the state two weeks ahead.

## A larger alphabet

With nine states the binary channels are no longer redundant and each
ordered pair carries its own value. The call below lists the strongest
flows of the `group_regulation` network.

``` r

regulation_te <- transfer_entropy(group_regulation, test = "none")
head(regulation_te, 8)
#>        from         to     te te_normalised     n
#> 1 consensus       plan 0.0586        0.0751 25533
#> 2   emotion   cohesion 0.0584        0.1648 25533
#> 3   discuss       plan 0.0458        0.0587 25533
#> 4   discuss  synthesis 0.0419        0.2352 25533
#> 5 consensus coregulate 0.0275        0.0677 25533
#> 6      plan coregulate 0.0190        0.0467 25533
#> 7 synthesis      adapt 0.0150        0.1035 25533
#> 8 consensus   cohesion 0.0149        0.0421 25533
```

The strongest flows run from `consensus` to `plan` (0.059 bits) and from
`emotion` to `cohesion` (0.058 bits). Among these eight flows,
`discuss -> synthesis` has the largest normalised value: it explains
23.5% of the leftover uncertainty about whether the next action is a
synthesis, against 7.5% for `consensus -> plan`. The surrogate test is
omitted here: its cost grows with `R` times the number of ordered pairs,
which is 72 for nine states.

The lag sequential analysis of the same data gives the signs.

``` r

transitions(lsa(group_regulation), significant = TRUE,
            sort = "strength") |> head(6)
#>        from         to lag count expected   prob prob_col adj_res         p yules_q
#> 1   emotion   cohesion   1   923    190.9 0.3253   0.5373    58.2  0.00e+00   0.860
#> 2   discuss  synthesis   1   557    106.8 0.1410   0.8072    48.0  0.00e+00   0.927
#> 3 synthesis      adapt   1   153     13.6 0.2347   0.2881    38.8  0.00e+00   0.904
#> 4   discuss       plan   1    46    961.6 0.0116   0.0074   -36.9 2.20e-298  -0.943
#> 5 consensus  consensus   1   519   1578.7 0.0820   0.0815   -35.5 5.32e-276  -0.661
#> 6 consensus coregulate   1  1188    519.3 0.1877   0.5671    35.3 3.36e-273   0.647
#>    kappa kappa_z   kappa_p    lift  sign significant
#> 1  0.439    55.0  0.00e+00  4.8353  over        TRUE
#> 2  0.721    46.1  0.00e+00  5.2168  over        TRUE
#> 3  0.256    37.0 2.44e-299 11.2837  over        TRUE
#> 4 -0.955   -38.2 2.88e-319  0.0478 under        TRUE
#> 5 -0.691   -37.6 7.65e-309  0.3287 under        TRUE
#> 6  0.412    34.6 5.06e-262  2.2877  over        TRUE
```

`emotion -> cohesion` is a strong flow in both analyses and
over-represented ($`z = 58.2`$): emotional expressions are followed by
cohesion. `discuss -> plan` is the third largest flow of information
(0.046 bits) but the most avoided transition of the table
($`z = -36.9`$, 46 observed against 962 expected): a discussion predicts
that the next action is *not* planning. Transfer entropy ranks the two
among the most informative relations of the network; only the residual
says that one is an attraction and the other an avoidance. Read transfer
entropy for the strength of directed dependence and the adjusted
residual for its direction.

## Two aligned series

When `y` is supplied,
[`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md)
computes the full-alphabet transfer entropy between two aligned series,
in both directions. The two series must have the same shape: the same
number of sequences and the same length within each. In `ai_long` every
coded behaviour carries both its `code` and its `cluster`, so the two
columns are aligned by construction. The log is ordered by session and
by position within session, and the columns are passed as one long
series each.

``` r

with(ai_long, transfer_entropy(code, cluster, R = 199, seed = 2026))
#>      from      to     te te_effective te_normalised     p    n
#> 1    code cluster 0.0715      0.06775        0.0642 0.005 8550
#> 2 cluster    code 0.0000     -0.00981        0.0000 1.000 8550
```

The coupling is one-directional. The code of the current behaviour
carries 0.0715 bits about the cluster of the next one beyond what the
current cluster says ($`p = 0.005`$). The reverse flow is exactly 0: the
cluster is a function of the code, so the current cluster cannot tell
anything about the next code that the current code does not already
tell. The negative effective value ($`-0.0098`$) is the bias correction
applied to a true zero. Passed as one long series, the 8,550 pairs
include the 427 that join the last behaviour of one session to the first
of the next, about 5% of the total; a series per session (one row per
session, `NA`-padded) avoids them.

## Transfer entropy and lag sequential analysis

|  | Lag sequential analysis | Transfer entropy |
|----|----|----|
| Question | Does $`i \to j`$ occur more or less than base rates predict? | How much does the source’s present tell about the target’s future beyond the target’s own past? |
| Sign | Signed: over- or under-represented | Unsigned: attraction and avoidance both count |
| Conditioning | On the margins of the transition table | On the target’s own history (`history`) |
| Unit | $`z`$, odds ratio, probability | Bits, or a share of leftover uncertainty |
| Test | Normal reference, bootstrap, permutation | Source-shuffling surrogates |
| Function | [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md) | [`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md) |

Use lag sequential analysis to name the transitions that occur more or
less than expected and to test them. Add transfer entropy when the
question is which states carry predictive information about others once
each state’s own persistence is accounted for, and when a process has
strong momentum, as the `engagement` data do. Report transfer entropy
with the surrogate p-value and next to the signed residual, never on its
own.

## References

Schreiber, T. (2000). Measuring information transfer. *Physical Review
Letters, 85*(2), 461-464.
