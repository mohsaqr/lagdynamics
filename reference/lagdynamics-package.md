# lagdynamics: Lag Sequential Analysis, Dynamics, and Lag Transition Networks

`lagdynamics` fits lag sequential analysis models, reads them as lag
transition networks, and assesses, compares and plots them. Lag
sequential analysis (LSA; Sackett, 1979; Bakeman & Gottman, 1997) tests
whether, in a stream of coded events, one behaviour follows another more
or less often than chance would predict.

The observed transitions at a given lag are counted into a contingency
table of antecedent by consequent states, an expected table is derived
under the hypothesis that the two are independent, and each cell is
tested with an adjusted residual (Haberman, 1979). The cells that depart
from independence form a lag transition network: a directed graph whose
edges are the over- and under-represented sequential contingencies.

The model is fitted with one constructor,
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), which
accepts event sequences, wide sequence tables, transition-count
matrices, long event logs and sequence objects from tna, Nestimate and
TraMineR. Results are read with verbs that return tidy data frames,
assessed with resampling and analytic inference, compared across groups,
and drawn as heatmaps, networks, chord diagrams and polar displays.

The documentation site, with tutorials and case studies, is at
<https://pak.dynasite.org/lagdynamics/>.

## Workflow

1.  **Fit** a model with
    [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md):
    choose the lag, the engine (the classical test of Bakeman & Quera,
    1995, or the two-cell, bidirectional and dominance tests of Sackett,
    1979, and Wampold, 1984), and any structurally forbidden
    transitions. A `group` argument fits one model per group on a shared
    state set.

2.  **Read** the fit with
    [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
    [`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
    [`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md),
    [`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md),
    [`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md)
    and [`summary()`](https://rdrr.io/r/base/summary.html).

3.  **Evaluate** how far the transitions can be trusted:
    [`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
    [`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md),
    [`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
    [`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md)
    and
    [`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md).

4.  **Compare** groups with
    [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
    (permutation test) or
    [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
    (Dirichlet-multinomial posterior).

5.  **Plot** the fit or a comparison:
    [`plot.lsa()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md),
    [`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md),
    [`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md),
    [`plot_polar()`](https://pak.dynasite.org/lagdynamics/reference/plot_polar.md),
    [`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)
    and
    [`plot.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md).

6.  **Go beyond LSA** with multi-lag profiles
    ([`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md),
    [`lag_profile()`](https://pak.dynasite.org/lagdynamics/reference/lag_profile.md)),
    directed transfer entropy
    ([`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md))
    and Transition Network Analysis through
    [`lsa_to_tna()`](https://pak.dynasite.org/lagdynamics/reference/lsa_to_tna.md).

## Functions by task

- Fitting:
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) and
  its engine wrappers
  [`lsa_classical()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
  [`lsa_two_cell()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
  [`lsa_bidirectional()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
  [`lsa_parallel_dominance()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
  [`lsa_nonparallel_dominance()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md);
  [`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md).

- Input:
  [`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md)
  (canonical sequences),
  [`lsa_transitions()`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md)
  (transition counts at a lag),
  [`lsa_ipf()`](https://pak.dynasite.org/lagdynamics/reference/lsa_ipf.md)
  (expected counts with structural zeros).

- Reading:
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
  [`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
  [`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md),
  [`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md),
  [`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md),
  [`lag_profile()`](https://pak.dynasite.org/lagdynamics/reference/lag_profile.md),
  and the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  in
  [lsa-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md).

- Inference:
  [`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
  [`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md),
  [`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
  [`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md),
  [`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md);
  every result is read with
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
  or [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- Group comparison:
  [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md),
  [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md).

- Plotting:
  [`plot.lsa()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md),
  [`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md),
  [`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md),
  [`plot_polar()`](https://pak.dynasite.org/lagdynamics/reference/plot_polar.md),
  [`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md),
  [`plot.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md).

- Engines:
  [`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md),
  [`get_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/get_lsa_engine.md),
  [`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md),
  [`unregister_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/unregister_lsa_engine.md).

- Other:
  [`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md),
  [`lsa_to_tna()`](https://pak.dynasite.org/lagdynamics/reference/lsa_to_tna.md).

## Data

[engagement](https://pak.dynasite.org/lagdynamics/reference/engagement.md)
(engagement states of 138 students over 15 time points),
[group_regulation](https://pak.dynasite.org/lagdynamics/reference/group_regulation.md)
and
[group_regulation_long](https://pak.dynasite.org/lagdynamics/reference/group_regulation_long.md)
(collaborative regulation actions, as a wide sequence table and as a
long event log), and
[ai_long](https://pak.dynasite.org/lagdynamics/reference/ai_long.md)
(coded AI behaviours in human-AI programming sessions).

## Formulas and sources

Each formula the package implements is listed with its primary source,
page and equation number in the file `REFERENCES.md` installed with the
package; its location is given by
`system.file("REFERENCES.md", package = "lagdynamics")`.

Fits are objects of class `c("lsa", "cograph_network")`, so the network
functions of cograph accept them directly.

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

Wickens, T. D. (1989). *Multiway contingency tables analysis for the
social sciences*. Lawrence Erlbaum.

## See also

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) to
start; <https://pak.dynasite.org/lagdynamics/> for the tutorials.

## Author

**Maintainer**: Mohammed Saqr <saqr@saqr.me>
([ORCID](https://orcid.org/0000-0001-5881-3109)) \[copyright holder\]

Authors:

- Mohammed Saqr <saqr@saqr.me>
  ([ORCID](https://orcid.org/0000-0001-5881-3109)) \[copyright holder\]

- Sonsoles López-Pernas ([ORCID](https://orcid.org/0000-0002-9621-1392))
