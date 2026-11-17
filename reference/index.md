# Package index

## Overview

- [`lagdynamics`](https://pak.dynasite.org/lagdynamics/reference/lagdynamics-package.md)
  [`lagdynamics-package`](https://pak.dynasite.org/lagdynamics/reference/lagdynamics-package.md)
  : lagdynamics: Lag Sequential Analysis, Dynamics, and Lag Transition
  Networks

## Fitting models

Lag sequential analysis of sequences or event logs at one lag or at
several, with a choice of engine.

- [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
  [`lsa_classical()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
  [`lsa_two_cell()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
  [`lsa_bidirectional()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
  [`lsa_parallel_dominance()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
  [`lsa_nonparallel_dominance()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
  : Fit a lag sequential analysis model
- [`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)
  [`print(`*`<lsa_lags>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)
  [`as.data.frame(`*`<lsa_lags>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)
  : Fit lag sequential analysis models at several lags

## Preparing data

Canonical sequence data, transition counts and expected counts under
structural zeros.

- [`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md)
  : Convert sequence input to the canonical form used by lag sequential
  analysis
- [`as.data.frame(`*`<lsa_data>`*`)`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_data.md)
  [`print(`*`<lsa_data>`*`)`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_data.md)
  : Read the events of canonical sequence data as a table
- [`lsa_transitions()`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md)
  [`print(`*`<lsa_transitions>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md)
  [`as.data.frame(`*`<lsa_transitions>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md)
  : Count transitions at a given lag
- [`lsa_ipf()`](https://pak.dynasite.org/lagdynamics/reference/lsa_ipf.md)
  : Fit expected counts with structural zeros by iterative proportional
  fitting

## Reading results

Tidy tables of a fitted model, one row per transition, state or test.

- [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
  : Read the transitions of a fit or an inference result
- [`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md) :
  Read the states of a lag sequential analysis fit
- [`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md) :
  Read the tablewise tests of independence of a fit
- [`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md)
  : Read the initial-state distribution of a fit
- [`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md)
  : Read the transition probabilities of a fit
- [`lag_profile()`](https://pak.dynasite.org/lagdynamics/reference/lag_profile.md)
  : Trace one transition across lags
- [`print(`*`<lsa_group>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
  [`print(`*`<lsa>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
  [`summary(`*`<lsa>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
  [`print(`*`<summary.lsa>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
  [`as.data.frame(`*`<summary.lsa>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
  [`as.data.frame(`*`<lsa>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
  [`as.data.frame(`*`<lsa_group>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
  [`summary(`*`<lsa_group>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
  [`nobs(`*`<lsa>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)
  : Print, summarise and tidy a lag sequential analysis fit

## Evaluating transitions

Uncertainty, stability, reliability and a permutation null for the
transitions of a fit.

- [`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md)
  : Estimate bootstrap confidence intervals for transitions
- [`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
  [`print(`*`<lsa_certainty>`*`)`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
  [`as.data.frame(`*`<lsa_certainty_group>`*`)`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
  [`print(`*`<lsa_certainty_group>`*`)`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
  : Compute the analytic certainty of transition probabilities
- [`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md)
  : Assess the case-drop stability of significant transitions
- [`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md)
  [`print(`*`<lsa_reliability_group>`*`)`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md)
  [`print(`*`<lsa_reliability>`*`)`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md)
  : Estimate the split-half reliability of a lag transition network
- [`as.data.frame(`*`<lsa_reliability>`*`)`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_reliability.md)
  [`as.data.frame(`*`<lsa_reliability_group>`*`)`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_reliability.md)
  : Read the tables of a split-half reliability result
- [`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)
  : Test transitions against no sequential structure by permutation
- [`print(`*`<lsa_bootstrap>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)
  [`print(`*`<lsa_permutation>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)
  [`print(`*`<lsa_stability>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)
  [`as.data.frame(`*`<lsa_bootstrap>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)
  [`as.data.frame(`*`<lsa_permutation>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)
  [`as.data.frame(`*`<lsa_stability>`*`)`](https://pak.dynasite.org/lagdynamics/reference/lsa-inference-methods.md)
  : Print and read resampling results

## Comparing groups

Transition-level differences between groups by permutation or by a
Bayesian Dirichlet-multinomial contrast.

- [`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
  : Compare transition structures between groups with a permutation test
- [`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
  [`print(`*`<lsa_bayes>`*`)`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
  [`print(`*`<lsa_bayes_pairwise>`*`)`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
  : Compare transition probabilities between groups with a Bayesian
  model
- [`as.data.frame(`*`<lsa_comparison>`*`)`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_comparison.md)
  [`as.data.frame(`*`<lsa_comparison_pairwise>`*`)`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_comparison.md)
  [`print(`*`<lsa_comparison>`*`)`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_comparison.md)
  [`print(`*`<lsa_comparison_pairwise>`*`)`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_comparison.md)
  : Read the tables of a group comparison

## Plots

Residual heatmaps, lag transition networks, chord diagrams, sunbursts,
uncertainty forests and comparison plots.

- [`plot(`*`<lsa>`*`)`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md)
  [`plot(`*`<lsa_group>`*`)`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md)
  : Plot a lag sequential analysis fit
- [`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md)
  : Plot a lag transition network
- [`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md)
  : Plot a chord diagram of transitions
- [`plot_polar()`](https://pak.dynasite.org/lagdynamics/reference/plot_polar.md)
  : Plot a polar sunburst of transitions
- [`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)
  [`plot(`*`<lsa_bootstrap>`*`)`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md)
  : Plot a circular forest of bootstrap confidence intervals
- [`plot(`*`<lsa_certainty>`*`)`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_certainty.md)
  : Plot credible intervals of transition probabilities
- [`plot(`*`<lsa_comparison>`*`)`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md)
  [`plot(`*`<lsa_comparison_pairwise>`*`)`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md)
  : Plot a comparison of transitions between groups

## Beyond LSA

Directed information transfer between states, and hand-off to transition
network analysis.

- [`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md)
  : Estimate directed transfer entropy between categorical sequences
- [`lsa_to_tna()`](https://pak.dynasite.org/lagdynamics/reference/lsa_to_tna.md)
  : Convert a lag sequential analysis fit to a transition network for
  tna

## Engines

The registry of test engines, open to user-defined engines.

- [`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md)
  : List the registered engines
- [`get_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/get_lsa_engine.md)
  : Retrieve a registered engine
- [`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md)
  : Register a lag sequential analysis engine
- [`unregister_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/unregister_lsa_engine.md)
  : Remove a registered engine

## Data

Engagement states, collaborative regulation and human-AI programming
sessions.

- [`engagement`](https://pak.dynasite.org/lagdynamics/reference/engagement.md)
  : Engagement states of students over time
- [`group_regulation`](https://pak.dynasite.org/lagdynamics/reference/group_regulation.md)
  : Regulation actions in collaborative learning, as sequences
- [`group_regulation_long`](https://pak.dynasite.org/lagdynamics/reference/group_regulation_long.md)
  : Regulation actions in collaborative learning, as an event log
- [`ai_long`](https://pak.dynasite.org/lagdynamics/reference/ai_long.md)
  : Coded AI behaviours in human-AI programming sessions
