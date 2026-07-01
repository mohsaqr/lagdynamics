# Package index

## Fit a model

Construct a lag-sequential model from sequences, an event log, or a
transition matrix.
[`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md) also
documents the convenience wrappers
([`lsa_classical()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md),
[`lsa_two_cell()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md),
…).

- [`lsa_classical()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)
  [`lsa_two_cell()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)
  [`lsa_bidirectional()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)
  [`lsa_parallel_dominance()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)
  [`lsa_nonparallel_dominance()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md) :
  Lag Sequential Analysis
- [`lsa_data()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_data.md)
  : Canonicalize Sequence Input for Lag Sequential Analysis
- [`lsa_transitions()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_transitions.md)
  : Tidy Transition Counts at a Given Lag
- [`lsa_ipf()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_ipf.md)
  : Iterative Proportional Fitting for Two-Way Tables with Structural
  Zeros

## Multi-lag analysis

- [`lsa_lags()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_lags.md)
  : Lag Sequential Analysis Across Several Lags
- [`lag_profile()`](https://mohsaqr.github.io/lagdynamics/reference/lag_profile.md)
  : Lag Profile of a Single Transition

## Read a fit

Tidy verbs that return one-row-per-observation data frames.

- [`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md)
  : Transitions of an LSA Fit (Tidy)
- [`nodes()`](https://mohsaqr.github.io/lagdynamics/reference/nodes.md)
  : Nodes of an LSA Fit (Tidy)
- [`tests()`](https://mohsaqr.github.io/lagdynamics/reference/tests.md)
  : Tablewise Independence Tests of an LSA Fit (Tidy)
- [`initial()`](https://mohsaqr.github.io/lagdynamics/reference/initial.md)
  : Initial-State Distribution of an LSA Fit (Tidy)
- [`transition_probabilities()`](https://mohsaqr.github.io/lagdynamics/reference/transition_probabilities.md)
  : Transition-Probability Matrix of an LSA Fit

## Inference

Uncertainty, robustness, and reliability of a fit.

- [`bootstrap_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/bootstrap_lsa.md)
  : Bootstrap Confidence Intervals for an LSA Fit
- [`certainty_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/certainty_lsa.md)
  : Analytic Certainty of Transition Edges (Dirichlet-Multinomial)
- [`permute_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/permute_lsa.md)
  : Permutation Test for an LSA Fit
- [`stability_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/stability_lsa.md)
  : Case-Drop Stability for an LSA Fit
- [`reliability_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/reliability_lsa.md)
  : Split-Half Reliability for an LSA Fit

## Group comparison

- [`compare_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/compare_lsa.md)
  : Compare Groups' Transition Structures
- [`bayes_compare_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/bayes_compare_lsa.md)
  : Bayesian Comparison of Group Transition Structures
  (Dirichlet-Multinomial)

## Plotting

- [`plot_transitions()`](https://mohsaqr.github.io/lagdynamics/reference/plot_transitions.md)
  : Plot the Transition Network
- [`plot(`*`<lsa>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/plot.lsa.md)
  [`plot(`*`<lsa_group>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/plot.lsa.md)
  : Plot an LSA Fit
- [`plot_chords()`](https://mohsaqr.github.io/lagdynamics/reference/plot_chords.md)
  : Circular (Chord) Diagram of an LSA Fit
- [`plot_forest()`](https://mohsaqr.github.io/lagdynamics/reference/plot_forest.md)
  [`plot(`*`<lsa_bootstrap>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/plot_forest.md)
  : Circular Bootstrap Forest of an LSA Fit
- [`plot_polar()`](https://mohsaqr.github.io/lagdynamics/reference/plot_polar.md)
  : Polar Sunburst of an LSA Fit
- [`plot(`*`<lsa_certainty>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/plot.lsa_certainty.md)
  : Plot an Analytic-Certainty Result
- [`plot(`*`<lsa_comparison>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/plot.lsa_comparison.md)
  [`plot(`*`<lsa_comparison_pairwise>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/plot.lsa_comparison.md)
  : Plot a Group Comparison

## Transfer entropy (experimental)

- [`transfer_entropy()`](https://mohsaqr.github.io/lagdynamics/reference/transfer_entropy.md)
  : Directed transfer entropy for categorical sequences (experimental)

## Engine registry

- [`register_lsa_engine()`](https://mohsaqr.github.io/lagdynamics/reference/register_lsa_engine.md)
  : Register a Lag Sequential Analysis Engine
- [`unregister_lsa_engine()`](https://mohsaqr.github.io/lagdynamics/reference/unregister_lsa_engine.md)
  : Remove a Registered LSA Engine
- [`get_lsa_engine()`](https://mohsaqr.github.io/lagdynamics/reference/get_lsa_engine.md)
  : Retrieve a Registered LSA Engine
- [`list_lsa_engines()`](https://mohsaqr.github.io/lagdynamics/reference/list_lsa_engines.md)
  : List All Registered LSA Engines

## Tidy result objects

- [`as.data.frame(`*`<lsa_data>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/as.data.frame.lsa_data.md)
  : Tidy the Canonical Sequence Object
- [`as.data.frame(`*`<lsa_comparison>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/as.data.frame.lsa_comparison.md)
  [`as.data.frame(`*`<lsa_comparison_pairwise>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/as.data.frame.lsa_comparison.md)
  : Tidy a Group Comparison
- [`as.data.frame(`*`<lsa_reliability>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/as.data.frame.lsa_reliability.md)
  [`as.data.frame(`*`<lsa_reliability_group>`*`)`](https://mohsaqr.github.io/lagdynamics/reference/as.data.frame.lsa_reliability.md)
  : Tidy the per-replicate split-half correlations

## Data

- [`ai_long`](https://mohsaqr.github.io/lagdynamics/reference/ai_long.md)
  : Human-AI Vibe Coding Interaction Events
- [`engagement`](https://mohsaqr.github.io/lagdynamics/reference/engagement.md)
  : Student Engagement Trajectories
- [`group_regulation`](https://mohsaqr.github.io/lagdynamics/reference/group_regulation.md)
  : Collaborative Learning Self-Regulation Sequences
- [`group_regulation_long`](https://mohsaqr.github.io/lagdynamics/reference/group_regulation_long.md)
  : Group Regulation Long Event Log

## Package

- [`lagdynamics`](https://mohsaqr.github.io/lagdynamics/reference/lagdynamics-package.md)
  [`lagdynamics-package`](https://mohsaqr.github.io/lagdynamics/reference/lagdynamics-package.md)
  : lagdynamics: Lag Sequential Analysis, Dynamics, and Lag Transition
  Networks
