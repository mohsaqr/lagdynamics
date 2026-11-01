# lagdynamics: Lag Sequential Analysis, Dynamics, and Lag Transition Networks

A unified, pipe-friendly interface for lag sequential analysis (LSA) of
categorical event sequences. Provides classical, two-cell,
bidirectional, parallel-dominance, and non-parallel-dominance engines
via a single
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
constructor with a pluggable engine registry, sequence-level bootstrap
and permutation inference, case-drop stability, and tidy edge tables
ready for transition-network visualization.

## Design principles

Numerical methods are implemented from primary literature and
cross-validated against base-R primitives and hand-formula identities
([`stats::loglin()`](https://rdrr.io/r/stats/loglin.html),
[`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html),
[`stats::pchisq()`](https://rdrr.io/r/stats/Chisquare.html)). See the
`inst/REFERENCES.md` file shipped with the package for the
formula-by-formula citation map.

The package returns S3 objects of class `c("lsa", "cograph_network")`,
which means fits inherit compatibility with the `cograph` package's
visualization layer when it is installed.

## Main functions

Construction:
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
(sequences, wide tables, transition matrices, long event logs and `tna`
/ `Nestimate` / `TraMineR` objects),
[`lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/lsa_data.md),
[`lsa_transitions()`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md),
and the multi-lag
[`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)
and
[`lag_profile()`](https://pak.dynasite.org/lagdynamics/reference/lag_profile.md).

Engine wrappers:
[`lsa_classical()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
[`lsa_two_cell()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
[`lsa_bidirectional()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
[`lsa_parallel_dominance()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
[`lsa_nonparallel_dominance()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
(documented with
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)).

Engine registry:
[`register_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/register_lsa_engine.md),
[`get_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/get_lsa_engine.md),
[`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md),
[`unregister_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/unregister_lsa_engine.md).

Reading a fit (tidy data frames):
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`nodes()`](https://pak.dynasite.org/lagdynamics/reference/nodes.md),
[`tests()`](https://pak.dynasite.org/lagdynamics/reference/tests.md),
[`initial()`](https://pak.dynasite.org/lagdynamics/reference/initial.md),
and [`summary()`](https://rdrr.io/r/base/summary.html) (see
[lsa-methods](https://pak.dynasite.org/lagdynamics/reference/lsa-methods.md)).
[`transition_probabilities()`](https://pak.dynasite.org/lagdynamics/reference/transition_probabilities.md)
returns the row-stochastic transition-probability matrix.

Inference:
[`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md),
[`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md),
and analytic certainty via
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md).
Every result is read with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md).

Group comparison:
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
(permutation) and
[`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
(Dirichlet-Multinomial) for between-group differences in transition
structure, on a grouped fit from `lsa(data, group = )`.

Plotting:
[`plot.lsa()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa.md)
(heatmap, network, chord, sunburst),
[`plot_transitions()`](https://pak.dynasite.org/lagdynamics/reference/plot_transitions.md),
[`plot_chords()`](https://pak.dynasite.org/lagdynamics/reference/plot_chords.md),
[`plot_polar()`](https://pak.dynasite.org/lagdynamics/reference/plot_polar.md),
[`plot_forest()`](https://pak.dynasite.org/lagdynamics/reference/plot_forest.md),
and
[`plot.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md).

Other:
[`transfer_entropy()`](https://pak.dynasite.org/lagdynamics/reference/transfer_entropy.md)
(directed information flow, experimental),
[`lsa_ipf()`](https://pak.dynasite.org/lagdynamics/reference/lsa_ipf.md)
(iterative proportional fitting for structural zeros), and
[`lsa_to_tna()`](https://pak.dynasite.org/lagdynamics/reference/lsa_to_tna.md)
(hand a fit to the `tna` package).

Data:
[engagement](https://pak.dynasite.org/lagdynamics/reference/engagement.md),
[group_regulation](https://pak.dynasite.org/lagdynamics/reference/group_regulation.md),
[group_regulation_long](https://pak.dynasite.org/lagdynamics/reference/group_regulation_long.md),
[ai_long](https://pak.dynasite.org/lagdynamics/reference/ai_long.md).

## Roadmap (not yet implemented)

Stationarity tests (`stationarity_lsa()`) are planned but not yet
exported.

## References

Bakeman, R., & Quera, V. (1995). *Analyzing interaction: Sequential
analysis*. Cambridge University Press.

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

Wampold, B. E. (1982). Sequential analysis: Stationarity, statistical
procedures, and applications. *Psychological Bulletin, 92*(2), 380-387.

Wampold, B. E. (1984). Tests of dominance in sequential categorical
data. *Psychological Bulletin, 96*(2), 424-429.

Wickens, T. D. (1989). *Multiway contingency tables analysis for the
social sciences*. Lawrence Erlbaum.

## See also

Useful links:

- <https://github.com/mohsaqr/lagdynamics>

- Report bugs at <https://github.com/mohsaqr/lagdynamics/issues>

## Author

**Maintainer**: Mohammed Saqr <saqr@saqr.me> \[copyright holder\]

Authors:

- Sonsoles López-Pernas
