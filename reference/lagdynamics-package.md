# lagdynamics: Lag Sequential Analysis, Dynamics, and Lag Transition Networks

A unified, pipe-friendly interface for lag sequential analysis (LSA) of
categorical event sequences. Provides classical, two-cell,
bidirectional, parallel-dominance, and non-parallel-dominance engines
via a single
[`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)
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
[`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md),
[`lsa_data()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_data.md),
[`lsa_transitions()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_transitions.md)

Convenience wrappers:
[`lsa_classical()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md),
[`lsa_two_cell()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md),
[`lsa_bidirectional()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md),
[`lsa_parallel_dominance()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md),
[`lsa_nonparallel_dominance()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)

Engine registry:
[`register_lsa_engine()`](https://mohsaqr.github.io/lagdynamics/reference/register_lsa_engine.md),
[`get_lsa_engine()`](https://mohsaqr.github.io/lagdynamics/reference/get_lsa_engine.md),
[`list_lsa_engines()`](https://mohsaqr.github.io/lagdynamics/reference/list_lsa_engines.md)

Inference:
[`bootstrap_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/permute_lsa.md),
[`stability_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/stability_lsa.md),
[`reliability_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/reliability_lsa.md),
and analytic certainty via
[`certainty_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/certainty_lsa.md).

Group comparison:
[`compare_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/compare_lsa.md)
and
[`bayes_compare_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/bayes_compare_lsa.md)
for between-group differences in transition structure.

Reading a fit (tidy data frames):
[`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md),
[`nodes()`](https://mohsaqr.github.io/lagdynamics/reference/nodes.md),
[`tests()`](https://mohsaqr.github.io/lagdynamics/reference/tests.md),
[`initial()`](https://mohsaqr.github.io/lagdynamics/reference/initial.md),
and [`summary()`](https://rdrr.io/r/base/summary.html).

Transition and initial probabilities:
[`transition_probabilities()`](https://mohsaqr.github.io/lagdynamics/reference/transition_probabilities.md)
returns the row-stochastic transition-probability matrix and
[`initial()`](https://mohsaqr.github.io/lagdynamics/reference/initial.md)
the initial-state probabilities.

## Roadmap (not yet implemented)

Stationarity tests (`stationarity_lsa()`) are planned but not yet
exported. Multi-group fits are built with `lsa(data, group = ...)`;
between-group differences are then tested with
[`compare_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/compare_lsa.md)
or
[`bayes_compare_lsa()`](https://mohsaqr.github.io/lagdynamics/reference/bayes_compare_lsa.md).

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

- <https://mohsaqr.github.io/lagdynamics/>

- Report bugs at <https://github.com/mohsaqr/lagdynamics/issues>

## Author

**Maintainer**: Mohammed Saqr <saqr@saqr.me> \[copyright holder\]
