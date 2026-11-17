# Register a lag sequential analysis engine

`register_lsa_engine()` adds a lag sequential analysis engine under a
name, after which `lsa(data, engine = name)` fits it; the fit is read,
resampled and plotted like any other. An engine is the part of a lag
sequential analysis that turns a table of observed transitions into
expected counts and cell-wise tests. The built-in engines implement the
classical analysis of Bakeman and Quera (1995), a two-cell test,
Sackett's (1979) bidirectional test, and the parallel and non-parallel
dominance tests (Wampold, 1984). Other tests of the same table, such as
a different residual or a model-based expectation, can be added as
engines and then used with every function of the package.

Re-registering a name replaces the engine for the session. The built-in
engines are registered when the package loads.

## Usage

``` r
register_lsa_engine(name, fn, description, requires = character())
```

## Arguments

- name:

  Name of the engine, used as `lsa(engine = name)`.

- fn:

  The engine function. It is called with `transitions` (the
  [`lsa_transitions()`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md)
  object), `structural_zeros` (`NULL` or a `K x K` 0/1 matrix),
  `alternative`, `n_events` (`NULL` for transition-matrix input) and any
  engine parameters given to
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md), so
  it should accept `transitions` and `...`. It must return a named list
  holding at least the `K x K` matrices `obs`, `exp`, `prob`, `adj_res`
  and `p`; further `K x K` matrices become extra columns of
  [`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md).

- description:

  One-line description shown by
  [`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md).

- requires:

  Names of packages the engine needs (default none);
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)
  checks that they are installed before fitting.

## Value

`name`, invisibly.

## Conditions

Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`
when `fn` has no `transitions` argument.

## References

Bakeman, R., & Quera, V. (1995). *Analyzing interaction: Sequential
analysis with SDIS and GSEQ*. Cambridge University Press.

Sackett, G. P. (1979). The lag sequential analysis of contingency and
cyclicity in behavioral interaction research. In J. D. Osofsky (Ed.),
*Handbook of infant development* (pp. 623-649). Wiley.

Wampold, B. E. (1984). Tests of dominance in sequential categorical
data. *Psychological Bulletin, 96*(2), 424-429.

## See also

[`list_lsa_engines()`](https://pak.dynasite.org/lagdynamics/reference/list_lsa_engines.md),
[`get_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/get_lsa_engine.md),
[`unregister_lsa_engine()`](https://pak.dynasite.org/lagdynamics/reference/unregister_lsa_engine.md),
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)

## Examples

``` r
# An engine that tests Pearson residuals under independence
pearson_engine <- function(transitions, ...) {
  obs <- transitions$obs
  expected <- outer(rowSums(obs), colSums(obs)) / sum(obs)
  residual <- (obs - expected) / sqrt(expected)
  list(obs = obs, exp = expected, prob = obs / rowSums(obs),
       adj_res = residual, p = 2 * stats::pnorm(-abs(residual)))
}
register_lsa_engine("pearson", pearson_engine,
                    description = "Pearson residuals under independence")
list_lsa_engines()
#>                    name
#> 1         bidirectional
#> 2             classical
#> 3 nonparallel_dominance
#> 4    parallel_dominance
#> 5               pearson
#> 6              two_cell
#>                                                            description requires
#> 1 Sackett's bidirectional / matched-pair test on the symmetrized table         
#> 2                    Bakeman & Quera classical lag sequential analysis         
#> 3       Sackett's non-parallel-dominance (observed-SE + binomial) test         
#> 4                      Sackett's parallel-dominance (expected-SE) test         
#> 5                                 Pearson residuals under independence         
#> 6                  2x2 cell test (odds ratio, log-OR Wald z, Yule's Q)         
fit <- lsa(engagement, engine = "pearson")
transitions(fit, significant = TRUE)
#>         from         to lag count  expected       prob   prob_col   adj_res
#> 1     Active     Active   1   459 247.03460 0.69756839 0.70506912 13.486091
#> 2    Average     Active   1   153 281.94983 0.20372836 0.23502304 -7.679534
#> 3 Disengaged     Active   1    39 122.01557 0.12000000 0.05990783 -7.515397
#> 4     Active    Average   1   176 289.53518 0.26747720 0.23066841 -6.672365
#> 5    Average    Average   1   458 330.45732 0.60985353 0.60026212  7.016133
#> 6     Active Disengaged   1    23 121.43022 0.03495441 0.07187500 -8.932336
#> 7 Disengaged Disengaged   1   157  59.97693 0.48307692 0.49062500 12.528033
#>              p yules_q kappa kappa_z kappa_p      lift  sign significant
#> 1 1.888519e-41      NA    NA      NA      NA 1.8580393  over        TRUE
#> 2 1.596688e-14      NA    NA      NA      NA 0.5426497 under        TRUE
#> 3 5.673846e-14      NA    NA      NA      NA 0.3196313 under        TRUE
#> 4 2.517140e-11      NA    NA      NA      NA 0.6078709 under        TRUE
#> 5 2.280924e-12      NA    NA      NA      NA 1.3859581  over        TRUE
#> 6 4.171055e-19      NA    NA      NA      NA 0.1894092 under        TRUE
#> 7 5.244751e-36      NA    NA      NA      NA 2.6176731  over        TRUE
unregister_lsa_engine("pearson")
```
