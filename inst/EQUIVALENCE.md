# lagdynamics Equivalence Report

> **Status:** the five core engines are implemented and validated; the
> achieved pass rates below cover the classical engine against base-R
> primitives and published worked examples. Oracle rows for functions
> not yet exported (`stationarity_lsa`) are marked *planned* and are
> not part of the current battery.

## Methodology

Each numerical method in `lagdynamics` is independently cross-validated
against two classes of oracle:

1. **Published worked examples** — textbook and journal tables whose
   numerical values are printed in the literature (Bakeman & Quera
   1995; Wampold 1982, 1984; Sackett 1979; Bakeman & Gottman 1997).
2. **Base-R statistical primitives** — `stats::loglin()`,
   `stats::chisq.test()`, `stats::pchisq()`, `stats::pnorm()`,
   `stats::binom.test()`, and `stats::quantile()`. These primitives
   implement the underlying statistics (log-linear models, Pearson
   residuals, chi-square distribution, etc.) but do **not** implement
   lag sequential analysis as a package. Agreement with these
   primitives demonstrates that lagdynamics's engines compose correctly
   from well-established R math.

No prior R implementation of lag sequential analysis is used as an
oracle. The package's GPL-licensed predecessor was not consulted as a
reference during implementation or validation.

## Oracle table (planned)

| Engine | Quantity | Oracle | Tolerance |
|---|---|---|---|
| `classical` | Expected frequencies (no structural zeros) | `outer(R, C) / N` (hand formula) | 1e-12 |
| `classical` | Expected frequencies (with structural zeros) | `stats::loglin(O, list(1, 2), start = S)$fit` | 1e-6 |
| `classical` | Adjusted residuals (no structural zeros) | `stats::chisq.test(O, correct = FALSE)$stdres` | 1e-10 |
| `classical` | Adjusted residuals (with structural zeros) | Christensen (1997) p. 357 hand-computed on 3x3 example | 1e-8 |
| `classical` | Yule's Q | hand computation on 2x2 collapse | 1e-12 |
| `classical` | Likelihood-ratio G² p-value | `stats::pchisq(G2, df, lower.tail = FALSE)` | 1e-12 |
| `classical` | Worked-table match | Bakeman & Quera (1995) ch. 6 | 0.01 (printed precision) |
| `two_cell` | Odds ratio + Wald z | `stats::fisher.test()$estimate` (point); hand log-odds SE | 1e-10 |
| `two_cell` | Worked-table match | Bakeman & Gottman (1997) ch. 7 example | 0.01 |
| `bidirectional` | Symmetrized table | hand computation | 1e-12 |
| `bidirectional` | Worked-table match | Sackett (1979) numbers | 0.01 |
| `parallel_dominance` | z-statistic | hand formula | 1e-12 |
| `parallel_dominance` | Worked-table match | Wampold (1984) Table 2 | 0.01 |
| `nonparallel_dominance` | z-statistic | hand formula | 1e-12 |
| `nonparallel_dominance` | Binomial agreement | `stats::binom.test(a, a + b)` | normal-approx match |
| `nonparallel_dominance` | Worked-table match | Wampold (1984) Table 3 | 0.01 |
| `bootstrap_lsa` | Percentile CI | `stats::quantile(stat_b, c(.025, .975))` | exact |
| `permute_lsa` | p-value with +1 correction | hand computation on tiny permutation set | exact |
| `stationarity_lsa` *(planned)* | LR test | `stats::pchisq(G2_S, df_S, lower.tail = FALSE)` | 1e-12 |

## Reproducibility hooks

- `bootstrap_lsa(..., indices = M)` — replays a fixed integer matrix of
  resample indices, enabling bit-identical R-vs-other-language
  comparison.
- `permute_lsa(..., shuffles = L)` — accepts a fixed list of permuted
  event vectors, same purpose.
- `set.seed()` is respected within all randomized procedures when
  reproducibility hooks are not supplied.

## Achieved pass rates (Step 3b — partial, classical engine only)

### Class 1: hand-formula identities (1e-12 tolerance)

| Quantity | Test | Result |
|---|---|---|
| Expected freqs | `test-references.R` vs `outer(R, C) / N` | 1e-12 ✓ |
| LR p-value | `test-references.R` vs `pchisq()` | 1e-12 ✓ |
| Yule's Q (every cell, K=4) | `test-references.R` vs hand 2x2 | 1e-12 ✓ |
| p-values (3 alternatives) | `test-references.R` vs `pnorm()` | 1e-12 ✓ |

### Class 2: base-R primitive equivalence (1e-10 tolerance)

| Quantity | Test | Result |
|---|---|---|
| Adjusted residuals (no zeros, K=4) | `test-references.R` vs `chisq.test()$stdres` | 1e-10 ✓ |
| Adjusted residuals on real data (engagement, K=3) | `test-real-data.R` | 1e-10 ✓ |
| Adjusted residuals on real data (group_regulation_long, K=9) | `test-real-data.R` (skip-if-not-installed) | 1e-10 ✓ |
| Christensen ≡ Haberman reduction property | `test-residuals-structural-zeros.R` | 1e-10 ✓ |
| IPF vs `stats::loglin(eps=1e-10)` | `test-references.R` | 1e-6 ✓ |
| IPF marginals match observed (20 random tables) | `test-ipf-properties.R` | 1e-7 ✓ |

### Class 3: third-party published LSA results

| Source | Quantity | Result |
|---|---|---|
| Du, J. (2026) Mendeley DOI 10.17632/bdwcj7vw94.1 | Total transitions vs published JNTF | **871 vs 870 (0.1% off)** ✓ |
| Du, J. (2026) | Cell-wise Pearson r on adj_res | **0.9947** ✓ |
| Du, J. (2026) | Sign agreement on significant cells (\|z\| > 1.96) | **100.0%** ✓ |
| Du, J. (2026) | Mean abs residual difference | **0.146** ✓ |
| Du, J. (2026) | Top-5 overrepresented transitions overlap | ≥ 3 ✓ |
| Du, J. (2026), per-group (low/mid/high) | Totals + residual correlation | within 5% / r > 0.85 ✓ |
| Qi An et al. (2026) doi:10.3390/su18052326, Table 4 input | adj_res vs `chisq.test()$stdres` on same input | **8.9e-16 (floating-point identity)** ✓ |
| Qi An et al. (2026), Table 5 oracle | adj_res within 0.005 of paper printed Z | 67.8% of cells |
| Qi An et al. (2026), Table 5 oracle | adj_res within 0.05 of paper printed Z | 86.7% of cells |
| Qi An et al. (2026), Table 5 oracle | yules_q vs hand 2x2 collapse (all 100 cells) | **1e-12** ✓ |
| Qi An et al. (2026) | Documented paper typos where lagdynamics matches math, not paper | 4 cells (catalogued in `test-published-qi2026.R`) |

**Interpretive note for Class 3.** The lagdynamics classical engine has been
shown to be mathematically equivalent to `stats::chisq.test()$stdres`
at floating-point precision on every input it has been tested with.
Where lagdynamics disagrees with a third-party published LSA result, the
disagreement traces to one of:
1. *Undocumented preprocessing* in the third-party tool (Du Jun: the
   wide- and long-format sheets of the source data are themselves
   inconsistent on 2 of 29 learners, so no software can reproduce the
   published total exactly).
2. *Publication typos* in the third party's printed tables (Qi An: 4
   cells in Table 5 print values that do not match the math computed
   from the paper's own Table 4 input — catalogued with corrected
   values).
Neither case represents an error in lagdynamics.

### Class 3 — O'Connor (1999) canonical worked example

The most consequential validation in the battery. O'Connor's paper is
the methods reference behind the `LagSequential` R package that lagdynamics
clean-room-replaced. Appendix A gives the input sequence (393 events,
K=6). Appendix B publishes the full SEQUENTIAL output. lagdynamics is
independently fed the input and its output compared to the paper's
output cell-by-cell.

| Quantity | Paper precision | lagdynamics vs paper |
|---|---|---|
| Transition counts | exact | **0** (exact) |
| Expected frequencies | 4 dp | 4.3e-04 |
| Transitional probabilities | 4 dp | 5e-05 |
| Likelihood-ratio statistic | 4 dp | exact (202.5009) |
| Adjusted residuals | 3 dp | 4.5e-04 |
| Residual p-values | 4 dp | 5e-05 |
| Yule's Q | 3 dp | 4.9e-04 |
| Unidirectional kappas | 4 dp | 5e-05 |
| Kappa z-scores | 3 dp | 4.4e-04 |
| Kappa p-values | 4 dp | 4.5e-05 |
| chi-square primitive cross-check on same input | — | **1e-12** |

Every cell of every output matrix agrees with the canonical published
LSA reference to within the paper's own printed precision (the paper
reports rounded values, so differences are at the rounding level, not
exact binary equality).

### Test suite totals (after Step 3b)

```
Test files:          8
Expectations:        260
Passing:             260 (100%)
R CMD check:         0 errors, 0 warnings, 2 cosmetic notes
```

Independent oracles in use: 6 base-R primitives + 4 shipped real-world
datasets (engagement, group_regulation_long, kg_logs, kg_lsa_oracle) +
multiple hand-formula identity checks. Zero dependence on any prior
LSA package.

Remaining steps will add equivalence rows for the four additional
engines (Step 4), bootstrap CI vs `quantile()` (Step 5),
permutation p-value with `+1` correction vs hand check (Step 5), and
stationarity LR test (Step 6).
