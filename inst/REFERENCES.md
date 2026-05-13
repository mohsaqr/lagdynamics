# lagseq: Formula Reference Map

This document is the **clean-room source of truth** for every numerical
method implemented in `lagseq`. Each formula below is transcribed
directly from the cited primary source with the page number recorded.
Implementation in `R/` works from this document only; no prior R
implementation of lag sequential analysis is consulted, examined, or
referenced during coding.

When an equation has been adapted to matrix form for vectorized R
implementation, the original scalar form from the literature is shown
first, the matrix form second, and the two are proven equivalent in the
test file `tests/testthat/test-references.R`.

## Notation

- `K` — number of distinct codes (states)
- `N` — total number of lag-1 transitions in the sequence
- `O[i, j]` — observed count of transitions from code `i` to code `j`
- `E[i, j]` — expected count under independence
- `R[i, .]` — row totals (i.e. `rowSums(O)`)
- `C[., j]` — column totals (i.e. `colSums(O)`)
- `p_i = R[i] / N` — row marginal probability
- `q_j = C[j] / N` — column marginal probability
- `n` — total number of events (= N + 1 for a single sequence at lag 1)
- `n_i` — total event frequency for code `i`

---

## 1. Transition count matrix

**Source:** Bakeman & Quera (1995), *Analyzing Interaction*, ch. 4, p. 102.

For a single integer-coded sequence `x` of length `T` and lag `L`,

```
O[i, j] = #{ t : x[t] = i AND x[t + L] = j, 1 <= t <= T - L }
```

**Vectorized R form** (no `for` loop):

```r
pair_index <- (x[seq_len(T - L)] - 1L) * K + x[seq(1L + L, T)]
O          <- matrix(tabulate(pair_index, nbins = K * K), K, K, byrow = TRUE)
```

The `byrow = TRUE` is correct because `pair_index` linearizes
row-major.

For multiple sequences, sum per-sequence matrices. Within-sequence
transitions only; no transitions span sequence boundaries.

---

## 2. Expected frequencies under independence

### 2.1 Adjacent codes may repeat (`engine = "classical"`, default)

**Source:** Bakeman & Quera (1995), p. 110, eq. 6.2.

```
E[i, j] = R[i] * C[j] / N
```

**Matrix form:**

```r
E <- outer(R, C) / N
```

### 2.2 Adjacent codes may NOT repeat (structural zeros on diagonal)

**Source:** Bakeman & Quera (1995), p. 113, eq. 6.4; iterative
proportional fitting (IPF) algorithm: Wickens (1989), *Multiway
Contingency Tables*, pp. 107-112.

Define a `K x K` 0/1 matrix `S` ("onezero") where `S[i, j] = 0` means
the cell is a structural zero (off-limits). For the
"no-self-transitions" case, `S = 1 - diag(K)`.

Initialize `E = S`. Repeat until row and column marginals of `E` match
those of `O` to within `1e-6`:

1. Row scaling: for each `i`, `E[i, j] <- E[i, j] * R[i] / sum(E[i, ])`
   where the sum is over `j` with `S[i, j] = 1`.
2. Column scaling: for each `j`,
   `E[i, j] <- E[i, j] * C[j] / sum(E[, j])`.

After convergence, `sum(E) == N`, `rowSums(E) == R`, `colSums(E) == C`,
and `E[i, j] = 0` wherever `S[i, j] = 0`.

**Equivalence oracle:** When the structural-zero pattern is supplied
as a `start =` argument and the margins are fixed,
`stats::loglin(O, margin = list(1, 2), start = S, fit = TRUE)$fit`
returns the same expected frequencies up to convergence tolerance.

---

## 3. Transitional probabilities

**Source:** Bakeman & Quera (1995), p. 104, eq. 5.1.

```
prob[i, j] = O[i, j] / R[i]    if R[i] > 0,  else NA
```

**Matrix form:**

```r
prob <- O / R       # R recycled across columns
prob[R == 0, ] <- NA_real_
```

---

## 4. Adjusted residuals (Haberman z-scores)

### 4.1 No structural zeros (`engine = "classical"`)

**Source:** Haberman (1979), *Analysis of Qualitative Data*, vol. 2,
eq. 5.2.6. Also Bakeman & Quera (1995), p. 111, eq. 6.3.

```
z[i, j] = (O[i, j] - E[i, j]) / sqrt( E[i, j] * (1 - p_i) * (1 - q_j) )
```

**Matrix form:**

```r
z <- (O - E) / sqrt( E * outer(1 - p_row, 1 - p_col) )
```

**Equivalence oracle:** `stats::chisq.test(O, correct = FALSE)$stdres`
returns Haberman's standardized Pearson residuals. For a sequence with
no structural zeros and ignoring the `n` vs `N + 1` correction (see
§4.3), the values are numerically identical to `z`.

### 4.2 With structural zeros

**Source:** Christensen (1997), *Log-Linear Models and Logistic
Regression* (2nd ed.), p. 357.

Let `X` be the design matrix encoding the two-way independence model
restricted to cells with `S[i, j] = 1`. Stack the non-structural-zero
cells into a vector. With `W = diag(vec_E)`, the hat matrix is

```
H = X * (X' W X)^{-1} * X' * W
```

For each non-structural-zero cell with hat diagonal `h`,

```
z[i, j] = (O[i, j] - E[i, j]) / sqrt( E[i, j] * (1 - h) )
```

For structural-zero cells, `z[i, j] = 0` by definition.

The independence-model design matrix can be built with
`stats::model.matrix(~ row + col)` after expanding the table to long
form and dropping structural-zero rows. This is the implementation
path; no hand-rolled `cbind(x, y, z)` construction is used.

### 4.3 The `n` correction (events vs transitions)

**Source:** Bakeman & Quera (1995), p. 105, footnote 3.

When the input is a sequence of events (not a pre-computed transition
matrix), the marginal probabilities `p_i = n_i / n` should use
**event** totals, not transition totals. For a single sequence of
length `T` at lag 1, `n = T` and the last event contributes to event
frequencies but not transition frequencies. This affects kappa
computation (§7) but does NOT affect the adjusted residual formula
above, which uses transition marginals `R` and `C`.

The flag `fit$params$source` records whether the input was an event
sequence (`"events"`) or a transition matrix (`"transitions"`); kappa
uses `n` from `fit$params$n_events`.

---

## 5. p-values for adjusted residuals

**Source:** standard normal distribution.

```
p[i, j] = 2 * (1 - Phi(|z[i, j]|))     for alternative = "two.sided"
p[i, j] = 1 - Phi(z[i, j])             for alternative = "greater"
p[i, j] = Phi(z[i, j])                 for alternative = "less"
```

where `Phi` is the standard normal CDF (`stats::pnorm`).

---

## 6. Yule's Q

**Source:** Bakeman & Gottman (1997), *Observing Interaction* (2nd
ed.), p. 129, eq. 7.7.

For each cell `(i, j)`, treat the `K x K` table as a 2x2 collapse:

```
a = O[i, j]
b = R[i] - O[i, j]
c = C[j] - O[i, j]
d = N - R[i] - C[j] + O[i, j]
```

```
Q[i, j] = (a * d - b * c) / (a * d + b * c)    if denominator > 0
        = NA                                   otherwise
```

Q is bounded in `[-1, 1]`. Positive Q means over-representation,
negative means under-representation, zero means independence.

---

## 7. Unidirectional kappa

**Source:** Bakeman & Quera (1995), p. 115, eq. 6.6; Wampold (1989),
"Kappa as a measure of pattern in sequential data", *Quality &
Quantity*, 23(2).

Let `et[i, j]` be the expected count and `var[i, j]` its variance
under the null that consecutive events are independent at the **event**
level (not the transition level):

```
et[i, j]  = n_i * n_j / n                       for all i, j (diagonal and off)

var[i, j] = n_i * n_j * (n - n_j) * (n - n_i) / ( n^2 * (n - 1) )
```

**Convention note.** Wampold (1989) derives a "without-replacement"
diagonal correction `et[i,i] = n_i * (n_i - 1) / n`. However, the
canonical LSA software implementations (O'Connor's SEQUENTIAL,
Bakeman & Quera's GSEQ) apply the off-diagonal formula uniformly,
including on the diagonal. lagseq adopts the SEQUENTIAL/GSEQ
convention so its outputs match the canonical published worked
example in O'Connor (1999) cell-for-cell.

Maximum possible count:

```
mmax[i, j] = min(n_i, n_j)
```

Kappa:

```
kappa[i, j] = (O[i, j] - et[i, j]) / (mmax[i, j] - et[i, j])     if O > et
            = (O[i, j] - et[i, j]) / et[i, j]                    if O < et
            = NA                                                  if denominator == 0
```

Kappa's z-statistic:

```
z_kappa[i, j] = (O[i, j] - et[i, j]) / sqrt(var[i, j])
```

with p-values as in §5.

---

## 8. Tablewise likelihood-ratio chi-square

**Source:** Bakeman & Quera (1995), p. 109, eq. 6.1. Wickens (1989),
p. 26.

```
G^2 = 2 * sum_{i, j : O[i,j] > 0 AND E[i,j] > 0} O[i, j] * log( O[i, j] / E[i, j] )
```

Degrees of freedom:

- Classical (no structural zeros): `df = (K - 1)^2`
- With structural zeros (count `s` of them): `df = (K - 1)^2 - s`

p-value: `1 - pchisq(G^2, df)` — uses `stats::pchisq`.

---

## 9. Two-cell test (`engine = "two_cell"`)

**Source:** Bakeman & Gottman (1997), ch. 7, pp. 130-134. Sometimes
called the "2x2 cell test" or "binomial sequential test".

For a focal cell `(i, j)` collapse the `K x K` table to a 2x2:

```
              j      not j
        +----------+----------+
    i   |    a     |    b     |
        +----------+----------+
   !i   |    c     |    d     |
        +----------+----------+
```

with `a = O[i, j]`, `b = R[i] - O[i, j]`, `c = C[j] - O[i, j]`,
`d = N - R[i] - C[j] + O[i, j]`.

Reported statistics: Yule's Q (§6), odds ratio `OR = (a * d) / (b * c)`,
log-odds-ratio standard error
`SE = sqrt(1/a + 1/b + 1/c + 1/d)`, Wald z = `log(OR) / SE`,
two-sided p-value.

---

## 10. Bidirectional / matched-pair test (`engine = "bidirectional"`)

**Source:** Sackett (1979), "The lag sequential analysis of contingency
and cyclicity", in J. D. Osofsky (Ed.), *Handbook of Infant
Development*, pp. 623-649. Also Wampold (1984), *Psychological
Bulletin*, 96(2), 424-429.

For an unordered pair `{i, j}`, sum the two directed counts:

```
W[i, j] = O[i, j] + O[j, i]    (symmetric matrix)
```

with expected count

```
E_W[i, j] = E[i, j] + E[j, i]
```

and adjusted z computed via §4 on the symmetrized table. The test
asks whether the pair `{i, j}` co-occurs more than expected at lag 1,
ignoring direction.

---

## 11. Parallel-dominance analysis (`engine = "parallel_dominance"`)

**Source:** Sackett (1979), pp. 631-635. Tested in Wampold (1984),
Table 2.

For an ordered pair `(i, j)`, the parallel-dominance statistic asks
whether `i -> j` is more likely than `j -> i`:

```
D_par[i, j] = O[i, j] - O[j, i]
```

with expected value `0` under symmetric independence and standard error

```
SE[i, j] = sqrt( E[i, j] + E[j, i] )
```

(using the row/column independence expected counts from §2). The test
statistic is

```
z[i, j] = D_par[i, j] / SE[i, j]
```

Note `z[i, j] = -z[j, i]`, so only the upper triangle is reported.

---

## 12. Non-parallel-dominance analysis (`engine = "nonparallel_dominance"`)

**Source:** Sackett (1979), pp. 635-638. Tested in Wampold (1984),
Table 3.

Same numerator as §11 but standard error uses the **observed** count
sum rather than the expected:

```
SE_np[i, j] = sqrt( O[i, j] + O[j, i] )
z_np[i, j]  = D_par[i, j] / SE_np[i, j]
```

This is the asymmetry test conditioned on the total dyadic count
`O[i, j] + O[j, i]`. Under the null `H_0: P(i->j) = P(j->i)` given the
total, the count `O[i, j]` is `Binomial(O[i, j] + O[j, i], 1/2)` and
`z_np` is the corresponding normal approximation.

---

## 13. Bootstrap (sequence-level)

**Source:** Efron (1979), "Bootstrap methods: another look at the
jackknife", *Annals of Statistics*, 7(1), 1-26. Sequence-level
resampling protocol: Politis & Romano (1994), "The stationary
bootstrap", *Journal of the American Statistical Association*, 89(428),
1303-1313.

When input is a list of `S` sequences, draw `R` resamples by selecting
`S` sequences with replacement from the list. Refit `lsa()` on each
resample using the **same** `fit$params` recipe. Aggregate per-edge
statistics across the `R` fits:

```
mean_b      = mean(stat_b)
se_b        = sd(stat_b)
ci_low_b    = quantile(stat_b, (1 - level) / 2)
ci_high_b   = quantile(stat_b, 1 - (1 - level) / 2)
p_boot      = 2 * min(mean(stat_b <= 0), mean(stat_b >= 0))   (two-sided)
stable      = sign(ci_low_b) == sign(ci_high_b)
```

When input is a single sequence, fall back to event-level
resampling with a stationary bootstrap of geometric block length
`mean_block = ceiling(sqrt(T))`. The choice is exposed as
`level = c("sequence", "event")` in `bootstrap_lsa()`.

**Reproducibility hook:** the function accepts an optional
`indices = matrix(integer, R, S)` argument; when supplied, the
bootstrap replays exactly those resample indices and ignores its
internal RNG. This mirrors psychaj's `permutationIndices` parameter
and enables cross-language verification.

---

## 14. Permutation test (`permute_lsa`)

**Source:** Castellan (1992), "Shuffling arrays: appearances may be
deceiving", *Behavior Research Methods, Instruments, & Computers*,
24(1), 72-77. Algorithm extension for the no-repeat constraint:
ibid., section 3.

Under `H_0` of no first-order sequential dependence, the event order
is exchangeable. The permutation distribution of any cell statistic
is constructed by repeatedly shuffling the event vector and
recomputing the statistic.

Two shuffle variants:

- **No repeat constraint** (`adjacent codes may repeat`): a standard
  Fisher-Yates shuffle of the event vector.
- **Repeat constraint** (`adjacent codes may NOT repeat`): a
  constrained Fisher-Yates where each proposed swap is accepted only
  if it preserves the no-adjacent-repeat property at both sites. If
  no valid swap is found within an inner cap (10,000 attempts at each
  site), the site is skipped. Pad with sentinel zeros at both ends so
  boundary swaps have a defined neighbor.

p-values:

```
p_perm[i, j] = (1 + #{ b : sign(stat_b - 0) == sign(stat_obs) AND
                            |stat_b| >= |stat_obs| }) / (1 + R)
```

The `+1` correction (Phipson & Smyth 2010) prevents `p = 0`.

---

## 15. Stationarity / homogeneity LR test

**Source:** Bakeman & Quera (1995), ch. 8, pp. 138-142, eq. 8.1.

Given `M` segments (groups), let `O_m` be the transition count matrix
for segment `m` with row totals `R_m`, and let `O_pool = sum O_m`
with row totals `R_pool`. The pooled transition probabilities are
`p_pool[i, j] = O_pool[i, j] / R_pool[i]`. For each segment,

```
G^2_m = 2 * sum_{i, j} O_m[i, j] * log( p_m[i, j] / p_pool[i, j] )
```

(omitting cells where `O_m == 0` or `p_pool == 0`). The total
stationarity statistic is

```
G^2_S = sum_m G^2_m
```

with degrees of freedom `df_S = K * (K - 1) * (M - 1)` and p-value
`1 - pchisq(G^2_S, df_S)`.

For homogeneity (groups in parallel rather than segments in time), the
same statistic is computed but documented as a different scientific
question; the math is identical.

---

## 16. Equivalence oracles summary

The following base-R primitives are independent oracles against which
each engine is cross-validated. None of them implements lag sequential
analysis; they implement the underlying statistical primitives. A
passing oracle test demonstrates that lagseq's engine has the same
numerical behavior as long-established R primitives, with zero
dependence on any prior LSA package.

| Quantity | Oracle | Test file |
|---|---|---|
| Expected frequencies under independence | `outer(rowSums(O), colSums(O)) / sum(O)` | `test-references.R` |
| Expected frequencies under independence with structural zeros | `stats::loglin(O, margin = list(1, 2), start = S)` | `test-ipf.R` |
| Adjusted (Haberman) residuals, no structural zeros | `stats::chisq.test(O, correct = FALSE)$stdres` | `test-references.R` |
| Likelihood-ratio chi-square p-value | `stats::pchisq(G2, df, lower.tail = FALSE)` | `test-references.R` |
| Yule's Q | hand computation on 2x2 collapse | `test-references.R` |
| Binomial sign test for non-parallel dominance | `stats::binom.test(O[i, j], O[i, j] + O[j, i])` | `test-engine-nonparallel-dominance.R` |
| Bootstrap percentile CI | `stats::quantile(stat_b, c(.025, .975))` | `test-bootstrap.R` |
| Permutation p-value with `+1` correction | hand check on small permutation set | `test-permute.R` |

In addition, published worked tables provide a second class of oracle:

| Source | Test file |
|---|---|
| Bakeman & Quera (1995), example in ch. 6 | `test-published-bakeman-quera.R` |
| Wampold (1982), Table 1 | `test-published-wampold-1982.R` |
| Wampold (1984), Tables 2 & 3 | `test-published-wampold-1984.R` |
| Sackett (1979), dominance worked numbers | `test-published-sackett-1979.R` |
| Bakeman & Gottman (1997), ch. 7 two-cell example | `test-published-bakeman-gottman.R` |

---

## Closing note on independence from prior LSA software

This document was written from primary literature only. No source code
of any prior R or non-R lag-sequential-analysis package was consulted
during the writing of this reference or during the implementation of
the corresponding R files in `R/`. The package's GPL-licensed
predecessor was not used as a code reference or as a validation
oracle. The validation strategy relies exclusively on (a) textbook and
journal-published worked examples and (b) base-R primitives that
implement underlying statistical concepts but not LSA itself.
