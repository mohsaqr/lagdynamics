# Directed transfer entropy for categorical sequences (experimental)

**Experimental.** Transfer entropy (Schreiber, 2000) measures *directed*
predictive coupling: how much knowing the source's past reduces
uncertainty about the target's future, **beyond what the target's own
past already explains**. Because it conditions on the target's own
history, transfer entropy is immune to the autocorrelation confound that
inflates plain lagged association when a process has strong momentum.

Unlike
[`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)'s
Yule's Q / adjusted residuals, transfer entropy is **sign-blind**: a
large value means "strong directed predictive structure", which may be
*facilitating* or *suppressing*. Read it alongside the signed measures
from
[`transitions()`](https://mohsaqr.github.io/lagdynamics/reference/transitions.md)
to interpret direction of effect.

Two modes:

- **State-flow network** (default, `y = NULL`): given categorical
  sequences, returns directed transfer entropy between every ordered
  pair of states, via a binary occupancy decomposition. Note: with very
  few states the non-target source channels become redundant; prefer the
  bivariate mode or a larger alphabet when that matters.

- **Bivariate** (`y` supplied): full-alphabet transfer entropy between
  two aligned categorical series, in both directions.

## Usage

``` r
transfer_entropy(
  x,
  y = NULL,
  lag = 1L,
  history = 1L,
  test = c("surrogate", "none"),
  R = 199L,
  normalize = TRUE,
  seed = NULL
)
```

## Arguments

- x:

  Categorical sequence data: a vector (one sequence), or a matrix /
  data.frame with one sequence per row and one time-step per column
  (`NA`-padded), exactly the shape
  [`lsa()`](https://mohsaqr.github.io/lagdynamics/reference/lsa.md)
  consumes.

- y:

  Optional second series, same shape as `x`, for bivariate transfer
  entropy. When `NULL`, the directed state-flow network of `x` is
  returned.

- lag:

  Integer \>= 1. Prediction horizon: the target's future is taken `lag`
  steps ahead. Default `1`.

- history:

  Integer \>= 1. Order of the target's own history conditioned on (and
  combined into a composite symbol). Default `1`.

- test:

  `"surrogate"` (default) runs a source-permutation null to give a
  p-value and the bias-corrected *effective* transfer entropy; `"none"`
  skips it.

- R:

  Integer. Number of surrogate permutations. Default `199`.

- normalize:

  Logical. Add `te_normalised`, transfer entropy as a share of the
  target's leftover uncertainty `H(future | history)`, in `[0, 1]`.
  Default `TRUE`.

- seed:

  Optional integer seed for the surrogate test.

## Value

A tidy `data.frame`, one row per ordered pair, with columns `from`,
`to`, `te` (bits), `te_effective` (surrogate-debiased), `te_normalised`
(0-1, if `normalize`), `p` (surrogate p-value), and `n` (pooled
transitions used). Rows are ordered by descending `te`.

## References

Schreiber, T. (2000). Measuring information transfer. *Physical Review
Letters*, 85(2), 461-464.

## Examples

``` r
# Directed information-flow network over engagement states
transfer_entropy(engagement, test = "none")
#>         from         to          te te_normalised    n
#> 1     Active Disengaged 0.036448215   0.060023389 1734
#> 2    Average Disengaged 0.036448215   0.060023389 1734
#> 3     Active    Average 0.006949150   0.007517345 1734
#> 4 Disengaged    Average 0.006949150   0.007517345 1734
#> 5    Average     Active 0.004793199   0.006345037 1734
#> 6 Disengaged     Active 0.004793199   0.006345037 1734

# Bivariate transfer entropy between two aligned series
a <- c("calm", "calm", "tense", "tense", "calm", "tense", "tense", "calm")
b <- c("low", "low", "low", "high", "high", "low", "high", "high")
transfer_entropy(a, b, test = "none")
#>   from to        te te_normalised n
#> 1    a  b 0.9649839     1.0000000 7
#> 2    b  a 0.6792696     0.7039181 7
```
