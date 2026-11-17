# Estimate directed transfer entropy between categorical sequences

`transfer_entropy()` (**experimental**) estimates directed transfer
entropy between categorical sequences. Transfer entropy (Schreiber,
2000) measures directed predictive information: how much knowing the
past of a source series reduces uncertainty about the future of a target
series, beyond what the target's own past already tells. Because it
conditions on the target's history, it is not inflated when a state
merely tends to persist, which inflates lagged association measures for
processes with strong momentum. Transfer entropy is unsigned: a large
value means strong directed dependence, which may make the target more
or less likely. Read it alongside the signed residuals and Yule's Q of
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md).

It works in two modes. With one set of sequences (`y = NULL`), it treats
each state as an indicator series (in the state or not) and returns the
transfer entropy from every state to every other state, a directed
information-flow network among the states. With very few states the
indicator series of the other states carry nearly the same information,
so several pairs can share a value; the bivariate mode avoids this. With
two aligned sets of sequences (`x` and `y`), it returns the transfer
entropy between the two series in both directions, using their full
alphabets.

Transfer entropy is in bits. The surrogate test (default) shuffles the
source to build a null distribution, which gives a p-value and an
effective transfer entropy with the average surrogate value subtracted.
For the classical test of sequential dependence use
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

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

  Categorical sequences: a vector (one sequence) or a matrix or data
  frame with one sequence per row and one time point per column, padded
  with `NA`, as accepted by
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

- y:

  Optional second set of sequences of the same shape as `x`, for the
  bivariate mode. Default `NULL` (state-flow network of `x`).

- lag:

  Prediction horizon, an integer of at least 1 (default `1`).

- history:

  Number of past steps of the target conditioned on, an integer of at
  least 1 (default `1`).

- test:

  `"surrogate"` (default) for the surrogate test, or `"none"`.

- R:

  Number of surrogates (default `199`).

- normalize:

  Add `te_normalised`, transfer entropy as a share of the target's
  remaining uncertainty given its own history, in `[0, 1]` (default
  `TRUE`).

- seed:

  Optional seed for the surrogate test (default `NULL`). The session's
  random number state is restored on exit.

## Value

A `data.frame` with one row per ordered pair (`from`, `to`) and the
columns `te` (bits), `te_effective` (with `test = "surrogate"`),
`te_normalised` (with `normalize = TRUE`), `p` (with
`test = "surrogate"`) and `n` (number of transitions used), sorted by
decreasing `te`. In the state-flow mode there is one row per ordered
pair of distinct states; in the bivariate mode two rows, named after the
arguments.

## Conditions

Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`
when `lag`, `history` or `R` is not a positive whole number, or when `x`
and `y` are not aligned (different numbers of sequences or different
lengths within a sequence).

## References

Schreiber, T. (2000). Measuring information transfer. *Physical Review
Letters, 85*(2), 461-464.

## See also

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md),
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md),
[`lsa_lags()`](https://pak.dynasite.org/lagdynamics/reference/lsa_lags.md)

## Examples

``` r
# Information flow among the engagement states
transfer_entropy(engagement, test = "none")
#>         from         to          te te_normalised    n
#> 1     Active Disengaged 0.036448215   0.060023389 1734
#> 2    Average Disengaged 0.036448215   0.060023389 1734
#> 3     Active    Average 0.006949150   0.007517345 1734
#> 4 Disengaged    Average 0.006949150   0.007517345 1734
#> 5    Average     Active 0.004793199   0.006345037 1734
#> 6 Disengaged     Active 0.004793199   0.006345037 1734

# With the surrogate test (an analysis would use more surrogates)
transfer_entropy(engagement, R = 19, seed = 1)
#>         from         to          te te_effective te_normalised    p    n
#> 1     Active Disengaged 0.036448215  0.035496615   0.060023389 0.05 1734
#> 2    Average Disengaged 0.036448215  0.035862001   0.060023389 0.05 1734
#> 3     Active    Average 0.006949150  0.005764179   0.007517345 0.05 1734
#> 4 Disengaged    Average 0.006949150  0.005940226   0.007517345 0.05 1734
#> 5    Average     Active 0.004793199  0.003806708   0.006345037 0.05 1734
#> 6 Disengaged     Active 0.004793199  0.004268635   0.006345037 0.05 1734

# Between two aligned series
mood <- c("calm", "calm", "tense", "tense", "calm", "tense", "tense", "calm")
effort <- c("low", "low", "low", "high", "high", "low", "high", "high")
transfer_entropy(mood, effort, test = "none")
#>     from     to        te te_normalised n
#> 1   mood effort 0.9649839     1.0000000 7
#> 2 effort   mood 0.6792696     0.7039181 7
```
