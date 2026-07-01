# Canonicalize Sequence Input for Lag Sequential Analysis

Coerces a wide variety of user input shapes into a single canonical
representation used by every downstream lagdynamics function (engines,
bootstrap, permutation, grouping, plotting).

## Usage

``` r
lsa_data(x, labels = NULL)
```

## Arguments

- x:

  Sequence input. See Details.

- labels:

  Optional character vector of label names for the states. When `NULL`,
  labels are extracted from the data: unique sorted values of character
  input, or `"Code 1", "Code 2", ...` for integer input.

## Value

An object of class `c("lsa_data", "list")` with elements:

- events:

  Integer vector of event codes (1-indexed), or `NULL` if input was a
  transition matrix.

- seq_id:

  Integer vector of sequence membership, same length as `events`, or
  `NULL` if input was a transition matrix.

- labels:

  Character vector of state labels.

- n_states:

  Number of distinct states (`K`).

- n_sequences:

  Integer count of independent sequences.

- n_events:

  Total number of events across all sequences.

- transitions_per_seq:

  Integer vector: number of transitions each sequence contributes at lag
  1.

- source:

  One of `"events"`, `"transitions"` — flags whether event-level data is
  available.

- obs_input:

  If `source = "transitions"`, the original `K x K` count matrix.
  Otherwise `NULL`.

## Details

Accepted input forms:

- An atomic vector of integer or character codes — treated as a single
  sequence.

- A list of atomic vectors — treated as multiple independent sequences;
  transitions are not counted across sequence boundaries.

- A wide matrix or data.frame with rows = sequences, columns = ordered
  time points. Missing values (`NA`) and empty strings are treated as
  missingness, not as a state: they are dropped wherever they occur in a
  row and the surrounding events close up, so no transition is counted
  into or out of a gap. To model missingness as its own state, recode it
  (e.g. `NA -> "missing"`) before calling
  [`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md).

- A square numeric matrix of pre-computed transition counts. Row `i`,
  column `j` is the count of `i -> j` transitions. In this case `events`
  and `seq_id` are not available and downstream resampling tools that
  need event-level data will error.

- A sequence-bearing object: a `tna` or `group_tna` (sequences read from
  its `$data` slot), a `tna_data` or `nestimate_data`
  (`$sequence_data`), or an `stslist`. The stored event sequences are
  recovered and analysed. A `tna` built from a bare matrix (no retained
  sequences) errors, because transition *counts* cannot be recovered
  from probability weights.

## See also

[`lsa_transitions()`](https://saqr.me/lagdynamics/reference/lsa_transitions.md),
[`lsa()`](https://saqr.me/lagdynamics/reference/lsa.md)

## Examples

``` r
# Single character sequence
d1 <- lsa_data(c("a", "b", "a", "c", "b"))
d1$n_events
#> [1] 5

# Multiple sequences
d2 <- lsa_data(list(c("a", "b", "a"), c("b", "c", "a", "b")))
d2$n_sequences
#> [1] 2

# Pre-computed transition matrix
tm <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3,
             dimnames = list(c("a","b","c"), c("a","b","c")))
d3 <- lsa_data(tm)
d3$source
#> [1] "transitions"
```
