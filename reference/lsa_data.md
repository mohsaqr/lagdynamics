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

  One of `"events"`, `"transitions"` – flags whether event-level data is
  available.

- obs_input:

  If `source = "transitions"`, the original `K x K` count matrix.
  Otherwise `NULL`.

## Details

Accepted input forms:

- An atomic vector of integer or character codes – treated as a single
  sequence.

- A list of atomic vectors – treated as multiple independent sequences;
  transitions are not counted across sequence boundaries.

- A wide matrix or data.frame with rows = sequences, columns = ordered
  time points. Missing values (`NA`) and empty strings are treated as
  missingness, not as a state: they are dropped wherever they occur in a
  row and the surrounding events close up, so no transition is counted
  into or out of a gap. To model missingness as its own state, recode it
  (e.g. `NA -> "missing"`) before calling
  [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

- A square numeric matrix of pre-computed transition counts. Row `i`,
  column `j` is the count of `i -> j` transitions. In this case `events`
  and `seq_id` are not available and downstream resampling tools that
  need event-level data will error.

- A sequence-bearing object: a `tna` or `group_tna` (sequences read from
  its `data` slot), a `tna_seq_data` code matrix, a `tna_data` or
  `nestimate_data` (its `sequence_data`), a `Nestimate` `netobject` (its
  prepared sequence data), or a `TraMineR` `stslist`. The stored event
  sequences are recovered and analysed. A `tna` built from a bare matrix
  (no retained sequences) errors, because transition *counts* cannot be
  recovered from probability weights.

## See also

[`lsa_transitions()`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md),
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md)

## Examples

``` r
# Single character sequence
lsa_data(c("a", "b", "a", "c", "b"))
#> <lsa_data>
#>   source:        events
#>   states (K):    3
#>   labels:        a, b, c
#>   sequences:     1
#>   events:        5
#>   transitions:   4 (lag 1, summed across sequences)

# Multiple sequences: no transition is counted across the boundary
lsa_data(list(c("a", "b", "a"), c("b", "c", "a", "b")))
#> <lsa_data>
#>   source:        events
#>   states (K):    3
#>   labels:        a, b, c
#>   sequences:     2
#>   events:        7
#>   transitions:   5 (lag 1, summed across sequences)

# Pre-computed transition matrix (source = "transitions")
tm <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3,
             dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
lsa_data(tm)
#> <lsa_data>
#>   source:        transitions
#>   states (K):    3
#>   labels:        a, b, c
#>   transitions:   16 (from input matrix)
```
