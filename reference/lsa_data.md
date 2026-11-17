# Convert sequence input to the canonical form used by lag sequential analysis

`lsa_data()` converts the forms in which sequence data usually arrive
into the canonical representation used by lag sequential analysis, and
it is what
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) calls
on its input. A lag sequential analysis needs, for every sequence, the
ordered list of coded events, and a common set of state labels across
sequences. Transitions are counted within sequences only, so the
boundaries between sequences must be kept.

Accepted forms are:

- an atomic vector of codes (one sequence);

- a list of vectors (several sequences);

- a wide matrix or data frame, one row per sequence and one column per
  ordered time point;

- a square numeric matrix of transition counts, row `i` and column `j`
  holding the count of `i -> j` (any square numeric matrix is read this
  way, so pass square integer-coded sequences as a data frame);

- a sequence object: a `tna` or `group_tna` model with its data, a
  `tna_seq_data`, `tna_data` or `nestimate_data` object, a Nestimate
  `netobject`, or a TraMineR `stslist`.

`NA` and empty strings are treated as missing observations, not as a
state: they are dropped wherever they occur and the surrounding events
close up. To model missingness as a state, recode it (for example to
`"missing"`) first. A transition-count matrix carries no sequences, so
the functions that resample sequences
([`bootstrap_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bootstrap_lsa.md),
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md),
[`stability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/stability_lsa.md),
[`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md),
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md))
cannot be used on a fit built from it. Long event logs are read by
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) through
its `actor` and `action` arguments.

## Usage

``` r
lsa_data(x, labels = NULL)
```

## Arguments

- x:

  Sequence input in one of the forms listed above. An `lsa_data` object
  is returned unchanged.

- labels:

  Optional character vector of state names. Default `NULL` takes the
  sorted distinct values of character input, the levels of a factor, or
  `"Code 1"`, `"Code 2"`, ... for integer codes. When supplied,
  character codes are matched to it by name and integer codes index into
  it.

## Value

An object of class `c("lsa_data", "list")`. Print it for the input
source, the number of states, sequences, events and lag-1 transitions;
read the events with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) (one row
per event, or per cell for a transition matrix). The object stores the
integer-coded events and their sequence membership (`NULL` for a
transition matrix), the state labels, the counts of states, sequences
and events, and the source (`"events"` or `"transitions"`).

## Conditions

Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`:
an unsupported input type, no usable events, invalid counts in a
transition matrix, invalid integer codes, or codes that do not match
`labels`.

## See also

[`as.data.frame.lsa_data()`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_data.md)
for the event table;
[`lsa_transitions()`](https://pak.dynasite.org/lagdynamics/reference/lsa_transitions.md)
for transition counts;
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

## Examples

``` r
lsa_data(engagement)
#> <lsa_data>
#>   source:        events
#>   states (K):    3
#>   labels:        Active, Average, Disengaged
#>   sequences:     136
#>   events:        1870
#>   transitions:   1734 (lag 1, summed across sequences)

# Several sequences: no transition is counted across the boundary
lsa_data(list(c("a", "b", "a"), c("b", "c", "a", "b")))
#> <lsa_data>
#>   source:        events
#>   states (K):    3
#>   labels:        a, b, c
#>   sequences:     2
#>   events:        7
#>   transitions:   5 (lag 1, summed across sequences)

# A transition-count matrix
counts <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3,
                 dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
lsa_data(counts)
#> <lsa_data>
#>   source:        transitions
#>   states (K):    3
#>   labels:        a, b, c
#>   transitions:   16 (from input matrix)
```
