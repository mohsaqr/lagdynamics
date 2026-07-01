# Student Engagement Trajectories

Wide-format categorical sequence data: 138 students observed over 15
weekly time points. Each row is one student; each column is a week.
Entries are the student's engagement state for that week, one of
`"Active"`, `"Average"`, `"Disengaged"`, or `NA` for missing weeks.

## Usage

``` r
engagement
```

## Format

A character matrix with 138 rows and 15 columns.

## Source

Derived without modification from the `trajectories` matrix in the
`Nestimate` package (<https://github.com/mohsaqr/Nestimate>), which is
MIT-licensed and produced by Saqr and collaborators as a synthetic
engagement trajectory example. Re-shipped here for convenience and
offline testing; both attribution and license are preserved.

## Details

This is a standard small-K, multi-sequence example for lag sequential
analysis: K = 3 states, S = 138 sequences, mean sequence length about
15. It exercises the wide-matrix input path of
[`lsa_data()`](https://mohsaqr.github.io/lagdynamics/reference/lsa_data.md)
and produces a stable transition pattern with clear adjusted-residual
signals.

## Examples

``` r
fit <- lsa(engagement, engine = "classical")
fit
#> Lag Sequential Analysis  —  classical  (lag 1, directed)
#>   3 states | 1734 transitions | 1870 events | 136 sequences
#>   states: Active, Average, Disengaged
#>   independence: G² = 618.3, df = 4, p <2e-16
#> 
#>   Significant transitions (p < 0.05): 7 of 9
#>   strongest over-represented (of 3):
#>     Active -> Active          z =  +21.7  ***
#>     Disengaged -> Disengaged  z =  +15.4  ***
#>     Average -> Average        z =  +12.5  ***
#> 
#>   Initial states:
#>     Active     0.382  ████████████████████████
#>     Average    0.368  ███████████████████████
#>     Disengaged 0.250  ████████████████
fit$adj_res
#>               Active    Average  Disengaged
#> Active      21.66281 -11.319132 -12.5569254
#> Average    -12.90604  12.452615   0.1757996
#> Disengaged -10.54949  -1.736461  15.3904661
```
