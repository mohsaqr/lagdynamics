# Engagement states of students over time

`engagement` is a coded dataset, in wide format, of engagement: 138
students, each observed at 15 ordered time points and coded at each as
`Active`, `Average` or `Disengaged`. Missing time points are `NA`. Two
students have no observation at all, so a fit reports 136 sequences; the
remaining sequences have on average about 14 observed states.

With three states and many short sequences, the data show the typical
pattern of a lag sequential analysis: each state tends to persist, and
moves between `Active` and `Disengaged` are rare. They are the main
example of the help pages.

## Usage

``` r
engagement
```

## Format

A data frame with 138 rows (students) and 15 character columns named `1`
to `15` (time points, in order). Values are `"Active"`, `"Average"`,
`"Disengaged"` or `NA`.

## Source

Derived from the synthetic `trajectories` example data of the Nestimate
package (<https://github.com/mohsaqr/Nestimate>, MIT license), with the
states kept and the matrix stored as a data frame.

## See also

[group_regulation](https://pak.dynasite.org/lagdynamics/reference/group_regulation.md)
for a larger set of sequences;
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

## Examples

``` r
fit <- lsa(engagement)
fit
#> Lag Sequential Analysis  -  classical  (lag 1, directed)
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
transitions(fit, significant = TRUE)
#>         from         to lag count  expected       prob   prob_col   adj_res
#> 1     Active     Active   1   459 247.03460 0.69756839 0.70506912  21.66281
#> 2    Average     Active   1   153 281.94983 0.20372836 0.23502304 -12.90604
#> 3 Disengaged     Active   1    39 122.01557 0.12000000 0.05990783 -10.54949
#> 4     Active    Average   1   176 289.53518 0.26747720 0.23066841 -11.31913
#> 5    Average    Average   1   458 330.45732 0.60985353 0.60026212  12.45261
#> 6     Active Disengaged   1    23 121.43022 0.03495441 0.07187500 -12.55693
#> 7 Disengaged Disengaged   1   157  59.97693 0.48307692 0.49062500  15.39047
#>               p    yules_q      kappa    kappa_z      kappa_p      lift  sign
#> 1 4.602765e-104  0.8278779  0.4438330  19.187756 4.684307e-82 1.8580393  over
#> 2  4.161672e-38 -0.6010580 -0.4994042 -14.696262 6.811736e-49 0.5426497 under
#> 3  5.107575e-26 -0.6983917 -0.7069460 -11.462654 2.031922e-30 0.3196313 under
#> 4  1.055115e-29 -0.5335259 -0.4241513 -12.481749 9.390092e-36 0.6078709 under
#> 5  1.353374e-35  0.5530510  0.2274894   9.834808 7.971952e-23 1.3859581  over
#> 6  3.642124e-36 -0.8264798 -0.8271733 -13.412058 5.139269e-41 0.1894092 under
#> 7  1.896612e-53  0.7544102  0.3135556  13.555614 7.339504e-42 2.6176731  over
#>   significant
#> 1        TRUE
#> 2        TRUE
#> 3        TRUE
#> 4        TRUE
#> 5        TRUE
#> 6        TRUE
#> 7        TRUE
```
