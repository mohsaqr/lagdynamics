# Regulation actions in collaborative learning, as sequences

`group_regulation` is a coded dataset, in wide format, of regulation
during collaborative learning: each row is one sequence of coded actions
and each column an ordered position. The nine actions are `adapt`,
`cohesion`, `consensus`, `coregulate`, `discuss`, `emotion`, `monitor`,
`plan` and `synthesis`. Sequences are 2 to 26 actions long; shorter
sequences are padded with `NA`.

With nine states and 2,000 sequences, the data give a dense transition
table in which most transitions depart from independence, and they are
used to show lag transition networks and group comparisons. The same
actions are available as an event log in
[group_regulation_long](https://pak.dynasite.org/lagdynamics/reference/group_regulation_long.md).

## Usage

``` r
group_regulation
```

## Format

A data frame with 2,000 rows (sequences) and 26 character columns `T1`
to `T26` (positions, in order).

## Source

Identical to
[`tna::group_regulation`](https://sonsoles.me/tna/reference/group_regulation.html)
(tna package, MIT license).

## See also

[group_regulation_long](https://pak.dynasite.org/lagdynamics/reference/group_regulation_long.md)
for the event log;
[engagement](https://pak.dynasite.org/lagdynamics/reference/engagement.md);
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md).

## Examples

``` r
fit <- lsa(group_regulation)
transitions(fit, direction = "over", sort = "strength") |> head()
#>        from         to lag count   expected      prob  prob_col  adj_res
#> 1   emotion   cohesion   1   923  190.88889 0.3253437 0.5372526 58.19539
#> 2   discuss  synthesis   1   557  106.77124 0.1409770 0.8072464 48.04631
#> 3 synthesis      adapt   1   153   13.55939 0.2346626 0.2881356 38.76588
#> 4 consensus coregulate   1  1188  519.29875 0.1877074 0.5670644 35.31581
#> 5 consensus       plan   1  2505 1540.29711 0.3957971 0.4031220 32.58405
#> 6      plan       plan   1  2304 1498.43724 0.3742082 0.3707757 27.46361
#>               p   yules_q     kappa  kappa_z       kappa_p      lift sign
#> 1  0.000000e+00 0.8599915 0.4392796 54.99496  0.000000e+00  4.835274 over
#> 2  0.000000e+00 0.9271844 0.7207889 46.05645  0.000000e+00  5.216761 over
#> 3  0.000000e+00 0.9041940 0.2564871 36.97961 2.435758e-299 11.283690 over
#> 4 3.358380e-273 0.6467473 0.4117397 34.57993 5.061358e-262  2.287700 over
#> 5 6.901601e-233 0.4647660 0.1744181 28.44885 5.034845e-178  1.626310 over
#> 6 4.779402e-166 0.4056979 0.1413267 23.45002 1.321402e-121  1.537602 over
#>   significant
#> 1        TRUE
#> 2        TRUE
#> 3        TRUE
#> 4        TRUE
#> 5        TRUE
#> 6        TRUE
```
