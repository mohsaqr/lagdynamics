# Regulation actions in collaborative learning, as an event log

`group_regulation_long` is a coded dataset, in long format, of
regulation during collaborative learning: each row is one coded action
of one student, with the student's collaboration group, course and
achievement level and the time of the action. It contains the nine
actions of
[group_regulation](https://pak.dynasite.org/lagdynamics/reference/group_regulation.md)
(`adapt`, `cohesion`, `consensus`, `coregulate`, `discuss`, `emotion`,
`monitor`, `plan`, `synthesis`) and serves to show how
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) reads
an event log and fits and compares groups.

The log holds 27,533 actions of 2,000 students in 200 collaboration
groups of 10 students each. Achievement (`High` or `Low`, 1,000 students
each) is constant within a student, and the course (`A`, `B` or `C`) is
constant within a group, so either can be passed to the `group` argument
of [`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) for
a between-group comparison.

## Usage

``` r
group_regulation_long
```

## Format

A data frame with 27,533 rows (actions) and 6 columns:

- `Actor`:

  Integer. Student identifier (2,000 students).

- `Achiever`:

  Character. Achievement level, `High` or `Low`.

- `Group`:

  Numeric. Collaboration group (200 groups).

- `Course`:

  Character. Course, `A`, `B` or `C`.

- `Time`:

  POSIXct. Time of the action.

- `Action`:

  Character. Regulation action (9 values).

## Source

Taken without modification from
[`tna::group_regulation_long`](https://sonsoles.me/tna/reference/group_regulation_long.html)
(tna package, MIT license).

## See also

[group_regulation](https://pak.dynasite.org/lagdynamics/reference/group_regulation.md)
for the wide sequence table of the same actions;
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md);
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md).

## Examples

``` r
fit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
           time = "Time")
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

# One model per achievement level
gfit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
            time = "Time", group = "Achiever")
summary(gfit)
#>   group    engine lag n_states n_sequences n_events n_transitions n_significant
#> 1  High classical   1        9        1000    13721         12721            66
#> 2   Low classical   1        9        1000    13812         12812            67
#>   alpha     lrx2 lrx2_df lrx2_p       x2 x2_df x2_p
#> 1  0.05 6146.793      64      0 6806.682    64    0
#> 2  0.05 7714.714      64      0 9197.143    64    0
```
