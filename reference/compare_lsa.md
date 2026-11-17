# Compare transition structures between groups with a permutation test

`compare_lsa()` tests whether groups of sequences (high and low
achievers, courses, conditions) differ in their transition structure,
transition by transition and as a whole, with a permutation test. Groups
can differ in which transitions they over- or under-produce. The
permutation test checks whether an observed difference between groups is
larger than would arise if group membership were unrelated to the
sequences: the group labels of whole sequences are shuffled many times,
the models are refitted, and the observed difference is set against the
differences produced by the shuffles. Because whole sequences are
permuted, the dependence between events of one sequence is kept and only
the link between group and transition structure is tested.

For each transition, `compare_lsa()` compares a per-group measure, by
default the log odds ratio of the transition against all other
transitions, which does not grow with the number of transitions
observed, so groups of different sizes can be compared. Each transition
gets the group difference, a permutation p-value (Phipson & Smyth, 2010)
and an adjusted p-value; transitions too rare to be compared are left
untested. An omnibus test of overall difference, the sum of squared
differences over the tested transitions, is computed from the same
permutations, in the spirit of the network comparison test of van
Borkulo et al. (2023). With more than two groups every pair is compared
and the `adjust` correction is applied once across the p-values of all
pairs.

The groups must share their states, engine, lag and structural zeros;
fitting them together with `lsa(data, group = )` ensures this. For a
Bayesian comparison of transition probabilities that reports the
plausible size of each difference without resampling, use
[`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md).
A difference between groups can also arise from a difference in their
composition; check subgroups when such a confound is plausible.

## Usage

``` r
compare_lsa(
  x,
  y = NULL,
  R = 1000L,
  measure = c("log_or", "adj_res", "yules_q", "prob", "count", "lift"),
  adjust = "none",
  min_count = 5L,
  parallel = FALSE,
  n_cores = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A grouped fit from `lsa(data, group = )` with two or more groups, or
  the `lsa` fit of the first group.

- y:

  The `lsa` fit of the second group when `x` is a single fit (the groups
  are then labelled `"a"` and `"b"`). Default `NULL`.

- R:

  Number of label permutations per comparison. Default `1000`.

- measure:

  Per-transition measure compared between groups: `"log_or"` (default;
  log odds ratio of the transition, with 0.5 added to every cell of its
  2 x 2 table when one cell is empty), `"yules_q"` (also independent of
  sample size, but bounded at -1 and 1), `"adj_res"` (adjusted residual;
  grows with the number of transitions, so a message is shown when one
  group has more than 1.5 times the transitions of another), `"prob"`
  (transition probability), `"count"` or `"lift"` (observed / expected
  count).

- adjust:

  Multiple-comparison correction, any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) such as
  `"BH"`, `"holm"` or `"bonferroni"`. Default `"none"`.

- min_count:

  Minimum combined count of a transition in the two groups for it to be
  tested. Default `5`; `0` tests every transition.

- parallel:

  Run the permutations on several cores. Default `FALSE`.

- n_cores:

  Number of cores when `parallel = TRUE`. Default `NULL`, all cores but
  one.

- verbose:

  Report progress. Default `FALSE`.

- ...:

  Not used.

## Value

With two groups, an object of class `c("lsa_comparison", "list")`; with
more than two, an object of class
`c("lsa_comparison_pairwise", "list")`. Read the per-transition table
with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md):
one row per transition (`from -> to` pair), and for more than two groups
one row per transition and pair of groups, with the columns

- `group_a`, `group_b` (more than two groups only): the pair compared;

- `from`, `to`: source and target state;

- `<measure>_a`, `<measure>_b` (for example `log_or_a`): the measure in
  each group;

- `diff`: the difference, group a minus group b;

- `p_perm`: permutation p-value (`NA` when not tested);

- `p_adj`: p-value after `adjust`;

- `significant`: `p_adj` below the `alpha` of the fits.

`transitions(x, significant = TRUE)` keeps the significant transitions;
`sort = "strength"` orders by the absolute difference. The omnibus test
is read with `as.data.frame(x, what = "global")` (one row per pair of
groups). [`print()`](https://rdrr.io/r/base/print.html) shows the
groups, the number of significant transitions and the omnibus test (see
[`as.data.frame.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_comparison.md));
[plot()](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md)
draws the differences. The object also stores the permuted differences,
the settings and the group fits.

## Conditions

Errors carry the class `lagdynamics_error` and one of:

- `lagdynamics_needs_events`: a fit was built from a transition matrix,
  so there are no sequences to permute.

- `lagdynamics_incompatible_fits`: the groups differ in their states,
  engine, lag or structural zeros.

- `lagdynamics_bad_input`: `adjust` is not a method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html); `x` is
  neither a grouped fit nor an `lsa` fit with a second fit in `y`; `y`
  is given together with a grouped fit; or the grouped fit has fewer
  than two groups.

## References

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: Calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology, 9*(1),
Article 39.

van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio,
P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2023). Comparing
network structures on three aspects: A permutation test. *Psychological
Methods, 28*(6), 1273-1285.

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
and
[`as.data.frame.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_comparison.md)
to read the result;
[`plot.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md)
to draw it;
[`bayes_compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/bayes_compare_lsa.md)
for a Bayesian comparison;
[`permute_lsa()`](https://pak.dynasite.org/lagdynamics/reference/permute_lsa.md)
to test one fit against no sequential structure;
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) with
`group =` to fit the groups.

## Examples

``` r
g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
         time = "Time", group = "Achiever")
cmp <- compare_lsa(g, R = 100, adjust = "BH")
cmp
#> <lsa_comparison>
#>   groups:   High vs Low
#>   measure:  log_or difference (High - Low)
#>   R:        100 label permutations
#>   edges:    36 significant of 78 tested (adjust = BH)
#>   omnibus:  statistic = 79.48, p = 0.009901
transitions(cmp, significant = TRUE, sort = "strength") |> head()
#>        from        to   log_or_a  log_or_b      diff     p_perm      p_adj
#> 1  cohesion     adapt -0.7935888 -3.920901  3.127313 0.00990099 0.02340234
#> 2   emotion synthesis -1.6186504 -4.521535  2.902885 0.00990099 0.02340234
#> 3      plan synthesis -2.1452241 -4.950504  2.805280 0.00990099 0.02340234
#> 4  cohesion synthesis -1.3553011 -3.944575  2.589274 0.00990099 0.02340234
#> 5 synthesis   monitor -3.2285363 -0.974550 -2.253986 0.00990099 0.02340234
#> 6 synthesis      plan -0.5805673 -2.680438  2.099871 0.00990099 0.02340234
#>   significant
#> 1        TRUE
#> 2        TRUE
#> 3        TRUE
#> 4        TRUE
#> 5        TRUE
#> 6        TRUE
as.data.frame(cmp, what = "global")
#>   group_a group_b statistic    p_value   R
#> 1    High     Low  79.47723 0.00990099 100
# \donttest{
# A full analysis uses more permutations
compare_lsa(g, R = 1000, adjust = "BH")
#> <lsa_comparison>
#>   groups:   High vs Low
#>   measure:  log_or difference (High - Low)
#>   R:        1000 label permutations
#>   edges:    39 significant of 78 tested (adjust = BH)
#>   omnibus:  statistic = 79.48, p = 0.000999

# Three courses: every pair, corrected across all pairs
courses <- lsa(group_regulation_long, actor = "Actor", action = "Action",
               time = "Time", group = "Course")
compare_lsa(courses, R = 100, adjust = "BH") |>
  as.data.frame(what = "global")
#>   group_a group_b statistic    p_value      p_adj
#> 1       A       B 57.753343 0.00990099 0.01485149
#> 2       A       C 55.297067 0.00990099 0.01485149
#> 3       B       C  8.364228 0.75247525 0.75247525
# }
```
