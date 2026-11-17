# Compare transition probabilities between groups with a Bayesian model

`bayes_compare_lsa()` compares the transition probabilities of groups
with a Bayesian model: given the counts, it estimates the plausible
range of the true difference between two groups' transition
probabilities and how precisely that difference is estimated. Under a
Dirichlet-Multinomial model the outgoing transitions of each state are a
multinomial draw whose probabilities have a Dirichlet posterior, so each
transition probability has a Beta posterior in each group. A transition
whose source state is rarely visited gets a wide interval, even when its
estimated probabilities look far apart. A permutation test
([`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md))
instead tests whether a group difference is larger than chance.

For every transition, `bayes_compare_lsa()` computes the posterior mean
difference of the transition probability between two groups (exact), its
credible interval and the probability of direction (the posterior
probability of the more probable sign of the difference), both from
Monte Carlo draws of the Beta posteriors, and a standardised effect
size. A transition is credibly different when its credible interval
excludes zero, the absolute mean difference exceeds `mean_threshold`,
and the interval end nearest zero lies further from zero than
`bound_threshold`; the thresholds keep detectable but negligible
differences out. With more than two groups every pair is compared. No
resampling of sequences is involved, so the comparison is fast.

The model treats transitions as independent draws and compares
transition probabilities only. For a test that respects the dependence
between events of one sequence, or that compares adjusted residuals or
odds ratios, use
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md).
For the posterior of one group's transition probabilities, use
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md).

## Usage

``` r
bayes_compare_lsa(
  x,
  y = NULL,
  prior = 0.5,
  draws = 10000L,
  ci = 0.95,
  mean_threshold = 0.01,
  bound_threshold = 0.001,
  adjust = "none",
  seed = NULL
)

# S3 method for class 'lsa_bayes'
print(x, ...)

# S3 method for class 'lsa_bayes_pairwise'
print(x, ...)
```

## Arguments

- x:

  A grouped fit from `lsa(data, group = )` with two or more groups, or
  the `lsa` fit of the first group. For the
  [`print()`](https://rdrr.io/r/base/print.html) methods, the comparison
  result.

- y:

  The `lsa` fit of the second group when `x` is a single fit (the groups
  are then labelled `"a"` and `"b"`). Default `NULL`.

- prior:

  Dirichlet concentration added to every cell. Default `0.5` (the
  Jeffreys prior); `1` gives a uniform prior.

- draws:

  Number of Monte Carlo draws for the credible intervals. Default
  `10000`.

- ci:

  Credible-interval level. Default `0.95`.

- mean_threshold:

  Minimum absolute posterior mean difference for a transition to be
  credibly different. Default `0.01`.

- bound_threshold:

  Minimum distance from zero of the interval end nearest zero. Default
  `0.001`.

- adjust:

  Multiple-comparison correction applied to the two-sided Bayesian
  p-value `2 * (1 - pd)` across transitions (and across all pairs for
  more than two groups); any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"none"`. It changes the `p_adj` column only: the `significant` flag
  is decided by the credible interval and the thresholds.

- seed:

  Optional seed for the Monte Carlo draws. Default `NULL`. The session's
  random number state is restored on exit.

- ...:

  Not used.

## Value

With two groups, an object of class
`c("lsa_bayes", "lsa_comparison", "list")`; with more than two, an
object of class
`c("lsa_bayes_pairwise", "lsa_comparison_pairwise", "list")`. Read the
per-transition table with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md):
one row per transition (`from -> to` pair), and for more than two groups
one row per transition and pair of groups, with the columns

- `group_a`, `group_b` (more than two groups only): the pair compared;

- `from`, `to`: source and target state;

- `prob_a`, `prob_b`: posterior mean transition probability in each
  group;

- `diff`: `prob_a - prob_b`; `ci_low`, `ci_high`: its credible interval;

- `pd`: probability of direction; `effect_size`: `diff` divided by its
  posterior standard deviation;

- `p_value`: `2 * (1 - pd)` (`NA` for transitions that cannot be
  estimated); `p_adj`: after `adjust`;

- `significant`: credibly different.

`transitions(x, significant = TRUE)` keeps the credibly different
transitions; `sort = "strength"` orders by the absolute difference.
`as.data.frame(x, what = "global")` gives the number of credibly
different transitions per pair of groups (see
[`as.data.frame.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_comparison.md));
[plot()](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md)
draws the differences. [`print()`](https://rdrr.io/r/base/print.html)
shows the settings and the number of credibly different transitions and
returns `x` invisibly.

## Conditions

Errors carry the class `lagdynamics_error` and one of:

- `lagdynamics_needs_events`: a fit was built from a transition matrix.

- `lagdynamics_incompatible_fits`: the groups differ in their states,
  engine, lag or structural zeros.

- `lagdynamics_bad_input`: `adjust` is not a method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html); `x` is
  neither a grouped fit nor an `lsa` fit with a second fit in `y`; `y`
  is given together with a grouped fit; or the grouped fit has fewer
  than two groups.

## References

Johnston, L., & Jendoubi, T. (2026). How delivery mode reshapes resource
engagement: A Bayesian differential network analysis. TNA Workshop 2026.

## See also

[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
and
[`as.data.frame.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/as.data.frame.lsa_comparison.md)
to read the result;
[`plot.lsa_comparison()`](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md)
to draw it;
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
for the permutation test;
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)
for one group's posterior.

## Examples

``` r
g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
         time = "Time", group = "Achiever")
bc <- bayes_compare_lsa(g, adjust = "BH", seed = 1)
bc
#> <lsa_bayes>  (Bayesian Dirichlet-Multinomial comparison)
#>   groups:    High vs Low
#>   prior:     Dirichlet(0.50)  |  draws: 10000  |  CI: 95%
#>   edges:     38 credibly different of 81 compared
transitions(bc, significant = TRUE, sort = "strength") |> head()
#>        from        to    prob_a     prob_b        diff      ci_low     ci_high
#> 1   discuss consensus 0.4241594 0.21434059  0.20981881  0.18162651  0.23797612
#> 2 synthesis consensus 0.5681416 0.38177015  0.18637145  0.10922167  0.26011889
#> 3 synthesis     adapt 0.1433628 0.29986790 -0.15650507 -0.21785545 -0.09399030
#> 4 synthesis      plan 0.1433628 0.02509908  0.11826376  0.07644129  0.16412292
#> 5 consensus   discuss 0.2324452 0.13649499  0.09595018  0.07709173  0.11529950
#> 6   discuss     adapt 0.0241594 0.12010243 -0.09594303 -0.11191197 -0.08030344
#>   pd effect_size p_value p_adj significant
#> 1  1   14.555641       0     0        TRUE
#> 2  1    4.863305       0     0        TRUE
#> 3  1   -4.965182       0     0        TRUE
#> 4  1    5.301216       0     0        TRUE
#> 5  1    9.918048       0     0        TRUE
#> 6  1  -11.770044       0     0        TRUE

# Three courses: every pair
courses <- lsa(group_regulation_long, actor = "Actor", action = "Action",
               time = "Time", group = "Course")
bayes_compare_lsa(courses, draws = 2000, seed = 1) |>
  as.data.frame(what = "global")
#>   group_a group_b n_credible
#> 1       A       B         31
#> 2       A       C         31
#> 3       B       C          4
```
