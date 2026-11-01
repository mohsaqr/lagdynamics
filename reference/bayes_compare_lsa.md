# Bayesian Comparison of Group Transition Structures (Dirichlet-Multinomial)

Closed-form Bayesian alternative to
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md)
for comparing the transition structures of two (or, pairwise, more)
groups. Each state's outgoing transitions are modelled as
Dirichlet-Multinomial with a Jeffreys prior, so each transition
probability is marginally Beta. The per-edge posterior mean difference
`prob_a - prob_b` is exact; a credible interval, the probability of
direction `pd`, and a two-sided Bayesian p-equivalent `2 * (1 - pd)`
come from a Monte Carlo draw on the Beta marginals.

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

  An `lsa_bayes` or `lsa_bayes_pairwise` object
  ([`print()`](https://rdrr.io/r/base/print.html) methods; for the
  comparison function itself `x` is the input described above).

- y:

  The second group's `lsa` fit when `x` is a single fit; otherwise
  `NULL`.

- prior:

  Numeric \> 0. Dirichlet prior concentration added to every cell.
  Default `0.5` (Jeffreys). Use `1` for a uniform (Laplace) prior.

- draws:

  Integer. Monte Carlo posterior draws for the credible intervals.
  Default `10000`.

- ci:

  Numeric in (0, 1). Credible-interval mass. Default `0.95`.

- mean_threshold, bound_threshold:

  An edge is flagged credibly different only if its credible interval
  excludes zero, `|posterior mean diff|` exceeds `mean_threshold`
  (default `0.01`), and the credible bound nearest zero exceeds
  `bound_threshold` (default `0.001`). The thresholds guard against
  differences that are detectable but negligibly small.

- adjust:

  Multiple-comparison correction applied to the two-sided Bayesian p
  across edges (and family-wide across pairs); any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"none"`.

- seed:

  Optional integer for reproducible credible intervals.

- ...:

  Ignored by the [`print()`](https://rdrr.io/r/base/print.html) methods.

## Value

For two groups, class `c("lsa_bayes", "lsa_comparison", "list")` with an
`edges` data frame (`from`, `to`, `prob_a`, `prob_b`, `diff`, `ci_low`,
`ci_high`, `pd`, `effect_size`, `p_value`, `p_adj`, `significant`), the
two `fits`, and the Bayesian settings. For more than two groups, an
all-pairwise
`c("lsa_bayes_pairwise", "lsa_comparison_pairwise", "list")` whose
`edges` table gains leading `group_a`, `group_b` columns. Read the
per-edge table with
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md);
[`print()`](https://rdrr.io/r/base/print.html) shows the number of
credibly different edges and returns `x` invisibly.

## Details

This complements
[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md):
the permutation test asks whether a difference is more extreme than
chance; the Bayesian comparison asks for the plausible range of the true
difference and how precisely it is estimated. An edge whose source state
is rarely visited gets a wide credible interval even when its
row-normalised probability looks decisive.

The result carries class `c("lsa_bayes", "lsa_comparison")` (and the
pairwise object `c("lsa_bayes_pairwise", "lsa_comparison_pairwise")`),
so
[plot()](https://pak.dynasite.org/lagdynamics/reference/plot.lsa_comparison.md)
(barrel / heatmap),
[`transitions()`](https://pak.dynasite.org/lagdynamics/reference/transitions.md)
and [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) work
as for a permutation comparison.

## References

Johnston, L. & Jendoubi, T. (2026). How Delivery Mode Reshapes Resource
Engagement: A Bayesian Differential Network Analysis. TNA Workshop 2026.

## See also

[`compare_lsa()`](https://pak.dynasite.org/lagdynamics/reference/compare_lsa.md),
[`certainty_lsa()`](https://pak.dynasite.org/lagdynamics/reference/certainty_lsa.md)

## Examples

``` r
g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
         time = "Time", group = "Achiever")
bc <- bayes_compare_lsa(g, seed = 1)
bc
#> <lsa_bayes>  (Bayesian Dirichlet-Multinomial comparison)
#>   groups:    High vs Low
#>   prior:     Dirichlet(0.50)  |  draws: 10000  |  CI: 95%
#>   edges:     38 credibly different of 81 compared
transitions(bc, significant = TRUE, sort = "strength")
#>          from         to      prob_a      prob_b        diff       ci_low
#> 1     discuss  consensus 0.424159402 0.214340589  0.20981881  0.181626512
#> 2   synthesis  consensus 0.568141593 0.381770145  0.18637145  0.109221672
#> 3   synthesis      adapt 0.143362832 0.299867900 -0.15650507 -0.217855451
#> 4   synthesis       plan 0.143362832 0.025099075  0.11826376  0.076441293
#> 5   consensus    discuss 0.232445164 0.136494986  0.09595018  0.077091731
#> 6     discuss      adapt 0.024159402 0.120102433 -0.09594303 -0.111911969
#> 7        plan       plan 0.327464181 0.414895244 -0.08743106 -0.111244164
#> 8    cohesion  consensus 0.534217507 0.448456993  0.08576051  0.037016381
#> 9    cohesion coregulate 0.081167109 0.166119501 -0.08495239 -0.117402669
#> 10 coregulate    discuss 0.234279354 0.304951411 -0.07067206 -0.110735644
#> 11    discuss  synthesis 0.106351183 0.176440461 -0.07008928 -0.091921063
#> 12  consensus       plan 0.364051229 0.431582526 -0.06753130 -0.091420405
#> 13       plan    emotion 0.182116347 0.115461080  0.06665527  0.048862514
#> 14  synthesis    discuss 0.030088496 0.088507266 -0.05841877 -0.093075092
#> 15 coregulate    emotion 0.203116305 0.145765849  0.05735046  0.024361375
#> 16    discuss    discuss 0.169115816 0.220998720 -0.05188290 -0.076976360
#> 17  synthesis coregulate 0.015929204 0.067371202 -0.05144200 -0.082091862
#> 18 coregulate       plan 0.265442404 0.216103656  0.04933875  0.010704872
#> 19 coregulate  consensus 0.108514190 0.155946321 -0.04743213 -0.077363441
#> 20    emotion    discuss 0.121904149 0.077816492  0.04408766  0.022264929
#> 21   cohesion    discuss 0.040848806 0.083388050 -0.04253924 -0.066575594
#> 22   cohesion   cohesion 0.044031830 0.007222587  0.03680924  0.022903389
#> 23  consensus coregulate 0.170911232 0.206867245 -0.03595601 -0.055349794
#> 24   cohesion    monitor 0.017506631 0.053184504 -0.03567787 -0.054266924
#> 25    emotion    emotion 0.062721132 0.094076655 -0.03135552 -0.051529084
#> 26    discuss   cohesion 0.062017435 0.033034571  0.02898286  0.015692336
#> 27      adapt coregulate 0.003436426 0.030872483 -0.02743606 -0.048510989
#> 28  consensus    monitor 0.035477698 0.059663437 -0.02418574 -0.034915324
#> 29    emotion coregulate 0.023480219 0.047619048 -0.02413883 -0.038226890
#> 30    discuss coregulate 0.072478207 0.096542894 -0.02406469 -0.041605435
#> 31  synthesis    monitor 0.001769912 0.022457067 -0.02068716 -0.038489431
#> 32  consensus    emotion 0.081407331 0.062723100  0.01868423  0.005764843
#> 33 coregulate coregulate 0.013912076 0.031929662 -0.01801759 -0.031037910
#> 34       plan    discuss 0.060245123 0.074782077 -0.01453695 -0.026946314
#> 35       plan coregulate 0.023994476 0.011469644  0.01252483  0.005853884
#> 36       plan   cohesion 0.031589850 0.019727787  0.01186206  0.004050320
#> 37    discuss    monitor 0.016687422 0.028425096 -0.01173767 -0.021247504
#> 38  consensus   cohesion 0.019873399 0.009348972  0.01052443  0.004667917
#>         ci_high     pd effect_size p_value  p_adj significant
#> 1   0.237976123 1.0000   14.555641  0.0000 0.0000        TRUE
#> 2   0.260118895 1.0000    4.863305  0.0000 0.0000        TRUE
#> 3  -0.093990303 1.0000   -4.965182  0.0000 0.0000        TRUE
#> 4   0.164122917 1.0000    5.301216  0.0000 0.0000        TRUE
#> 5   0.115299504 1.0000    9.918048  0.0000 0.0000        TRUE
#> 6  -0.080303439 1.0000  -11.770044  0.0000 0.0000        TRUE
#> 7  -0.063876098 1.0000   -7.205805  0.0000 0.0000        TRUE
#> 8   0.132777794 0.9996    3.493701  0.0008 0.0008        TRUE
#> 9  -0.053961240 1.0000   -5.232337  0.0000 0.0000        TRUE
#> 10 -0.032544723 0.9998   -3.562048  0.0004 0.0004        TRUE
#> 11 -0.048648906 1.0000   -6.346308  0.0000 0.0000        TRUE
#> 12 -0.043574916 1.0000   -5.498556  0.0000 0.0000        TRUE
#> 13  0.084388132 1.0000    7.330156  0.0000 0.0000        TRUE
#> 14 -0.023975759 0.9995   -3.323002  0.0010 0.0010        TRUE
#> 15  0.090668818 0.9999    3.370883  0.0002 0.0002        TRUE
#> 16 -0.027176250 1.0000   -4.109155  0.0000 0.0000        TRUE
#> 17 -0.023199110 0.9998   -3.448102  0.0004 0.0004        TRUE
#> 18  0.087584857 0.9933    2.529063  0.0134 0.0134        TRUE
#> 19 -0.017476324 0.9996   -3.114592  0.0008 0.0008        TRUE
#> 20  0.065536818 1.0000    3.993001  0.0000 0.0000        TRUE
#> 21 -0.019744211 0.9999   -3.567142  0.0002 0.0002        TRUE
#> 22  0.051834469 1.0000    4.984778  0.0000 0.0000        TRUE
#> 23 -0.016338973 0.9998   -3.610483  0.0004 0.0004        TRUE
#> 24 -0.018722298 1.0000   -3.873508  0.0000 0.0000        TRUE
#> 25 -0.011232325 0.9988   -3.070456  0.0024 0.0024        TRUE
#> 26  0.042077508 1.0000    4.322233  0.0000 0.0000        TRUE
#> 27 -0.007643052 0.9922   -2.658680  0.0156 0.0156        TRUE
#> 28 -0.013811984 1.0000   -4.491818  0.0000 0.0000        TRUE
#> 29 -0.010385772 0.9999   -3.384110  0.0002 0.0002        TRUE
#> 30 -0.006945683 0.9972   -2.703052  0.0056 0.0056        TRUE
#> 31 -0.007071379 0.9978   -2.568484  0.0044 0.0044        TRUE
#> 32  0.031246248 0.9970    2.869829  0.0060 0.0060        TRUE
#> 33 -0.004977067 0.9961   -2.742499  0.0078 0.0078        TRUE
#> 34 -0.001795777 0.9870   -2.279701  0.0260 0.0260        TRUE
#> 35  0.019372937 0.9999    3.637226  0.0002 0.0002        TRUE
#> 36  0.020008912 0.9979    2.914653  0.0042 0.0042        TRUE
#> 37 -0.002685564 0.9930   -2.472852  0.0140 0.0140        TRUE
#> 38  0.016330192 1.0000    3.538363  0.0000 0.0000        TRUE
```
