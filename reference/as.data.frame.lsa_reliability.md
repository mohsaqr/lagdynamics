# Tidy the Per-Replicate Split-Half Correlations

Returns the split-half correlation of every replicate from
[`reliability_lsa()`](https://pak.dynasite.org/lagdynamics/reference/reliability_lsa.md)
as a tidy data frame.

## Usage

``` r
# S3 method for class 'lsa_reliability'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_reliability_group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `lsa_reliability` (or `lsa_reliability_group`) object.

- row.names, optional, ...:

  Ignored (method signature compatibility).

## Value

A `data.frame`, one row per replicate, with columns `replicate` and
`correlation` (a grouped object gains a leading `group` column). `NA`
correlations from degenerate splits are kept.

## Examples

``` r
fit <- lsa(engagement)
rel <- reliability_lsa(fit, R = 20)
rel
#> <lsa_reliability>
#>   engine:        classical
#>   replicates:    20
#>   weights:       prob
#>   method:        pearson
#>   n sequences:   136
#>   split-half r:  0.965  (sd = 0.033)
#>   95% CI:        [0.889, 0.994]
as.data.frame(rel)
#>    replicate correlation
#> 1          1   0.9813133
#> 2          2   0.9683861
#> 3          3   0.9885036
#> 4          4   0.9873595
#> 5          5   0.9729948
#> 6          6   0.9741100
#> 7          7   0.9725761
#> 8          8   0.9394135
#> 9          9   0.9720272
#> 10        10   0.9462399
#> 11        11   0.9605240
#> 12        12   0.9386932
#> 13        13   0.8443619
#> 14        14   0.9863998
#> 15        15   0.9939222
#> 16        16   0.9939688
#> 17        17   0.9655112
#> 18        18   0.9625765
#> 19        19   0.9930543
#> 20        20   0.9619175

g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
         time = "Time", group = "Achiever")
grel <- reliability_lsa(g, R = 10)
grel
#> <lsa_reliability_group>
#>   groups: 2
#> 
#> [High]
#> <lsa_reliability>
#>   engine:        classical
#>   replicates:    10
#>   weights:       prob
#>   method:        pearson
#>   n sequences:   1000
#>   split-half r:  0.988  (sd = 0.003)
#>   95% CI:        [0.984, 0.992]
#> 
#> [Low]
#> <lsa_reliability>
#>   engine:        classical
#>   replicates:    10
#>   weights:       prob
#>   method:        pearson
#>   n sequences:   1000
#>   split-half r:  0.986  (sd = 0.005)
#>   95% CI:        [0.978, 0.992]
head(as.data.frame(grel))
#>   group replicate correlation
#> 1  High         1   0.9880121
#> 2  High         2   0.9914672
#> 3  High         3   0.9926938
#> 4  High         4   0.9853287
#> 5  High         5   0.9903779
#> 6  High         6   0.9862367
```
