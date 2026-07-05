# Tidy a Group Comparison

Returns the per-edge comparison table (the same data frame as `x$edges`)
so a comparison can be read with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) like the
other result objects, without reaching into the object.

## Usage

``` r
# S3 method for class 'lsa_comparison'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'lsa_comparison_pairwise'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `lsa_comparison` or `lsa_comparison_pairwise` object.

- row.names, optional, ...:

  Standard
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  arguments (unused; present for method consistency).

## Value

The tidy per-edge data frame.

## Examples

``` r
# \donttest{
grp <- ifelse(group_regulation$T1 == "plan", "plan", "other")
g <- lsa(group_regulation, group = grp)
cmp <- compare_lsa(g, R = 200)
head(as.data.frame(cmp))
#>         from    to  log_or_a    log_or_b        diff    p_perm     p_adj
#> 1      adapt adapt -2.934854 -1.22410616 -1.71074827        NA        NA
#> 2   cohesion adapt -1.878259 -2.41906008  0.54080059 0.5970149 0.5970149
#> 3  consensus adapt -1.773051 -1.53109175 -0.24195921 0.6368159 0.6368159
#> 4 coregulate adapt -0.309429 -0.09747994 -0.21194909 0.6766169 0.6766169
#> 5    discuss adapt  1.882091  1.89742114 -0.01533041 0.9353234 0.9353234
#> 6    emotion adapt -2.247226 -2.33556543  0.08833958 1.0000000 1.0000000
#>   significant
#> 1       FALSE
#> 2       FALSE
#> 3       FALSE
#> 4       FALSE
#> 5       FALSE
#> 6       FALSE
# }
```
