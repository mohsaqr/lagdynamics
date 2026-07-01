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
