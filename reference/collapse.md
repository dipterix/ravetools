# Collapse array

Collapse array

## Usage

``` r
collapse(x, keep, ...)

# S3 method for class 'array'
collapse(
  x,
  keep,
  average = TRUE,
  transform = c("asis", "10log10", "square", "sqrt"),
  ...
)
```

## Arguments

- x:

  A numeric multi-mode tensor (array), without `NA`

- keep:

  Which dimension to keep

- ...:

  passed to other methods

- average:

  collapse to sum or mean

- transform:

  transform on the data before applying collapsing; choices are `'asis'`
  (no change), `'10log10'` (used to calculate decibel), `'square'`
  (sum-squared), `'sqrt'` (square-root and collapse)

## Value

a collapsed array with values to be mean or summation along collapsing
dimensions

## Examples

``` r

# Set ncores = 2 to comply to CRAN policy. Please don't run this line
ravetools_threads(n_threads = 2L)

# Example 1
x = matrix(1:16, 4)

# Keep the first dimension and calculate sums along the rest
collapse(x, keep = 1)
#> [1]  7  8  9 10
rowMeans(x)  # Should yield the same result
#> [1]  7  8  9 10

# Example 2
x = array(1:120, dim = c(2,3,4,5))
result = collapse(x, keep = c(3,2))
compare = apply(x, c(3,2), mean)
sum(abs(result - compare)) # The same, yield 0 or very small number (1e-10)
#> [1] 5.684342e-14


# \donttest{

ravetools_threads(n_threads = -1)

# Example 3 (performance)

# Small data, no big difference
x = array(rnorm(240), dim = c(4,5,6,2))
microbenchmark::microbenchmark(
  result = collapse(x, keep = c(3,2)),
  compare = apply(x, c(3,2), mean),
  times = 1L, check = function(v) {
    max(abs(range(do.call('-', v)))) < 1e-10
  }
)
#> Unit: microseconds
#>     expr     min      lq    mean  median      uq     max neval
#>   result 406.459 406.459 406.459 406.459 406.459 406.459     1
#>  compare 312.003 312.003 312.003 312.003 312.003 312.003     1

# large data big difference
x = array(rnorm(prod(300,200,105)), c(300,200,105,1))
microbenchmark::microbenchmark(
  result = collapse(x, keep = c(3,2)),
  compare = apply(x, c(3,2), mean),
  times = 1L , check = function(v) {
    max(abs(range(do.call('-', v)))) < 1e-10
  })
#> Unit: milliseconds
#>     expr       min        lq      mean    median        uq       max neval
#>   result  20.45195  20.45195  20.45195  20.45195  20.45195  20.45195     1
#>  compare 178.92030 178.92030 178.92030 178.92030 178.92030 178.92030     1

# }
```
