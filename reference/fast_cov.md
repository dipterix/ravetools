# Calculate massive covariance matrix in parallel

Speed up covariance calculation for large matrices. The default behavior
is the same as [`cov`](https://rdrr.io/r/stats/cor.html) (`'pearson'`,
no `NA` handling).

## Usage

``` r
fast_cov(x, y = NULL, col_x = NULL, col_y = NULL, df = NA)
```

## Arguments

- x:

  a numeric vector, matrix or data frame; a matrix is highly recommended
  to maximize the performance

- y:

  NULL (default) or a vector, matrix or data frame with compatible
  dimensions to x; the default is equivalent to `y = x`

- col_x:

  integers indicating the subset indices (columns) of `x` to calculate
  the covariance, or `NULL` to include all the columns; default is
  `NULL`

- col_y:

  integers indicating the subset indices (columns) of `y` to calculate
  the covariance, or `NULL` to include all the columns; default is
  `NULL`

- df:

  a scalar indicating the degrees of freedom; default is `nrow(x)-1`

## Value

A covariance matrix of `x` and `y`. Note that there is no `NA` handling.
Any missing values will lead to `NA` in the resulting covariance
matrices.

## Examples

``` r
# Set ncores = 2 to comply to CRAN policy. Please don't run this line
ravetools_threads(n_threads = 2L)

x <- matrix(rnorm(400), nrow = 100)

# Call `cov(x)` to compare
fast_cov(x)
#>             [,1]         [,2]         [,3]         [,4]
#> [1,]  1.11343789 -0.018191052  0.111347455 -0.016289440
#> [2,] -0.01819105  1.110481918  0.007465818 -0.009387727
#> [3,]  0.11134745  0.007465818  0.734855745 -0.137734199
#> [4,] -0.01628944 -0.009387727 -0.137734199  1.281456246

# Calculate covariance of subsets
fast_cov(x, col_x = 1, col_y = 1:2)
#>          [,1]        [,2]
#> [1,] 1.113438 -0.01819105

# \donttest{

# Speed comparison, better to use multiple cores (4, 8, or more)
# to show the differences.

ravetools_threads(n_threads = -1)
x <- matrix(rnorm(100000), nrow = 1000)
microbenchmark::microbenchmark(
  fast_cov = {
    fast_cov(x, col_x = 1:50, col_y = 51:100)
  },
  cov = {
    cov(x[,1:50], x[,51:100])
  },
  unit = 'ms', times = 10
)
#> Unit: milliseconds
#>      expr      min       lq     mean   median       uq      max neval
#>  fast_cov 1.571815 1.597952 1.656024 1.643054 1.700178 1.776818    10
#>       cov 5.448677 5.485155 5.732695 5.546323 5.930236 6.719019    10

# }

```
