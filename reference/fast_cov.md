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
#>             [,1]        [,2]        [,3]        [,4]
#> [1,]  0.89287065  0.04864466 -0.02881969  0.01654759
#> [2,]  0.04864466  0.94437950 -0.12577546  0.01717759
#> [3,] -0.02881969 -0.12577546  0.87766960 -0.25643174
#> [4,]  0.01654759  0.01717759 -0.25643174  0.91763758

# Calculate covariance of subsets
fast_cov(x, col_x = 1, col_y = 1:2)
#>           [,1]       [,2]
#> [1,] 0.8928707 0.04864466

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
#>  fast_cov 1.469737 1.485350 1.647569 1.503451 1.528674 2.951971    10
#>       cov 6.576223 6.582833 6.595625 6.590279 6.600439 6.653628    10

# }

```
