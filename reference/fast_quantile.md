# Compute quantiles

Compute quantiles

## Usage

``` r
fast_quantile(x, prob = 0.5, na.rm = FALSE, ...)

fast_median(x, na.rm = FALSE, ...)

fast_mvquantile(x, prob = 0.5, na.rm = FALSE, ...)

fast_mvmedian(x, na.rm = FALSE, ...)
```

## Arguments

- x:

  numerical-value vector for `fast_quantile` and `fast_median`, and
  column-major matrix for `fast_mvquantile` and `fast_mvmedian`

- prob:

  a probability with value from 0 to 1

- na.rm:

  logical; if true, any `NA` are removed from `x` before the quantiles
  are computed

- ...:

  reserved for future use

## Value

`fast_quantile` and `fast_median` calculate univariate quantiles
(single-value return); `fast_mvquantile` and `fast_mvmedian` calculate
multivariate quantiles (for each column, result lengths equal to the
number of columns).

## Examples

``` r

fast_quantile(runif(1000), 0.1)
#> [1] 0.09972501
fast_median(1:100)
#> [1] 50.5

x <- matrix(rnorm(100), ncol = 2)
fast_mvquantile(x, 0.2)
#> [1] -0.7926140 -0.7404126
fast_mvmedian(x)
#> [1] -0.14570401 -0.09950553

# Compare speed for vectors (usually 30% faster)
x <- rnorm(10000)
microbenchmark::microbenchmark(
  fast_median = fast_median(x),
  base_median = median(x),
  # bioc_median = Biobase::rowMedians(matrix(x, nrow = 1)),
  times = 100, unit = "milliseconds"
)
#> Unit: milliseconds
#>         expr      min        lq      mean    median       uq      max neval
#>  fast_median 0.099739 0.1444555 0.1594043 0.1604035 0.174565 0.230834   100
#>  base_median 0.108441 0.1214410 0.1458746 0.1380050 0.164310 0.296261   100

# Multivariate cases
# (5~7x faster than base R)
# (3~5x faster than Biobase rowMedians)
x <- matrix(rnorm(100000), ncol = 20)
microbenchmark::microbenchmark(
  fast_median = fast_mvmedian(x),
  base_median = apply(x, 2, median),
  # bioc_median = Biobase::rowMedians(t(x)),
  times = 10, unit = "milliseconds"
)
#> Unit: milliseconds
#>         expr      min       lq      mean    median       uq      max neval
#>  fast_median 0.784359 0.803307 0.8515959 0.8273725 0.880853 0.978999    10
#>  base_median 2.920581 2.968462 3.0051758 2.9988375 3.026859 3.158565    10
```
