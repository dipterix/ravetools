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
#>         expr      min        lq      mean   median        uq      max neval
#>  fast_median 0.084538 0.1181705 0.1337705 0.137166 0.1482115 0.177431   100
#>  base_median 0.115976 0.1386385 0.1485625 0.147044 0.1562810 0.264012   100

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
#>  fast_median 0.752325 0.768042 0.8252348 0.7938805 0.812413 1.174230    10
#>  base_median 2.755808 2.796774 2.8535154 2.8009730 2.823394 3.332814    10
```
