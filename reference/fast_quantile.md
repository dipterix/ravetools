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
#> [1] 0.1157489
fast_median(1:100)
#> [1] 50.5

x <- matrix(rnorm(100), ncol = 2)
fast_mvquantile(x, 0.2)
#> [1] -0.8980551 -1.1179596
fast_mvmedian(x)
#> [1] -0.3705184 -0.2642868

# Compare speed for vectors (usually 30% faster)
x <- rnorm(10000)
microbenchmark::microbenchmark(
  fast_median = fast_median(x),
  base_median = median(x),
  # bioc_median = Biobase::rowMedians(matrix(x, nrow = 1)),
  times = 100, unit = "milliseconds"
)
#> Unit: milliseconds
#>         expr      min        lq      mean    median        uq      max neval
#>  fast_median 0.095067 0.1304335 0.1406933 0.1419950 0.1514325 0.193831   100
#>  base_median 0.135382 0.1613965 0.1728805 0.1696065 0.1804415 0.309908   100

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
#>         expr      min       lq      mean   median       uq      max neval
#>  fast_median 0.732255 0.747714 0.8307049 0.817494 0.855996 1.056590    10
#>  base_median 2.764945 2.799209 2.8659879 2.821942 2.874389 3.235012    10
```
