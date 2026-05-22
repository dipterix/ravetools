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
#>         expr      min       lq      mean    median        uq      max neval
#>  fast_median 0.101541 0.145736 0.1595494 0.1606135 0.1741735 0.209882   100
#>  base_median 0.104726 0.119688 0.1529826 0.1394170 0.1789760 0.310281   100

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
#>  fast_median 0.712599 0.738888 0.7673243 0.7614815 0.778597 0.862953    10
#>  base_median 2.921902 2.945938 2.9816570 2.9768790 2.988672 3.141378    10
```
