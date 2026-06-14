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
#> [1] 0.112295
fast_median(1:100)
#> [1] 50.5

x <- matrix(rnorm(100), ncol = 2)
fast_mvquantile(x, 0.2)
#> [1] -0.9341112 -0.8614320
fast_mvmedian(x)
#> [1] -0.09334316 -0.15043863

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
#>  fast_median 0.087754 0.1296470 0.1430562 0.1456665 0.155976 0.181409   100
#>  base_median 0.132858 0.1476355 0.1604058 0.1585905 0.167272 0.338601   100

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
#>         expr      min       lq     mean    median       uq      max neval
#>  fast_median 0.685179 0.770308 0.824381 0.8072055 0.828075 1.132993    10
#>  base_median 2.639795 2.708413 2.784651 2.7689810 2.803570 3.002913    10
```
