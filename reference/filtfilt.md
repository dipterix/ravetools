# Forward and reverse filter a one-dimensional signal

The result has been tested against 'Matlab' `filtfilt` function.
Currently this function only supports one filter at a time.

## Usage

``` r
filtfilt(b, a = 1, x)
```

## Arguments

- b:

  one-dimensional real numerical vector, the moving-average coefficients
  of an `ARMA` filter; alternatively, a `Sos` (second-order sections)
  object from `gsignal`, in which case `a` is ignored and
  [`gsignal::filtfilt`](https://rdrr.io/pkg/gsignal/man/filtfilt.html)
  is used internally

- a:

  the auto-regressive (recursive) coefficients of an `ARMA` filter;
  ignored when `b` is a `Sos` object

- x:

  numerical vector or matrix input (real value)

## Value

The filtered signal, normally the same length as the input signal `x`.

## Examples

``` r

t <- seq(0, 1, by = 0.01)
x <- sin(2 * pi * t * 2.3)
bf <- gsignal::butter(2, c(0.15, 0.3))

res <- filtfilt(bf$b, bf$a, x)

## Matlab (2022a) equivalent:
# t = [0:0.01:1];
# x = sin(2 * pi * t * 2.3);
# [b,a] = butter(2,[.15,.3]);
# res = filtfilt(b, a, x)
```
